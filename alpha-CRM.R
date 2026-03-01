# setwd('/rsrch8/home/biostatistics/syang10/IPDE-Project')
# .libPaths("/rsrch8/home/biostatistics/syang10/R/x86_64-pc-linux-gnu-library/4.4")

library(coda)
library(rjags)
library(mvtnorm)
build_temp_rows <- function(patient_df, t_now, window) {
  
  # only administrations that already occurred
  tmp <- patient_df[patient_df$arrival_time <= t_now, , drop=FALSE]
  
  # 1) has the DLT already been observed?
  tmp$y_obs <- as.integer(tmp$DLT_time <= t_now)
  
  # 2) is the toxicity outcome now known?
  # known if either DLT occurred OR window finished
  tmp$complete <- as.integer(tmp$y_obs == 1L | tmp$eval_time <= t_now)
  
  # keep only rows whose toxicity status is known
  tmp <- tmp[tmp$complete == 1, , drop=FALSE]
  
  return(tmp)
}
estimate_MTD_alphaCRM_JAGS <- function(tmp,
                                       d_grid, s_grid,
                                       TARGET = 0.3,
                                       cutoff = 0.95,
                                       T = 28,
                                       model_file = "alpha-crm.bug",
                                       n.chains = 3,
                                       n.adapt  = 1000,
                                       n.burn   = 2000,
                                       n.iter   = 5000,
                                       thin     = 2,
                                       monitor = c("theta","alpha")) {
  
  # If no data, return prior-based estimate (same as your estimate_MTD_prior logic)
  if (is.null(tmp) || nrow(tmp) == 0L) {
    theta <- rnorm(n.iter, 0, sqrt(2))  # match dnorm(0,0.5)
    posttox <- vapply(seq_along(s_grid), \(j) mean(s_grid[j]^exp(theta)), 0.0)
    p1_draws <- s_grid[1]^exp(theta)
    prob_overtox <- mean(p1_draws > TARGET)
    stop_flag <- as.integer(prob_overtox > cutoff)
    mtd <- which.min(abs(posttox - TARGET))
    return(list(MTD=mtd, posttox=posttox, prob_overtox=prob_overtox, stop=stop_flag,
                theta_draws=theta, alpha_draws=NA))
  }
  
  data_jags <- make_alpha_crm_jags_data(tmp, d_grid, s_grid, T = T)
  
  # JAGS model
  jags <- rjags::jags.model(
    file     = model_file,
    data     = data_jags,
    n.chains = n.chains,
    n.adapt  = n.adapt,
    quiet    = TRUE
  )
  
  update(jags, n.burn, progress.bar = "none")
  
  smp <- rjags::coda.samples(
    model          = jags,
    variable.names = monitor,
    n.iter         = n.iter,
    thin           = thin,
    progress.bar   = "none"
  )
  
  draws <- as.matrix(smp)
  
  if (!("theta" %in% colnames(draws))) stop("theta not found in posterior draws.")
  theta <- draws[, "theta"]
  
  # α is optional
  alpha <- if ("alpha" %in% colnames(draws)) draws[, "alpha"] else NA
  
  # ---- Decision rule (paper): use theta draws with ORIGINAL skeleton s_grid
  K <- length(s_grid)
  posttox <- vapply(seq_len(K), function(k) mean(s_grid[k]^exp(theta)), numeric(1))
  
  pi1_draws <- s_grid[1]^exp(theta)
  prob_overtox <- mean(pi1_draws > TARGET)
  stop_flag <- as.integer(prob_overtox > cutoff)
  
  diff <- abs(posttox - TARGET)
  dose.best <- which(diff == min(diff))[1]
  
  list(
    MTD = dose.best,
    posttox = posttox,
    prob_overtox = prob_overtox,
    stop = stop_flag,
    theta_draws = theta,
    alpha_draws = alpha
  )
}
make_alpha_crm_jags_data <- function(tmp, dgrid, sgrid, T) {
  # tmp must contain only rows you want in the likelihood (completed)
  stopifnot(all(c("id","cycle","dose","arrival_time","y_obs") %in% names(tmp)))
  
  ids <- sort(unique(tmp$id))
  I <- length(ids)
  
  # cycles per patient
  ncyc <- vapply(ids, function(pid) max(tmp$cycle[tmp$id == pid]), integer(1))
  Jmax <- max(ncyc)
  
  dose_mat <- matrix(0, nrow = I, ncol = Jmax)
  time_mat <- matrix(0, nrow = I, ncol = Jmax)
  y_mat    <- matrix(0, nrow = I, ncol = Jmax)
  use_mat  <- matrix(0, nrow = I, ncol = Jmax)
  
  for (ii in seq_along(ids)) {
    pid <- ids[ii]
    sub <- tmp[tmp$id == pid, ]
    sub <- sub[order(sub$cycle), ]
    
    for (r in seq_len(nrow(sub))) {
      j <- sub$cycle[r]
      dose_mat[ii, j] <- sub$dose[r]           # physical dose
      time_mat[ii, j] <- sub$arrival_time[r]  # admin time
      y_mat[ii, j]    <- sub$y_obs[r]
      use_mat[ii, j]  <- 1L
    }
  }
  
  list(
    I = I,
    Jmax = Jmax,
    ncyc = as.integer(ncyc),
    dose = dose_mat,
    time = time_mat,
    y = y_mat,
    use = use_mat,
    K = length(dgrid),
    dgrid = as.numeric(dgrid),
    sgrid = as.numeric(sgrid),
    T = as.numeric(T)
  )
}
# ============================================================
# Deterministic alpha-CRM via alpha-grid + integrate() over theta
# Includes overtoxicity stopping rule:
#   stop if  Pr( p1(theta) > TARGET | data ) > cutoff
# where p1(theta) = s1^{exp(theta)}  (since S(d1)=s1)
# ============================================================

# ---- 1) Utilities: S(D) interpolation from (d_grid, s_grid) ----
interp_S_vec <- function(D, d_grid, s_grid) {
  K <- length(d_grid)
  out <- numeric(length(D))
  
  idx_lo <- D <= d_grid[1]
  out[idx_lo] <- s_grid[1]
  
  idx_hi <- D >= d_grid[K]
  if (any(idx_hi)) {
    t <- (D[idx_hi] - d_grid[K-1]) / (d_grid[K] - dgrid[K-1])
    out[idx_hi] <- pmin(1, s_grid[K-1] + t * (s_grid[K] - s_grid[K-1]))
  }
  
  idx_mid <- !(idx_lo | idx_hi)
  if (any(idx_mid)) {
    Dm <- D[idx_mid]
    k <- findInterval(Dm, d_grid)
    k <- pmax(1, pmin(k, K-1))
    t <- (Dm - d_grid[k]) / (d_grid[k+1] - d_grid[k])
    out[idx_mid] <- s_grid[k] + t * (s_grid[k+1] - s_grid[k])
  }
  out
}

compute_Dij_alpha <- function(tmp, alpha, T,  d_grid) {
  tmp <- tmp[order(tmp$id, tmp$cycle), ]
  Dij <- numeric(nrow(tmp))
  idx_by_id <- split(seq_len(nrow(tmp)), tmp$id)
  
  for (idx in idx_by_id) {
    sub <- tmp[idx, , drop = FALSE]
    for (r in seq_len(nrow(sub))) {
      t_j <- sub$arrival_time[r]
      s <- 0
      for (rp in seq_len(r)) {
        s <- s + d_grid[sub$dose[rp]] * alpha^((t_j - sub$arrival_time[rp]) / T)
      }
      Dij[idx[r]] <- s
    }
  }
  Dij
}

# vector-safe log posterior in theta (theta may be vector)
logpost_theta_given_alpha_vec <- function(theta, y, S_at_D, eps = 1e-12) {
  lp <- dnorm(theta, 0, sqrt(2), log = TRUE)
  et <- exp(theta)
  ll <- vapply(et, function(et1) {
    p <- S_at_D^et1
    p <- pmin(1 - eps, pmax(eps, p))
    sum(y * log(p) + (1 - y) * log1p(-p))
  }, numeric(1))
  lp + ll
}

estimate_MTD_alphaCRM_integrate <- function(
    tmp,
    d_grid, s_grid,
    TARGET = 0.30,
    cutoff = 0.95,
    T = 28,
    eps = 1e-12,
    alpha_grid = seq(0.01, 0.99, length.out = 61),
    L = 8,
    rel.tol = 1e-8,
    n_draw_prior = 5000
) {
  
  if (is.null(tmp) || nrow(tmp) == 0L) {
    theta <- rnorm(n_draw_prior, 0, sqrt(2))
    posttox <- vapply(seq_along(s_grid), \(k) mean(s_grid[k]^exp(theta)), 0.0)
    prob_overtox <- mean(s_grid[1]^exp(theta) > TARGET)
    stop_flag <- as.integer(prob_overtox > cutoff)
    mtd <- which.min(abs(posttox - TARGET))
    return(list(MTD=mtd, posttox=posttox, prob_overtox=prob_overtox, stop=stop_flag))
  }
  
  tmp <- tmp[order(tmp$id, tmp$cycle), , drop = FALSE]
  y <- as.integer(tmp$y_obs)
  
  sigma <- sqrt(2)
  lower <- -L * sigma
  upper <-  L * sigma
  
  theta_star <- log(log(TARGET) / log(s_grid[1]))
  
  K <- length(s_grid)
  A <- length(alpha_grid)
  
  m_alpha <- numeric(A)
  m_lt_alpha <- numeric(A)
  cache_S <- vector("list", A)
  
  for (a in seq_len(A)) {
    alpha <- alpha_grid[a]
    Dij <- compute_Dij_alpha(tmp, alpha, T, d_grid)
    S_at_D <- interp_S_vec(Dij, d_grid, s_grid)
    S_at_D <- pmin(1 - eps, pmax(eps, S_at_D))
    cache_S[[a]] <- S_at_D
    
    f_den <- function(theta) exp(logpost_theta_given_alpha_vec(theta, y, S_at_D, eps))
    m_alpha[a] <- integrate(f_den, lower, upper, rel.tol = rel.tol)$value
    
    up2 <- min(theta_star, upper)
    m_lt_alpha[a] <- if (up2 <= lower) 0 else integrate(f_den, lower, up2, rel.tol = rel.tol)$value
  }
  
  w <- m_alpha / sum(m_alpha)
  prob_overtox <- sum(w * (m_lt_alpha / m_alpha))
  stop_flag <- as.integer(prob_overtox > cutoff)
  
  posttox <- numeric(K)
  for (k in seq_len(K)) {
    Ek_alpha <- numeric(A)
    for (a in seq_len(A)) {
      S_at_D <- cache_S[[a]]
      f_den <- function(theta) exp(logpost_theta_given_alpha_vec(theta, y, S_at_D, eps))
      f_num <- function(theta) {
        base <- exp(logpost_theta_given_alpha_vec(theta, y, S_at_D, eps))
        base * (s_grid[k]^exp(theta))
      }
      den <- m_alpha[a]
      num <- integrate(f_num, lower, upper, rel.tol = rel.tol)$value
      Ek_alpha[a] <- num / den
    }
    posttox[k] <- sum(w * Ek_alpha)
  }
  
  mtd <- which.min(abs(posttox - TARGET))
  
  list(
    MTD = mtd,
    posttox = posttox,
    prob_overtox = prob_overtox,
    stop = stop_flag
  )
}

# ============================================================
# Example usage:
#
# tmp_complete must contain:
#   id, cycle, dose (physical amounts), arrival_time, y_obs (0/1)
#
# d_grid <- c(15,20,30,35,45)
# s_grid <- c(0.05,0.10,0.20,0.30,0.45)
#
# est <- estimate_MTD_alphaCRM_integrate(tmp_complete, d_grid, s_grid,
#                                       TARGET=0.30, cutoff=0.95, T=28)
# est$MTD
# est$posttox
# est$prob_overtox
# est$stop
# ============================================================

simulate_aCRM_trial <- function(
    PI, PI_ipde, J = length(PI),
    COHORTSIZE = 3,
    ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
    Kmax = 3,                # max IPDE cycles per participant
    TARGET = 0.30,
    cutoff = 0.95,
    window = 28,
    arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
    seed = 1,
    verbose = FALSE,
    ordering,                # ordering for IPDE-POCRM
    model_file = 'alpha-crm.bug'
) {
  set.seed(seed)
  Nmax <- COHORTSIZE * ncohort
  ndose <- length(PI)
  # administration-level table; each row = one dosing (cycle) for some patient
  patient <- data.frame(
    id = integer(0),
    cycle = integer(0),
    dose = integer(0),
    arrival_time = numeric(0),  # this row's admin time (use arrival_time name to match your block)
    eval_time = numeric(0),
    DLT_time = numeric(0),                 
    ipde_true = integer(0),
    ipde_ok = integer(0),
    y = integer(0),           # latent truth (1 if will have DLT in this window)
    MTD = integer(0)                   # Current MTD level when patient enrolled
  )
  
  # participant-level state
  active <- rep(FALSE, Nmax)
  current_dose <- rep(NA_integer_, Nmax)
  current_cycle <- rep(0L, Nmax)      # number of doses administered so far
  stop_reason <- rep(NA_character_, Nmax)
  
  # accrual bookkeeping
  t_last <- 0
  next_pid <- 0L
  j_recent <- 1L
  MTD_hat <- 1L
  tmp_last  <- patient[0, ]   # empty
  t_mtd     <- -Inf           # time of last MTD update
  # j_S_prev <- 1L
  #stop the trial due to overtoxicity
  trial_stop = 0
  #stop IPDE due to patient cap
  stop_ipde = 0
  # event loop needs next cohort arrival time. We'll generate cohort-by-cohort using your block.
  # But evaluation events can occur between cohort arrivals, so we manage them explicitly.
  
  for (coh in 1:ncohort) {
    
    # -------- 0) process ALL evaluation events that occur before next cohort arrival (unknown yet)
    # We'll first generate the next cohort arrival times (for this cohort), and treat the FIRST arrival
    # in the cohort as the "next cohort arrival event time" for flushing pending evals.
    # This matches your cohort-level generation while still being event-driven enough.
    #
    # If you want per-participant arrivals (not batched cohorts), switch COHORTSIZE=1.
    
    # -------- 1) Generate this cohort's arrival times using YOUR code (Poisson process) ----------
    interarrival <- rexp(COHORTSIZE, rate = arrival_rate)
    arrival_time <- t_last + cumsum(interarrival) # vector length COHORTSIZE
    t_last <- tail(arrival_time, 1)
    
    # Before enrolling them at their respective arrival times,
    # we must process eval events up to each arrival in order.
    for (k in 1:COHORTSIZE) {
      
      t_now <- arrival_time[k]
      
      # -------- 2) process evaluation completions up to time t_now (can cascade with IPDE) ------
      repeat {
        # find the earliest unevaluated administration whose eval_time <= t_now
        pending_idx <- which(patient$eval_time <= t_now & patient$ipde_ok)  # evaluated by calendar time
        if (length(pending_idx) == 0) break
        # process them in chronological order, one by one
        r <- pending_idx[which.min(patient$eval_time[pending_idx])]
        t_admin <- patient$eval_time[r]
        pid <- patient$id[r]
        d_cur <- patient$dose[r]
        cyc <- patient$cycle[r]
        
        # if participant already inactive, nothing to do (but row is "evaluated" conceptually)
        if (!active[pid]) {
          # drop row from further consideration by setting eval_time to Inf (or mark a flag)
          patient$ipde_ok[r] = 0
          next
        }
        
        # build temp dataset at evaluation time (this is the "temporary dataset" for estimation)
        tmp <- build_temp_rows(patient, t_admin, window)  # your complete-only builder
        
        est <- estimate_MTD_alphaCRM_integrate(
          tmp = tmp,
          d_grid = c(15,20,30,35,45),
          s_grid = skeleton,
          TARGET = TARGET,
          cutoff = cutoff,
          T = window
        )
        cat('Post Tox', est$posttox, '\n')
        cat('MTD', est$MTD, '\n')
        MTD_hat <- est$MTD
        tmp_last <- tmp
        MTD_hat  <- est$MTD
        t_mtd    <- t_admin
        stop <- est$stop
        if(stop){
          trial_stop = TRUE
          break
        }
        # observe DLT by time t_now for this row
        y_obs <- as.integer(patient$DLT_time[r] <= t_admin)
        
        if (verbose) {
          cat(sprintf("[EVAL] t=%.2f pid=%d cycle=%d dose=%d y_obs=%d MTD=%d n_eval_d=%d\n",
                      t_now, pid, cyc, d_cur, y_obs, MTD_hat,
                      evaluated_count(patient, d_cur, t_admin)))
        }
        
        # "consume" this evaluated row so we don't process it again
        
        # stop on DLT
        if (y_obs == 1L) {
          active[pid] <- FALSE
          if (is.na(stop_reason[pid])) stop_reason[pid] <- "DLT"
          next
        }
        
        # stop if reached max cycles
        if (cyc >= Kmax) {
          active[pid] <- FALSE
          if (is.na(stop_reason[pid])) stop_reason[pid] <- "maxK"
          next
        }
        
        n_complete_d <- sum(tmp$dose == d_cur & tmp$complete == 1L)
        
        can_escalate <- (y_obs == 0L) &&
          (MTD_hat > d_cur) &&
          (n_complete_d >= 3L) &&
          (d_cur < J)
        
        if (!can_escalate) {
          active[pid] <- FALSE
          if (is.na(stop_reason[pid])) stop_reason[pid] <- "no_escalation"
          next
        }
        
        # -------- 3) Escalate immediately: add a NEW administration row at time t_now -----------
        new_dose <- d_cur + 1L
        new_cycle <- cyc + 1L
        # -------- 3) Escalate immediately: add a NEW administration row -----------
        if (nrow(patient) >= Nmax) {
          stop_ipde <- 1L
          break
        }
        # generate latent truth for the new administration (your mechanism)
        y_new <- rbinom(1L, 1L, PI_ipde[new_dose])
        eval_new <- t_admin + window 
        DLT_new <- Inf
        if (y_new == 1L) DLT_new <- t_admin + sample.int(window, size = 1, replace = TRUE)
        patient$ipde_ok[r] = 0
        patient <- rbind(
          patient,
          data.frame(
            id = pid,
            cycle = new_cycle,
            dose = new_dose,
            arrival_time = t_admin,
            eval_time = eval_new,
            DLT_time = DLT_new,
            ipde_true = 1,
            ipde_ok = 1,
            y = y_new,
            MTD = MTD_hat
          )
        )
        
        current_dose[pid] <- new_dose
        current_cycle[pid] <- new_cycle
        # j_recent <- new_dose
        
        if (verbose) {
          cat(sprintf("      -> ESCALATE pid=%d to dose=%d (cycle=%d), new eval=%.2f\n",
                      pid, new_dose, new_cycle, eval_new))
        }
      } # end repeat eval processing
      #if the trial stops in inner intra-patient escalation stage, jump out the loop
      if(trial_stop | stop_ipde){
        break
      }
      # -------- 4) Decide starting dose for this arriving participant --------------------------
      # build temp dataset at arrival time and estimate MTD (placeholder)
      # tmp <- build_temp_rows(patient, t_now, window)
      # est <- estimate_MTD_POCRM(ordering,
      #                           tmp,
      #                           skeleton,
      #                           ndose,
      #                           TARGET,
      #                           cutoff)
      # j_MTD <- est$MTD
      # #if trial stops when new patient arrives, jump out the loop
      # if(est$stop){
      #   trial_stop = 1
      #   break
      # }
      n_complete_recent <- sum(tmp_last$dose == j_recent & tmp_last$complete == 1L)
      
      if (MTD_hat > j_recent && n_complete_recent >= 3L) {
        j_S_curr <- min(j_recent + 1L, J)
      } else if (MTD_hat > j_recent && n_complete_recent < 3L) {
        j_S_curr <- j_recent
      } else if (MTD_hat < j_recent) {
        j_S_curr <- max(j_recent - 1L, 1L)
      } else {
        j_S_curr <- MTD_hat
      }
      j_recent <- j_S_curr
      # enroll ONE participant at this arrival time (k-th in the cohort)
      if (nrow(patient) >= Nmax) {break}
      next_pid <- next_pid + 1L
      pid_new <- next_pid
      
      active[pid_new] <- TRUE
      current_dose[pid_new] <- j_S_curr
      current_cycle[pid_new] <- 1L
      
      # -------- 5) Generate this participant's first administration row USING YOUR BLOCK -------
      y0 <- rbinom(1L, 1L, PI[j_S_curr])
      eval0 <- t_now + window
      DLT0 <- Inf
      if (y0 == 1L) DLT0 <- t_now + sample.int(window, size = 1, replace = TRUE)
      new_pat <- data.frame(
        id = pid_new,
        cycle = 1L,
        dose = j_S_curr,
        arrival_time = t_now,  # admin time
        eval_time = eval0,
        DLT_time = DLT0,
        y = y0,
        ipde_true = 0,
        ipde_ok = 1,
        MTD = MTD_hat
      )
      patient <- rbind(patient, new_pat)
      
      if (verbose) {
        cat(sprintf("[ARR] t=%.2f pid=%d startdose=%d eval=%.2f y(latent)=%d\n",
                    t_now, pid_new, j_S_curr, eval0, y0))
      }
      #trial stops with total of N assignments
      if (nrow(patient) >= Nmax) break
    } # end within-cohort arrivals
    
    if (nrow(patient) >= Nmax) break
  } # end cohorts
  if(nrow(patient) == 0){t_final = 0}
  else{
    t_start <- min(patient$arrival_time[patient$ipde_true == 0L])
    t_end   <- max(patient$eval_time)
    t_final <- t_end - t_start
  }
  
  if (trial_stop) {
    final_MTD = 99
    est_final <- list(MTD=NA, posttox=rep(NA, ndose))
  }
  # After accrual ends, you can optionally flush all remaining evals to finalize MTD
  # Here we just compute final estimate at t_now = max arrival + window
  else{
    
    # build temp dataset at evaluation time (this is the "temporary dataset" for estimation)
    tmp_final <- build_temp_rows(patient, t_final, window)  # your complete-only builder
    
    est_final <- estimate_MTD_alphaCRM_integrate(
      tmp = tmp_final,
      d_grid = c(15,20,30,35,45),
      s_grid = skeleton,
      TARGET = TARGET,
      cutoff = cutoff,
      T = window
    )
    # p_real <- est_final$posttox
    # jmax   <- max(tmp_final$dose) 
    # 
    # cand   <- 1:jmax
    # final_MTD <- cand[ which.min(abs(p_real[cand] - TARGET)) ]
    final_MTD <- est_final$MTD
    if(est_final$stop){
      final_MTD = 99
      est_final <- list(MTD=NA, posttox=rep(NA, ndose))
    }
  }
  list(
    patient = patient,
    stop_reason = stop_reason[1:next_pid],
    final_MTD = final_MTD,
    posttox = est_final$posttox,
    trial_time = t_final,
    trial_stop = trial_stop
  )
}