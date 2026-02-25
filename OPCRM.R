# setwd('/rsrch8/home/biostatistics/syang10/IPDE-Project')
# .libPaths("/rsrch8/home/biostatistics/syang10/R/x86_64-pc-linux-gnu-library/4.4")

library(coda)
library(rjags)
library(mvtnorm)

estimate_MTD_prior <- function(skeleton, TARGET=0.3, cutoff=0.96,
                               n_draw=5000, sigma_theta=sqrt(2)) {
  theta <- rnorm(n_draw, 0, sigma_theta)   # match your JAGS prior precision=0.5
  posttox <- vapply(seq_along(skeleton), \(j) mean(skeleton[j]^exp(theta)), 0.0)
  p1_draws <- skeleton[1]^exp(theta)
  prob_overtox <- mean(p1_draws > TARGET)
  stop <- as.integer(prob_overtox > cutoff)
  mtd <- which.min(abs(posttox - TARGET))
  list(MTD=mtd, posttox=posttox, prob_overtox=prob_overtox, stop=stop)
}
# build the "temporary dataset" at time t_now (evaluated only)
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
assemble_ipde_data <- function(tmp, ordering) {
  stopifnot(all(c("y_obs","complete","dose","ipde_true") %in% names(tmp)))
  if (nrow(tmp) == 0) {
    return(list(N=0, y=integer(0), d=integer(0), complete=integer(0),
                real_pos=integer(0), pseudo_pos=integer(0)))
  }
  
  ord <- as.character(ordering)
  L <- length(ord)
  
  # Parse labels like "p3" (real) and "p3p" (pseudo)
  ord_j <- as.integer(gsub("[^0-9]", "", ord))
  ord_is_pseudo <- grepl("p$", ord)   # TRUE for "p3p"
  
  K <- max(tmp$dose, na.rm = TRUE)
  
  # Safety: IPDE indicator cannot occur at dose 1
  if (any(tmp$ipde_true == 1L & tmp$dose == 1L, na.rm = TRUE)) {
    stop("Invalid data: ipde_true==1 for dose 1. IPDE only exists for doses 2..K.")
  }
  
  # Build maps: dose j -> position in ordering for real; and for pseudo (only j>=2)
  pos_real <- rep(NA_integer_, K)
  pos_pseudo <- rep(NA_integer_, K)  # pos_pseudo[1] will remain NA by design
  
  for (j in 1:K) {
    r <- which(ord_j == j & !ord_is_pseudo)
    if (length(r) != 1) {
      stop(sprintf("Need exactly one real label for dose %d (e.g., 'p%d')", j, j))
    }
    pos_real[j] <- r
    
    if (j >= 2) {
      p <- which(ord_j == j & ord_is_pseudo)
      if (length(p) != 1) {
        stop(sprintf("Need exactly one pseudo label for dose %d (e.g., 'p%dp')", j, j))
      }
      pos_pseudo[j] <- p
    }
  }
  
  # Map observed doses to ordering indices
  d_mapped <- ifelse(tmp$ipde_true == 1L,
                     pos_pseudo[tmp$dose],   # safe because ipde_true==1 implies dose>=2
                     pos_real[tmp$dose])
  
  # More safety: make sure we didn't produce NAs
  if (anyNA(d_mapped)) {
    bad <- which(is.na(d_mapped))
    stop(sprintf("Dose mapping produced NA at rows: %s", paste(bad, collapse=", ")))
  }
  
  list(
    N          = nrow(tmp),
    y          = as.integer(tmp$y_obs),
    d          = as.integer(d_mapped),
    complete   = as.integer(tmp$complete),
    real_pos   = which(!ord_is_pseudo),  # positions of p1..pK inside ordering
    pseudo_pos = which(ord_is_pseudo)    # positions of p2p..pKp inside ordering
  )
}
# helper: count evaluated at dose d at time t_now (includes the just-evaluated row)
evaluated_count <- function(patient_df, dose_d, t_now) {
  sum(patient_df$dose == dose_d & patient_df$eval_time <= t_now)
}
#Build all sequences for OPCRM
# Return all distinct simple orderings used by IPDE-POCRM:
# orderings of {d1, d2, c*d2, ..., dK, c*dK} for c in [1,2]
# Labels:
#   d1..dK   = doses without IPDE
#   d2p..dKp = doses after IPDE (modeled as c * d_k), k>=2

ipde_pocrm_orderings <- function(d, c_min = 1, c_max = 2, tol = 1e-12) {
  stopifnot(is.numeric(d), length(d) >= 2, all(is.finite(d)))
  K <- length(d)
  if (is.unsorted(d, strictly = TRUE)) stop("d must be strictly increasing (d1 < ... < dK).")
  
  # critical c values where relative order can change: c * d_k == d_j  -> c = d_j / d_k
  crit <- c(c_min, c_max)
  for (k in 2:K) {
    for (j in 1:K) {
      r <- d[j] / d[k]
      if (r > c_min + tol && r < c_max - tol) crit <- c(crit, r)
    }
  }
  crit <- sort(unique(round(crit, 14)))  # mild rounding to merge near-duplicates
  
  # helper: compute ordering at a specific c
  ord_at_c <- function(cval) {
    vals <- c(d, cval * d[2:K])
    labs <- c(paste0("d", 1:K), paste0("d", 2:K, "p"))
    o <- order(vals, decreasing = FALSE) # strictly increasing ordering
    labs[o]
  }
  
  # evaluate at endpoints and midpoints of each open interval
  cand_orders <- list()
  cand_orders[[length(cand_orders) + 1]] <- ord_at_c(c_min + 0*tol)
  for (i in seq_len(length(crit) - 1)) {
    a <- crit[i]; b <- crit[i + 1]
    if (b - a > 10 * tol) {
      mid <- (a + b) / 2
      cand_orders[[length(cand_orders) + 1]] <- ord_at_c(mid)
    }
  }
  cand_orders[[length(cand_orders) + 1]] <- ord_at_c(c_max - 0*tol)
  
  # deduplicate orderings
  keys <- vapply(cand_orders, paste, collapse = " < ", FUN.VALUE = character(1))
  uniq_idx <- !duplicated(keys)
  out <- cand_orders[uniq_idx]
  
  # return both the orderings and the c-interval breakpoints we used
  list(
    M = length(out),
    orderings = out,
    ordering_strings = keys[uniq_idx],
    critical_c = crit
  )
}
estimate_MTD_JAGS <- function(y, d_level, skeleton,
                              TARGET = 0.3,
                              cutoff = 0.96,
                              model_file = "crm.bug",
                              n.chains = 3,
                              n.adapt  = 1000,
                              n.burn   = 2000,
                              n.iter   = 5000,
                              thin     = 2) {
  
  N <- length(y)
  data_jags <- list(
    N = N,
    y = as.integer(y),
    d = as.integer(d_level),
    ske = as.numeric(skeleton)
  )
  jags <- rjags::jags.model(
    file   = model_file,
    data   = data_jags,
    n.chains = n.chains,
    n.adapt  = n.adapt,
    quiet = TRUE
  )
  
  update(jags, n.burn, progress.bar = "none")
  
  var.names <- c("theta")
  smp <- coda.samples(
    model          = jags,
    variable.names = var.names,
    n.iter         = n.iter,
    thin           = thin,
    progress.bar   = "none"
  )
  
  draws <- as.matrix(smp)
  theta <- draws[, "theta"]
  
  J <- length(skeleton)
  posttox <- vapply(seq_len(J), function(j) {
    mean(skeleton[j]^exp(theta))
  }, numeric(1))
  pi1_draws <- skeleton[1]^exp(theta)
  prob_overtox <- mean(pi1_draws > TARGET)
  stop_flag <- as.integer(prob_overtox > cutoff)
  
  # cat('posttox', posttox, '\n')
  # cat('prob_overtox', prob_overtox, '\n')
  # MTD = closest to TARGET
  diff <- abs(posttox - TARGET)
  dose.best <- which(diff == min(diff))[1]
  
  
  
  
  list(
    MTD = dose.best,
    posttox = posttox,
    prob_overtox = prob_overtox,
    stop = stop_flag
  )
}
log_marginal_lik <- function(y, d, skeleton, sigma = sqrt(2),  # Falke uses Var=2
                             L = 8, rel.tol = 1e-8) {
  
  stopifnot(length(y) == length(d))
  s_i <- skeleton[d]
  if (any(s_i <= 0 | s_i >= 1, na.rm = TRUE)) stop("skeleton entries must be in (0,1).")
  log_integrand <- function(beta) {
    
    # clamp away from 0/1 to avoid log(0)
    eps <- 1e-12
    ll <- 0
    for(i in 1:length(y)){
      p <- skeleton[d[i]]^exp(beta)
      p <- pmin(1-eps, pmax(eps, p))
      ll <- ll + y[i]*log(p) + (1-y[i])*log(1-p)
    }
    
    prior <- dnorm(beta, 0, sigma, log = TRUE)
    ll + prior
  }
  
  # integrate over ±L*sigma (practically all prior mass)
  lower <- -L * sigma
  upper <-  L * sigma
  
  f <- function(beta) exp(log_integrand(beta))
  val <- integrate(f, lower = lower, upper = upper, rel.tol = rel.tol)$value
  if (!is.finite(val) || val <= 0) stop("Integration failed / returned non-positive value.")
  log(val)
}
compute_model_weights <- function(logml){
  a <- max(logml)
  w <- exp(logml - a)
  w / sum(w)
}

# Example:
# d <- c(15, 20, 30, 35, 45)
# res <- ipde_pocrm_orderings(d)
# res$M
# res$orderings[[1]]
# cat(res$ordering_strings, sep = "\n\n")

estimate_MTD_POCRM <- function(orderings,
                               tmp,
                               skeleton,
                               ndose = 5,
                               TARGET = 0.3,
                               cutoff = 0.96,
                               model_file = "crm.bug", n.chains = 3,
                               n.adapt  = 1000,
                               n.burn   = 2000,
                               n.iter   = 5000,
                               thin     = 2 ){
  if (nrow(tmp)==0) return(estimate_MTD_prior(skeleton, TARGET, cutoff))
  posttox_order = matrix(nrow = length(orderings), ncol = ndose)
  p_overtox_order = rep(0,length(orderings))
  marginal = rep(0,length(orderings))
  for(i in 1:length(orderings)){
    combined = assemble_ipde_data(tmp, orderings[[i]])
    posttox = estimate_MTD_JAGS(combined$y, 
                                combined$d, 
                                skeleton,
                                TARGET,
                                cutoff,
                                model_file,
                                n.chains,
                                n.adapt,
                                n.burn,
                                n.iter,
                                thin)
    posttox_order[i,] = posttox$posttox[combined$real_pos]
    p_overtox_order[i] = posttox$prob_overtox
    marginal[i] = log_marginal_lik(combined$y,combined$d,skeleton, sigma = sqrt(2))
  }
  weights = compute_model_weights(marginal)
  weighted_p <- as.numeric(t(weights) %*% posttox_order)
  weighted_overtox = sum(p_overtox_order * weights)
  stop_flag = as.integer(weighted_overtox > cutoff)
  
  cat('posttox', weighted_p, '\n')
  cat('prob_overtox', weighted_overtox, '\n')
  # MTD = closest to TARGET
  diff <- abs(weighted_p - TARGET)
  dose.best <- which(diff == min(diff))[1]
  
  list(
    MTD = dose.best,
    posttox = weighted_p,
    prob_overtox = weighted_overtox,
    stop = stop_flag
  )
}

# MAIN: one trial
# OPCRM
simulate_IPDE_trial <- function(
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
    ordering                 # ordering for IPDE-POCRM
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
        tmp <- build_temp_rows(patient, t_admin, window)
        est <- estimate_MTD_POCRM(ordering,
                                  tmp,
                                  skeleton,
                                  ndose,
                                  TARGET,
                                  cutoff)
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
        if (y_new == 1L) DLT_new <- runif(1L, t_admin, eval_new) 
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
        j_recent <- new_dose
        
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
      if (y0 == 1L) DLT0 <- runif(1L, min = t_now, max = eval0)
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
  t_final <- if (nrow(patient) == 0) 0 else (max(patient$arrival_time) + window)
  if (trial_stop) {
    final_MTD = 99
    est_final <- list(MTD=NA, posttox=rep(NA, ndose))
  }
  # After accrual ends, you can optionally flush all remaining evals to finalize MTD
  # Here we just compute final estimate at t_now = max arrival + window
  else{
    tmp_final <- build_temp_rows(patient, t_final, window)
    est_final <- estimate_MTD_POCRM(ordering,
                              tmp_final,
                              skeleton,
                              ndose,
                              TARGET,
                              cutoff)
    if(est_final$stop){
      final_MTD = 99
      est_final <- list(MTD=NA, posttox=rep(NA, ndose))
    }
  }
  list(
    patient = patient,
    stop_reason = stop_reason[1:next_pid],
    final_MTD = est_final$MTD,
    posttox = est_final$posttox,
    trial_time = t_final,
    trial_stop = trial_stop
  )
}
