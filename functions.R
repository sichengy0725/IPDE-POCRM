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

# ---- 1) Utilities: S(D) interpolation from (d_grid, s_grid) ----
interp_S_vec <- function(D, d_grid, s_grid) {
  K <- length(d_grid)
  out <- numeric(length(D))
  
  idx_lo <- D <= d_grid[1]
  out[idx_lo] <- s_grid[1]
  
  idx_hi <- D >= d_grid[K]
  if (any(idx_hi)) {
    t <- (D[idx_hi] - d_grid[K-1]) / (d_grid[K] - d_grid[K-1])
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


p_base_scenario <- function(c, j = 1:5, phi = 0.3) {
  1 / (1 + exp((c - j)/2 - log(phi/(1 - phi))))
}
interp_tox <- function(D, d_grid, p_base) {
  
  K <- length(d_grid)
  
  # D inside grid
  if (D < d_grid[K]) {
    k <- findInterval(D, d_grid)
    
    if (k == 0) k <- 1
    
    d_k  <- d_grid[k]
    d_k1 <- d_grid[k+1]
    
    p_k  <- p_base[k]
    p_k1 <- p_base[k+1]
    
    return( p_k + (D - d_k)/(d_k1 - d_k)*(p_k1 - p_k) )
  }
  
  # D beyond largest dose → extrapolate & cap at 1
  d_Km1 <- d_grid[K-1]
  d_K   <- d_grid[K]
  
  p_Km1 <- p_base[K-1]
  p_K   <- p_base[K]
  
  p <- p_Km1 + (D - d_Km1)/(d_K - d_Km1)*(p_K - p_Km1)
  
  return(min(1, p))
}
ipde_probabilities <- function(c, alpha,
                               d_grid = c(15,20,30,35,45),
                               phi = 0.3){
  
  # baseline p_k
  p_base <- baseline_curve(c, 1:length(d_grid), phi)
  
  K <- length(d_grid)
  p_ipde <- numeric(K)
  
  # first dose has no IPDE
  p_ipde[1] <- NA
  
  # compute p'_k
  for(k in 2:K){
    
    d_prev <- d_grid[k-1]
    d_curr <- d_grid[k]
    
    D_eff <- alpha*d_prev + d_curr
    
    p_ipde[k] <- interp_tox(D_eff, d_grid, p_base)
  }
  
  return(list(
    p_base = p_base,   # p_k
    p_ipde = p_ipde    # p'_k
  ))
}
baseline_curve <- function(c, j = 1:5, phi = 0.3){
  1/(1 + exp((c-j)/2 - log(phi/(1-phi))))
}
build_temp_rows <- function(patient_df, t_now, window, p = rep(0,5)) {
  
  # only administrations that already occurred
  tmp <- patient_df[patient_df$arrival_time <= t_now, , drop = FALSE]
  
  if (nrow(tmp) == 0) return(tmp)
  
  # make sure rows are ordered within patient by cycle / time
  tmp <- tmp[order(tmp$id, tmp$cycle, tmp$arrival_time), , drop = FALSE]
  
  # raw dose amount corresponding to dose index
  tmp$dose_raw <- p[tmp$dose]
  
  # cumulative dose PRIOR to current administration, using raw dose
  tmp$cumu_dose <- ave(
    tmp$dose_raw,
    tmp$id,
    FUN = function(x) c(0, cumsum(x)[-length(x)])
  )
  
  # 1) has the DLT already been observed?
  tmp$y_obs <- as.integer(tmp$DLT_time <= t_now)
  
  # 2) is the toxicity outcome now known?
  # known if either DLT occurred OR window finished
  tmp$complete <- as.integer(tmp$y_obs == 1L | tmp$eval_time <= t_now)
  
  # keep only rows whose toxicity status is known
  tmp <- tmp[tmp$complete == 1, , drop = FALSE]
  
  return(tmp)
}
# ------------------------------------------------------------
# Posterior MTD estimator (handles all 4 model types)
# ------------------------------------------------------------
estimate_MTD_JAGS <- function(y, d_level,
                              p = NULL,
                              Tcum = NULL,
                              pid = NULL,
                              skeleton = NULL,
                              TARGET = 0.3,
                              cutoff = 0.96,
                              model_file = "logit.bug",
                              n.chains = 3,
                              n.adapt  = 1000,
                              n.burn   = 2000,
                              n.iter   = 5000,
                              thin     = 2,
                              T_ref_for_curve = 0) {
  
  # model flags
  is_CRM <- model_file %in% c("crm.bug")
  is_CO  <- model_file %in% c("logit_CO.bug", "logit_CORF.bug")
  is_RF  <- model_file %in% c("logit_RF.bug", "logit_CORF.bug")
  
  if (is_CRM) {
    if (is.null(skeleton)) stop("For crm.bug, 'skeleton' must be provided.")
    stopifnot(length(y) == length(d_level))
    
    N <- length(y)
    data_jags <- list(
      N   = N,
      y   = as.integer(y),
      d   = as.integer(d_level),
      ske = as.numeric(skeleton)
    )
    
    jags <- rjags::jags.model(
      file     = model_file,
      data     = data_jags,
      n.chains = n.chains,
      n.adapt  = n.adapt,
      quiet    = TRUE
    )
    
    update(jags, n.burn, progress.bar = "none")
    
    smp <- coda.samples(
      model          = jags,
      variable.names = "theta",
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
    
  } else {
    if (is.null(p)) stop("For logit-type models, 'p' must be provided.")
    if (is.null(Tcum)) Tcum <- rep(0, length(y))
    if (is.null(pid))  pid  <- seq_along(y)
    
    stopifnot(length(y) == length(d_level),
              length(y) == length(Tcum),
              length(y) == length(pid))
    
    Nobs <- length(y)
    nPat <- if (Nobs == 0) 0L else max(pid)
    
    d_numeric <- as.numeric(p[d_level])
    T_numeric <- as.numeric(Tcum)
    pid_int   <- as.integer(pid)
    
    data_jags <- list(
      Nobs  = Nobs,
      nPat  = nPat,
      y_bin = as.integer(y),
      d     = d_numeric,
      T     = T_numeric,
      pid   = pid_int
    )
    
    jags <- rjags::jags.model(
      file     = model_file,
      data     = data_jags,
      n.chains = n.chains,
      n.adapt  = n.adapt,
      quiet    = TRUE
    )
    
    update(jags, n.burn, progress.bar = "none")
    
    var.names <- c("beta0", "beta1")
    if (is_CO) var.names <- c(var.names, "beta2")
    if (is_RF) var.names <- c(var.names, "sigma")
    
    smp <- coda.samples(
      model          = jags,
      variable.names = var.names,
      n.iter         = n.iter,
      thin           = thin,
      progress.bar   = "none"
    )
    
    draws <- as.matrix(smp)
    b0 <- draws[, "beta0"]
    b1 <- draws[, "beta1"]
    b2 <- if (is_CO) draws[, "beta2"] else 0
    sigma <- if (is_RF) draws[, "sigma"] else 0
    
    J <- length(p)
    
    # non-CO, non-RF: plain logit
    if (!is_CO && !is_RF) {
      p_j_draws <- sapply(seq_len(J), function(j) plogis(b0 + b1 * p[j]))
      posttox <- colMeans(p_j_draws)
      p1_draws <- p_j_draws[, 1]
      prob_overtox <- mean(p1_draws > TARGET)
      stop_flag <- as.integer(prob_overtox > cutoff)
    }
    
    # RF only
    if (!is_CO && is_RF) {
      M <- 200
      U <- matrix(rnorm(length(b0) * M), nrow = length(b0), ncol = M) * sigma
      eta <- outer(b0, rep(1, J)) + outer(b1, p)
      
      p_j_draws <- sapply(seq_len(J), function(j) {
        rowMeans(plogis(eta[, j] + U))
      })
      
      posttox <- colMeans(p_j_draws)
      p1_draws <- p_j_draws[, 1]
      prob_overtox <- mean(p1_draws > TARGET)
      stop_flag <- as.integer(prob_overtox > cutoff)
    }
    
    # CO without RF
    if (is_CO && !is_RF) {
      posttox <- vapply(seq_len(J), function(j) {
        mean(plogis(b0 + b1 * p[j] + b2 * T_ref_for_curve))
      }, numeric(1))
      
      pi1_draws <- plogis(b0 + b1 * p[1] + b2 * T_ref_for_curve)
      prob_overtox <- mean(pi1_draws > TARGET)
      stop_flag <- as.integer(prob_overtox > cutoff)
    }
    
    # CO + RF
    if (is_CO && is_RF) {
      M <- 200
      U <- matrix(rnorm(length(b0) * M), nrow = length(b0), ncol = M) * sigma
      
      eta_base <- outer(b0, rep(1, J)) + outer(b1, p) + outer(b2, rep(T_ref_for_curve, J))
      
      p_j_draws <- sapply(seq_len(J), function(j) {
        rowMeans(plogis(eta_base[, j] + U))
      })
      
      posttox <- colMeans(p_j_draws)
      p1_draws <- p_j_draws[, 1]
      prob_overtox <- mean(p1_draws > TARGET)
      stop_flag <- as.integer(prob_overtox > cutoff)
    }
  }
  
  diff <- abs(posttox - TARGET)
  dose.best <- which(diff == min(diff))[1]
  
  list(
    MTD = dose.best,
    posttox = posttox,
    prob_overtox = prob_overtox,
    stop = stop_flag
  )
}

#POCRM Helpers
log_marginal_lik <- function(y, d, skeleton, sigma = sqrt(2),  # Falke uses Var=2
                             L = 8, rel.tol = 1e-8) {
  
  stopifnot(length(y) == length(d))
  s_i <- skeleton[d]
  if (any(s_i <= 0 | s_i >= 1, na.rm = TRUE)) stop("skeleton entries must be in (0,1).")
  # cat('y', y, '\n')
  # cat('d', d, '\n')
  # cat('ske', s_i, '\n')
  log_integrand <- function(beta) {
    
    # clamp away from 0/1 to avoid log(0)
    eps <- 1e-12
    ll <- 0
    for(i in 1:length(y)){
      p <- pmin(1 - eps, pmax(eps, skeleton[d[i]]^exp(beta)))
      ll <- ll + (y[i] * log(p))+ (1 - y[i]) * log(1-p)
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

estimate_MTD_POCRM <- function(orderings,
                               tmp,
                               skeleton,
                               ndose = 5,
                               TARGET = 0.3,
                               cutoff = 0.96,
                               model_file = "crm.bug",
                               n.chains = 3,
                               n.adapt  = 1000,
                               n.burn   = 2000,
                               n.iter   = 5000,
                               thin     = 2) {
  if (nrow(tmp) == 0) return(estimate_MTD_prior(skeleton, TARGET, cutoff))
  
  posttox_order   <- matrix(NA_real_, nrow = length(orderings), ncol = ndose)
  p_overtox_order <- numeric(length(orderings))
  marginal        <- numeric(length(orderings))
  
  for (i in seq_along(orderings)) {
    combined <- assemble_ipde_data(tmp, orderings[[i]])
    
    posttox <- estimate_MTD_JAGS(
      y         = combined$y,
      d_level   = combined$d,
      skeleton  = skeleton,
      TARGET    = TARGET,
      cutoff    = cutoff,
      model_file = model_file,
      n.chains  = n.chains,
      n.adapt   = n.adapt,
      n.burn    = n.burn,
      n.iter    = n.iter,
      thin      = thin
    )
    
    posttox_order[i, ]  <- posttox$posttox[combined$real_pos]
    p_overtox_order[i]  <- posttox$prob_overtox
    marginal[i]         <- log_marginal_lik(
      combined$y,
      combined$d,
      skeleton,
      sigma = sqrt(2)
    )
  }
  
  weights <- compute_model_weights(marginal)
  weighted_p <- as.numeric(t(weights) %*% posttox_order)
  weighted_overtox <- sum(p_overtox_order * weights)
  stop_flag <- as.integer(weighted_overtox > cutoff)
  
  diff <- abs(weighted_p - TARGET)
  dose.best <- which(diff == min(diff))[1]
  
  list(
    MTD = dose.best,
    posttox = weighted_p,
    prob_overtox = weighted_overtox,
    stop = stop_flag
  )
}