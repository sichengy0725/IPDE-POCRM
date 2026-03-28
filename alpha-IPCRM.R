setwd('/rsrch8/home/biostatistics/syang10/IPDE-POCRM-main/IPDE-POCRM-main')
.libPaths("/rsrch8/home/biostatistics/syang10/R/x86_64-pc-linux-gnu-library/4.4")
library(coda)
library(rjags)
library(mvtnorm)
source('functions.R')
simulate_IPCRM_trial <- function(
    PI, PI_ipde, 
    J = length(PI),
    COHORTSIZE = 3,
    ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
    Kmax = 3,                # max IPDE cycles per participant
    TARGET = 0.30,
    window = 28,
    arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
    seed = 1,
    verbose = FALSE,             
    model = 'IPCRM',
    parameters,
    escalation_rule_ipde = 1,
    escalation_rule_new = 1     #1: alphacrm with n >= 3, 2: no >= 3
) {
  set.seed(seed)
  Nmax <- COHORTSIZE * ncohort
  ndose <- length(PI)
  is_CO <- model %in% c("IPCRM-CO")
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
        if(model %in% c('IPCRM-CO', 'IPCRM')){
          tmp <- build_temp_rows(patient, t_admin, window, parameters$d)  # your complete-only builder
        }
        else{
          tmp <- build_temp_rows(patient, t_admin, window)
        }
        if(model %in% c('IPCRM', 'IPCRM-CO')){
          est <- estimate_MTD_JAGS(tmp$y, tmp$dose, parameters$d, tmp$cumu_dose, tmp$id,
                                   TARGET = TARGET,
                                   cutoff = parameters$cutoff,
                                   model_file = parameters$model_file)
        }
        if(model %in% c('Alpha-CRM')){
          est <- estimate_MTD_alphaCRM_integrate(
            tmp = tmp,
            d_grid = parameters$d_grid,
            s_grid = parameters$s_grid,
            TARGET = TARGET,
            cutoff = parameters$cutoff,
            T = window
          )
        }
        if(model %in% c("POCRM")){
          est <- estimate_MTD_POCRM(parameters$ordering,
                                    tmp,
                                    parameters$skeleton,
                                    ndose,
                                    TARGET,
                                    parameters$cutoff,
                                    model_file = parameters$model_file)
        }
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
        if(escalation_rule_ipde == 1){
          can_escalate <- (y_obs == 0L) &&
            (MTD_hat > d_cur) &&
            (n_complete_d >= 3L) &&
            (d_cur < J)
        }
        if(escalation_rule_ipde == 2){
          can_escalate <- (y_obs == 0L) &&
            (MTD_hat > d_cur) &&
            (d_cur < J)
        }
        #IPCRM
        if(escalation_rule_ipde == 3){
          can_escalate <- (y_obs == 0L) &&
            (MTD_hat >= d_cur) &&
            (d_cur < J)
        }
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
        if(window == 0){
          if (y_new == 1L) DLT_new <- t_admin 
        }
        else{
          if (y_new == 1L) DLT_new <- t_admin + sample.int(window, size = 1, replace = TRUE)
        }
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
      n_complete_recent <- sum(tmp_last$dose == j_recent & tmp_last$complete == 1L)
      if(escalation_rule_new == 1){
        if (MTD_hat > j_recent && n_complete_recent >= 3L) {
          j_S_curr <- min(j_recent + 1L, J)
        } else if (MTD_hat > j_recent && n_complete_recent < 3L) {
          j_S_curr <- j_recent
        } else if (MTD_hat < j_recent) {
          j_S_curr <- max(j_recent - 1L, 1L)
        } else {
          j_S_curr <- MTD_hat
        }
      }
      if(escalation_rule_new == 2){
        if (MTD_hat > j_recent) {
          j_S_curr <- min(j_recent + 1L, J)
        } else if (MTD_hat < j_recent) {
          j_S_curr <- max(j_recent - 1L, 1L)
        } else {
          j_S_curr <- MTD_hat
        }
      }
      if(escalation_rule_new == 3) {
        if(nrow(patient) == 0){
          j_S_curr = 1
        } else {
          j_H <- max(patient$dose)
          
          if (MTD_hat > j_H) {
            j_S_curr <- min(j_H + 1L, J)
          } else if (MTD_hat == j_H) {
            j_S_curr <- j_H
          } else {
            # MTD_hat < j_H
            if (MTD_hat >= j_recent) {
              j_S_curr <- MTD_hat
            } else {
              j_S_curr <- max(j_recent - 1L, 1L)
            }
          }  
        }

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
      if(window == 0){
        if (y0 == 1L) DLT0 <- t_now 
      }
      else{
        if (y0 == 1L) DLT0 <- t_now + sample.int(window, size = 1, replace = TRUE)
      }
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
    # tmp_final <- build_temp_rows(patient, t_end, window, p_model)  # your complete-only builder
    if(model %in% c('IPCRM-CO', 'IPCRM')){
      tmp_final <- build_temp_rows(patient, t_end, window, parameters$d)  # your complete-only builder
    }
    else{
      tmp_final <- build_temp_rows(patient, t_end, window)
    }
    if(model %in% c('IPCRM', 'IPCRM-CO')){
      est_final <- estimate_MTD_JAGS(tmp_final$y, 
                                     tmp_final$dose, 
                                     parameters$d, 
                                     tmp_final$cumu_dose, 
                                     tmp_final$id,
                                     TARGET = TARGET,
                                     cutoff = parameters$cutoff,
                                     model_file = parameters$model_file)
      
    }
   
    if(model %in% c('Alpha-CRM')){
      est_final <- estimate_MTD_alphaCRM_integrate(
        tmp = tmp_final,
        d_grid = parameters$d_grid,
        s_grid = parameters$s_grid,
        TARGET = TARGET,
        cutoff = parameters$cutoff,
        T = window
      )
    }
    if(model %in% c("POCRM")){
      est_final <- estimate_MTD_POCRM(parameters$ordering,
                                tmp_final,
                                parameters$skeleton,
                                ndose,
                                TARGET,
                                parameters$cutoff,
                                model_file = parameters$model_file)
    }
    final_MTD <- est_final$MTD
    if(est_final$stop){
      final_MTD = 99
      est_final <- list(MTD=NA, posttox=rep(NA, ndose))
    }
    else{
      p_real <- est_final$posttox
      jmax   <- max(tmp_final$dose) 
      # 
      cand   <- 1:jmax
      final_MTD <- cand[ which.min(abs(p_real[cand] - TARGET))]
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

backsol <- function(ske, mu_beta0, mu_beta1){
  dose = (log(ske/(1-ske)) - mu_beta0)/mu_beta1
}
ske1 = c(0.02, 0.12, 0.3, 0.5, 0.65)
d_grid <- c(15, 20, 30, 35, 45)
skeleton <- c(0.12, 0.20, 0.30, 0.40, 0.50)
#backsolve for d_j for IPCRM
dose = backsol(ske1, mu_beta0 = 3, mu_beta1 = 1)
#dose for IPCRM-CO
co_dose <- c(.1, .3, .5, .7, .9)
co_dose <- co_dose / (2 * sd(co_dose))
model <- 'Alpha-CRM'

parameters <- switch(
  model,
  "IPCRM" = list(cutoff = 0.96, d = dose, model_file = 'logit.bug'),
  "IPCRM-CO" = list(cutoff = 0.73, d = co_dose, model_file = 'logit_CO.bug'),
  "Alpha-CRM" = list(cutoff = 0.95, 
                     s_grid = c(0.12, 0.20, 0.30, 0.40, 0.50),
                     d_grid = c(15, 20, 30, 35, 45)),
  "POCRM" = list(cutoff = 0.95, 
                 skeleton = c(0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                 ordering = ipde_pocrm_orderings(d_grid)$orderings,
                 model_file = 'crm.bug'),
  stop("Unknown model")
)

Kmax = 2
window = 28
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide the job index as the first argument.")
job_i <- as.integer(args[1])
cat('Seed', job_i, '\n')
alpha_list <- c(0, 0.3, 0.6, 0.9)
for(j in 1){
  for(k in 1){
    for (alpha in alpha_list) {
      cat('Alpha', alpha, '\n')
      for (sc in 0:5) {
        cat('sce', sc, '\n')
        sce <- ipde_probabilities(c = sc, alpha)
        
        res <- simulate_IPCRM_trial(
          PI = sce$p_base,
          PI_ipde = sce$p_ipde,
          J = length(ske1),
          COHORTSIZE = 3,
          ncohort = 10,
          Kmax = Kmax,
          TARGET = 0.30,
          window = window,
          arrival_rate = 1/14,
          seed = job_i,
          verbose = FALSE,
          model = model,
          parameters = parameters,
          escalation_rule_ipde = j,
          escalation_rule_new = k
        )
        
        foldername <- paste0(
          "results/3-27-2026/Alpha-IPCRM/",
          "res", sc + 1,
          "-alpha-", alpha,
          "-Kmax-", Kmax,
          "-Window-", window,
          "-model-", model,
          "-escalation_rule_ipde-", j,
          "-escalation_rule_new-", k
        )
        
        if (!dir.exists(foldername)) {
          dir.create(foldername, recursive = TRUE)
        }
        
        saveRDS(res, paste0(foldername, "/trial-", job_i))
      }
    }
  }
}
