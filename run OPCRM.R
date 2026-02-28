# baseline p_j from Scenario c (exactly as in the paper)
# Eq (3) carryover data-generation:
#   D = alpha*d1 + d2
#   P(Y2=1|d1,d2,alpha) = linear interpolation of baseline p_k at dose grid d_k
#   with upper-tail extrapolation capped at 1.
source('OPCRM.R')
source('alpha-CRM.R')
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

d_grid <- c(15, 20, 30, 35, 45)  # paper's dose grid
p_base <- p_base_scenario(c = 3) # pick scenario 3 baseline p1..p5
# skeleton = c(0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
skeleton = c(0.12,0.20,0.3,0.4,0.5)
alpha <- 0.3
model_file <- 'alpha-crm.bug'
sce1 <- ipde_probabilities(c = 0, alpha)
ordering <- ipde_pocrm_orderings(d_grid)
sce1_res <- simulate_aCRM_trial(
                                sce1$p_base, 
                                sce1$p_ipde, 
                                J = length(sce1$p_base),
                                COHORTSIZE = 3,
                                ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
                                Kmax = 2,                # max IPDE cycles per participant
                                TARGET = 0.30,
                                cutoff = 0.95,
                                window = 28,
                                arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
                                seed = 1,
                                verbose = TRUE,
                                ordering$orderings,     # ordering for IPDE-POCRM
                                model_file = model_file
                            )


sce2 <- ipde_probabilities(c = 1, alpha)
ordering <- ipde_pocrm_orderings(d_grid)
sce2_res <- simulate_aCRM_trial(
  sce2$p_base, 
  sce2$p_ipde, 
  J = length(sce1$p_base),
  COHORTSIZE = 3,
  ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
  Kmax = 2,                # max IPDE cycles per participant
  TARGET = 0.30,
  cutoff = 0.95,
  window = 28,
  arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
  seed = 1,
  verbose = TRUE,
  ordering$orderings,     # ordering for IPDE-POCRM
  model_file = model_file                # ordering for IPDE-POCRM
)

sce3 <- ipde_probabilities(c = 2, alpha)
ordering <- ipde_pocrm_orderings(d_grid)
sce3_res <- simulate_aCRM_trial(
  sce3$p_base, 
  sce3$p_ipde, 
  J = length(sce3$p_base),
  COHORTSIZE = 3,
  ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
  Kmax = 2,                # max IPDE cycles per participant
  TARGET = 0.30,
  cutoff = 0.95,
  window = 28,
  arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
  seed = 1,
  verbose = TRUE,
  ordering$orderings                # ordering for IPDE-POCRM
)
sce4 <- ipde_probabilities(c = 3, alpha)
ordering <- ipde_pocrm_orderings(d_grid)
sce4_res <- simulate_aCRM_trial(
  sce4$p_base, 
  sce4$p_ipde, 
  J = length(sce4$p_base),
  COHORTSIZE = 3,
  ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
  Kmax = 2,                # max IPDE cycles per participant
  TARGET = 0.30,
  cutoff = 0.95,
  window = 28,
  arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
  seed = 1,
  verbose = TRUE,
  ordering$orderings,     # ordering for IPDE-POCRM
  model_file = model_file                # ordering for IPDE-POCRM
)
sce5 <- ipde_probabilities(c = 4, alpha)
ordering <- ipde_pocrm_orderings(d_grid)
sce5_res <- simulate_aCRM_trial(
  sce5$p_base, 
  sce5$p_ipde, 
  J = length(sce5$p_base),
  COHORTSIZE = 3,
  ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
  Kmax = 2,                # max IPDE cycles per participant
  TARGET = 0.30,
  cutoff = 0.95,
  window = 28,
  arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
  seed = 1,
  verbose = TRUE,
  ordering$orderings,     # ordering for IPDE-POCRM
  model_file = model_file                # ordering for IPDE-POCRM
)
sce6 <- ipde_probabilities(c = 5, alpha)
ordering <- ipde_pocrm_orderings(d_grid)
sce6_res <- simulate_aCRM_trial(
  sce6$p_base, 
  sce6$p_ipde, 
  J = length(sce6$p_base),
  COHORTSIZE = 3,
  ncohort = 10,            # cap by cohorts (Nmax = COHORTSIZE*ncohort)
  Kmax = 2,                # max IPDE cycles per participant
  TARGET = 0.30,
  cutoff = 0.95,
  window = 28,
  arrival_rate = 1/14,     # Poisson: mean inter-arrival 14 days
  seed = 2,
  verbose = TRUE,
  ordering$orderings,     # ordering for IPDE-POCRM
  model_file = model_file                # ordering for IPDE-POCRM
)