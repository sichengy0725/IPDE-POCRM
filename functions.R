#returns Dij = \sum j' dij'*alpha^((tj-tj')/T)
effective_dose <- function(d,tau,alpha,T){
 D = d * alpha^(tau/T)
}
#takes in dosage, time of each dose, parameter alpha,t
#returns Dij = \sum j' dij'*alpha^((tj-tj')/T)
#d - an array indicating dosage amount received
#t - an array indicating date dose received
#T - evaluation window
#alpha - decay parameter
effective_cum_dose <- function(d, t, alpha, T) {
  stopifnot(length(d) == length(t))
  o <- order(t)
  d <- d[o]; t <- t[o]
  
  m <- length(d)
  D <- numeric(m)
  D[1] <- d[1]
  
  for (j in 2:m) {
    decay <- alpha^((t[j] - t[j-1]) / T)
    D[j] <- d[j] + decay * D[j-1]
  }
  D
}

#S(Dij)
s_dij <- function(dij, dose_level, skeleton){
  if (dij < dose_level[1]) stop("dij below minimum dose_level; check units / dosing grid")
  K = length(skeleton)
  if(dij > max(dose_level)){
    d_new =  skeleton[K-1] + 
      ((dij - dose_level[K-1])/(dose_level[K] - dose_level[K-1])) * 
      (skeleton[K] - skeleton[K-1])
    sdij = min(1,d_new)
  }
  else{
    for(i in 2:length(skeleton)){
      #if Dij is in a certain interval
      if(dij <= dose_level[i] & dij >= dose_level[i-1]){
        sdij = skeleton[i-1] + 
          ((dij - dose_level[i-1])/(dose_level[i] - dose_level[i-1])) * 
          (skeleton[i] - skeleton[i-1])
      } 
    }
  }
  sdij
  
}


draws_JAGS <- function(y, d, p, e,
                       model_file = "eff.bug",
                       n.chains = 3,
                       n.adapt  = 1000,
                       n.burn   = 2000,
                       n.iter   = 5000,
                       thin     = 2) {
  
  data_jags <- list(
    N     = length(y),
    y_bin = as.integer(y),
    e_bin = as.integer(e),
    dose_index = as.integer(d),
    D = length(p),
    d_j     = as.numeric(p[d])  # p is effective dose vector (dtilde)
  )
  
  jags <- rjags::jags.model(
    file    = model_file,
    data    = data_jags,
    n.chains = n.chains,
    n.adapt = n.adapt,
    quiet   = TRUE
  )
  
  update(jags, n.burn, progress.bar = "none")
  
  smp <- coda.samples(
    model          = jags,
    variable.names = c("beta0", "beta1", "gamma0", "gamma"),
    n.iter         = n.iter,
    thin           = thin,
    progress.bar   = "none"
  )
  draws <- as.matrix(smp)  # works for mcmc.list; stacks chains
  beta0 <- draws[, "beta0"]
  beta1 <- draws[, "beta1"]
  gamma0 <- draws[, "gamma0"]
  gamma <- draws[, 3:(ncol(draws)-1)]
  list(
    beta0 = beta0,
    beta1 = beta1,
    gamma0 = gamma0,
    gamma = gamma
  )
}
