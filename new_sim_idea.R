mu1=c(-3,-3,-3,-1.5,-1.5,-1.5,0,0,0,0)
mu2=rep(0,10)
mu3=abs(mu1)

# Outlier
Y[sample(nrow(Y),nrow(Y)*0.05),]=runif(nrow(Y)*0.05*ncol(Y),-100,100)



sim_data_stud_t_FWJM <- function(seed = 123,
                            TT,
                            P,
                            Pcat = NULL,
                            Ktrue = 3,
                            mu = 3,
                            rho = 0,
                            nu = 4,
                            phi = .8,
                            pers = .95,
                            W = NULL,
                            outlier_frac = 0,
                            Out_bound = 100) {
  
  if (is.null(W)) {
    W <- matrix(1, Ktrue, P)
  }
  
  set.seed(seed)
  
  # state means base levels
  MU <- seq(-mu, mu, length.out = Ktrue)
  
  # Markov chain
  x <- numeric(TT)
  
  Q <- matrix((1 - pers) / (Ktrue - 1),
              nrow = Ktrue,
              ncol = Ktrue)
  
  diag(Q) <- pers
  
  init <- rep(1 / Ktrue, Ktrue)
  
  x[1] <- sample(1:Ktrue, 1, prob = init)
  
  for (i in 2:TT) {
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Covariance
  Sigma <- matrix(rho, ncol = P, nrow = P)
  diag(Sigma) <- 1
  
  # scale adjustment so Var = Sigma
  Sigma_t <- (nu - 2) / nu * Sigma
  
  SimData <- matrix(0, TT, P)
  
  # Simulate by state
  for (k in 1:Ktrue) {
    
    # state-specific mean vector
    delta_k <- MU[k] * W[k, ]
    
    u_k <- mvtnorm::rmvt(
      n = TT,
      sigma = Sigma_t,
      df = nu,
      delta = delta_k
    )
    
    idx <- which(x == k)
    
    if (length(idx) > 0)
      SimData[idx, ] <- u_k[idx, ]
  }
  
  # OUTLIER CONTAMINATION
  if (outlier_frac > 0) {
    
    n_out <- ceiling(TT * outlier_frac)
    
    idx_out <- sample(seq_len(TT), n_out)
    
    SimData[idx_out, ] <-
      matrix(runif(n_out * P, -Out_bound, Out_bound),
             nrow = n_out,
             ncol = P)
  }
  
  SimData <- as.data.frame(SimData)
  
  # Optional categorical
  if (!is.null(Pcat)) {
    
    for (j in 1:Pcat) {
      SimData[, j] <- get_cat_t(
        SimData[, j],
        x,
        MU,
        phi = phi,
        df = nu
      )
      SimData[, j] <- factor(SimData[, j], levels = 1:Ktrue)
    }
  }
  
  return(list(
    SimData = SimData,
    mchain = x,
    TT = TT,
    P = P,
    K = Ktrue,
    Ktrue = Ktrue,
    pers = pers,
    seed = seed,
    W = W
  ))
}
