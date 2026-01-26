
sim_hmm_subspace <- function(seed = 123,
                             TT = 200,
                             P = 5,
                             Ktrue = 3,
                             a_list = NULL,            # list di lunghezza Ktrue, ciascuno vettore length P
                             mu_tilde_list = NULL,     # list di base means
                             Sigma_tilde_list = NULL,  # list di base covariance matrices
                             pers = 0.95,
                             eps = 1e-6) {
  
  set.seed(seed)
  
  #persistenza
  Q <- matrix((1 - pers) / (Ktrue - 1), nrow = Ktrue, ncol = Ktrue)
  diag(Q) <- rep(pers, Ktrue)
  
  x <- integer(TT)
  init_prob <- rep(1 / Ktrue, Ktrue)
  x[1] <- sample.int(Ktrue, size = 1, prob = init_prob)
  for (t in 2:TT) {
    x[t] <- sample.int(Ktrue, size = 1, prob = Q[x[t - 1], ])
  }
  
  mu_k_list <- vector("list", Ktrue)
  Sigma_k_list <- vector("list", Ktrue)
  for (k in seq_len(Ktrue)) {
    a_k <- pmax(a_list[[k]], eps)                    
    A_k <- diag(a_k, nrow = P, ncol = P)
    Ainv_k <- diag(1 / a_k, nrow = P, ncol = P)
    mu_k_list[[k]] <- as.numeric(A_k %*% as.numeric(mu_tilde_list[[k]]))
    Sigma_k_list[[k]] <- Ainv_k %*% Sigma_tilde_list[[k]] %*% Ainv_k
    # simmetry
    Sigma_k_list[[k]] <- (Sigma_k_list[[k]] + t(Sigma_k_list[[k]])) / 2
  }
  
  # Simulazione osservazioni
  Y <- matrix(NA_real_, nrow = TT, ncol = P)
  for (k in seq_len(Ktrue)) {
    Y_k <- mvtnorm::rmvnorm(n = TT, mean = mu_k_list[[k]], sigma = Sigma_k_list[[k]])
    idx_k <- which(x == k)
    if (length(idx_k) > 0) {
      Y[idx_k, ] <- Y_k[idx_k, , drop = FALSE]
    }
  }
  
  colnames(Y) <- paste0("y", seq_len(P))
  SimData <- as.data.frame(Y)
  
  return(list(
    SimData = SimData,
    mchain = x,
    TT = TT,
    P = P,
    K = Ktrue,
    a_list = a_list,
    mu_tilde_list = mu_tilde_list,
    Sigma_tilde_list = Sigma_tilde_list,
    mu_k_list = mu_k_list,
    Sigma_k_list = Sigma_k_list,
    pers = pers,
    seed = seed
  ))
}

P <- 5; Ktrue <- 3; TT <- 1000
a1 <- c(3, 1, 0.1, 0.1, 0.1)
a2 <- c(0.1, 0.1, 3, 1, 0.1)
a3 <- c(.1, 0.1, 1, 0.1, 3)
a_list <- list(a1, a2,a3)
mu_tilde_list <- list(rep(1, P), rep(1, P), rep(1, P))
Sigma_tilde_list <- list(diag(1, P), diag(1, P), diag(1, P))

res <- sim_hmm_subspace(seed = 123,
                        TT = TT,
                        P = P,
                        Ktrue = Ktrue,
                        a_list = a_list,
                        mu_tilde_list = mu_tilde_list,
                        Sigma_tilde_list = Sigma_tilde_list,
                        pers = 0.95)

x11()
par(mfrow=c(2,3))
plot(res$SimData$y1,col=res$mchain)
plot(res$SimData$y2,col=res$mchain)
plot(res$SimData$y3,col=res$mchain)
plot(res$SimData$y4,col=res$mchain)
plot(res$SimData$y5,col=res$mchain)
par(mfrow=c(1,1))

source("Utils_sparse_robust_2.R")
 

fit=robust_sparse_jump(Y=res$SimData,
                                   zeta0=.2,
                                   lambda=.5,
                                   K=3,
                                   tol     = 1e-16,
                                   n_init  = 5,
                                   n_outer = 20,
                                   n_inner = 10,
                                   alpha   = 0.1,
                                   verbose = F,
                                   knn     = 10,
                                   M       = NULL,
                                   qt      = NULL,
                                   c       = NULL,
                                   mif     = 3,
                                   hd      = FALSE,
                                   n_hd    = 500,
                                   outlier = F,
                                   truth   = NULL,
                                   ncores=6)

round(fit$W,3)
