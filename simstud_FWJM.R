library(scales)
source("Utils_feat_weight_robust.R")


P=10

zeta0=c(0,0.05,.1,.15,.20,.25,.3,.4,.5)
zeta0[1]=.01
alpha=.1
lambda=seq(0,1,.25)

K=3

nseed=50

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  lambda=lambda,
  P=P,
  K=K
)

ncores=9

Ktrue <- K; TT <- 1000
a1 <- c(2, 1, 0.1, 0.1, 1,rep(.1,5))
a2 <- c(0.1, 0.1, 2, 1, 1,rep(.1,5))
a3 <- c(.1, 0.1, .1, 0.1, 1,rep(.1,5))
a_list <- list(a1, a2,a3)
mu_tilde_list <- list(rep(1, P), rep(2, P), rep(3, P))
Sigma_tilde_list <- list(diag(1, P), diag(1, P), diag(1, P))



start_full <- Sys.time()


res_list_K3 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    K =hp$K[i]
    
    simDat <- sim_hmm_SNR(seed = seed,
                          TT = TT,
                          P = P,
                          Ktrue = K,
                          a_list = a_list,            # list di lunghezza Ktrue, ciascuno vettore length P
                          mu_tilde_list = mu_tilde_list,     # list di base means
                          Sigma_tilde_list = Sigma_tilde_list,  # list di base covariance matrices
                          pers = 0.95,
                          eps = 1e-6)
    
    
    
    truth <- simDat$truth
    
    Y = as.matrix(simDat$SimData)
    
    fit <- feat_weight_jump(
      Y = Y,
      zeta0 = zeta0,
      lambda = lambda,
      K = K,
      tol = NULL,
      n_init = 1,
      n_outer = 25,
      n_inner=5,
      alpha = 0.1,
      verbose = T,
      hd=F,
      n_hd=NULL,
      mif=5,
      truth=truth
    )
    
    fit <- feat_weight_jump(
      Y = Y,
      zeta0 = zeta0,
      lambda = lambda,
      K = K,
      tol = 1e-2,
      n_init = 1,
      n_outer = 25,
      n_inner=5,
      alpha = 0.1,
      verbose = F,
      hd=F,
      n_hd=NULL,
      mif=5,
      truth=truth
    )
    
    end=Sys.time()
    
    est_s <- fit$s
    
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      W = fit$W,
      s = est_s,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss_vec,
      elapsed=fit$elapsed
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_K3_,elapsed_full, hp,file='simstud_FWJM_K3.Rdata')

# Results -----------------------------------------------------------------

library(scales)
source("Utils_sparse_robust_2.R")


# K=3 ---------------------------------------------------------------------

#load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K3_OUT(2).Rdata")


