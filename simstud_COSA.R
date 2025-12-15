library(scales)
source("Utils_sparse_robust_2.R")


TT=2000

P=10

zeta0=c(0,0.05,.1,.15,.20,.25,.3,.4,.5)
zeta0[1]=.01
alpha=.1


perc_out=c(.05,.01)

nseed=50

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  P=P,
  perc_out=perc_out
)

ncores=20

start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

K=3

res_list_K3_OUT_COSA <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- 0
    P <- hp$P[i]
    perc_out <- hp$perc_out[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P,
                              Pcat=NULL,
                              Ktrue=K,
                              mu=3,
                              rho=0,
                              nu=10,
                              pers=.95)
    
    simDat_sparse <- simulate_sparse_hmm(
      Y = simDat$SimData,
      rel_ = rel_,
      true_stat = simDat$mchain,
      perc_out = perc_out,
      out_sigma = 5,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    
    fit <- robust_sparse_jump(
      Y = as.matrix(simDat_sparse$Y),
      zeta0 = zeta0,
      lambda = lambda,
      K = K,
      tol = 1e-2,
      n_init = 1,
      n_outer = 25,
      n_inner=5,
      alpha = 0.1,
      verbose = F,
      knn = 10,
      qt=NULL,
      c = NULL,
      M = NULL,
      hd=T,
      n_hd=500,
      outlier=F,
      mif=5,
      truth=truth
    )
    
    end=Sys.time()
    
    est_s <- fit$s
    est_s[fit$v==0]=0
    
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      c = c,
      W = fit$W,
      s = est_s,
      v = fit$v,
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

save(res_list_K3_OUT_COSA,elapsed_full, hp,file='simstud_SRJM_K3_OUT_COSA.Rdata')


start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,5,6)
rel_[[2]]=c(2,5,6)
rel_[[3]]=c(3,5,6)
rel_[[4]]=c(4,5,6)

K=4

res_list_K4_OUT_COSA <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    P <- hp$P[i]
    perc_out <- hp$perc_out[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P,
                              Pcat=NULL,
                              Ktrue=K,
                              mu=3,
                              rho=0,
                              nu=10,
                              pers=.95)
    
    simDat_sparse <- simulate_sparse_hmm(
      Y = simDat$SimData,
      rel_ = rel_,
      true_stat = simDat$mchain,
      perc_out = perc_out,
      out_sigma = 5,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    
    fit <- robust_sparse_jump(
      Y = as.matrix(simDat_sparse$Y),
      zeta0 = zeta0,
      lambda = 0,
      K = K,
      tol = 1e-2,
      n_init = 1,
      n_outer = 25,
      n_inner=5,
      alpha = 0.1,
      verbose = F,
      knn = 10,
      qt=NULL,
      c = NULL,
      M = NULL,
      hd=T,
      n_hd=500,
      outlier=T,
      mif=5,
      truth=truth
    )
    
    end=Sys.time()
    
    est_s <- fit$s
    est_s[fit$v==0]=0
    
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      TT = TT,
      P = P,
      W = fit$W,
      s = est_s,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      elapsed=fit$elapsed
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_K4_OUT_COSA,elapsed_full, hp,file='simstud_SRJM_K4_OUT_COSA.Rdata')

# Results COSA ------------------------------------------------------------

source("Utils_sparse_robust_2.R")

# K=3 ---------------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K3_OUT_COSA.Rdata")

# Clean results
outK3_COSA <- summarize_COSA(res_list_K3_OUT_COSA, hp)
outK3_COSA

# K=4 ---------------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K4_OUT_COSA.Rdata")

outK4_COSA <- summarize_COSA(res_list_K4_OUT_COSA, hp)
outK4_COSA
