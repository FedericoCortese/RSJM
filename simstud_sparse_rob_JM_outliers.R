library(scales)
source("Utils_sparse_robust_2.R")


TT=2000

P=10

zeta0=c(0,0.05,.1,.15,.20,.25,.3,.4,.5)
zeta0[1]=.01
alpha=.1
lambda=seq(0,1,.25)

# Fare prima delle prove a mano per capire i valori di c
qt=c(.95,.99)
perc_out=c(.05,.01)

nseed=50

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  lambda=lambda,
  P=P,
  qt=qt,
  perc_out=perc_out
)

ncores=20

start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

K=3

res_list_K3_OUT <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    qt <- hp$qt[i]
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
      qt=qt,
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

save(res_list_K3_OUT,elapsed_full, hp,file='simstud_SRJM_K3_OUT.Rdata')


start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,5,6)
rel_[[2]]=c(2,5,6)
rel_[[3]]=c(3,5,6)
rel_[[4]]=c(4,5,6)

K=4

res_list_K4_OUT <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    qt <- hp$qt[i]
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
      qt=qt,
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

save(res_list_K4_noOUT,elapsed_full, hp,file='simstud_SRJM_K4_OUT.Rdata')


# Results -----------------------------------------------------------------

library(scales)
source("Utils_sparse_robust_2.R")


# K=3 ---------------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K3_OUT.Rdata")

## Caso 1: perc_out = 0.05 e qt = 0.95
idx_95 <- which(hp$perc_out == 0.05 & hp$qt == 0.95)
res_list_95 <- res_list_K3_OUT[idx_95]
hp_95 <- hp[idx_95, ]

## Caso 2: perc_out = 0.01 e qt = 0.99
idx_99 <- which(hp$perc_out == 0.01 & hp$qt == 0.99)
res_list_99 <- res_list_K3_OUT[idx_99]
hp_99 <- hp[idx_99, ]

out_95 <- analyze_results(res_list_95, hp_95, P = 10, K = 3)
out_99 <- analyze_results(res_list_99, hp_99, P = 10, K = 3)

#1000x300
out_95$BAC_plot
out_99$BAC_plot

#1500X700
out_95$heatmap_plot
out_99$heatmap_plot

out_95$best_row
out_99$best_row

out_95$time_summary
out_99$time_summary


out_v_95 <- analyze_v_truth_boxgrid(res_list_95, hp_95)
#1000x800
out_v_95$grid_plot


out_v_99 <- analyze_v_truth_boxgrid(res_list_99, hp_99)
out_v_99$grid_plot

# K=4 ---------------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K4_OUT.Rdata")

## Caso 1: perc_out = 0.05 e qt = 0.95
idx_95 <- which(hp$perc_out == 0.05 & hp$qt == 0.95)
res_list_95 <- res_list_K4_OUT[idx_95]
hp_95 <- hp[idx_95, ]

## Caso 2: perc_out = 0.01 e qt = 0.99
idx_99 <- which(hp$perc_out == 0.01 & hp$qt == 0.99)
res_list_99 <- res_list_K4_OUT[idx_99]
hp_99 <- hp[idx_99, ]

out_95 <- analyze_results(res_list_95, hp_95, P = 10, K = 3)
out_99 <- analyze_results(res_list_99, hp_99, P = 10, K = 3)

#1000x300
out_95$BAC_plot
out_99$BAC_plot

#1500X700
out_95$heatmap_plot
out_99$heatmap_plot

out_95$time_summary
out_99$time_summary

out_v_95 <- analyze_v_truth_boxgrid(res_list_95, hp_95)
#1000x800
out_v_95$grid_plot


out_v_99 <- analyze_v_truth_boxgrid(res_list_99, hp_99)
out_v_99$grid_plot

save(res_list_K4_OUT,elapsed_full, hp,file='simstud_SRJM_K4_OUT.Rdata')
