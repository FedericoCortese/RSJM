# Function to get results -------------------------------------------------

get_best_results <- function(res_list) {
  # Namespaces diretti per evitare conflitti
  dplyr  <- asNamespace("dplyr")
  purrr  <- asNamespace("purrr")
  tidyr  <- asNamespace("tidyr")
  ggplot2 <- asNamespace("ggplot2")
  scales <- asNamespace("scales")
  
  # --- Estrai solo i campi scalari e i pesi
  df <- purrr$map_dfr(res_list, function(x) {
    data.frame(
      lambda  = x$lambda,
      kappa   = x$kappa,
      K       = x$K,
      BAC_s   = x$BAC_s,
      elapsed = as.numeric(x$elapsed),
      W       = I(list(x$W))
    )
  })
  
  # --- Riassunto generale per (kappa, lambda, K)
  res_summary <- dplyr$group_by(df, kappa, lambda, K)
  res_summary <- dplyr$summarise(
    res_summary,
    median_BAC  = median(BAC_s, na.rm = TRUE),
    q025        = quantile(BAC_s, 0.025, na.rm = TRUE),
    q975        = quantile(BAC_s, 0.975, na.rm = TRUE),
    median_time = median(elapsed, na.rm = TRUE),
    median_W    = list(apply(do.call(rbind, W), 2, median, na.rm = TRUE)),
    .groups = "drop"
  )
  
  # --- Tabella 1: risultati SJM / JM / k-means
  best_SJM <- dplyr$slice_max(res_summary, median_BAC, n = 1, with_ties = FALSE)
  best_SJM$Method <- "SJM"
  
  max_kappa <- max(df$kappa, na.rm = TRUE)
  
  best_JM <- dplyr$filter(res_summary, kappa == max_kappa)
  best_JM <- dplyr$slice_max(best_JM, median_BAC, n = 1, with_ties = FALSE)
  best_JM$Method <- "JM"
  
  best_KMEANS <- dplyr$filter(res_summary, kappa == max_kappa & lambda == 0)
  best_KMEANS <- dplyr$slice_max(best_KMEANS, median_BAC, n = 1, with_ties = FALSE)
  best_KMEANS$Method <- "k-means"
  
  tab_summary <- dplyr$bind_rows(best_SJM, best_JM, best_KMEANS)
  tab_summary <- tab_summary[, c("Method", "kappa", "lambda", "K",
                                 "median_BAC", "q025", "q975", "median_time")]
  
  # --- Tabella 2: pesi SJM
  best_per_kappa <- dplyr$group_by(res_summary, kappa)
  best_per_kappa <- dplyr$slice_max(best_per_kappa, median_BAC, n = 1, with_ties = FALSE)
  best_per_kappa <- dplyr$ungroup(best_per_kappa)
  
  best_per_kappa$median_W <- lapply(best_per_kappa$median_W, function(w) {
    if (sum(w, na.rm = TRUE) == 0 || all(!is.finite(w))) rep(NA_real_, length(w))
    else w / sum(w, na.rm = TRUE)
  })
  
  # --- Long format per heatmap
  W_long <- dplyr$mutate(best_per_kappa, id = seq_len(nrow(best_per_kappa)))
  W_long <- W_long[, c("id", "kappa", "lambda", "median_W")]
  W_long <- tidyr$unnest_longer(W_long, median_W, values_to = "weight")
  W_long <- dplyr$group_by(W_long, kappa)
  W_long <- dplyr$mutate(W_long, variable = seq_len(dplyr$n()))
  W_long <- dplyr$ungroup(W_long)
  
  # --- Arrotonda κ, rimuovi κ = 1
  W_long$kappa <- round(W_long$kappa, 2)
  W_long <- dplyr$filter(W_long, kappa != 1)
  
  # --- Pulizia: rimpiazza NA e clamp nel range [min,max]
  eps <- 1e-8
  W_long$weight[is.na(W_long$weight)] <- 0
  min_w <- round(min(W_long$weight, na.rm = TRUE) - eps, 2)
  max_w <- round(max(W_long$weight, na.rm = TRUE) + eps, 2)
  W_long$weight <- pmin(pmax(W_long$weight, min_w), max_w)
  
  # --- Heatmap: stessa palette di analyze_results
  heatmap_plot <- ggplot2$ggplot(W_long, ggplot2$aes(
    x = factor(variable),
    y = factor(kappa),
    fill = weight
  )) +
    ggplot2$geom_tile(color = "grey80") +
    ggplot2$geom_text(
      ggplot2$aes(label = ifelse(is.na(weight), "", sprintf("%.2f", weight))),
      color = "black",
      size = 3
    ) +
    ggplot2$scale_fill_gradientn(
      colours = c("white", "yellow", "orange", "red2"),
      values = scales$rescale(c(min_w, 0.05, 0.10, max_w)),
      limits = c(min_w, max_w),
      breaks = c(min_w, 0.05, 0.10, max_w),
      labels = scales$number_format(accuracy = 0.02),
      guide = "none"
    ) +
    ggplot2$labs(
      x = "Variable",
      y = expression(kappa)
      # ,
      # title = "Median normalized weights per κ"
    ) +
    ggplot2$theme_minimal(base_size = 14) +
    ggplot2$theme(
      legend.position = "none",
      plot.title = ggplot2$element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2$element_text(
        size = 11,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = ggplot2$element_text(size = 11)
    )
  
  # --- Output finale
  list(
    summary = tab_summary,
    weights_long = W_long,
    heatmap_plot = heatmap_plot
  )
}


library(parallel)
library(Rcpp)
#library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)
library(MixSim)

Rcpp::sourceCpp("sjm.cpp")
source("Utils_sparse_robust_2.R")


# No outliers ------------------------------------------------------------------

TT=2000

P=c(10,25)

kappa10=seq(1,sqrt(P[1]),length.out=9)
kappa25=seq(1,sqrt(P[2]),length.out=9)

lambda_grid=seq(0,50,length.out=30)

nseed=50

ncores=20

K=3
rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

start_full <- Sys.time()
hp=expand.grid(
  seed=1:nseed,
  kappa=kappa10,
  lambda=lambda_grid
  #,c=c
)

res_list_SJM_K3_P10 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    kappa <- hp$kappa[i]
    lambda <- hp$lambda[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P[1],
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
      perc_out = 0,
      out_sigma = 0,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    st=Sys.time()
    fit <- sparse_jump(Y=as.matrix(simDat_sparse$Y),
                       n_states=K,
                       max_features=kappa,
                       jump_penalty=lambda)
    end=Sys.time()
    
    est_s <- fit$states
    feat_w <- fit$feat_w/sum(fit$feat_w)
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      kappa = kappa,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = feat_w,
      s = est_s,
      #v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss,
      elapsed=end-st
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_SJM_K3_P10,elapsed_full, hp,file='res_list_SJM_K3_P10.Rdata')

start_full <- Sys.time()
hp=expand.grid(
  seed=1:nseed,
  kappa=kappa25,
  lambda=lambda_grid
  #,c=c
)

res_list_SJM_K3_P25 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    kappa <- hp$kappa[i]
    lambda <- hp$lambda[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P[2],
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
      perc_out = 0,
      out_sigma = 0,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    
    st=Sys.time()
    fit <- sparse_jump(Y=as.matrix(simDat_sparse$Y),
                       n_states=K,
                       max_features=kappa,
                       jump_penalty=lambda)
    end=Sys.time()
    
    est_s <- fit$states
    feat_w <- fit$feat_w/sum(fit$feat_w)
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      kappa = kappa,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = feat_w,
      s = est_s,
      #v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss,
      elapsed=end-st
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_SJM_K3_P25,elapsed_full, hp,file='res_list_SJM_K3_P25.Rdata')


rel_=list()
rel_[[1]]=c(1,5,6)
rel_[[2]]=c(2,5,6)
rel_[[3]]=c(3,5,6)
rel_[[4]]=c(4,5,6)

K=4

start_full <- Sys.time()

hp=expand.grid(
  seed=1:nseed,
  kappa=kappa10,
  lambda=lambda_grid
  #,c=c
)

res_list_SJM_K4_P10 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    kappa <- hp$kappa[i]
    lambda <- hp$lambda[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P[1],
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
      perc_out = 0,
      out_sigma = 0,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    
    st=Sys.time()
    fit <- sparse_jump(Y=as.matrix(simDat_sparse$Y),
                       n_states=K,
                       max_features=kappa,
                       jump_penalty=lambda)
    end=Sys.time()
    
    est_s <- fit$states
    feat_w <- fit$feat_w/sum(fit$feat_w)
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      kappa = kappa,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = feat_w,
      s = est_s,
      #v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss,
      elapsed=end-st
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_SJM_K4_P10,elapsed_full, hp,file='res_list_SJM_K4_P10.Rdata')

start_full <- Sys.time()
hp=expand.grid(
  seed=1:nseed,
  kappa=kappa25,
  lambda=lambda_grid
  #,c=c
)

res_list_SJM_K4_P25 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    kappa <- hp$kappa[i]
    lambda <- hp$lambda[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P[2],
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
      perc_out = 0,
      out_sigma = 0,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    
    st=Sys.time()
    fit <- sparse_jump(Y=as.matrix(simDat_sparse$Y),
                       n_states=K,
                       max_features=kappa,
                       jump_penalty=lambda)
    end=Sys.time()
    
    est_s <- fit$states
    feat_w <- fit$feat_w/sum(fit$feat_w)
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      kappa = kappa,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = feat_w,
      s = est_s,
      #v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss,
      elapsed=end-st
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_SJM_K4_P25,elapsed_full, hp,file='res_list_SJM_K4_P25.Rdata')


end_full <- Sys.time()


# Results -----------------------------------------------------------------

# K3 P10

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/res_list_SJM_K3_P10.Rdata")
best_results_K3_P10 <- get_best_results(res_list_SJM_K3_P10)
best_results_K3_P10$summary

#900X400
best_results_K3_P10$heatmap_plot

# K3 P25
load("D:/CNR/OneDrive - CNR/simres_robusts_parse/res_list_SJM_K3_P25.Rdata")
best_results_K3_P25 <- get_best_results(res_list_SJM_K3_P25)
best_results_K3_P25$summary
#1100X400
best_results_K3_P25$heatmap_plot

# K4 P10
load("D:/CNR/OneDrive - CNR/simres_robusts_parse/res_list_SJM_K4_P10.Rdata")
best_results_K4_P10 <- get_best_results(res_list_SJM_K4_P10)
best_results_K4_P10$summary

#900X400
best_results_K4_P10$heatmap_plot

# K4 P25
load("D:/CNR/OneDrive - CNR/simres_robusts_parse/res_list_SJM_K4_P25.Rdata")
best_results_K4_P25 <- get_best_results(res_list_SJM_K4_P25)
best_results_K4_P25$summary

#1100X400
best_results_K4_P25$heatmap_plot



# Outliers ----------------------------------------------------------------

TT=2000

ncores=20
P=10

kappa10=seq(1,sqrt(P),length.out=9)

lambda_grid=seq(0,50,length.out=30)

qt=c(.95,.99)
nseed=50

ncores=20

K_true=3
rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

start_full <- Sys.time()
K=c(3,4)
hp=expand.grid(
  seed=1:nseed,
  kappa=kappa10,
  lambda=lambda_grid,
  qt=qt,
  K=K
)

res_list_SJM_K3_out <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    kappa <- hp$kappa[i]
    lambda <- hp$lambda[i]
    K=hp$K[i]
    qt=hp$qt[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P,
                              Pcat=NULL,
                              Ktrue=K_true,
                              mu=3,
                              rho=0,
                              nu=10,
                              pers=.95)
    
    simDat_sparse <- simulate_sparse_hmm(
      Y = simDat$SimData,
      rel_ = rel_,
      true_stat = simDat$mchain,
      perc_out = 1-qt,
      out_sigma = 5,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    st=Sys.time()
    fit <- sparse_jump(Y=as.matrix(simDat_sparse$Y),
                       n_states=K,
                       max_features=kappa,
                       jump_penalty=lambda)
    end=Sys.time()
    
    est_s <- fit$states
    feat_w <- fit$feat_w/sum(fit$feat_w)
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      kappa = kappa,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = feat_w,
      s = est_s,
      #v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss,
      elapsed=end-st
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()
elapsed_full <- end_full - start_full
save(res_list_SJM_K3_out,elapsed_full, hp,file='res_list_SJM_K3_out.Rdata')

rel_=list()
rel_[[1]]=c(1,5,6)
rel_[[2]]=c(2,5,6)
rel_[[3]]=c(3,5,6)
rel_[[4]]=c(4,5,6)

K_true=4

start_full <- Sys.time()

K=c(4,5)
hp=expand.grid(
  seed=1:nseed,
  kappa=kappa10,
  lambda=lambda_grid,
  qt=qt,
  K=K
)

res_list_SJM_K4_out <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    kappa <- hp$kappa[i]
    lambda <- hp$lambda[i]
    K=hp$K[i]
    qt=hp$qt[i]
    
    simDat <- sim_data_stud_t(seed=seed,
                              TT=TT,
                              P=P,
                              Pcat=NULL,
                              Ktrue=K_true,
                              mu=3,
                              rho=0,
                              nu=10,
                              pers=.95)
    
    simDat_sparse <- simulate_sparse_hmm(
      Y = simDat$SimData,
      rel_ = rel_,
      true_stat = simDat$mchain,
      perc_out = 1-qt,
      out_sigma = 5,
      seed = seed
    )
    
    truth <- simDat_sparse$truth
    
    st=Sys.time()
    fit <- sparse_jump(Y=as.matrix(simDat_sparse$Y),
                       n_states=K,
                       max_features=kappa,
                       jump_penalty=lambda)
    end=Sys.time()
    
    est_s <- fit$states
    feat_w <- fit$feat_w/sum(fit$feat_w)
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      kappa = kappa,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = feat_w,
      s = est_s,
      #v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      BAC_s=BAC_s,
      loss=fit$loss,
      elapsed=end-st
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end_full <- Sys.time()


elapsed_full=end_full-start_full

save(res_list_SJM_K4_out,elapsed_full, hp,file='res_list_SJM_K4_out.Rdata')


# Results -----------------------------------------------------------------


# Ktrue=3 -----------------------------------------------------------------
load("D:/CNR/OneDrive - CNR/simres_robusts_parse/res_list_SJM_K3_out.Rdata")

# Ktrue=3 K=3 qt=95
idx_3_3_95 <- which(hp$K == 3 & hp$qt == 0.95)
res_list_SJM_K3_3_95 <- res_list_SJM_K3_out[idx_3_3_95]
res_SJM_K3_3_95=get_best_results(res_list_SJM_K3_3_95)
res_SJM_K3_3_95$summary
res_SJM_K3_3_95$heatmap_plot


# Ktrue=3 K=3 qt=99
idx_3_3_99 <- which(hp$K == 3 & hp$qt == 0.99)
res_list_SJM_K3_3_99 <- res_list_SJM_K3_out[idx_3_3_99]
res_SJM_K3_3_99=get_best_results(res_list_SJM_K3_3_99)
res_SJM_K3_3_99$summary
 
# NON NE VALE LA PENA DI FARE IL SEGUENTE I RISULTATI SONO IDENTICI
# Ktrue=3 K=4 qt=95
idx_3_4_95 <- which(hp$K == 4 & hp$qt == 0.95)
res_list_SJM_K3_4_95 <- res_list_SJM_K3_out[idx_3_4_95]
res_SJM_K3_4_95=get_best_results(res_list_SJM_K3_4_95)
res_SJM_K3_4_95$summary
res_SJM_K3_4_95$heatmap_plot

# Ktrue=3 K=4 qt=99
idx_3_4_99 <- which(hp$K == 4 & hp$qt == 0.99)
res_list_SJM_K3_4_99 <- res_list_SJM_K3_out[idx_3_4_99]
res_SJM_K3_4_99=get_best_results(res_list_SJM_K3_4_99)
res_SJM_K3_4_99$summary


# Ktrue=4 -----------------------------------------------------------------

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/res_list_SJM_K4_out.Rdata")

# Ktrue=4 K=4 qt=95
idx_4_4_95 <- which(hp$K == 4 & hp$qt == 0.95)
res_list_SJM_K4_4_95 <- res_list_SJM_K4_out[idx_4_4_95]
res_SJM_K4_4_95=get_best_results(res_list_SJM_K4_4_95)
res_SJM_K4_4_95$summary
res_SJM_K4_4_95$heatmap_plot


# Ktrue=4 K=4 qt=99
idx_4_4_99 <- which(hp$K == 4 & hp$qt == 0.99)
res_list_SJM_K4_4_99 <- res_list_SJM_K4_out[idx_4_4_99]
res_SJM_K4_4_99=get_best_results(res_list_SJM_K4_4_99)
res_SJM_K4_4_99$summary

# NON NE VALE LA PENA DI FARE IL SEGUENTE I RISULTATI SONO IDENTICI
# Ktrue=4 K=5 qt=95
idx_4_5_95 <- which(hp$K == 5 & hp$qt == 0.95)
res_list_SJM_K4_5_95 <- res_list_SJM_K4_out[idx_4_5_95]
res_SJM_K4_5_95=get_best_results(res_list_SJM_K4_5_95)
res_SJM_K4_5_95$summary

# Ktrue=4 K=5 qt=99
idx_4_5_99 <- which(hp$K == 5 & hp$qt == 0.99)
res_list_SJM_K4_5_99 <- res_list_SJM_K4_out[idx_4_5_99]
res_SJM_K4_5_99=get_best_results(res_list_SJM_K4_5_99)
res_SJM_K4_5_99$summary
