library(scales)
source("Utils_sparse_robust_2.R")


# T=2000 ------------------------------------------------------------------

TT=2000

P=c(10,25)

zeta0=c(0,0.05,.1,.15,.20,.25,.3,.4,.5)
zeta0[1]=.01
alpha=.1
lambda=seq(0,1,.25)

nseed=50

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  lambda=lambda,
  P=P
  #,c=c
)

ncores=parallel::detectCores()-1

start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

K=3

res_list_K3_noOUT <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    
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
      perc_out = 0,
      out_sigma = 0,
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
      knn = NULL,
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
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = fit$W,
      s = est_s,
      #v = fit$v,
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

save(res_list_K3_noOUT,elapsed_full, hp,file='simstud_SRJM_K3_noOUT.Rdata')


start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,5,6)
rel_[[2]]=c(2,5,6)
rel_[[3]]=c(3,5,6)
rel_[[4]]=c(4,5,6)

K=4

res_list_K4_noOUT <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    
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
      perc_out = 0,
      out_sigma = 0,
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
      knn = NULL,
      c = NULL,
      M = NULL,
      hd=T,
      n_hd=500,
      outlier=F,
      mif=5,
      truth=truth
    )

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
      #c = c,
      W = fit$W,
      s = est_s,
      #v = fit$v,
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

save(res_list_K4_noOUT,elapsed_full, hp,file='simstud_SRJM_K4_noOUT.Rdata')


# heatmap - Truth ---------------------------------------------------------

# ---- Libraries ----
# install.packages(c("ggplot2","reshape2"))
library(ggplot2)
library(reshape2)

# ============== K = 3 ==============

# pattern: rows are the three states, columns 1..10
row1 <- c("W","W","O","O","O","W","W","W","W","W")   # "WWOOOWWWW..."
row2 <- c("W","O","W","O","O","W","W","W","W","W")   # "WOWOOWWWW..."
row3 <- c("O","W","W","O","O","W","W","W","W","W")   # "OWWOOWWWW..."

mat  <- rbind(row1,row2,row3)
df <- melt(mat)
names(df) <- c("State","Variable","val")
df$State <- factor(df$State, levels = rev(unique(df$State)))  # top-to-bottom
df$Variable <- as.integer(df$Variable)

# colors: orange for O, white for W
fill_map <- c(W = "white", O = "#f39c12")

p3 <- ggplot(df, aes(x = Variable, y = State, fill = val)) +
  geom_tile(color = "grey90", size = 0.3) +
  scale_fill_manual(values = fill_map, guide = "none") +
  coord_fixed() +
  # X labels: 1..9 then blank for the 10th tick
  scale_x_continuous(breaks = 1:10, labels = c(1:9, " ")) +
  scale_y_discrete(labels = c("1","2","3")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(x = "Variable", y = "State")

# add "..." text inside the 10th column tiles
p3 <- p3 + geom_text(
  data = subset(df, Variable == 10),
  aes(label = "..."),
  size = 4
)

p3


# ============== K = 4 ==============

# pattern: rows are the four states, columns 1..10
row1 <- c("W","W","W","O","O","O","W","W","W","W")   # "WWWOOOWWWW..."
row2 <- c("W","W","O","W","O","O","W","W","W","W")   # "WWOWOOWWWW..."
row3 <- c("W","O","W","W","O","O","W","W","W","W")   # "WOWWOOWWWW..."
row4 <- c("O","W","W","W","O","O","W","W","W","W")   # "OWWWOOWWWW..."

mat <- rbind(row1,row2,row3,row4)
df4 <- melt(mat)
names(df4) <- c("State","Variable","val")
df4$State <- factor(df4$State, levels = rev(unique(df4$State)))  # top-to-bottom
df4$Variable <- as.integer(df4$Variable)

p4 <- ggplot(df4, aes(x = Variable, y = State, fill = val)) +
  geom_tile(color = "grey90", size = 0.3) +
  scale_fill_manual(values = fill_map, guide = "none") +
  coord_fixed() +
  # X labels: 1..9 then blank for the 10th tick
  scale_x_continuous(breaks = 1:10, labels = c(1:9, " ")) +
  scale_y_discrete(labels = 1:4) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(x = "Variable", y = "State")

# add "..." text inside the 10th column tiles
p4 <- p4 + geom_text(
  data = subset(df4, Variable == 10),
  aes(label = "..."),
  size = 4
)

p4


# Results -----------------------------------------------------------------
source("Utils_sparse_robust_2.R")

# K=3 ----------------------------------------------------------------

# P=10 
load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K3_noOUT.Rdata")

res_K3_P10=analyze_results(res_list_K3_noOUT, hp, P=10, K=3, 
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12
                           )


# CONSIGLIO: salva le heatmap in jpeg con dimensione 1500X700 per P=10
# e 2000X1000 per P=25

res_K3_P10$best_row
res_K3_P10$time_summary

#1000X300
x11()
res_K3_P10$BAC_plot

#1500X700
x11()
res_K3_P10$heatmap_plot

# P=25 
res_K3_P25=analyze_results(res_list_K3_noOUT, hp, P=25, K=3,
                           #label_size = 2,
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12
                           # ,
                           # x_axis_font_size = 10
                           )

res_K3_P25$best_row
res_K3_P25$time_summary

x11()
res_K3_P25$BAC_plot

#2000X1000
#x11()
jpeg("hmK3P25.jpeg", height = 21, width = 39.7, units = "cm", res = 300)
res_K3_P25$heatmap_plot
dev.off()

# K=4  ---------------------------------------------------------------

# P=10

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K4_noOUT.Rdata")

res_K4_P10=analyze_results(res_list_K4_noOUT, hp, P=10, K=4, 
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12
)

res_K4_P10$best_row
res_K4_P10$time_summary

x11()
res_K4_P10$BAC_plot

#1500X700
x11()
res_K4_P10$heatmap_plot

res_K4_P10$time_summary

#P=25

res_K4_P25=analyze_results(res_list_K4_noOUT, hp, P=25, K=4,
                           #label_size = 2,
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12
                           # ,
                           # x_axis_font_size = 10
                           # 
                           )


res_K4_P25$best_row
res_K4_P25$time_summary

x11()
res_K4_P25$BAC_plot
x11()
jpeg("hmK4P25.jpeg", height = 21, width = 39.7, units = "cm", res = 300)
res_K4_P25$heatmap_plot
dev.off()

res_K4_P25$time_summary


#Cairo::CairoPDF("W_heatmap_K3_P10_lambda1_c10.pdf", width = 12, height = 12)


#dev.off()


# T=4000 ------------------------------------------------------------------

# Tutto identico a prima, ora vogliamo vedere se il subsampling impatta sui risultati
TT=4000

P=c(10,25)

zeta0=c(0,0.05,.1,.15,.20,.25,.3,.4,.5)
zeta0[1]=.01
alpha=.1
lambda=seq(0,1,.25)

nseed=50

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  lambda=lambda,
  P=P
  #,c=c
)

ncores=20

start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

K=3

res_list_K3_noOUT_T4000 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    
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
      perc_out = 0,
      out_sigma = 0,
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
      knn = NULL,
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
    truth <- simDat_sparse$truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    BAC_s = balanced_accuracy(est_s, truth)
    
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      #c = c,
      W = fit$W,
      s = est_s,
      #v = fit$v,
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

save(res_list_K3_noOUT_T4000,elapsed_full, hp,file='simstud_SRJM_K3_noOUT_T4000.Rdata')


start_full <- Sys.time()
rel_=list()
rel_[[1]]=c(1,5,6)
rel_[[2]]=c(2,5,6)
rel_[[3]]=c(3,5,6)
rel_[[4]]=c(4,5,6)

K=4

res_list_K4_noOUT_T4000 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    lambda <- hp$lambda[i]
    P <- hp$P[i]
    
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
      perc_out = 0,
      out_sigma = 0,
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
      knn = NULL,
      c = NULL,
      M = NULL,
      hd=T,
      n_hd=500,
      outlier=F,
      mif=5,
      truth=truth
    )
    
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
      #c = c,
      W = fit$W,
      s = est_s,
      #v = fit$v,
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

save(res_list_K4_noOUT_T4000,elapsed_full, hp,file='simstud_SRJM_K4_noOUT_T4000.Rdata')

# Results -----------------------------------------------------------------
source("Utils_sparse_robust_2.R")

# K=3 ----------------------------------------------------------------

# P=10 
load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K3_noOUT_T4000.Rdata")

res_K3_P10=analyze_results(res_list_K3_noOUT_T4000, hp, P=10, K=3, 
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12
)


# CONSIGLIO: salva le heatmap in jpeg con dimensione 1500X700 per P=10
# e 2000X1000 per P=25

res_K3_P10$best_row
res_K3_P10$time_summary

#1000X300
x11()
res_K3_P10$BAC_plot

#1500X700
x11()
res_K3_P10$heatmap_plot

# P=25 
res_K3_P25=analyze_results(res_list_K3_noOUT_T4000, hp, P=25, K=3,
                           label_size = 2,
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12,
                           x_axis_font_size = 10)

res_K3_P25$best_row
res_K3_P25$time_summary

x11()
res_K3_P25$BAC_plot

#2000X1000
#x11()
jpeg("hmK3P25_T4000.jpeg", height = 21, width = 39.7, units = "cm", res = 300)
res_K3_P25$heatmap_plot
dev.off()

# K=4  ---------------------------------------------------------------

# P=10

load("D:/CNR/OneDrive - CNR/simres_robusts_parse/simstud_SRJM_K4_noOUT_T4000.Rdata")

res_K4_P10=analyze_results(res_list_K4_noOUT_T4000, hp, P=10, K=4, 
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12
)

res_K4_P10$best_row
res_K4_P10$time_summary

x11()
res_K4_P10$BAC_plot

#1500X700
x11()
res_K4_P10$heatmap_plot

res_K4_P10$time_summary

#P=25

res_K4_P25=analyze_results(res_list_K4_noOUT_T4000, hp, P=25, K=4,
                           label_size = 2,
                           ylim_BAC = c(0.25, 1), 
                           show_legend = F,
                           facet_font_size = 12,
                           x_axis_font_size = 10)


res_K4_P25$best_row
res_K4_P25$time_summary

x11()
res_K4_P25$BAC_plot
x11()
jpeg("hmK4P25_T4000.jpeg", height = 21, width = 39.7, units = "cm", res = 300)
res_K4_P25$heatmap_plot
dev.off()

res_K4_P25$time_summary