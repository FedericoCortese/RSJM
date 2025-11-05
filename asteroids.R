
source("Utils_asteroids.R")

# propagation_2001GO2_new_v2

df_2001GO2=read.table("data_asteroids/propagation_2001GO2_new_v2.txt",header = T)
unique(df_2001GO2$type)

# propagation_2002AA29_new_v2
df_2002AA29=read.table("data_asteroids/propagation_2002AA29_new_v2.txt",header = T)
df_2002AA29=trans_theta(df_2002AA29)
unique(df_2002AA29$type)
plot(df_2002AA29$theta,type='p',col=df_2002AA29$type+1)


# 2014OL339 ---------------------------------------------------------------

# propagation_2014OL339_new_v2

source("Utils_asteroids.R")
df_2014OL339=read.table("data_asteroids/propagation_2014OL339_new_v2.txt",header = T)
df_2014OL339=trans_theta(df_2014OL339)

# Select relevant variables
df_2014OL339=df_2014OL339[,c("t","a","theta","type")]
#df_2014OL339$type[which(df_2014OL339$type==100)]=10

# Recode states
old_vals <- sort(unique(df_2014OL339$type))
df_2014OL339$type <- match(df_2014OL339$type, old_vals)

plot(df_2014OL339$theta,col=df_2014OL339$type)

# Compute features

df_2014OL339_maxmins_theta=max_min_feat(df_2014OL339,
                                        "theta",tt_thres_maxmin=2.5,
                                        l=5)

df_2014OL339_stat_theta=compute_feat(df_2014OL339,"theta",l=50,ma_flag   = F)

df_2014OL339_stat_a=compute_feat(df_2014OL339,"a",l=50,sd_flag   = F)
df_2014OL339_stat_a=df_2014OL339_stat_a[,c("t","ma_a")]


# Merge by t
features_2014OL339=merge(df_2014OL339,
                         df_2014OL339_maxmins_theta,by="t")
features_2014OL339=merge(features_2014OL339,
                         df_2014OL339_stat_theta,by="t")
features_2014OL339=merge(features_2014OL339,
                         df_2014OL339_stat_a,by="t")

features_2014OL339=features_2014OL339[complete.cases(features_2014OL339), ]

summary(features_2014OL339)

sel_features_2014OL339=features_2014OL339[,c("value_max_theta", "value_min_theta",
                                             "dtheta","sd_theta","sd_dtheta",
                                             #"I_TP",
                                             "I_HS",
                                             "I_QS","I_CP",
                                             "ma_a")]

# Extract ground truth and time
gt_2014OL339=features_2014OL339$type
time_2014OL339=features_2014OL339$t
a_2014OL339=features_2014OL339$a
theta_2014OL339=features_2014OL339$theta

# Artificially replicate data
n_times=3
sel_features_2014OL339=sel_features_2014OL339[rep(1:nrow(sel_features_2014OL339),
                                                  times=n_times),]

dataplot_2014OL339=data.frame(t=1:nrow(sel_features_2014OL339),
                              a=rep(a_2014OL339,n_times),
                              theta=rep(theta_2014OL339,n_times),
                              ground_truth=rep(gt_2014OL339,n_times))

plot(dataplot_2014OL339$theta,col=dataplot_2014OL339$ground_truth+1)
plot(dataplot_2014OL339$a,col=dataplot_2014OL339$ground_truth+1)
plot(scale(sel_features_2014OL339$ma_a),col=dataplot_2014OL339$ground_truth)
plot(sel_features_2014OL339$value_min_theta,col=dataplot_2014OL339$ground_truth)

sel_features_2014OL339=subset(sel_features_2014OL339,select=-c(dtheta,I_HS))

source("Utils_sparse_robust_source.R")

# Select lambda, K, and c
st=Sys.time()
cv_2014OL339=cv_robust_sparse_jump(
  Y=sel_features_2014OL339,
  true_states=dataplot_2014OL339$ground_truth,
  K_grid=2:4,
  zeta0_grid=seq(0.05,0.5,.05),
  lambda_grid=seq(0,1,.1),
  n_folds = 5,
  parallel=T,
  n_cores=detectCores()-1,
  cv_method="blocked-cv",
  knn=10,
  c_grid=10,
  M=NULL,
  hd=T,
  n_hd=1000,
  outlier=F
)
en=Sys.time()
en-st

# Select zeta0 based on the best lambda, K and c
# st=Sys.time()
# gap_2014OL339=gap_robust_sparse_jump(
#   Y=sel_features_2014OL339,
#   K_grid=NULL,
#   zeta0_grid=seq(0.05,.5,0.05),
#   lambda=0,
#   B=10,
#   parallel=F,
#   n_cores=NULL,
#   knn=10,
#   c=10,
#   M=NULL
# )
# en=Sys.time()
# en-st

# Final fit

fit_2014OL339=robust_sparse_jump(
  Y=sel_features_2014OL339,
  K=3,
  zeta0=.05,
  lambda=.7,
  c=10,
  knn=10,
  M=NULL,
  n_init=3,
  verbose=T,
  tol=1e-6,
  hd=T,
  n_hd=1000
)

est_s_2014OL339=fit_2014OL339$s
est_s_2014OL339[fit_2014OL339$v==0]=0

dataplot_2014OL339$s=est_s_2014OL339

plot(dataplot_2014OL339$theta,col=dataplot_2014OL339$s+1)
plot(dataplot_2014OL339$a,col=dataplot_2014OL339$s+1)

est_W_2014OL339= data.frame(round(fit_2014OL339$W,2))
colnames(est_W_2014OL339)=names(sel_features_2014OL339)

est_W_2014OL339


# 2015XX169 ---------------------------------------------------------------

source("Utils_asteroids.R")
df_2015XX169=read.table("data_asteroids/propagation_2015XX169_new_v2.txt",header = T)
df_2015XX169=trans_theta(df_2015XX169)

# Select relevant variables
df_2015XX169=df_2015XX169[,c("t","a","theta","type")]
#df_2015XX169$type[which(df_2015XX169$type==100)]=10

# Recode states
old_vals <- sort(unique(df_2015XX169$type))
df_2015XX169$type <- match(df_2015XX169$type, old_vals)

plot(df_2015XX169$theta,col=df_2015XX169$type)

# Compute features

df_2015XX169_maxmins_theta=max_min_feat(df_2015XX169,
                                        "theta",tt_thres_maxmin=2.5,
                                        l=5)

df_2015XX169_stat_theta=compute_feat(df_2015XX169,"theta",l=50,ma_flag   = F)

df_2015XX169_stat_a=compute_feat(df_2015XX169,"a",l=50,sd_flag   = F)
df_2015XX169_stat_a=df_2015XX169_stat_a[,c("t","ma_a")]


# Merge by t
features_2015XX169=merge(df_2015XX169,
                         df_2015XX169_maxmins_theta,by="t")
features_2015XX169=merge(features_2015XX169,
                         df_2015XX169_stat_theta,by="t")
features_2015XX169=merge(features_2015XX169,
                         df_2015XX169_stat_a,by="t")

features_2015XX169=features_2015XX169[complete.cases(features_2015XX169), ]

summary(features_2015XX169)

sel_features_2015XX169=features_2015XX169[,c("value_max_theta", "value_min_theta",
                                             "dtheta","sd_theta","sd_dtheta",
                                             #"I_TP",
                                             "I_HS",
                                             "I_QS","I_CP",
                                             "ma_a")]

# Extract ground truth and time
gt_2015XX169=features_2015XX169$type
time_2015XX169=features_2015XX169$t
a_2015XX169=features_2015XX169$a
theta_2015XX169=features_2015XX169$theta


T_origin=dim(sel_features_2015XX169)[1]

# Artificially replicate data
n_times=3
sel_features_2015XX169=sel_features_2015XX169[rep(1:nrow(sel_features_2015XX169),
                                                  times=n_times),]

dataplot_2015XX169=data.frame(t=1:nrow(sel_features_2015XX169),
                              a=rep(a_2015XX169,n_times),
                              theta=rep(theta_2015XX169,n_times),
                              ground_truth=rep(gt_2015XX169,n_times))

plot(dataplot_2015XX169$theta,col=dataplot_2015XX169$ground_truth+1)

source("Utils_sparse_robust_2.R")

# Select lambda, K, and c
st=Sys.time()
cv_2015XX169=cv_robust_sparse_jump(
  Y=sel_features_2015XX169,
  true_states=dataplot_2015XX169$ground_truth,
  K_grid=2:4,
  zeta0_grid=seq(0.05,0.5,.05),
  lambda_grid=seq(0,1,.1),
  n_folds = 5,
  parallel=T,
  n_cores=detectCores()-1,
  cv_method="blocked-cv",
  knn=10,
  c_grid=10,
  M=NULL,
  hd=T,
  n_hd=1000,
  outlier=F
)
en=Sys.time()
en-st

# Select zeta0 based on the best lambda, K and c
# st=Sys.time()
# gap_2015XX169=gap_robust_sparse_jump(
#   Y=sel_features_2015XX169,
#   K_grid=NULL,
#   zeta0_grid=seq(0.05,.5,0.05),
#   lambda=0,
#   B=10,
#   parallel=F,
#   n_cores=NULL,
#   knn=10,
#   c=10,
#   M=NULL
# )
# en=Sys.time()
# en-st

# Final fit

fit_2015XX169=robust_sparse_jump(
  Y=sel_features_2015XX169,
  K=3,
  zeta0=.20,
  lambda=.5,
  c=10,
  knn=10,
  M=NULL,
  n_init=3,
  verbose=T,
  tol=1e-6,
  hd=T,
  n_hd=1000,
  outlier=F
)

est_s_2015XX169=fit_2015XX169$s
est_s_2015XX169[fit_2015XX169$v==0]=0

dataplot_2015XX169$s=est_s_2015XX169

library(ggplot2)

dat <- dataplot_2015XX169[1:T_origin,]
dat$s_lab <- factor(dat$s, levels = c(1, 2, 3), labels = c("QS", "CP", "HS"))

library(ggplot2)

ggplot(dat, aes(x = seq_along(theta), y = theta, color = s_lab)) +
  geom_point(shape = 16, size = 2.5) + # filled points
  scale_color_manual(
    values = c(CP = "green", HS = "black", QS = "red"),
    name = NULL
  ) +
  labs(x = NULL, y = expression(theta)) +
  theme_minimal(base_size = 14) + # bigger default text size
  theme(
    axis.title = element_text(size = 16), # bigger axis labels
    axis.text = element_text(size = 14),  # bigger tick labels
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10) # add top space for legend
  )



plot(dataplot_2015XX169$a,col=dataplot_2015XX169$s+1)

est_W_2015XX169= data.frame(round(fit_2015XX169$W,2))
colnames(est_W_2015XX169)=names(sel_features_2015XX169)

est_W_2015XX169



# 2016CA138 ---------------------------------------------------------------

source("Utils_asteroids.R")
df_2016CA138=read.table("data_asteroids/propagation_2016CA138_new_v2.txt",header = T)
df_2016CA138=trans_theta(df_2016CA138)

# Select relevant variables
df_2016CA138=df_2016CA138[,c("t","a","theta","type")]
#df_2016CA138$type[which(df_2016CA138$type==100)]=10

# Recode states
old_vals <- sort(unique(df_2016CA138$type))
df_2016CA138$type <- match(df_2016CA138$type, old_vals)

plot(df_2016CA138$theta,col=df_2016CA138$type)

# Compute features

df_2016CA138_maxmins_theta=max_min_feat(df_2016CA138,
                                        "theta",tt_thres_maxmin=2.5,
                                        l=5)

df_2016CA138_stat_theta=compute_feat(df_2016CA138,"theta",l=50,ma_flag   = F)

df_2016CA138_stat_a=compute_feat(df_2016CA138,"a",l=100,sd_flag   = F)
df_2016CA138_stat_a=df_2016CA138_stat_a[,c("t","ma_a")]


# Merge by t
features_2016CA138=merge(df_2016CA138,
                         df_2016CA138_maxmins_theta,by="t")
features_2016CA138=merge(features_2016CA138,
                         df_2016CA138_stat_theta,by="t")
features_2016CA138=merge(features_2016CA138,
                         df_2016CA138_stat_a,by="t")

features_2016CA138=features_2016CA138[complete.cases(features_2016CA138), ]

summary(features_2016CA138)

sel_features_2016CA138=features_2016CA138[,c("value_max_theta", "value_min_theta",
                                             "dtheta","sd_theta","sd_dtheta",
                                             #"I_TP",
                                             "I_HS",
                                             "I_QS","I_CP",
                                             "ma_a")]

# Extract ground truth and time
gt_2016CA138=features_2016CA138$type
time_2016CA138=features_2016CA138$t
a_2016CA138=features_2016CA138$a
theta_2016CA138=features_2016CA138$theta

# Artificially replicate data
n_times=2
sel_features_2016CA138=sel_features_2016CA138[rep(1:nrow(sel_features_2016CA138),
                                                  times=n_times),]

dataplot_2016CA138=data.frame(t=1:nrow(sel_features_2016CA138),
                              a=rep(a_2016CA138,n_times),
                              theta=rep(theta_2016CA138,n_times),
                              ground_truth=rep(gt_2016CA138,n_times))

plot(dataplot_2016CA138$theta,col=dataplot_2016CA138$ground_truth+1)
plot(sel_features_2016CA138$ma_a,col=dataplot_2016CA138$ground_truth+1)


source("Utils_sparse_robust_2.R")

# Select lambda, K, and c
st=Sys.time()
cv_2016CA138=cv_robust_sparse_jump(
  Y=sel_features_2016CA138,
  true_states=dataplot_2016CA138$ground_truth,
  K_grid=2:4,
  zeta0_grid=seq(0.05,0.5,.05),
  lambda_grid=seq(0,1,.1),
  n_folds = 5,
  parallel=T,
  n_cores=detectCores()-1,
  cv_method="blocked-cv",
  knn=10,
  c_grid=10,
  M=NULL,
  hd=T,
  n_hd=1000,
  outlier=F
)
en=Sys.time()
en-st

# Select zeta0 based on the best lambda, K and c
# st=Sys.time()
# gap_2016CA138=gap_robust_sparse_jump(
#   Y=sel_features_2016CA138,
#   K_grid=NULL,
#   zeta0_grid=seq(0.05,.5,0.05),
#   lambda=0,
#   B=10,
#   parallel=F,
#   n_cores=NULL,
#   knn=10,
#   c=10,
#   M=NULL
# )
# en=Sys.time()
# en-st

# Final fit

fit_2016CA138=robust_sparse_jump(
  Y=sel_features_2016CA138,
  K=4,
  zeta0=.2,
  lambda=.8,
  c=10,
  knn=10,
  M=NULL,
  n_init=5,
  verbose=T,
  tol=1e-6,
  hd=T,
  n_hd=1000,
  outlier=F
)

est_s_2016CA138=fit_2016CA138$s
est_s_2016CA138[fit_2016CA138$v==0]=0

dataplot_2016CA138$s=est_s_2016CA138

plot(dataplot_2016CA138$theta,col=dataplot_2016CA138$s+1)
plot(dataplot_2016CA138$a,col=dataplot_2016CA138$s+1)

table(fit_2016CA138$s,dataplot_2016CA138$ground_truth)

est_W_2016CA138= data.frame(round(fit_2016CA138$W,2))
colnames(est_W_2016CA138)=names(sel_features_2016CA138)

est_W_2016CA138


# Unlabelled --------------------------------------------------------------


# 2005UH6 -----------------------------------------------------------------

source("Utils_asteroids.R")
df_2005UH6=read.table("data_asteroids/unlabelled/propagation_2005UH6_new_v2.txt",header = T)
df_2005UH6=trans_theta(df_2005UH6)

# Select relevant variables
df_2005UH6=df_2005UH6[,c("t","a","theta")]
#df_2005UH6$type[which(df_2005UH6$type==100)]=10

plot(df_2005UH6$theta)

# Compute features

df_2005UH6_maxmins_theta=max_min_feat(df_2005UH6,
                                      "theta",tt_thres_maxmin=2.5,
                                      l=5)

df_2005UH6_stat_theta=compute_feat(df_2005UH6,"theta",l=50,ma_flag   = F)

df_2005UH6_stat_a=compute_feat(df_2005UH6,"a",l=50,sd_flag   = F)
df_2005UH6_stat_a=df_2005UH6_stat_a[,c("t","ma_a")]


# Merge by t
features_2005UH6=merge(df_2005UH6,
                       df_2005UH6_maxmins_theta,by="t")
features_2005UH6=merge(features_2005UH6,
                       df_2005UH6_stat_theta,by="t")
features_2005UH6=merge(features_2005UH6,
                       df_2005UH6_stat_a,by="t")

features_2005UH6=features_2005UH6[complete.cases(features_2005UH6), ]

summary(features_2005UH6)

sel_features_2005UH6=features_2005UH6[,c("value_max_theta", "value_min_theta",
                                         "dtheta","sd_theta","sd_dtheta",
                                         #"I_TP",
                                         "I_HS",
                                         "I_QS","I_CP",
                                         "ma_a")]

# Extract ground truth and time
gt_2005UH6=features_2005UH6$type
time_2005UH6=features_2005UH6$t
a_2005UH6=features_2005UH6$a
theta_2005UH6=features_2005UH6$theta

# Artificially replicate data
n_times=3
sel_features_2005UH6=sel_features_2005UH6[rep(1:nrow(sel_features_2005UH6),
                                              times=n_times),]

dataplot_2005UH6=data.frame(t=1:nrow(sel_features_2005UH6),
                            a=rep(a_2005UH6,n_times),
                            theta=rep(theta_2005UH6,n_times),
                            ground_truth=rep(gt_2005UH6,n_times))

plot(dataplot_2005UH6$theta,col=dataplot_2005UH6$ground_truth+1)


# 100000045 ---------------------------------------------------------------

source("Utils_asteroids.R")
df_100000045=read.table("data_asteroids/unlabelled/100000045.txt",header = T)


# 100000241 ---------------------------------------------------------------

source("Utils_asteroids.R")
df_100000241=read.table("data_asteroids/unlabelled/100000241.txt",header = T)


