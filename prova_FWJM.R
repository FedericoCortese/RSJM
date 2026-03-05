library(scales)
source("Utils_feat_weight_robust.R")


# P=10 K=3--------------------------------------------------------------------

P=10

K=3

W=matrix(c(1,1,.5,0,0,0,0,0,0,0,
           1,0,1,1,0,0,0,0,0,0,
           .5,.5,0,.5,.5,0,0,0,0,0),byrow = T,nrow=3)

TT=300

P=10
outlier_frac=0.05

simDat <- sim_data_stud_t_FWJM(seed = 123,
                               TT = TT,
                               P = P,
                               Ktrue = K,
                               mu = 3,
                               rho = 0.2,
                               nu = 4,
                               pers = .95,
                               W = W,
                               outlier_frac = outlier_frac,
                               Out_bound = 100)


truth=simDat$mchain
Y = simDat$SimData

zeta0=.1
lambda=0.1

# Fit model
fit4 <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = 4,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  tukey=T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=truth
)
plot(fit4$loss_vec$loss)
table(fit4$s,simDat$mchain)

fit3 <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  tukey=T,
  verbose = T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=truth
)
plot(fit3$loss_vec$loss)
table(fit3$s,simDat$mchain)

fit2 <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = 2,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  tukey=T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=truth
)
plot(fit2$loss_vec$loss)

c(fit4$loss,fit3$loss,fit2$loss)

fit1 <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = 1,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  tukey=T,
  verbose = T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=NULL
)

plot(fit1$loss_vec$loss)


# COSA

# Fit model
fit4COSA <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = 0,
  K = 4,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  tukey=T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=truth
)
plot(fit4COSA$loss_vec$loss)
table(fit4COSA$s,simDat$mchain)

fit3COSA <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = 0,
  K = K,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  tukey=T,
  verbose = T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=truth
)
plot(fit3COSA$loss_vec$loss)
table(fit3COSA$s,simDat$mchain)

fit2COSA <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = 0,
  K = 2,
  tol = NULL,
  n_init = 1,
  n_outer = 100,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  tukey=T,
  mif=1, # mif=1 ordina gli stati in base alle mediane condizionate della feat numero 1
  truth=truth
)
plot(fit2COSA$loss_vec$loss)

c(fit4COSA$loss,fit3COSA$loss,fit2COSA$loss)

# IGNORE FROM HERE --------------------------------------------------------


# Fit with no Tukey adj

fit_no_rob <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 3,
  n_outer = 25,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  mif=5, # mif=5 ordina gli stati in base alle mediane condizionate della feat numero 5
  truth=truth,
  ncores = 3, # parallel only for Mac and Linux
  tukey=F
)

table(fit_no_rob$s,truth)
balanced_accuracy(fit_no_rob$s,truth)
# 0.71 vs 0.93 of the previous

# Fit with high zeta0

fit_zeta1 <- feat_weight_jump(
  Y = Y,
  zeta0 = 1,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 3,
  n_outer = 25,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  hd=F,
  n_hd=NULL,
  mif=5, # mif=5 ordina gli stati in base alle mediane condizionate della feat numero 5
  truth=truth,
  ncores = 3, # parallel only for Mac and Linux
  tukey=T
)

table(fit_zeta1$s,truth)
balanced_accuracy(fit_zeta1$s,truth)
#0.83

round(fit_zeta1$W,2)

# P=25 K=3 ----------------------------------------------------------------


P=25

zeta0=.125
lambda=.3

K=3


Ktrue <- K; TT <- 1000
a1 <- c(2, 1.5, 0.1, 0.1, 1.5,2,rep(0.1,19))
a2 <- c(0.1, 0.1, 2, 1.5, 1.5,1.5,rep(0.1,19))
a3 <- c(.5, 0.1, .1, 0.1, 2,1.5,rep(0.1,19))
a_list <- list(a1, a2,a3)
mu_tilde_list <- list(rep(1, P), rep(2, P), rep(3, P))
Sigma_tilde_list <- list(diag(1, P), diag(1, P), diag(1, P))

# Sim data
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

Y = simDat$SimData


# Fit model
fit_25 <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 3,
  n_outer = 25,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  hd=F,
  n_hd=NULL,
  mif=5,
  truth=truth,
  ncores = 3
)


# Results
table(fit_25$s,truth)
balanced_accuracy(fit_25$s,truth)
adjustedRandIndex(fit_25$s,truth)


round(fit_25$W,2)


# K=4 ---------------------------------------------------------------------

P=10

zeta0=.15
lambda=.3

K=4


Ktrue <- K; TT <- 1000
a1 <- c(2, 1.5, 0.1, 0.1, 1.5,2,rep(0.1,4))
a2 <- c(0.1, 0.1, 2, 1.5, 1.5,1.5,rep(0.1,4))
a3 <- c(.5, 0.1, .1, 0.1, 2,1.5,rep(0.1,4))
a4 <- c(0.1, 2, 0.1, 1.5, 0.1,1.5,rep(0.1,4))

a_list <- list(a1, a2,a3,a4)
mu_tilde_list <- list(rep(1, P), rep(2, P), rep(-1, P), rep(-2, P))
Sigma_tilde_list <- list(diag(1, P), diag(1, P), diag(1, P), diag(1, P))

# Sim data
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

Y = simDat$SimData

# Fit model
fit_4 <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 3,
  n_outer = 25,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  mif=6,
  truth=truth
)


# Results
table(fit_4$s,truth)

tab <- table(fit_4$s, truth)
perm <- apply(tab, 1, which.max)
s_relab <- perm[fit_4$s]

balanced_accuracy(s_relab,truth)
adjustedRandIndex(fit_4$s,truth)

W_reordered <- fit_4$W[ order(perm), , drop = FALSE ]

round(W_reordered,2)
