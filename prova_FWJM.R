library(scales)
source("Utils_feat_weight_robust.R")


# P=10 K=3--------------------------------------------------------------------

P=10

zeta0=.2
lambda=.3

K=3


Ktrue <- K; TT <- 1000
a1 <- c(2, 1.5, 0.1, 0.1, 1.5,2,rep(0.1,4))
a2 <- c(0.1, 0.1, 2, 1.5, 1.5,1.5,rep(0.1,4))
a3 <- c(.5, 0.1, .1, 0.1, 2,1.5,rep(0.1,4))
a_list <- list(a1, a2,a3)
mu_tilde_list <- list(rep(1, P), rep(2, P), rep(3, P))
Sigma_tilde_list <- list(diag(1, P), diag(1, P), diag(1, P))

seed=123

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

Y = as.matrix(simDat$SimData)


# Fit model
fit <- feat_weight_jump(
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
table(fit$s,truth)
balanced_accuracy(fit$s,truth)
adjustedRandIndex(fit$s,truth)

round(fit$W,2)


# P=25 K=3 ----------------------------------------------------------------


P=25

zeta0=.1
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

Y = as.matrix(simDat$SimData)


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

zeta0=.2
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

Y = as.matrix(simDat$SimData)

x11()
pairs(Y,col=truth)

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
  hd=F,
  n_hd=NULL,
  mif=5,
  truth=truth,
  ncores = 3
)


# Results
table(fit_4$s,truth)

library(clue)

tab <- table(fit_4$s, truth)
perm <- apply(tab, 1, which.max)
s_relab <- perm[fit_4$s]

balanced_accuracy(s_relab,truth)
adjustedRandIndex(fit_4$s,truth)

W_reordered <- fit_4$W[ order(perm), , drop = FALSE ]

round(W_reordered,2)
