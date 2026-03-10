source("Utils_feat_weight_robust.R")
# Common params -----------------------------------------------------------
rho=.2
nu=3

zeta0_grid_lP_sT=c(.025,.05,.075,.1,.125,.15,.175,.2,.3)
zeta0_grid_sP_lT=c(.05,.1,.15,.2,.3)
lambda_grid=c(0,.25,.5,.75,1)

seed_grid=1:200
q=50
out_frac=.05
tukey_grid=c(TRUE,FALSE)

n_init=1
n_outer=20
n_inner=10
tol_val=NULL


# K=2, P=5 ----------------------------------------------------------------
K=2
P=5
mu=1
MUs=seq(-mu,mu,length.out=K)
TT <- 1000


# More informative --------------------------------------------------------
P_true=3
P_false=2

a1 <- c(1,.5,0.1,0.1,0.1)
a2 <- c(0.1,.5,1,0.1,0.1)

a_list <- list(a1, a2)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})

#lapply(Sigma_tilde_list, isSymmetric.matrix)
seed=seed_grid[1]

simDat <- sim_hmm_SNR(seed = seed,
                      TT = TT,
                      P = P,
                      Ktrue = K,
                      a_list = a_list,            # list di lunghezza Ktrue, ciascuno vettore length P
                      mu_tilde_list = mu_tilde_list,     # list di base means
                      Sigma_tilde_list = Sigma_tilde_list,  # list di base covariance matrices
                      pers = 0.95,
                      eps = 1e-6,
                      nu=nu,
                      Out_bound = q,
                      outlier_frac = out_frac)

x11()
pairs(simDat$Y,col=simDat$truth)
summary(simDat$Y)
tapply(simDat$Y[,1],simDat$truth,summary)
tapply(simDat$Y[,2],simDat$truth,summary)
tapply(simDat$Y[,3],simDat$truth,summary)
tapply(simDat$Y[,4],simDat$truth,summary)
tapply(simDat$Y[,5],simDat$truth,summary)
Y=simDat$Y
x=simDat$truth

tukey_val=tukey_grid[1]
zeta0_val=zeta0_grid_sP_lT[3]
lambda_val=lambda_grid[1]

est=feat_weight_jump(Y=Y,
                 zeta0=zeta0_val,
                 lambda=lambda_val,
                 K=K,
                 tol     = tol_val,
                 n_init  = n_init,
                 n_outer = n_outer,
                 n_inner = n_inner,
                 alpha   = 0.1,
                 verbose = T,
                 mif     = 2,
                 truth   = x,
                 ncores=NULL,
                 tukey=tukey_val)

table(est$s,x)
round(est$W,2)

# Less informative --------------------------------------------------------
P_true=2
P_false=3

a1 <- c(1,0,0,0,0)
a2 <- c(0,1,0,0,0)

a_list <- list(a1, a2)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})

simDat <- sim_hmm_SNR(seed = seed,
                      TT = TT,
                      P = P,
                      Ktrue = K,
                      a_list = a_list,            # list di lunghezza Ktrue, ciascuno vettore length P
                      mu_tilde_list = mu_tilde_list,     # list di base means
                      Sigma_tilde_list = Sigma_tilde_list,  # list di base covariance matrices
                      pers = 0.95,
                      eps = 1e-6,
                      Out_bound = q,
                      outlier_frac = out_frac)




# K=2, P=50 ---------------------------------------------------------------
K=2
P=50
mu=1
#seq(-mu,mu,length.out=K)
MUs=seq(-mu,mu,length.out=K)
TT <- 50

# (A)
P_true=30
P_false=20

a1 <- c(rep(1,10),rep(.5,10),rep(0,30))
a2 <- c(rep(0,10),rep(.5,10),rep(1,10),rep(0,20))

a_list <- list(a1, a2)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})

# (B)
P_true=20
P_false=30

a1 <- c(rep(1,5),rep(.5,5),rep(0,40))
a2 <- c(rep(0,10),rep(.5,5),rep(1,5),rep(0,30))

a_list <- list(a1, a2)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})


# K=4, P=5 ----------------------------------------------------------------
K=4
P=5

mu=3
MUs=seq(-mu,mu,length.out=K)
TT=1000

# (A)
P_true=3
P_false=2

a1 <- c(1,.5,0,0,0)
a2 <- c(1,.5,0,0,0)
a3 <- c(0,.5,1,0,0)
a3 <- c(0,.5,1,0,0)

a_list <- list(a1, a2,a3,a4)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P), rep(MUs[3], P), rep(MUs[4], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})

# (B)
P_true=2
P_false=3

a1 <- c(1,0,0,0,0)
a2 <- c(1,0,0,0,0)
a3 <- c(0,1,0,0,0)
a3 <- c(0,1,0,0,0)

a_list <- list(a1, a2,a3,a4)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P), rep(MUs[3], P), rep(MUs[4], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})



# K=4, P=50 ---------------------------------------------------------------
K=4
P=50
mu=3
MUs=seq(-mu,mu,length.out=K)
TT=50

# (A)
P_true=30
P_false=20

a1 <- c(rep(1,10),rep(.5,10),rep(0,30))
a2 <- c(rep(1,10),rep(.5,10),rep(0,30))
a3 <- c(rep(0,10),rep(.5,10),rep(1,10),rep(0,20))
a4 <- c(rep(0,10),rep(.5,10),rep(1,10),rep(0,20))

a_list <- list(a1, a2,a3,a4)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P), rep(MUs[3], P), rep(MUs[4], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})

# (B)
P_true=20
P_false=30

a1 <- c(rep(1,5),rep(.5,5),rep(0,40))
a2 <- c(rep(1,5),rep(.5,5),rep(0,40))
a3 <- c(rep(0,10),rep(1,5),rep(.5,5),rep(0,30))
a4 <- c(rep(0,10),rep(1,5),rep(.5,5),rep(0,30))

a_list <- list(a1, a2,a3,a4)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P), rep(MUs[3], P), rep(MUs[4], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})