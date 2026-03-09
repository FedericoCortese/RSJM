rho=.2
# K=2, P=5 ----------------------------------------------------------------
K=2
P=5
mu=1
MUs=seq(-mu,mu,length.out=K)
TT <- 1000

# (A)
P_true=3
P_false=2

a1 <- c(1,1,0,0,0)
a2 <- c(0,1,1,0,0)

a_list <- list(a1, a2)

mu_tilde_list <- list(rep(MUs[1], P), rep(MUs[2], P))
Sigma_tilde_list <- list(matrix(rho,P,P), matrix(rho,P,P))
Sigma_tilde_list <- lapply(Sigma_tilde_list, function(x) {
  diag(x) <- rep(1,P)
  x
})

#lapply(Sigma_tilde_list, isSymmetric.matrix)

# (B)
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
a2 <- c(rep(0,10),rep(1,5),rep(.5,5),rep(0,30))

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

a1 <- c(1,1,0,0,0)
a2 <- c(1,1,0,0,0)
a3 <- c(0,1,1,0,0)
a3 <- c(0,1,1,0,0)

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