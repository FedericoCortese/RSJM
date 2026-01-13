gdp=readxl::read_excel("GDP2_cleaned.xlsx")

library(dplyr)
library(tidyr)
gdp_long <- gdp %>%
  pivot_longer(
    cols = -State,
    names_to  = "time",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = State,
    values_from = value
  ) %>%
  arrange(time)

matplot(gdp_long[,-1],type='l')

summary(gdp_long)


source("Utils_sparse_robust_2.R")

zeta0=.3
lambda=1
K=3
qt=99

Y=as.matrix(gdp_long[,-1])
Y=apply(Y,2,diff)
Y=scale(Y)

matplot(Y,type='l')

fit <- robust_sparse_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = 1e-16,
  n_init = 5,
  n_outer = 100,
  n_inner=5,
  alpha = 0.1,
  verbose = T,
  knn = 10,
  qt=qt,
  c = NULL,
  M = NULL,
  hd=F,
  n_hd=500,
  outlier=F,
  mif=5
)

plot(fit$loss_vec$loss,type='l')


library(zoo)

time_q <- as.yearqtr(gdp_long$time, format = "%Y-Q%q")

plot(
  x = time_q[-1],
  y = fit$s,
  type = "l",
  xlab = "Time",
  ylab = "State",
  main = "Latent state over time"
)

est_W=round(fit$W,3)
colnames(est_W)=colnames(gdp_long)[-1]

est_W*100
