load("KMB.RData")

head(counts)
matplot(x=ti,y=counts,type='l')

summary(counts)

counts_2=data.frame(counts[,1:2])
counts_2$X3=counts[,3]+counts[,4]

matplot(counts_2,type='l')

library(zoo)

w <- 5  

ma <- rollapply(counts_2, 
                width = w, 
                FUN = mean, 
                align = "right", 
                fill = NA)
colnames(ma) <- paste0("MA_", colnames(counts_2))

msd <- rollapply(counts_2,
                 width = w,
                 FUN = sd,
                 align = "right",
                 fill = NA)
colnames(msd) <- paste0("MSD_", colnames(counts_2))

matplot(ma,type='l')
matplot(msd,type='l')

data_input=data.frame(ma,msd)
data_input=data_input[complete.cases(data_input),]

source("Utils_feat_weight_robust.R")

zeta0=.2
lambda=.3
K=3

fit_kmb <- feat_weight_jump(
  Y = data_input,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 1,
  n_outer = 25,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  mif=1)


plot(fit_kmb$s,type='l')

res_=data.frame(counts_2[-(1:(w-1)),],cluster=fit_kmb$s)

tapply(res_$X1,res_$cluster,mean)
tapply(res_$X2,res_$cluster,mean)
tapply(res_$X3,res_$cluster,mean)

est_feat=fit_kmb$W
colnames(est_feat)=colnames(data_input)
round(est_feat,2)
