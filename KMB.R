load("KMB.RData")

head(counts)
matplot(x=ti,y=counts,type='l')

summary(counts)

counts_2=data.frame(counts[,1:2])

# Aggrego colonne 3 e 4 
counts_2$X3=counts[,3]+counts[,4]

matplot(counts_2,type='l')

library(zoo)

w <- 5  

# Medie mobili
ma <- rollapply(counts_2, 
                width = w, 
                FUN = mean, 
                align = "right", 
                fill = NA)
colnames(ma) <- paste0("MA_", colnames(counts_2))

# Volatilità mobili
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

zeta0=.15
lambda=.5
K=3

dim(data_input)
# 629x6

fit_kmb <- feat_weight_jump(
  Y = data_input,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = 1e-16,
  n_init = 3,
  n_outer = 30,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  mif=1)


res_=data.frame(counts_2[-(1:(w-1)),],cluster=fit_kmb$s)

tapply(res_$X1,res_$cluster,mean)
tapply(res_$X2,res_$cluster,mean)
tapply(res_$X3,res_$cluster,mean)

tapply(res_$X1,res_$cluster,sd)
tapply(res_$X2,res_$cluster,sd)
tapply(res_$X3,res_$cluster,sd)

par(mfrow=c(1,3))
plot(res_$X1,col=res_$cluster,type='p',pch=16)
plot(res_$X2,col=res_$cluster,type='p',pch=16)
plot(res_$X3,col=res_$cluster,type='p',pch=16)

est_feat=fit_kmb$W
colnames(est_feat)=colnames(data_input)
round(est_feat,2)
