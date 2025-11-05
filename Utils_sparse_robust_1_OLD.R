order_states_condMed=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional median for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMed=sort(tapply(y,s,median,na.rm=T))
  
  states_temp=match(s,names(condMed))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

# OLD -------------------------------------------------------------------------

weight_inv_exp_dist_old <- function(Y,
                                    s, 
                                    W, zeta) {
  TT <- nrow(Y)
  P <- ncol(Y)
  
  sk=apply(Y,2,function(x)IQR(x)/1.35)
  
  # 2. Genera indici delle coppie (i < j)
  pairs <- combn(TT, 2)
  i_idx <- pairs[1, ]
  j_idx <- pairs[2, ]
  n_pairs <- ncol(pairs)
  
  # 3. Estrai le righe corrispondenti
  Yi <- Y[i_idx, , drop = FALSE]
  Yj <- Y[j_idx, , drop = FALSE]
  diff <- abs(Yi - Yj)
  #sk=apply(diff,2,function(x)sum(x)/TT^2)
  
  diff=sweep(diff, 2, sk, FUN = "/")
  #diff=sweep(diff, 2, range_Y, FUN = "/")
  
  
  # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
  W_si <- W[s[i_idx], , drop = FALSE]
  W_sj <- W[s[j_idx], , drop = FALSE]
  # W_si <- W[i_idx, , drop = FALSE]
  # W_sj <- W[j_idx, , drop = FALSE]
  max_w <- pmax(W_si, W_sj)
  
  # 5. Calcola la distanza finale
  weighted_exp <- exp(-diff / zeta) * max_w
  dist_vals <- -zeta * log(rowSums(weighted_exp))
  
  # 6. Ricostruzione della matrice simmetrica
  mat <- matrix(0, TT, TT)
  mat[cbind(i_idx, j_idx)] <- dist_vals
  mat[cbind(j_idx, i_idx)] <- dist_vals
  
  return(mat)
}

# weight_inv_exp_dist <- function(Y,
#                                 s, 
#                                 W, zeta) {
#   
#   indx_num=sapply(Y,is.numeric)
#   
#   TT <- nrow(Y)
#   P <- ncol(Y)
#   Y_orig=Y
#   
#   # continuous vars
#   Y=Y[,indx_num]
#   sk=apply(Y,2,function(x)IQR(x)/1.35)
#   
#   # 2. Genera indici delle coppie (i < j)
#   pairs <- combn(TT, 2)
#   i_idx <- pairs[1, ]
#   j_idx <- pairs[2, ]
#   n_pairs <- ncol(pairs)
#   
#   diff=matrix(0,ncol=P,nrow=n_pairs)
#   
#   # 3. Estrai le righe corrispondenti
#   Yi <- Y[i_idx, , drop = FALSE]
#   Yj <- Y[j_idx, , drop = FALSE]
#   
#   diff_cont <- abs(Yi - Yj)
#   #sk=apply(diff,2,function(x)sum(x)/TT^2)
#   
#   diff_cont=sweep(diff_cont, 2, sk, FUN = "/")
#   diff[,indx_num]=as.matrix(diff_cont[, names(indx_num)[indx_num]])
#   #diff=sweep(diff, 2, range_Y, FUN = "/")
#   
#   # Categorical vars
#   if(sum(indx_num)!=P){
#     
#     # FUNZIONA SOLO SE LE VARIABILI CATEGORIALI SONO ALMENO DUE
#     
#     Y=Y_orig[,!indx_num]
#     
#     # 3. Estrai le righe corrispondenti
#     Yi <- Y[i_idx, , drop = FALSE]
#     Yj <- Y[j_idx, , drop = FALSE]
#     diff_cat <- apply(Yi != Yj, 2, as.numeric)
#     #sk=apply(diff,2,function(x)sum(x)/TT^2)
#     
#     diff[,!indx_num]=as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
#     
#     
#   }
#   
#   
#   # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
#   W_si <- W[s[i_idx], , drop = FALSE]
#   W_sj <- W[s[j_idx], , drop = FALSE]
#   max_w <- pmax(W_si, W_sj)
#   
#   # 5. Calcola la distanza finale
#   weighted_exp <- exp(-diff / zeta) * max_w
#   dist_vals <- -zeta * log(rowSums(weighted_exp))
#   
#   # 6. Ricostruzione della matrice simmetrica
#   mat <- matrix(0, TT, TT)
#   mat[cbind(i_idx, j_idx)] <- dist_vals
#   mat[cbind(j_idx, i_idx)] <- dist_vals
#   
#   return(mat)
# }
# 
# 
# weight_inv_exp_dist_medoids <- function(Y, Ymedoids, s, W, zeta) {
#   indx_num <- sapply(Y, is.numeric)
#   
#   TT <- nrow(Y)
#   K <- nrow(Ymedoids)
#   P <- ncol(Y)
#   
#   Y_orig <- Y
#   Ymedoids_orig <- Ymedoids
#   
#   # Continuous variables
#   Y <- Y[, indx_num, drop = FALSE]
#   Ymedoids_cont <- Ymedoids[, indx_num, drop = FALSE]
#   sk <- apply(Y, 2, function(x) IQR(x) / 1.35)
#   
#   # Pre-allocate output matrix
#   mat <- matrix(0, nrow = TT, ncol = K)
#   
#   # Ciclo su tutti i medoids
#   for (k in 1:K) {
#     diff_cont <- abs(sweep(Y, 2, as.numeric(Ymedoids_cont[k,]), FUN = "-"))
#     
#     diff_cont <- sweep(diff_cont, 2, sk, FUN = "/") # normalizzazione
#     diff <- matrix(0, nrow = TT, ncol = P)
#     diff[, indx_num] <- as.matrix(diff_cont[, names(indx_num)[indx_num]])
#     
#     # Categorical variables
#     if (sum(indx_num) != P) {
#       Y_cat <- Y_orig[, !indx_num, drop = FALSE]
#       medoid_cat <- Ymedoids_orig[k, !indx_num, drop = FALSE]
#       diff_cat <- sweep(as.matrix(Y_cat), 2, as.matrix(medoid_cat), FUN = "!=") * 1
#       diff[, !indx_num] <- as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
#     }
#     
#     # Pesatura con pesi W
#     W_si <- W[s, , drop = FALSE]
#     W_sj <- matrix(W[k, ], nrow = TT, ncol = ncol(W), byrow = TRUE)
#     max_w <- pmax(W_si, W_sj)
#     
#     weighted_exp <- exp(-diff / zeta) * max_w
#     dist_vals <- -zeta * log(rowSums(weighted_exp))
#     
#     mat[, k] <- dist_vals
#   }
#   
#   return(mat)  # matrice T x K
# }

WCD_old=function(s,Y,K){
  #TT <- nrow(Y)
  P <- ncol(Y)
  
  wcd=matrix(0,nrow=K,ncol=P)
  
  sk=apply(Y,2,function(x)IQR(x)/1.35)
  
  for(i in 1:K){
    Ys=Y[s==i,]
    TTk <- nrow(Ys)
    pairs <- combn(TTk, 2)
    i_idx <- pairs[1, ]
    j_idx <- pairs[2, ]
    n_pairs <- ncol(pairs)
    
    Yi <- Ys[i_idx, , drop = FALSE]
    Yj <- Ys[j_idx, , drop = FALSE]
    diff <- abs(Yi - Yj)
    #sk=apply(diff,2,function(x)sum(x)/TTk^2)
    
    diff=sweep(diff, 2, sk, FUN = "/")
    # diff=sweep(diff, 2, range_Y, FUN = "/")
    
    for(p in 1:P){
      mat <- matrix(0, TTk, TTk)
      mat[cbind(i_idx, j_idx)] <- diff[,p]
      mat[cbind(j_idx, i_idx)] <- diff[,p]
      wcd[i,p]=mean(apply(mat,1,median))
    }
    # 
    #wcd[i,]=colSums(diff)/TTk^2
    #wcd[i,]=apply(diff,2,median)/TTk
  }
  return(wcd)
  
}

# WCD=function(s,Y,K){
#   #TT <- nrow(Y)
#   P <- ncol(Y)
#   
#   Y_orig=Y
#   
#   wcd=matrix(0,nrow=K,ncol=P)
#   
#   indx_num=sapply(Y,is.numeric)
#   
#   # 1. Normalizzazione Gower
#   # range_Y <- apply(Y, 2, function(col) {
#   #   r <- max(col) - min(col)
#   #   if (r == 0) 1 else r
#   # })
#   sk=apply(Y[,indx_num],2,function(x)IQR(x)/1.35)
#   
#   for(i in 1:K){
#     Y=Y_orig[,indx_num]
#     Ys=Y[s==i,]
#     TTk <- nrow(Ys)
#     pairs <- combn(TTk, 2)
#     i_idx <- pairs[1, ]
#     j_idx <- pairs[2, ]
#     n_pairs <- ncol(pairs)
#     
#     diff=matrix(0,ncol=P,nrow=n_pairs)
#     
#     Yi <- Ys[i_idx, , drop = FALSE]
#     Yj <- Ys[j_idx, , drop = FALSE]
#     
#     diff_cont <- abs(Yi - Yj)
#     #sk=apply(diff,2,function(x)sum(x)/TT^2)
#     
#     diff_cont=sweep(diff_cont, 2, sk, FUN = "/")
#     diff[,indx_num]=as.matrix(diff_cont[, names(indx_num)[indx_num]])
#     
#     
#     if(sum(indx_num)!=P){
#       
#       # FUNZIONA SOLO SE LE VARIABILI CATEGORIALI SONO ALMENO DUE
#       
#       Y=Y_orig[,!indx_num]
#       
#       # 3. Estrai le righe corrispondenti
#       Yi <- Y[i_idx, , drop = FALSE]
#       Yj <- Y[j_idx, , drop = FALSE]
#       diff_cat <- apply(Yi != Yj, 2, as.numeric)
#       #sk=apply(diff,2,function(x)sum(x)/TT^2)
#       
#       diff[,!indx_num]=as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
#       
#       
#     }
#     
#     
#     for(p in 1:P){
#       mat <- matrix(0, TTk, TTk)
#       mat[cbind(i_idx, j_idx)] <- diff[,p]
#       mat[cbind(j_idx, i_idx)] <- diff[,p]
#       wcd[i,p]=mean(apply(mat,1,median))
#     }
#     # 
#     #wcd[i,]=colSums(diff)/TTk^2
#     #wcd[i,]=apply(diff,2,median)/TTk
#   }
#   return(wcd)
#   
# }


# Sim ---------------------------------------------------------------------

sim_data_stud_t=function(seed=123,
                         TT,
                         P,
                         Pcat,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=4,
                         phi=.8,
                         pers=.95){
  
  
  MU=seq(mu, -mu, length.out=Ktrue)
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    # u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    u = mvtnorm::rmvt(TT, sigma = (nu-2)*Sigma/nu, df = nu, delta = rep(MU[k],P))
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  SimData=data.frame(SimData)
  
  if(!is.null(Pcat)){
    for (j in 1:Pcat) {
      SimData[, j] <- get_cat_t(SimData[, j], x, MU, phi=phi, df = nu)
      SimData[, j]=factor(SimData[, j],levels=1:Ktrue)
    }  
  }
  
  
  
  return(list(
    SimData=SimData,
    mchain=x,
    TT=TT,
    P=P,
    K=Ktrue,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}



# COSA=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F){
#   
#   # Simple version of COSA (no outlier detection)
#   
#   P=ncol(Y)
#   TT=nrow(Y)
#   
#   # best_loss <- NULL
#   # best_s <- NULL
#   # best_W = NULL
#   
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   s=initialize_states(Y,K)
#   # Ymedoids=cluster::pam(Y,k=K)
#   # s=Ymedoids$clustering
#   # Ymedoids=Ymedoids$medoids
#   #s=sample(1:K,TT,replace = T)
#   
#   for (outer in 1:n_outer){
#     
#     ## Clustering
#     #for(inner in 1:n_inner){
#     
#     #Compute distances
#     DW=weight_inv_exp_dist(Y,
#                            s,
#                            W,zeta)
#     medoids=cluster::pam(x=DW,k=K,diss=TRUE)
#     #Ymedoids=Y[medoids$medoids,]
#     s=medoids$clustering
#     
#     # Compute weights
#     
#     Spk=WCD(s,Y,K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     w_loss=sum(W*Spk)  
#     eps_W=mean((W-W_old)^2)
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
#     
#     # print(W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
#     }
#     
#   }
#   
#   
#   return(list(W=W,s=s,medoids=medoids,w_loss=w_loss))
# }

library(DescTools)

lof_star=function(x,knn){
  lof_x=DescTools::LOF(x, knn)
  mean_lof=mean(lof_x)
  sd_lof=sd(lof_x)
  lof_st=(lof_x-mean_lof)/sd_lof
  return(lof_st)
}

# v_1=function(x,knn=10,c=2,M=NULL){
#   
#   lof_st=lof_star(x,knn)
#   
#   if(is.null(M)){
#     M=median(lof_st)+mad(lof_st)
#   }
#   
#   v=rep(1,dim(x)[1])
#   v[lof_st>=c]=0
#   indx=which(M<lof_st&lof_st<c)
#   v[indx]=(1-((lof_st[indx]-M)/(c-M))^2)^2
#   return(v)
# }

# robust_COSA=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F,knn=10,c=2,M=NULL){
#   
#   # Robust version of COSA
# 
#   
#   P=ncol(Y)
#   TT=nrow(Y)
#   
#   # best_loss <- NULL
#   # best_s <- NULL
#   # best_W = NULL
#   
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   s=initialize_states(Y,K)
#   
#   for (outer in 1:n_outer){
#     
#     ## Clustering
#     #for(inner in 1:n_inner){
#     
#     v1=v_1(W[s,]*Y,knn=knn,c=c,M=M)
#     v2=v_1(Y,knn=knn,c=c,M=M)
#     
#     v=apply(cbind(v1,v2),1,min)
#     
#     #Compute distances
#     DW=weight_inv_exp_dist(Y * v,
#                            s,
#                            W,zeta)
#     medoids=cluster::pam(x=DW,k=K,diss=TRUE)
#     #Ymedoids=Y[medoids$medoids,]
#     s=medoids$clustering
#     
#     # Compute weights
#     
#     Spk=WCD(s,Y * v,K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     w_loss=sum(W*Spk)  
#     eps_W=mean((W-W_old)^2)
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
#     
#     # print(W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
#     }
#     
#   }
#   
#   
#   return(list(W=W,s=s,medoids=medoids,w_loss=w_loss,v=v))
# }



JM_COSA=function(Y,zeta0,lambda,K,tol,n_outer=20,alpha=.1,verbose=F,Ts=NULL){
  P=ncol(Y)
  TT=nrow(Y)
  
  # if(is.null(Ts)){
  #   Ts=round(TT/2)
  # }
  
  Gamma <- lambda * (1 - diag(K))
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  # Ottimizzare questo che segue??
  s=initialize_states(Y,K)
  
  loss_old=1e10
  
  for (outer in 1:n_outer){
    
    # subsample=sample(1:TT,Ts,replace = F)
    # Ys=Y[subsample,]
    # ss=s[subsample]
    
    # Questo qui sotto lentissimo se T>1000
    
    DW=weight_inv_exp_dist(Y,
                           s,
                           W,zeta)
    
    # PAM only on subsample
    # DW_1=weight_inv_exp_dist(Ys,
    #                        ss,
    #                        W,zeta)
    
    medoids=cluster::pam(x=DW,k=K,diss=TRUE)
    #Ymedoids=Ys[medoids$medoids,]
    #Ymedoids=Y[medoids$medoids,]
    
    # Questo se lavoro su tutto il dato
    loss_by_state=DW[,medoids$id.med]
    
    # Questo se lavoro su sottocampioni per ottimizzare i tempi
    #loss_by_state=weight_inv_exp_dist_medoids(Y, Ymedoids, s, W, zeta)
    
    V <- loss_by_state
    for (t in (TT-1):1) {
      V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
    }
    s_old=s
    s[1] <- which.min(V[1,])
    for (t in 2:TT) {
      s[t] <- which.min(V[t,] + Gamma[s[t-1],])
    }
    loss <- min(V[1,])
    if (length(unique(s)) < K) {
      s=s_old
      break
    }
    if (!is.null(tol)) {
      epsilon <- loss_old - loss
      if (epsilon < tol) {
        break
      }
    } else if (all(s == s_old)) {
      break
    }
    loss_old <- loss
    #}
    
    # Compute weights
    
    Spk=WCD(s,Y,K)
    wcd=exp(-Spk/zeta0)
    W=wcd/rowSums(wcd)
    
    #}
    
    eps_W=mean((W-W_old)^2)
    if (!is.null(tol)) {
      if (eps_W < tol) {
        break
      }
    }
    
    W_old=W
    zeta=zeta+alpha*zeta0
    
    #print(W)
    # print(zeta)
    # print(Spk)
    # print(zeta0)
    # print(range(DW))
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
      #cat(sprintf('Out iteration %d (# inn iterations %d): %.6e\n', outer, inner, eps_W))
    }
    
  }
  return(list(W=W,s=s,medoids=medoids))
}

invert_rel <- function(rel_, P) {
  # rel_: list of length K, each rel_[[k]] is a vector of indices in 1:P
  # P    : total number of elements
  #
  # returns a list inv of length P, where
  #   inv[[i]] = all k such that i %in% rel_[[k]]
  
  # initialize with empty integer vectors
  inv <- vector("list", P)
  for(i in seq_len(P)) inv[[i]] <- integer(0)
  
  # for each group k, append k to every member i in rel_[[k]]
  for(k in seq_along(rel_)) {
    members <- rel_[[k]]
    for(i in members) {
      inv[[i]] <- c(inv[[i]], k)
    }
  }
  
  inv
}

simulate_sparse_hmm <- function(Y,
                                rel_,
                                true_stat,
                                perc_out   = 0.02,
                                out_sigma  = 100,
                                seed       = NULL) {
  # Y         : T x P data matrix or data.frame
  # rel_      : list of length K; rel_[[k]] is a vector of features in state k
  # true_stat : integer vector length T with values in 1:K
  # perc_out  : fraction of rows to turn into outliers
  # out_sigma : sd of Gaussian noise added for outliers
  # seed      : optional RNG seed for reproducibility
  
  if(!is.null(seed)) set.seed(seed)
  Y <- as.matrix(Y)
  TT  <- nrow(Y)
  P  <- ncol(Y)
  
  # 1) invert the rel_ list: for each feature p, which states mention p?
  inv_rel <- invert_rel(rel_, P)
  
  # 2) Irrelevant features = those never mentioned in rel_
  irrelevant <- which(vapply(inv_rel, length, integer(1)) == 0)
  if(length(irrelevant) > 0) {
    # permute their rows globally
    Y[, irrelevant] <- Y[sample(TT), irrelevant]
  }
  
  # 3) Relevant features = those that appear in at least one state
  relevant <- which(vapply(inv_rel, length, integer(1)) > 0)
  for(p in relevant) {
    relevant_states <- inv_rel[[p]]
    # indices of rows belonging to any of those states
    idx_in_state <- which(true_stat %in% relevant_states)
    # rows not in those states:
    idx_out_state <- setdiff(seq_len(TT), idx_in_state)
    if(length(idx_out_state) > 1) {
      # permute only those rows of column p
      Y[idx_out_state, p] <- sample(Y[idx_out_state, p])
    }
  }
  
  # 4) Introduce outliers
  N_out <- ceiling(TT * perc_out)
  t_out <- sort(sample(seq_len(TT), N_out))
  # add Gaussian noise to all features of those rows
  Y[t_out, ] <- Y[t_out, ] + matrix(rnorm(N_out * P, 0, out_sigma),
                                    nrow = N_out, ncol = P)
  
  # 5) update truth: set outlier rows to 0
  new_truth <- true_stat
  new_truth[t_out] <- 0L
  
  W_truth <- matrix(FALSE, nrow = K, ncol = P)
  
  for (k in seq_len(K)) {
    W_truth[k, rel_[[k]]] <- TRUE
  }
  
  list(
    Y          = Y,
    truth      = new_truth,
    out_indices = t_out,
    W_truth = W_truth
  )
}

library(Rcpp)

Rcpp::sourceCpp("robJM.cpp")


robust_sparse_jump <- function(Y,
                               zeta0,
                               lambda,
                               K,
                               tol     = 1e-16,
                               n_init  = 5,
                               n_outer = 20,
                               alpha   = 0.1,
                               verbose = FALSE,
                               knn     = 10,
                               c       = 10,
                               M       = NULL,
                               mif=NULL,
                               hd=F,
                               n_hd=NULL,
                               outlier=T) {
  
  P  <- ncol(Y)
  TT <- nrow(Y)
  
  library(Rcpp)
  
  Rcpp::sourceCpp("robJM.cpp")
  
  if(!is.matrix(Y)){
    Y=as.matrix(
      data.frame(
        lapply(Y, function(col) {
          if (is.factor(col)||is.character(col)) as.integer(as.character(col))
          else              as.numeric(col)
        })
      )
    )
  }
  
  # Check if which features are constant
  constant_features <- apply(Y, 2, function(col) length(unique(col)) == 1)
  
  # Transform into noise (?)
  Y[,constant_features]=rnorm(TT,sd=10)
  
  # Standardize
  Y=scale(Y)
  
  Gamma <- lambda * (1 - diag(K))
  
  run_one <- function(init_id) {
    # 1) initial W, zeta, s
    W        <- matrix(1/P, nrow=K, ncol=P)
    W_old    <- W
    zeta     <- zeta0
    loss_old <- Inf
    #s = sample(1:K,TT,replace=T)
    s        <- initialize_states(Y, K)
    
    for (outer in seq_len(n_outer)) {
      
      if(!outlier){
        v=rep(1,TT)
      }
      else{
        v1 <- v_1(W[s, , drop=FALSE] * Y, knn=knn, c=c, M=M)
        v2 <- v_1(Y,                     knn=knn, c=c, M=M)
        v  <- pmin(v1, v2)
      }
      
      
      if(hd){
        if(is.null(n_hd)){
          n_hd=500
        }
        
        sel_idx=sort(sample(1:TT,n_hd,replace=F))
        Y_search=Y[sel_idx,]
        
        Y_search=as.matrix(Y_search*v[sel_idx])
        
      }
      
      else{
        Y_search=as.matrix(Y * v)
        sel_idx=1:TT
      }
      
      DW      <- weight_inv_exp_dist(Y_search, s[sel_idx], W, zeta)
      pam_out <- cluster::pam(DW, k=K, diss=TRUE)
      #medoids <- pam_out$id.med
      medoids=sel_idx[pam_out$id.med]
      
      # 4) build loss-by-state
      if(!hd){
        loss_by_state <- DW[, medoids, drop=FALSE]  # TT x K
      }
      
      else{
        loss_by_state <- weight_inv_exp_dist(Y=as.matrix(Y * v),
                                             s=s,W=W,zeta=zeta,
                                             medoids=medoids)
      }
      
      
      
      # 5) DP forward: V[t,j] = loss[t,j] + min_i( V[t+1,i] + Gamma[i,j] )
      # V <- loss_by_state
      # 
      # # Input: V is TT x K, Gamma is KxK, s and s_old have length TT
      # for (t in (TT-1):1) {
      #   for (j in seq_len(K)) {
      #     # look at row t+1 of V plus column j of Gamma:
      #     V[t, j] <- loss_by_state[t, j] +
      #       min( V[t+1, ] + Gamma[, j] )
      #   }
      # }
      # 
      # # 6) backtrack to get s
      # s_old <- s
      # # first time‐point
      # s[1] <- which.min(V[1, ])
      # loss  <- V[1, s[1]]
      # # subsequent
      # for (t in 2:TT) {
      #   prev <- s[t-1]
      #   # pick state j minimizing V[t,j] + penalty from prev
      #   scores <- V[t, ] + Gamma[prev, ]
      #   s[t] <- which.min(scores)
      # }
      
      # Cpp alternative
      s=E_step(loss_by_state,
               Gamma)
      loss  <- loss_by_state[1, s[1]]
      
      # 7) must have all K states or revert
      if (length(unique(s)) < K) {
        s <- s_old
        break
      }
      
      # 8) loss‐convergence
      if (!is.null(tol) && (loss_old - loss) < tol) break
      loss_old <- loss
      
      # 9) update W via WCD + exp
      if(hd){
        Spk <- WCD(s[sel_idx], as.matrix(Y_search * v[sel_idx]), K)
      }
      else{
        Spk <- WCD(s, as.matrix(Y * v), K)
      }
      
      wcd <- exp(-Spk / zeta0)
      W   <- wcd / rowSums(wcd)
      
      # 10) W‐convergence
      epsW <- mean((W - W_old)^2)
      if (!is.null(tol) && epsW < tol) break
      W_old <- W
      
      # 11) bump zeta
      zeta <- zeta + alpha * zeta0
      
      if (verbose) {
        cat(sprintf("init %2d, outer %2d → loss=%.4e, epsW=%.4e, zeta=%.3f\n",
                    init_id, outer, loss, epsW, zeta))
      }
    }
    
    list(W      = W,
         s      = s,
         medoids= medoids,
         v      = v,
         loss   = loss)
  }
  
  # run n_init times, pick the one with smallest loss
  res_list <- lapply(seq_len(n_init), run_one)
  losses   <- vapply(res_list, `[[`, numeric(1), "loss")
  best_run <- res_list[[ which.min(losses) ]]
  
  best_s   <- best_run$s
  best_loss<- best_run$loss
  best_W = best_run$W
  best_medoids  <- Y[best_run$medoids,]
  best_v <- best_run$v
  
  # Most important features (mif)
  if(is.null(mif)){
    mif=which.max(apply(best_W,2,sum))
  }
  
  # Re‐order states based on most important feature state-conditional median
  new_best_s <- order_states_condMed(Y[, mif], best_s)
  
  tab <- table(best_s, new_best_s)
  new_order <- apply(tab, 1, which.max)
  
  best_W <- best_W[new_order,]
  
  ret_list=list(W = best_W,
                s = new_best_s,
                medoids = best_medoids,
                v = best_v,
                loss = best_loss,
                zeta0 = zeta0,
                lambda = lambda,
                c = c,
                knn=knn,
                M = M)
  
  return(ret_list)
  
}

cv_robust_sparse_jump <- function(
    Y,
    true_states,
    K_grid=NULL,
    zeta0_grid=NULL,
    lambda_grid=NULL,
    n_folds = 5,
    parallel=F,
    n_cores=NULL,
    cv_method="blocked-cv",
    knn=10,
    c_grid=NULL,
    M=NULL,
    n_init =5,
    n_outer=5,
    tol=1e-8,
    hd=F,
    n_hd=NULL,
    outlier=T
) {
  
  # cv_sparse_jump: Cross-validate Sparse Jump Model parameters (K and lambda)
  
  # Arguments:
  #   Y           - data matrix (N × P)
  #   true_states - vector of true states for ARI computation
  #   K_grid      - vector of candidate numbers of states
  #   zeta0       - sparsity hyperparameter (if NULL it is set to 0.2)
  #   lambda_grid - vector of candidate lambdas
  #   n_folds     - number of folds for cross-validation (default: 5)
  #   parallel    - logical; TRUE for parallel execution (default: FALSE)
  #   n_cores     - number of cores to use for parallel execution (default:  NULL)
  #   cv_method   - method for cross-validation: "blocked-cv" or "forward-chain"
  #   knn         - number of nearest neighbors for LOF (default: 10)
  #   c_grid           - lower threshold for LOF (default: NULL)
  #   M           - upper threshold for LOF (default: NULL, uses median + mad)
  
  
  # Value:
  #   A data.frame with one row per (K, zeta0, lambda) combination, containing:
  #     K      - number of states tested
  #     zeta0  - sparsity hyperparameter value
  #     lambda - jump penalty value
  #     ARI    - mean Adjusted Rand Index across folds
  
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  
  if(is.null(lambda_grid)) {
    lambda_grid <- seq(0,1,.1)  # Default range for lambda
  }
  
  if(is.null(zeta0_grid)) {
    zeta0_grid <- 0.2  # Default sparsity hyperparameter
  }
  
  if(is.null(c_grid)) {
    c_grid <- c(10,15)
  }
  
  # Libreria per ARI
  library(mclust)
  
  TT <- nrow(Y)
  P <- ncol(Y)
  
  # Suddivido gli N campioni in n_folds blocchi contigui
  if(cv_method=="blocked-cv"){
    fold_indices <- split(
      1:TT,
      rep(1:n_folds, each = ceiling(TT / n_folds), length.out = TT)
    )
    
  }
  else if(cv_method=="forward-chain"){
    fold_indices <- lapply(seq_len(n_folds), function(k) {
      idx_end <- TT - (k - 1)
      1:(idx_end-1)
    })
    names(fold_indices) <- as.character(seq_len(n_folds))
  }
  else{
    stop("cv_method must be either 'blocked-cv' or 'forward-chain'")
  }
  
  # Funzione che, per una tripla (K, kappa, lambda) e un fold (train_idx, val_idx),
  # calcola l’ARI sui punti di validazione
  fold_ari <- function(Y,K, zeta0, lambda,c, train_idx, val_idx,true_states) {
    # 1) Fit del modello sparse_jump su soli dati di TRAIN
    res <- robust_sparse_jump(Y=as.matrix(Y[train_idx, , drop = FALSE]),
                              zeta0=zeta0,
                              lambda=lambda,
                              K=K,
                              tol        = tol,
                              n_init     = n_init,
                              n_outer    = n_outer,
                              alpha      = 0.1,
                              verbose    = F,
                              knn        = knn,
                              c          = c,
                              M          = M,
                              hd=hd,
                              n_hd=n_hd,
                              outlier=outlier)
    states_train <- res$s
    feat_idx     <- which(colSums(res$W) > 0.025)
    
    # Se non vengono selezionate feature, restituisco ARI = 0
    if (length(feat_idx) == 0) {
      return(0)
    }
    
    # 2) Extract medoids
    medoids=res$medoids[,feat_idx]
    
    # 3) Assegno ciascun punto in VAL a uno stato: 
    #    lo stato k che minimizza la distanza gower su feature selezionate
    Y_val_feats     <- as.matrix(Y[val_idx, feat_idx, drop = FALSE])
    pred_val_states <- integer(length(val_idx))
    
    dists<- gower_dist(Y_val_feats,medoids)
    
    pred_val_states<- apply(dists,1,which.min)
    
    
    # 4) Calcolo ARI tra etichette vere e quelle predette sul blocco VAL
    return(
      adjustedRandIndex(true_states[val_idx], pred_val_states)
    )
  }
  
  # Costruisco la griglia di tutte le combinazioni di (K, kappa, lambda)
  grid <- expand.grid(
    K      = K_grid,
    # For kappa we take a representative value, we will select kappa later based on GAP stat
    zeta0  = zeta0_grid,
    lambda = lambda_grid,
    c =c_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Data.frame in cui raccogliere (K, kappa, lambda, media_ARI)
  results <- data.frame(
    K      = integer(0),
    zeta0  = integer(0),
    lambda = numeric(0),
    c=numeric(0),
    ARI    = numeric(0)
  )
  
  # Loop su ciascuna riga della griglia
  if (!parallel) {
    # applichiamo una funzione su ogni riga di 'grid'
    results_list <- lapply(seq_len(nrow(grid)), function(row) {
      K_val     <- grid$K[row]
      zeta0_val <- grid$zeta0[row]
      lambda_val<- grid$lambda[row]
      c_val     <- grid$c[row]
      
      # calcolo ARI su ciascun fold
      ari_vals <- numeric(n_folds)
      for (f in seq_len(n_folds)) {
        val_idx   <- fold_indices[[f]]
        train_idx <- setdiff(seq_len(TT), val_idx)
        ari_vals[f] <- fold_ari(Y,K_val, zeta0_val, lambda_val, c_val,
                                train_idx, val_idx,true_states)
      }
      mean_ari <- mean(ari_vals)
      
      # restituisco un data.frame di una sola riga
      data.frame(
        K      = K_val,
        zeta0  = zeta0_val,
        lambda = lambda_val,
        c= c_val,
        ARI    = mean_ari,
        stringsAsFactors = FALSE
      )
    })
    
    # combino tutti i data.frame in un unico data.frame
    results <- do.call(rbind, results_list)
  }
  
  
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    
    results_list <- parallel::mclapply(
      seq_len(nrow(grid)),
      function(row) {
        # === re‐load any packages / source code you need ===
        library(mclust)             # for adjustedRandIndex()
        library(cluster)            # for pam() or any cluster functions
        # if gower_dist() comes from a package, load it here;
        # otherwise, ensure that you have defined gower_dist() in the global env
        # and that it’s visible to the children.
        
        # If your robust_sparse_jump() relies on an Rcpp sourceCpp,
        # you’ll also need to re‐source it inside each worker:
        Rcpp::sourceCpp("robJM.cpp")
        
        # now pull out the grid parameters
        K_val     <- as.integer(grid$K[row])
        zeta0_val <- grid$zeta0[row]
        lambda_val<- grid$lambda[row]
        c_val     <- grid$c[row]
        
        # a quick local copy of true_states so it’s in scope
        ts <- true_states
        
        # compute ARI over folds just like before
        ari_vals <- vapply(seq_len(n_folds), function(f) {
          val_idx   <- fold_indices[[f]]
          train_idx <- setdiff(seq_len(TT), val_idx)
          fold_ari(
            Y, K_val, zeta0_val, lambda_val, c_val,
            train_idx, val_idx, ts
          )
        }, numeric(1))
        
        data.frame(
          K      = K_val,
          zeta0  = zeta0_val,
          lambda = lambda_val,
          c      = c_val,
          ARI    = mean(ari_vals),
          stringsAsFactors = FALSE
        )
      },
      mc.cores = n_cores
    )
    
    results <- do.call(rbind, results_list)
  }
  
  
  
  return(results)
}

permute_gap <- function(Y) {
  if (!is.matrix(Y)) {
    stop("Input Y must be a matrix.")
  }
  TT <- nrow(Y)
  P <- ncol(Y)
  
  Y_perm <- apply(Y, 2, function(col) {
    sample(col, size = TT, replace = FALSE)
  })
  
  if (P == 1) {
    Y_perm <- matrix(Y_perm, nrow = T, ncol = 1)
  }
  
  rownames(Y_perm) <- rownames(Y)
  colnames(Y_perm) <- colnames(Y)
  
  return(Y_perm)
}

gap_robust_sparse_jump=function(
    Y,
    K_grid=NULL,
    zeta0_grid=NULL,
    lambda=0,
    B=10,
    parallel=F,
    n_cores=NULL,
    knn=10,
    c=10,
    M=NULL
){
  
  # gap_robust_sparse_jump: Compute the Gap Statistic for Robust Sparse Jump Model
  
  # Arguments:
  #   Y           - data matrix (N × P)
  #   K_grid      - vector of candidate numbers of states (default: seq(2, 4, by = 1))
  #   zeta0_grid  - vector of candidate kappa values (default: seq(0.05, 0.4, by = 0.05))
  #   lambda      - jump penalty value (default: 0.2)
  #   B           - number of bootstrap samples (default: 10)
  #   parallel    - logical; TRUE for parallel execution (default: FALSE)
  #   n_cores     - number of cores to use for parallel execution (default: NULL)
  #   knn         - number of nearest neighbors for LOF (default: 10)
  #   c           - lower threshold for LOF (default: 5)
  #   M           - upper threshold for LOF (default: NULL, uses median + mad)
  
  # Value:
  #   A list containing:
  #     gap_stats - data.frame with Gap Statistic results
  #     plot_res  - ggplot object of Gap Statistic vs zeta0
  #     meta_df   - data.frame with detailed results for each (K, zeta0, lambda, b)
  
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  
  if(is.null(zeta0_grid)) {
    zeta0_grid <- seq(0.05,.4,.05)  # Default range for lambda
  }
  
  if(is.null(lambda)){
    lambda=0.2
  }
  
  # Libreria per ARI
  library(mclust)
  
  N <- nrow(Y)
  P <- ncol(Y)
  
  grid <- expand.grid(zeta0 = zeta0_grid, 
                      lambda = lambda, 
                      K = K_grid, 
                      b = 0:B)
  
  gap_one_run=function(zeta0,lambda,K,b){
    if (b == 0) {
      Y_input <- Y
      permuted <- FALSE
    } else {
      # Permute features for kappa
      Y_input <- permute_gap(Y)
      permuted <- TRUE
    }
    # Fit the model
    
    fit=robust_sparse_jump(Y,
                           zeta0=zeta0,
                           lambda=lambda,
                           K=K,
                           tol     = NULL,
                           n_init  = 5,
                           n_outer = 20,
                           alpha   = 0.1,
                           verbose = FALSE,
                           knn     = knn,
                           c       = c,
                           M       = M)
    
    return(list(loss=fit$loss,
                K=K,
                permuted=permuted,
                zeta0=zeta0,
                lambda=lambda
    ))
    
  }
  
  if(parallel){
    if(is.null(n_cores)){
      n_cores <- parallel::detectCores() - 1
    }
    results <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
      params <- grid[i, ]
      gap_one_run(
        zeta0 = params$zeta0,
        lambda= params$lambda,
        K     = params$K,
        b     = params$b
      )
    }, mc.cores = mc_cores)
  }
  else{
    results <- lapply(seq_len(nrow(grid)), function(i) {
      params <- grid[i, ]
      gap_one_run(
        zeta0 = params$zeta0,
        lambda= params$lambda,
        K     = params$K,
        b     = params$b
      )
    })
  }
  
  meta_df <- do.call(rbind.data.frame, c(results, make.row.names = FALSE))
  
  library(dplyr)
  gap_stats <- meta_df %>%
    group_by(K, zeta0) %>%
    summarise(
      log_O = log(loss[!permuted]),
      log_O_star_mean = mean(log(loss[permuted])),
      se_log_O_star=sd(log(loss[permuted])),
      GAP =  log_O - log_O_star_mean,
      .groups = 'drop'
    )
  
  library(ggplot2)
  
  plot_res=ggplot(gap_stats, aes(x = zeta0, y = GAP, color = factor(K))) +
    geom_line() +
    geom_point() +
    scale_color_discrete(name = "Number of clusters\n(K)") +
    labs(
      x = expression(kappa),
      y = "Gap Statistic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  return(list(gap_stats=gap_stats,
              plot_res=plot_res,
              meta_df=meta_df))
  
  
}


# robust_JM_COSA=function(Y,zeta0,lambda,K,tol,n_outer=20,alpha=.1,
#                         verbose=F,knn=10,c=2,M=NULL){
#   library(Rcpp)
#   Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
#   Rcpp::sourceCpp("wcd.cpp")
#   P=ncol(Y)
#   TT=nrow(Y)
# 
#   Gamma <- lambda * (1 - diag(K))
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   # Multiple initialization, keep best one (lower loss)
#   s=initialize_states(Y,K)
#   
#   loss_old=1e10
#   
#   for (outer in 1:n_outer){
#     
#     # subsample=sample(1:TT,Ts,replace = F)
#     # Ys=Y[subsample,]
#     # ss=s[subsample]
#     
#     v1=v_1(W[s,]*Y,knn=knn,c=c,M=M)
#     v2=v_1(Y,knn=knn,c=c,M=M)
#     
#     v=apply(cbind(v1,v2),1,min)
#     
#     #Compute distances
#     DW=weight_inv_exp_dist(as.matrix(Y * v),
#                            s,
#                            W,zeta)
#     medoids=cluster::pam(x=DW,k=K,diss=TRUE)
#     
#     # Questo se lavoro su tutto il dato
#     loss_by_state=DW[,medoids$id.med]
#     
#     # Questo se lavoro su sottocampioni per ottimizzare i tempi
#     #loss_by_state=weight_inv_exp_dist_medoids(Y, Ymedoids, s, W, zeta)
#     
#     V <- loss_by_state
#     for (t in (TT-1):1) {
#       V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
#     }
#     s_old=s
#     s[1] <- which.min(V[1,])
#     for (t in 2:TT) {
#       s[t] <- which.min(V[t,] + Gamma[s[t-1],])
#     }
#     loss <- min(V[1,])
#     if (length(unique(s)) < K) {
#       s=s_old
#       break
#     }
#     
#     epsilon <- loss_old - loss
#     if (!is.null(tol)) {
#       if (epsilon < tol) {
#         break
#       }
#     } 
#     # else if (all(s == s_old)) {
#     #   break
#     # }
#     loss_old <- loss
#     #}
#     
#     # Compute weights
#     
#     Spk=WCD(s,as.matrix(Y * v),K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     #}
#     
#     eps_W=mean((W-W_old)^2)
#     
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
# 
#     # print(W)
#     # print(epsilon)
#     # print(eps_W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Iteration %d: %.6e\n', outer, loss))
#       #cat(sprintf('Out iteration %d (# inn iterations %d): %.6e\n', outer, inner, eps_W))
#     }
#     
#   }
#   return(list(W=W,s=s,medoids=medoids,v=v,loss=loss))
# }

# RJM_COSA_gap=function(Y,
#                       zeta_grid=seq(0.1,.7,.1),
#                       lambda_grid=seq(0,1,.1),
#                       K_grid=2:6,
#                       tol=NULL,n_outer=20,alpha=.1,verbose=F,n_cores=NULL,
#                       B=10, knn=10,c=2,M=NULL){
#   
#   # B is the number of permutations
#   
#   grid <- expand.grid(zeta0 = zeta_grid, 
#                       lambda = lambda_grid, 
#                       K = K_grid, b = 0:B)
#   
#   library(foreach)
#   library(doParallel)
#   
#   if(is.null(n_cores)){
#     n_cores <- parallel::detectCores() - 1
#   } 
#   
#   # Set up cluster
#   cl <- makeCluster(n_cores)
#   registerDoParallel(cl)
#   
#   results_list <- foreach(i = 1:nrow(grid), .combine = 'list',
#                           .packages = c("cluster","Rcpp","DescTools"),
#                           .multicombine = TRUE,
#                           .export = c("Y", "robust_JM_COSA", 
#                                       #"WCD", "weight_inv_exp_dist",
#                                       "initialize_states",
#                                       "v_1","lof_star",
#                                       "grid", "tol", "n_outer", "alpha",
#                                       "knn","c","M")) %dopar% {
#                                         K_val <- grid$K[i]
#                                         zeta_val <- grid$zeta0[i]
#                                         lambda_val=grid$lambda[i]
#                                         b <- grid$b[i]
#                                         
#                                         set.seed(b + 1000 * i)
#                                         
#                                         if (b == 0) {
#                                           Y_input <- Y
#                                           permuted <- FALSE
#                                         } else {
#                                           # Permute features for zeta0
#                                           Y_input <- apply(Y, 2, sample)
#                                           # Permute rows for lambda
#                                           Y_input <- Y_input[ sample(nrow(Y_input),
#                                                                      size = nrow(Y_input),
#                                                                      replace = FALSE), ]
#                                           permuted <- TRUE
#                                         }
#                                         
#                                         res <- robust_JM_COSA(Y_input, zeta0 = zeta_val, 
#                                                               lambda = lambda_val, K = K_val, tol = tol,
#                                                               n_outer = n_outer, alpha = alpha, verbose = FALSE,
#                                                               knn=knn,c=c,M=M)
#                                         
#                                         list(
#                                           meta = data.frame(zeta0 = zeta_val, lambda=lambda_val,K = K_val, 
#                                                             loss = res$loss, permuted = permuted),
#                                           cosa = if (!permuted) list(zeta0 = zeta_val, lambda=lambda_val,K = K_val, 
#                                                                      W = res$W, s = res$s, 
#                                                                      medoids = res$medoids$medoids,
#                                                                      v=res$v) else NULL
#                                         )
#                                       }
#   
#   
#   stopCluster(cl)
#   
#   # Flatten results
#   meta_df <- do.call(rbind, lapply(results_list, `[[`, "meta"))
#   cosa_results <- Filter(Negate(is.null), lapply(results_list, `[[`, "cosa"))
#   
#   
#   # Compute GAP
#   library(dplyr)
#   gap_stats <- meta_df %>%
#     group_by(K, zeta0) %>%
#     summarise(
#       log_O = log(loss[!permuted]),
#       log_O_star_mean = mean(log(loss[permuted])),
#       se_log_O_star=sd(log(loss[permuted])),
#       GAP = log_O_star_mean - log_O,
#       .groups = 'drop'
#     )
#   
#   return(list(
#     gap_stats = gap_stats,
#     cosa_results = cosa_results
#   ))
#   
# }

plot_W=function(W){
  library(reshape)
  df <- as.data.frame(W)
  df$Cluster <- factor(paste0("Cluster_", 1:nrow(df)))
  
  # Riorganizziamo in formato lungo
  df_long <- melt(df, id.vars = "Cluster", variable.name = "Feature", value.name = "Weight")
  
  # Converti Feature in fattore per ordinare le colonne
  #df_long$Feature <- factor(df_long$Feature, levels = paste0("V", 1:ncol(fit$W)))
  
  # Bar plot
  library(ggplot2)
  p=ggplot2::ggplot(df_long, aes(x = variable, y = value, fill = Cluster)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Feature Weights by Cluster",
         x = "Feature",
         y = "Weight") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}


