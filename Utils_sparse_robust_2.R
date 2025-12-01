library(parallel)
library(Rcpp)
#library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)
library(MixSim)

order_states_condMed <- function(y, s, decreasing = FALSE) {
  
  # This function organizes states by assigning 1 to the state with the smallest (or largest, if decreasing=TRUE)
  # conditional median for vector y, and sequentially numbering each new state.
  
  condMed <- sort(tapply(y, s, median, na.rm = TRUE), decreasing = decreasing)
  
  states_temp <- match(s, names(condMed))
  
  return(states_temp)
}

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

balanced_accuracy <- function(pred, truth) {
  
  conf=caret::confusionMatrix(as.factor(pred), as.factor(truth))
  b=conf$byClass[, "Balanced Accuracy"]
  return(mean(b))
}



# Simulation ---------------------------------------------------------------------

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
  
  
  #MU=seq(mu, -mu, length.out=Ktrue)
  MU=seq(-mu, mu, length.out=Ktrue)
  
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

sim_data_stud_t_overlap <- function(seed = 123,
                                    TT,
                                    P,
                                    Ktrue = 3,
                                    omega = 0.05,   # overlap target
                                    nu = 4,
                                    pers = 0.95,
                                    mif=5) {
  
  
  
  # Generate mixture parameters with target overlap
  set.seed(Ktrue*P)
  simpars <- MixSim::MixSim(K = Ktrue, p = P, BarOmega = omega)
  MU <- simpars$Mu            # K × P matrix
  Sigma <- simpars$S           # P × P × K array
  
  # Order states by decreasing mean on variable mif
  ord <- order(MU[, mif], decreasing = TRUE)
  MU <- MU[ord, , drop = FALSE]
  Sigma <- Sigma[, , ord, drop = FALSE]
  
  # --- Markov chain simulation ---
  set.seed(seed)
  x <- numeric(TT)
  Q <- matrix(rep((1 - pers) / (Ktrue - 1), Ktrue * Ktrue), 
              ncol = Ktrue, byrow = TRUE)
  diag(Q) <- rep(pers, Ktrue)
  init <- rep(1 / Ktrue, Ktrue)
  
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for (i in 2:TT) {
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # --- Continuous variables simulation ---
  SimData <- matrix(0, TT, P)
  set.seed(seed)
  for (k in 1:Ktrue) {
    u <- mvtnorm::rmvt(
      TT,
      sigma = (nu - 2) * Sigma[,,k] / nu,
      df = nu,
      delta = MU[k, ]
    )
    idx <- which(x == k)
    SimData[idx, ] <- u[idx, , drop = FALSE]
  }
  
  SimData <- data.frame(SimData)
  
  
  return(list(
    SimData = SimData,
    mchain  = x,
    TT      = TT,
    P       = P,
    K       = Ktrue,
    pers    = pers,
    seed    = seed,
    omega   = omega,
    Mu      = MU,
    Sigma   = Sigma
  ))
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
                                out_sigma  = 5,
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
    # Y[, irrelevant] <- Y[sample(TT), irrelevant]
    
    # Uniform noise
    Y[, irrelevant] <- matrix(
      runif(length(irrelevant) * TT,
            min = min(Y),
            max = max(Y)),
      # rnorm(length(irrelevant) * TT, 
      #       mean = mean(Y), 
      #       sd = sd(Y)),
      nrow = TT, ncol = length(irrelevant)
    )
    
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
      # Y[idx_out_state, p] <- sample(Y[idx_out_state, p])
      
      # Uniform noise
      Y[idx_out_state, p] =
        runif(length(idx_out_state),
              min = min(Y),
              max = max(Y))
      # rnorm(length(idx_out_state), 
      #       mean = mean(Y[idx_out_state,]), 
      #       sd = sd(Y[idx_out_state,]))
      
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


# Estimation --------------------------------------------------------------



library(Rcpp)
Rcpp::sourceCpp("robJM.cpp")


robust_sparse_jump <- function(Y,
                               zeta0,
                               lambda,
                               K,
                               tol     = 1e-16,
                               n_init  = 5,
                               n_outer = 20,
                               n_inner=10,
                               alpha   = 0.1,
                               verbose = FALSE,
                               knn     = 10,
                               M       = NULL,
                               qt=NULL,
                               c=NULL,
                               mif=NULL,
                               hd=F,
                               n_hd=500,
                               outlier=T,
                               truth=NULL) {
  
  P  <- ncol(Y)
  TT <- nrow(Y)
  
  library(Rcpp)
  
  Rcpp::sourceCpp("robJM.cpp")
  
  start=Sys.time()
  Y=scale(Y)
  
  Gamma <- lambda * (1 - diag(K))
  if(outlier){
    v2 <- v_1(Y,knn=knn, M=M,qt=qt)
  }
  
  run_one <- function(init_id) {
    # 1) initial W, zeta, s
    W        <- matrix(1/P, nrow=K, ncol=P)
    W_old    <- W
    zeta     <- zeta0
    loss_old <- Inf
    #ARI=NA
    BAC=NA
    loss_vec=data.frame(iter=0,loss=loss_old,
                        #ARI=ARI,
                        BAC=BAC,
                        zeta=zeta0)
    
    
    #s = sample(1:K,TT,replace=T)
    
    s  <- initialize_states(Y, K)
    
    for (outer in seq_len(n_outer)) {
      
      for(inner in seq_len(n_inner)){
        if(!outlier){
          v=rep(1,TT)
        }
        else{
          v1 <- v_1(W[s, , drop=FALSE] * Y, knn=knn, c=c, M=M,qt=qt)
          v  <- pmin(v1, v2)
        }
        
        
        if(hd){
          
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
        
        medoids=sel_idx[pam_out$id.med]
        
        # 4) build loss-by-state
        if(!hd){
          loss_by_state <- DW[, medoids, drop=FALSE]  # TT x K
        }
        
        else{
          loss_by_state <- weight_inv_exp_dist(Y=as.matrix(Y * v),s=s,W=W,zeta=zeta,medoids=medoids)
          
        }
        
        s_old <- s
        
        Estep=E_step(loss_by_state,
                     Gamma)
        V=Estep$V
        s=Estep$s
        
        # 7) must have all K states or revert
        if (length(unique(s)) < K) {
          s <- s_old
          break
        }
        
        # 9) update W via WCD + exp
        
        Spk <- WCD(s[sel_idx], as.matrix(Y_search 
                                         * v[sel_idx]
        ), K)
        
        
        wcd <- exp(-Spk / zeta0)
        W   <- wcd / rowSums(wcd)
        
        # Loss computation
        loss=sum(DW[upper.tri(DW)])+
          zeta0*sum(W*log(W))+
          lambda*sum(s[-1]!=s[-TT])
        
        # Loss convergence
        if (!is.null(tol) && abs(loss - loss_old) < tol) break
        loss_old <- loss
        
        
        # 10) W‐convergence as in Kampert 2017
        epsW <- sum(abs(W - W_old))
        if (!is.null(tol) && epsW < tol) break
        W_old <- W
      }
      if(!is.null(truth)){
        ARI=mclust::adjustedRandIndex(truth,s)
        #BAC=balanced_accuracy(s,truth)
      }
      #loss_vec=rbind(loss_vec,c(outer,loss,ARI,zeta))
      loss_vec=rbind(loss_vec,c(outer,loss,BAC,zeta))
      
      # 11) bump zeta
      zeta <- zeta + alpha * zeta0
      
      if (verbose) {
        cat(sprintf("init %2d, outer %2d → loss=%.4e, epsW=%.4e, zeta=%.3f, ARI=%.4f\n",
                    init_id, outer, loss, epsW, zeta, ARI))
        # cat(sprintf("init %2d, outer %2d → loss=%.4e, epsW=%.4e, zeta=%.3f, BAC=%.4f\n",
        #             init_id, outer, loss, epsW, zeta, BAC))
      }
    }
    
    list(W      = W,
         s      = s,
         medoids= medoids,
         v      = v,
         loss   = loss,
         loss_vec=loss_vec[-1,],
         ARI=ARI
         #BAC=BAC
         )
  }
  
  # run n_init times, pick the one with smallest loss
  res_list <- lapply(seq_len(n_init), run_one)
  losses   <- vapply(res_list, `[[`, numeric(1), "loss")
  best_run <- res_list[[ which.min(losses) ]]
  
  best_s   <- best_run$s
  best_loss<- best_run$loss
  loss_vec=best_run$loss_vec
  
  best_W = best_run$W
  best_medoids  <- Y[best_run$medoids,]
  best_v <- best_run$v
  
  # Most important features (mif)
  if(is.null(mif)){
    mif=which.max(apply(best_W,2,sum))
  }
  
  # Re‐order states based on most important feature state-conditional median
  new_best_s <- order_states_condMed(Y[, mif], best_s,decreasing=F)
  
  # tab <- table(best_s, new_best_s)
  # new_order <- apply(tab, 1, which.max)
  # 
  # best_W <- best_W[new_order,]
  tab <- table(factor(best_s, levels = 1:K),
               factor(new_best_s, levels = 1:K))
  
  perm <- apply(tab, 2, which.max)  # per ogni nuova etichetta (colonna), la vecchia (riga)
  best_W <- best_W[perm, , drop = FALSE]
  end=Sys.time()
  
  ret_list=list(W = best_W,
                s = new_best_s,
                medoids = best_medoids,
                v = best_v,
                loss = best_loss,
                loss_vec=loss_vec,
                zeta0 = zeta0,
                lambda = lambda,
                c = c,
                knn=knn,
                M = M,
                elapsed=end-start)
  
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



# Function to analyze sim results -----------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)
#library(Cairo)
analyze_results <- function(res_list, hp, P, K, 
                            label_size = 3.2,
                            ylim_BAC = c(0, 1), 
                            show_legend = F,
                            facet_font_size = 22,
                            x_axis_font_size = 13) {
  
  # --- Filtra per P ---
  idx_P <- hp$P == P
  res_list <- res_list[idx_P]
  hp <- hp[idx_P, ]
  
  # --- Validate results ---
  valid_idx <- sapply(res_list, function(x) {
    is.list(x) && !is.null(x$BAC_s) && is.finite(x$BAC_s)
  })
  
  res_valid <- res_list[valid_idx]
  hp_valid  <- hp[valid_idx, ]
  
  # --- Data frame con BAC e tempi (convertiti in secondi) ---
  BAC_vals <- sapply(res_valid, function(x) x$BAC_s)
  elapsed_vals <- sapply(res_valid, function(x) as.numeric(x$elapsed, units = "secs"))
  
  df <- cbind(hp_valid, BAC_s = BAC_vals, elapsed = elapsed_vals)
  
  # --- Aggregation ---
  summary_df <- df %>%
    group_by(lambda, zeta0) %>%
    summarise(
      median_BAC = median(BAC_s, na.rm = TRUE),
      q025 = quantile(BAC_s, 0.025, na.rm = TRUE),
      q975 = quantile(BAC_s, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(lambda, zeta0), round, 2))
  
  # --- Tempo medio per P (aggregato su lambda e zeta0) ---
  time_summary <- df %>%
    group_by(P) %>%
    summarise(
      mean_elapsed = mean(elapsed, na.rm = TRUE),
      sd_elapsed   = sd(elapsed, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # --- Miglior combinazione ---
  best_indx <- which.max(summary_df$median_BAC)
  best_row  <- summary_df[best_indx, ]
  print(best_row)
  
  # --- Plot per ogni λ ---
  make_plot_lambda <- function(lam, show_y = TRUE) {
    ggplot(dplyr::filter(summary_df, lambda == lam),
           aes(x = zeta0, y = median_BAC, group = 1)) +
      geom_line(color = "darkorange3", linewidth = 1.2) +
      geom_point(color = "darkorange3", size = 2.8) +
      geom_errorbar(aes(ymin = q025, ymax = q975),
                    width = 0.015, color = "darkorange3", alpha = 0.6) +
      coord_cartesian(ylim = ylim_BAC) +
      scale_y_continuous(
        limits = ylim_BAC,
        breaks = seq(ylim_BAC[1], ylim_BAC[2], length.out = 5),
        labels = function(x) sprintf("%.2f", x)
      ) +
      labs(
        x = expression(zeta[0]),
        y = if (show_y) "Median BAC" else NULL,
        title = bquote(lambda == .(lam))
      ) +
      theme_minimal(base_size = 20) +
      theme(
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.spacing = unit(1, "lines"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text.y = if (show_y) element_text(size = 14) else element_blank()
      )
  }
  
  # --- Plot BAC ---
  lambda_vals <- sort(unique(summary_df$lambda))
  plot_list <- lapply(seq_along(lambda_vals), function(i) {
    show_y <- (i == 1)
    make_plot_lambda(lambda_vals[i], show_y = show_y)
  })
  
  BAC_plot <- wrap_plots(plot_list, ncol = length(plot_list)) &
    theme(
      text = element_text(size = 20),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 18),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # --- Compute median weights per (λ, ζ0) ---
  W_median_list <- list()
  param_combos <- hp_valid %>% distinct(lambda, zeta0)
  
  for (i in seq_len(nrow(param_combos))) {
    lam <- param_combos$lambda[i]
    zet <- param_combos$zeta0[i]
    idx <- which(hp_valid$lambda == lam & hp_valid$zeta0 == zet)
    
    Ws <- lapply(res_valid[idx], function(x) x$W)
    arrW <- simplify2array(Ws)
    W_med <- apply(arrW, c(1, 2), median, na.rm = TRUE)
    W_med <- W_med / rowSums(W_med)
    
    name_key <- paste0("lambda=", lam, "_zeta0=", zet)
    W_median_list[[name_key]] <- W_med
  }
  
  # --- Long format per heatmap ---
  W_long <- lapply(names(W_median_list), function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    lambda <- sub("lambda=", "", parts[1])
    zeta0  <- sub("zeta0=", "", parts[2])
    mat <- W_median_list[[nm]]
    df <- as.data.frame(mat)
    df$row <- seq_len(nrow(mat))
    
    df %>%
      tidyr::pivot_longer(-row, names_to = "col", values_to = "weight") %>%
      dplyr::mutate(
        lambda = as.numeric(lambda),
        zeta0  = as.numeric(zeta0),
        state  = row,
        variable = as.integer(gsub("V", "", col))
      ) %>%
      dplyr::select(lambda, zeta0, state, variable, weight)
  }) %>% dplyr::bind_rows()
  
  # --- Heatmap per best parametri ---
  lab_zeta <- labeller(zeta0 = function(x) paste0("\u03B6", "0 = ", x))
  W_sub <- dplyr::filter(W_long, lambda == best_row$lambda)
  
  # Se P > 10, comprimiamo a 10 colonne e usiamo la penultima per colorare la 10ª
  if (P > 10) {
    P_curr <- max(W_sub$variable, na.rm = TRUE)
    
    # Dati delle prime 9 variabili
    W_keep <- dplyr::filter(W_sub, variable <= 9)
    
    # Dati della penultima variabile (per colore della colonna 10)
    W_penult <- dplyr::filter(W_sub, variable == (P_curr - 1)) %>%
      dplyr::mutate(variable = 10L) %>%      # nuova colonna 10
      dplyr::mutate(lbl = "...")             # etichetta "..." nella 10ª
    
    # Etichette numeriche per le prime 9
    W_keep <- W_keep %>%
      dplyr::mutate(lbl = sprintf("%.2f", weight))
    
    # Unisco e fisso l'ordine 1..10
    W_plot <- dplyr::bind_rows(W_keep, W_penult) %>%
      dplyr::mutate(variable = factor(variable, levels = 1:10),
                    state    = factor(state))
    
  } else {
    # Caso P <= 10: tutto come prima (etichette numeriche ovunque)
    W_plot <- W_sub %>%
      dplyr::mutate(
        lbl      = sprintf("%.2f", weight),
        variable = factor(variable, levels = sort(unique(variable))),
        state    = factor(state)
      )
  }
  
  # Range per la scala dei colori
  min_w <- round(min(W_plot$weight, na.rm = TRUE), 2)
  max_w <- round(max(W_plot$weight, na.rm = TRUE), 2)
  
  # Etichette sull'asse x
  if (P > 10) {
    x_labels <- c(1:9, " ")
  } else {
    x_labels <- levels(W_plot$variable)
  }
  
  heatmap_plot <- ggplot(W_plot, aes(x = variable, y = state, fill = weight)) +
    geom_tile() +
    geom_text(aes(label = lbl), color = "black", size = label_size) +
    facet_wrap(~ zeta0, labeller = lab_zeta) +
    scale_x_discrete(labels = x_labels) +  # <-- qui aggiungiamo "..." sull’asse x
    scale_fill_gradientn(
      colours = c("white", "yellow", "orange", "red2"),
      values  = scales::rescale(c(min_w, 0.05, 0.10, max_w)),
      limits  = c(min_w, max_w),
      breaks  = c(min_w, 0.05, 0.10, max_w),
      labels  = scales::number_format(accuracy = 0.02),
      guide   = if (show_legend) "colourbar" else "none"
    ) +
    labs(
      x = "Variable",
      y = "State",
      fill = "Weight"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      strip.text   = element_text(size = facet_font_size),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = x_axis_font_size),
      axis.text.y  = element_text(size = 13)
    )
  
  
  # --- Return outputs ---
  list(
    summary_df   = summary_df,
    best_row     = best_row,
    best_indx    = best_indx,
    BAC_plot     = BAC_plot,
    heatmap_plot = heatmap_plot,
    time_summary = time_summary
  )
}


analyze_v_truth_boxgrid <- function(res_list, hp) {
  
  # --- raccogli tutti i v e truth con lambda e zeta0 ---
  run_df <- lapply(seq_along(res_list), function(i) {
    r <- res_list[[i]]
    tibble(
      v = r$v,
      truth = factor(ifelse(r$truth == 0, "Outlier", "Inlier"),
                     levels = c("Outlier", "Inlier")),
      lambda = hp$lambda[i],
      zeta0  = hp$zeta0[i],
      seed   = hp$seed[i]
    )
  }) %>% bind_rows()
  
  # --- trova best combo per BAC ---
  best_row <- lapply(seq_along(res_list), function(i){
    tibble(
      lambda = hp$lambda[i],
      zeta0  = hp$zeta0[i],
      BAC_s  = ifelse(!is.null(res_list[[i]]$BAC_s), res_list[[i]]$BAC_s, NA_real_)
    )
  }) %>% bind_rows() %>%
    group_by(lambda, zeta0) %>%
    summarise(med_BAC = median(BAC_s, na.rm = TRUE), .groups = "drop") %>%
    filter(med_BAC == max(med_BAC, na.rm = TRUE)) %>% 
    slice(1)
  
  cat("\n>> Best λ, ζ₀ =", best_row$lambda, best_row$zeta0,
      "(median BAC =", round(best_row$med_BAC, 3), ")\n")
  
  # --- boxplot grid v vs truth per lambda × zeta0 con cornici ---
  p <- ggplot(run_df, aes(x = truth, y = v, fill = truth)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6, color = "black", linewidth = 0.3) +
    facet_grid(rows = vars(lambda), cols = vars(zeta0),
               labeller = label_bquote(rows = λ == .(lambda), cols = ζ[0] == .(zeta0))) +
    scale_fill_manual(values = c("Outlier" = "firebrick3", "Inlier" = "steelblue")) +
    labs(x = NULL, y = "v") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 11, face = "bold"),
      strip.text.y = element_text(size = 11, face = "bold", angle = 0),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 13),
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.background = element_rect(fill = "white", colour = "black"),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "black")
    )
  
  list(run_df = run_df, best_row = best_row, grid_plot = p)
}

clean_results <- function(res_list, hp) {
  # Keep only valid elements (no $error field)
  keep <- vapply(res_list, function(x) {
    !(is.list(x) && !is.null(x$error))
  }, logical(1))
  
  # Filter both lists accordingly
  res_list <- res_list[keep]
  hp <- hp[keep, , drop = FALSE]
  
  # Recompute BAC for each valid element
  for (i in seq_along(res_list)) {
    res_list[[i]]$BAC_s <- balanced_accuracy(
      res_list[[i]]$truth,
      res_list[[i]]$s
    )
  }
  
  # Return both cleaned objects
  list(res_list = res_list, hp = hp)
}

