library(scales)
source("Utils_sparse_robust_2.R")

TT=1000

P=10

K_grid=c(3,4)
zeta0=c(0,.05,.1,.2)
zeta0[1]=.01
alpha=.1
lambda=seq(0,1,.25)

rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

K=3

simDat <- sim_data_stud_t(seed=123,
                          TT=TT,
                          P=P,
                          Pcat=NULL,
                          Ktrue=K,
                          mu=3,
                          rho=0,
                          nu=10,
                          pers=.95)

simDat_sparse <- simulate_sparse_hmm(
  Y = simDat$SimData,
  rel_ = rel_,
  true_stat = simDat$mchain,
  perc_out = 0,
  out_sigma = 0,
  seed = 123
)

truth <- simDat_sparse$truth

# Y = as.matrix(simDat_sparse$Y)
# K_grid        = K_grid
# zeta0_grid    = zeta0
# lambda_grid   = lambda
# B             = 10
# parallel      = TRUE
# n_cores       = 9
# tol           = 1e-4
# n_init        = 1
# n_outer       = 50
# n_inner       = 5
# alpha         = 0.1
# verbose       = FALSE
# knn           = 10
# M             = NULL
# qt            = NULL
# c             = NULL
# mif           = NULL
# hd            = TRUE
# n_hd          = 500
# outlier       = FALSE
# truth         = NULL
# row_block     = NULL
# seed          = NULL
# 

gap_temp=gap_robust_sparse_jump(
    Y=as.matrix(simDat_sparse$Y),
    K_grid        = K_grid,
    zeta0_grid    = zeta0,
    lambda_grid   = lambda,
    B             = 10,
    parallel      = TRUE,
    n_cores       = 9,
    tol           = 1e-4,
    n_init        = 1,
    n_outer       = 50,
    n_inner       = 5,
    alpha         = 0.1,
    verbose       = FALSE,
    knn           = 10,
    M             = NULL,
    qt            = NULL,
    c             = NULL,
    mif           = NULL,
    hd            = TRUE,
    n_hd          = 500,
    outlier       = FALSE,
    truth         = NULL,
    row_block     = NULL,
    seed          = NULL
) 

library(dplyr)

gaps <- gap_temp$gap_stats %>%
  mutate(
    K      = as.integer(K),
    zeta0  = as.numeric(zeta0),
    lambda = as.numeric(lambda),
    GAP    = as.numeric(GAP),
    s      = as.numeric(s)   # 1-SE uses s at the maximizer, per Witten & Tibshirani
  )

# Witten & Tibshirani rule (sparse k-means style):
# For a 1D tuning curve: pick the MOST PARSIMONIOUS tuning value whose GAP is within 1 SE of the maximum.
# 'prefer' controls "parsimonious": "smallest" or "largest" along the tuning axis.
wt_pick_1se <- function(df, tune_col, prefer = c("smallest", "largest")) {
  prefer <- match.arg(prefer)
  
  df <- df %>% filter(!is.na(GAP), !is.na(s), !is.na(.data[[tune_col]]))
  if (nrow(df) == 0) return(NULL)
  
  i_max <- which.max(df$GAP)
  gap_max <- df$GAP[i_max]
  se_max  <- df$s[i_max]
  thr <- gap_max - se_max
  
  cand <- df %>% filter(GAP >= thr)
  if (nrow(cand) == 0) cand <- df[i_max, , drop = FALSE]
  
  if (prefer == "smallest") {
    cand <- cand %>% arrange(.data[[tune_col]], desc(GAP))
  } else {
    cand <- cand %>% arrange(desc(.data[[tune_col]]), desc(GAP))
  }
  
  cand[1, , drop = FALSE]
}

# ---- Procedure (WT-style) adapted to (K, lambda, zeta0) ----
# Step 1: for each (K, lambda), pick zeta0* within 1-SE of max over zeta0 (parsimonious: smallest zeta0)
zeta_star <- gaps %>%
  group_by(K, lambda) %>%
  group_modify(~{
    out <- wt_pick_1se(.x, "zeta0", prefer = "largest")
    if (is.null(out)) return(tibble())
    tibble(
      K = out$K,
      lambda = out$lambda,
      zeta0_star = out$zeta0,
      GAP_prof   = out$GAP,
      s_prof     = out$s
    )
  }) %>%
  ungroup()

# Step 2: for each K, pick lambda* within 1-SE of max over lambda on profiled GAP (parsimonious: largest lambda)
lambda_star <- zeta_star %>%
  rename(GAP = GAP_prof, s = s_prof) %>%
  group_by(K) %>%
  group_modify(~{
    out <- wt_pick_1se(.x, "lambda", prefer = "largest")
    if (is.null(out)) return(tibble())
    tibble(
      K = out$K,
      lambda_star = out$lambda,
      zeta0_star  = out$zeta0_star,
      GAP_star    = out$GAP,
      s_star      = out$s
    )
  }) %>%
  ungroup()

# Step 3: choose K that maximizes GAP_star (ties -> smallest K), as in Witten & Tibshirani
choice <- lambda_star %>%
  arrange(desc(GAP_star), K) %>%
  slice(1)

list(
  per_K = as.data.frame(lambda_star),
  choice = as.data.frame(choice)
)




