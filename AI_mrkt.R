#
# Data downloading and cleaning -------------------------------------------

library(quantmod)   # dati da Yahoo Finance
library(TTR)        # runMean / runSD per rolling stats
library(xts)
library(dplyr)
library(lubridate)
library(zoo)

start_date <- as.Date("2023-01-01") # Ossia da quando ChatGPT è esplosa, approx
end_date   <- Sys.Date()

vol_window <- 5L  # ampiezza finestra rolling (es. 5 giorni)

# --- Tickers (con descrizione azienda) ---------------------------------------
# NVDA = NVIDIA 
# AMD  = Advanced Micro Devices, Inc. (CPU/GPU)
# TSM  = Taiwan Semiconductor Manufacturing Company (TSMC) 
# ASML = ASML Holding N.V. 
# SMCI = Super Micro Computer, Inc. 
# MSFT = Microsoft 
# GOOGL= Alphabet 
# META = Meta 
# AAPL = Apple 
# MU   = Micron Technology, Inc. (memorie DRAM/NAND)
# EQIX = Equinix, Inc. (data center / digital infrastructure) 

ai_core      <- c("NVDA","AMD","TSM","ASML","SMCI")   # 5
ai_platform  <- c("MSFT","GOOGL","META","AAPL")       # 4
supply_chain <- c("MU","EQIX")                        # 2

ai_assets    <- c(ai_core, ai_platform, supply_chain) # 11

get_adj_prices <- function(symbols, from, to) {
  out_list <- lapply(symbols, function(sym) {
    ohlc <- getSymbols(
      Symbols     = sym,
      src         = "yahoo",
      from        = from,
      to          = to,
      auto.assign = FALSE
    )
    Ad(ohlc)  # Adjusted close
  })
  names(out_list) <- symbols
  prices_xts <- do.call(merge, out_list)
  colnames(prices_xts) <- symbols
  prices_xts
}

prices_ai   <- get_adj_prices(ai_assets, start_date, end_date)

# Log-rendimenti
returns_all <- na.omit(diff(log(prices_ai)))

# MA
ret_ma <- rollapply(
  data      = returns_all,
  width     = vol_window,
  FUN       = function(x) colMeans(x, na.rm = TRUE),
  by.column = FALSE,
  align     = "right",
  fill      = NA
)

# Moving sd
ret_sd <- rollapply(
  data      = returns_all,
  width     = vol_window,
  FUN       = function(x) apply(x, 2, sd, na.rm = TRUE),
  by.column = FALSE,
  align     = "right",
  fill      = NA
)


colnames(ret_ma) <- paste0(colnames(ret_ma), "_ma", vol_window)
colnames(ret_sd) <- paste0(colnames(ret_sd), "_sd", vol_window)

rolling_stats <- merge(ret_ma, ret_sd)

features_ai_dataset=na.omit(rolling_stats)
dim(features_ai_dataset)

matplot(features_ai_dataset[,1:11], type='l')

matplot(features_ai_dataset[,12:22], type='l')

returns_all=returns_all[-(1:(vol_window-1)),]

save(features_ai_dataset,returns_all, file = "features_ai_dataset.RData")

# --------------------------------------------

load("features_ai_dataset.RData")

library(xts)
library(zoo)
library(ggplot2)
library(corrplot)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

df_df=data.frame(date=index(features_ai_dataset), coredata(features_ai_dataset))

apply(df_df[,-1],2,mean)

cor_mat <- cor(df_df[,2:12], use = "pairwise.complete.obs")

corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "red",   # <<< colore numerini
  tl.cex = 0.7,
  tl.col = "red"         # colore dei nomi delle variabili
)


## 3. Time Series Plots (6 × 4 grid)

cols <- colnames(df_df)[-1]   
n    <- length(cols)

target_panels <- 22

if (n < target_panels) {
  cols <- c(cols, rep(NA, target_panels - n))
} else if (n > target_panels) {
  cols <- cols[1:target_panels]
}

plot_list <- list()

for (i in seq_along(cols)) {
  if (is.na(cols[i])) {
    # empty plot
    p <- ggplot() + theme_void()
  } else {
    p <- ggplot(df_df, aes_string(x = "date", y = cols[i])) +
      geom_line(color = "steelblue") +
      ggtitle(cols[i]) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(size = 8),   # titolo più piccolo
        axis.title.x = element_blank(),        # rimuove label asse X
        axis.title.y = element_blank()         # rimuove label asse Y
      )
  }
  plot_list[[i]] <- p
}

# 6x4 grid
x11()
grid.arrange(grobs = plot_list, ncol = 4, nrow = 6)


## 4. Scaling all features → matrix Y

Y <- scale(df_df %>% dplyr::select(-date))


# Fit FWJM ----------------------------------------------------------------

source("Utils_feat_weight_robust.R")


zeta0=.1
lambda=.5
K=2

fit_ai <- feat_weight_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = NULL,
  n_init = 1,
  n_outer = 30,
  n_inner=10,
  alpha = 0.1,
  verbose = T,
  mif=1)

plot(fit_ai$s,type='l')

plot(fit_ai$loss_vec$loss)

est_W=fit_ai$W
colnames(est_W)=colnames(Y)
round(est_W,3)

res_ai=data.frame(date=df_df$date,returns_all, cluster=fit_ai$s)

tapply(res_ai$NVDA, res_ai$cluster, mean)*100
tapply(res_ai$NVDA, res_ai$cluster, sd)*100

tapply(res_ai$META, res_ai$cluster, mean)*100
tapply(res_ai$META, res_ai$cluster, sd)*100

