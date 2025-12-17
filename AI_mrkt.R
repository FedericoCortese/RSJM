
# Data downloading and cleaning -------------------------------------------

## 0. Pacchetti
library(quantmod)   # dati da Yahoo Finance
library(TTR)        # runSD per la volatilità rolling
library(xts)
library(dplyr)
library(lubridate)
library(gtrendsR)   # Google Trends per sentiment
library(zoo) 
## 1. Parametri generali

# Circa 2000 giorni lavorativi indietro da oggi ≈ aprile 2018
start_date <- as.Date("2018-04-18")
end_date   <- Sys.Date()

vol_window <- 5L  # finestra per la volatilità rolling

# Tickers
ai_core      <- c("NVDA","AMD","TSM","ASML","SMCI")           # 5
ai_platform  <- c("MSFT","GOOGL","META","AAPL")               # 4
supply_chain <- c("MU","EQIX")                                # 2

ai_assets    <- c(ai_core, ai_platform, supply_chain)         # 11 in totale

mkt_symbol   <- "^NDX"    # indice di mercato per AIglob excess (Nasdaq 100)

## 2. Funzione helper per scaricare prezzi (Adjusted Close)

get_adj_prices <- function(symbols, from, to) {
  out_list <- lapply(symbols, function(sym) {
    ohlc <- getSymbols(
      Symbols     = sym,
      src         = "yahoo",
      from        = from,
      to          = to,
      auto.assign = FALSE
    )
    Ad(ohlc)
  })
  names(out_list) <- symbols
  prices_xts <- do.call(merge, out_list)
  colnames(prices_xts) <- symbols
  prices_xts
}

## 3. Scarico prezzi AI e mercato

# Prezzi per gli 11 asset AI
prices_ai  <- get_adj_prices(ai_assets, start_date, end_date)

# Prezzi per l'indice di mercato (Nasdaq 100)
mkt_ohlc   <- getSymbols(
  Symbols     = mkt_symbol,
  src         = "yahoo",
  from        = start_date,
  to          = end_date,
  auto.assign = FALSE
)
prices_mkt <- Ad(mkt_ohlc)
colnames(prices_mkt) <- "Mkt"

## 4. Log-rendimenti

# Merge ai + mercato su base comune prima di differenziare
prices_all <- merge(prices_ai, prices_mkt, join = "inner")

returns_all <- na.omit(diff(log(prices_all)))

# Split AI asset vs mercato
returns_ai  <- returns_all[, ai_assets]
returns_mkt <- returns_all[, "Mkt", drop = FALSE]

## 5. Volatilità rolling per ciascun asset AI

## 5. Volatilità rolling per ciascun asset AI (con rollapplyr) ----

vol_ai <- returns_ai  # stessa struttura (date + colonne), poi sovrascriviamo

for (sym in ai_assets) {
  vol_ai[, sym] <- rollapplyr(
    data   = returns_ai[, sym],
    width  = vol_window,
    FUN    = sd,
    fill   = NA,
    na.rm  = TRUE
  )
}

# Allineo per sicurezza
vol_ai <- vol_ai[index(returns_ai)]

## 6. AI globale e excess return

# AI globale equal-weighted sui 11 asset
r_AIglob <- xts(
  rowMeans(returns_ai, na.rm = TRUE),
  order.by = index(returns_ai)
)
colnames(r_AIglob) <- "AIglob"

# Mercato allineato
returns_mkt_aligned <- returns_mkt[index(returns_ai)]

# Excess return AIglob - Mkt
r_excess_AIglob <- r_AIglob - returns_mkt_aligned
colnames(r_excess_AIglob) <- "AIglob_Mkt_excess"

# Volatilità rolling dell'excess
vol_excess_AIglob <- runSD(r_excess_AIglob, n = vol_window, sample = TRUE)
colnames(vol_excess_AIglob) <- "AIglob_Mkt_excess_vol"


## 7. Sentiment AI con gtrendsR ----


# NA


## 8. Costruzione del dataset di feature

# Base temporale comune: giorni in cui abbiamo returns_ai
common_index <- index(returns_ai)

returns_ai_aligned         <- returns_ai[common_index]
vol_ai_aligned             <- vol_ai[common_index]
r_excess_AIglob_aligned    <- r_excess_AIglob[common_index]
vol_excess_AIglob_aligned  <- vol_excess_AIglob[common_index]

# Rinomino colonne di returns/vol per chiarezza (r_<TICKER>, vol_<TICKER>)
colnames(returns_ai_aligned) <- paste0("r_", colnames(returns_ai_aligned))
colnames(vol_ai_aligned)     <- paste0("vol_", colnames(vol_ai_aligned))

# Combino tutte le feature:
# - 11 log-rendimenti asset  (r_*)
# - 11 volatilitá asset      (vol_*)
# - 1 log-rend AIglob excess (AIglob_Mkt_excess)
# - 1 vol AIglob excess      (AIglob_Mkt_excess_vol)

features_xts <- merge(
  returns_ai_aligned,
  vol_ai_aligned,
  r_excess_AIglob_aligned,
  vol_excess_AIglob_aligned,
  join = "inner"
)

# Rimuovo righe con NA dovuti alla finestra della volatilità
features_xts <- na.omit(features_xts)

## 9. Check finale

print(head(features_xts))
print(dim(features_xts))
print(colnames(features_xts))

save(features_xts, file = "features_ai_dataset.RData")


# Data loading and description --------------------------------------------

load("features_ai_dataset.RData")

## 0. Preliminaries

library(xts)
library(zoo)
library(ggplot2)
library(corrplot)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# Rename object if needed
df <- features_xts   # ensure consistent naming

# Convert to data.frame while preserving dates
df_df <- data.frame(date = index(df), coredata(df))

## 1. Summary Statistics

summary_stats <- df_df %>%
  dplyr::select(-date) %>%
  summarise(
    across(
      everything(),
      list(
        mean  = ~ mean(.x, na.rm = TRUE),
        sd    = ~ sd(.x, na.rm = TRUE),
        min   = ~ min(.x, na.rm = TRUE),
        q25   = ~ quantile(.x, 0.25, na.rm = TRUE),
        median= ~ median(.x, na.rm = TRUE),
        q75   = ~ quantile(.x, 0.75, na.rm = TRUE),
        max   = ~ max(.x, na.rm = TRUE)
      )
    )
  )

print("Summary statistics:")
print(t(summary_stats))


## 2. Correlation Matrix + Corrplot

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

cols <- colnames(df_df)[-1]   # all feature names
n    <- length(cols)

# If >24 features, truncate or adjust; if <24, fill empty panels
target_panels <- 24

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

print("Scaled matrix Y created:")
print(dim(Y))
print(colnames(Y))
class(Y)


# Fit RSJM ----------------------------------------------------------------

source("Utils_sparse_robust_2.R")

zeta0=.15
lambda=1
K=3
qt=99

fit <- robust_sparse_jump(
  Y = Y,
  zeta0 = zeta0,
  lambda = lambda,
  K = K,
  tol = 1e-2,
  n_init = 1,
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
  mif=1
)

plot(fit$loss_vec$loss,type='l')

estW=round(fit$W,2)
colnames(estW)=colnames(Y)
estW

df_df$State=fit$s

plot_list <- list()

for (i in seq_along(cols)) {
  if (is.na(cols[i])) {
    # empty plot
    p <- ggplot() + theme_void()
  } else {
    p <- ggplot(df_df, aes(x = date, y = .data[[cols[i]]], color = factor(State))) +
      geom_point(size = 0.6) +
      scale_color_manual(
        values = c(
          "1" = "red",
          "2" = "yellow",
          "3" = "green"
        )
      ) +
      ggtitle(cols[i]) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"   # opzionale: nasconde la legenda
      )
  }
  plot_list[[i]] <- p
}


x11()
grid.arrange(grobs = plot_list, ncol = 4, nrow = 6)

round(tapply(df_df$r_NVDA,df_df$State,mean)*100,2)
round(tapply(df_df$r_NVDA,df_df$State,sd)*100,2)
