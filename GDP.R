euro_countries <- c(
  "Austria",
  "Belgium",
  "Estonia",
  "Finland",
  "France",
  "Greece",
  "Ireland",
  "Italy",
  "Malta",
  "Netherlands",
  "Portugal",
  "Spain"
)

load("gmd.RData")

names(dati)
str(dati)
unique(dati$countryname)
unique(dati$year)

cor(dati[,1:5])

dati_euro <- dati[dati$countryname %in% euro_countries, ]

vars <- c("cons_GDP","inv_GDP","finv_GDP","exports_GDP","imports_GDP")

countries <- sort(unique(dati_euro$countryname))

x11()
par(mfrow = c(3,2), mar = c(4,4,2,1))
for (v in vars) {
  wide <- xtabs(as.formula(paste(v, "~ year + countryname")), data = dati_euro)
  years <- as.numeric(rownames(wide))
  
  matplot(years, wide,
          type = "l",
          lty  = 1:ncol(wide),
          xlab = "Year",
          ylab = v,
          main = v,
          col  = 1:ncol(wide))
}

plot.new()
legend("center", legend = countries, 
       col = 1:length(countries), lty = 1:ncol(wide), cex = 0.8, ncol = 2)

#

dati_euro_wide <- dati_euro %>%
  pivot_wider(names_from = countryname,
              values_from = -c(countryname, year)) %>%
  arrange(year)

# fit ---------------------------------------------------------------------



Y_input=dati_euro_wide[,-1]
dim(Y_input)  # 49x60

source("Utils_feat_weight_robust.R")
Y_input <- as.data.frame(Y_input)

fit=feat_weight_jump(Y_input,zeta0=.05,lambda=0.5,K=3,
                     tol     = 1e-16,
                     n_init  = 5,
                     n_outer = 50,
                     n_inner = 10,
                     verbose=T)

est_W=fit$W
colnames(est_W)=colnames(Y_input)
rownames(est_W) <- c("1","2","3")
round(est_W,2)

res_fit=data.frame(year=dati_euro_wide$year,Y_input,s=fit$s)


years <- res_fit$year

state_cols <- c("magenta", "cyan", "orange")

# helper
plot_panel <- function(prefix, main_title = prefix) {
  cols <- grep(paste0("^", prefix, "_"), names(res_fit), value = TRUE)
  if (length(cols) == 0) stop("No columns for: ", prefix)
  
  mat <- as.matrix(res_fit[, cols, drop = FALSE])
  yrng <- range(mat, na.rm = TRUE)
  
  plot(years, mat[,1], type = "n", ylim = yrng,
       xlab = "Year", ylab = prefix, main = main_title)
  
  s <- as.integer(res_fit$s)
  for (k in seq_len(ncol(mat))) {
    y <- mat[, k]
    lines(years, y, col = "grey80", lwd = 1)
    # punti colorati dal vettore s (uno per anno)
    points(years, y, pch = 16, col = state_cols[s], cex = 0.6)
  }
}

x11()
par(mfrow = c(3,2), mar = c(4,4,2,1))

plot_panel("cons_GDP",    "cons_GDP (all countries)")
plot_panel("inv_GDP",     "inv_GDP (all countries)")
plot_panel("finv_GDP",    "finv_GDP (all countries)")
plot_panel("exports_GDP", "exports_GDP (all countries)")
plot_panel("imports_GDP", "imports_GDP (all countries)")
plot.new()
legend("topright", legend = c("s=1","s=2","s=3"), col = state_cols, pch = 16, bty = "n")

# Pesi
vars_all <- colnames(est_W)
vars_ic  <- setdiff(vars_all, "cross_corr")

split_parts <- strsplit(vars_ic, "_")
indicator <- vapply(split_parts, function(x) paste(x[1:2], collapse = "_"), character(1))

country <- vapply(split_parts, function(x) paste(x[-c(1,2)], collapse = "_"), character(1))
ind_levels <- c("cons_GDP","inv_GDP","finv_GDP","exports_GDP","imports_GDP")
cty_levels <- sort(unique(country))  # 12 paesi

make_grid <- function(state_k) {
  M <- matrix(NA_real_, nrow = length(ind_levels), ncol = length(cty_levels),
              dimnames = list(ind_levels, cty_levels))
  for (j in seq_along(vars_ic)) {
    M[indicator[j], country[j]] <- est_W[state_k, vars_ic[j]]
  }
  M
}

grid_list <- lapply(1:3, make_grid)

plot_weight_grid <- function(M, main, digits = 3) {
  
  Mr <- M[nrow(M):1, , drop = FALSE]
  
  image(x = seq_len(ncol(Mr)), y = seq_len(nrow(Mr)), z = t(Mr),
        axes = FALSE, xlab = "", ylab = "", main = main)
  
  axis(1, at = seq_len(ncol(M)), labels = colnames(M), las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(nrow(M)), labels = rev(rownames(M)), las = 2, cex.axis = 0.9)
  box()
  
  # numeri dentro le celle
  lab <- format(round(Mr, digits), nsmall = digits)
  for (i in seq_len(nrow(Mr))) {
    for (j in seq_len(ncol(Mr))) {
      if (is.finite(Mr[i, j])) text(j, i, labels = lab[i, j], cex = 0.65)
    }
  }
}

x11()
par(mfrow = c(1,3), mar = c(8,4,3,1))

for (k in 1:3) {
  plot_weight_grid(grid_list[[k]],
                   main = paste0("Weights grid (s=", k, ")"),
                   digits = 3)
  mtext("Indicator × Country", side = 3, line = 0.2, cex = 0.8, col = state_cols[k])
}


