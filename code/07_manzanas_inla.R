## ---------- INLA ----------##

### 1. Set up ###

# Libraries
rm(list = ls())
packages <- c(
  "sf", "dplyr", "tsModel", "INLA", "dlnm", "lubridate", "Hmisc",
  "dichromat", "RColorBrewer", "png"
)
lapply(packages, require, character.only = TRUE)
inla.setOption(inla.mode = "experimental")

# Import
manzanas <- st_read("data/manzanas.shp") %>% st_drop_geometry()
manzanas <- manzanas[order(
  manzanas$year,
  manzanas$month,
  manzanas$manzana
), ]

# Removing geometries without population (therefore no demographic data)
manzanas <- manzanas[which(manzanas$population != 0), ]

# Renaming manzanas (consecutive)
n_manzanas <- nrow(manzanas) / (12 * 6)
new_manzanas <- rep(1:n_manzanas, 12 * 6)
manzanas$manzana <- new_manzanas

### 2. Lagged variables (environmental) ###
# Number of lags
nlag <- 6

# Lag knots
lagknot <- equalknots(0:nlag, 2)

# Lag weather variables (up to 6 months)
lag_mean_temp <- tsModel::Lag(manzanas$mean_temp,
  group = manzanas$manzana,
  k = 0:nlag
)
lag_mean_temp <- lag_mean_temp[as.Date(paste(manzanas$year,
  manzanas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_tot_precip <- tsModel::Lag(manzanas$tot_precip,
  group = manzanas$manzana,
  k = 0:nlag
)
lag_tot_precip <- lag_tot_precip[as.Date(paste(manzanas$year,
  manzanas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ] # Over june 2013


lag_ndvi <- tsModel::Lag(manzanas$ndvi,
  group = manzanas$manzana,
  k = 0:nlag
)
lag_ndvi <- lag_ndvi[as.Date(paste(manzanas$year,
  manzanas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_wet_days <- tsModel::Lag(manzanas$wet_days,
  group = manzanas$manzana,
  k = 0:nlag
)
lag_wet_days <- lag_wet_days[as.Date(paste(manzanas$year,
  manzanas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_days_ov_32 <- tsModel::Lag(manzanas$days_ov_32,
  group = manzanas$manzana,
  k = 0:nlag
)
lag_days_ov_32 <- lag_days_ov_32[as.Date(paste(manzanas$year,
  manzanas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

# Redefining time window for the rest of the data
manzanas <- manzanas[as.Date(paste(manzanas$year,
  manzanas$month, "1",
  sep = "-"
)) >
  "2013-6-30", ] # Over jun-2013

## 2.1. Crossbasis matrix ##
basis_mean_temp <- crossbasis(lag_mean_temp,
  argvar = list(fun = "ns", knots = equalknots(manzanas$mean_temp, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_mean_temp) <- paste0(
  "basis_mean_temp.",
  colnames(basis_mean_temp)
)

basis_tot_precip <- crossbasis(lag_tot_precip,
  argvar = list(fun = "ns", knots = equalknots(manzanas$tot_precip, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_tot_precip) <- paste0(
  "basis_tot_precip.",
  colnames(basis_tot_precip)
)

basis_ndvi <- crossbasis(lag_ndvi,
  argvar = list(fun = "ns", knots = equalknots(manzanas$ndvi, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_ndvi) <- paste0(
  "basis_ndvi.",
  colnames(basis_ndvi)
)

basis_wet_days <- crossbasis(lag_wet_days,
  argvar = list(fun = "ns", knots = equalknots(manzanas$wet_days, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_wet_days) <- paste0(
  "basis_wet_days.",
  colnames(basis_wet_days)
)
basis_days_ov_32 <- crossbasis(lag_days_ov_32,
  argvar = list(fun = "ns", knots = equalknots(manzanas$days_ov_32, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_days_ov_32) <- paste0(
  "basis_days_ov_32.",
  colnames(basis_days_ov_32)
)

### 3. INLA ###

## 3.1 General data ##
manzanas_inla <- manzanas %>%
  summarise(
    manzana = manzana,
    year_index = year - 2012,
    month = month,
    density = scale(density, scale = T, center = T),
    population = population,
    houses = houses,
    sewage = scale(sewage / houses, scale = T, center = T),
    gas = scale(gas / houses, scale = T, center = T),
    garbage = scale(garbage / houses, scale = T, center = T),
    higher_ed = scale(higher_ed / population, scale = T, center = T),
    kids = scale(kids / population, scale = T, center = T),
    women = scale(women / population, scale = T, center = T),
    cases = cases
  ) %>%
  replace(is.na(.), 0)

manzanas_inla <- manzanas_inla[order(
  manzanas_inla$year_index,
  manzanas_inla$month,
  manzanas_inla$manzana
), ]


### 3.2 Prior definition ##
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

## 3.4 Spatial data import ##

nb_manzanas_inla <- inla.read.graph("maps/manzanas.adj") # INLA geometries

## 3.5 Formulas ##
# Base formula demographic (random effects)
formula_manzanas.0 <- cases ~ 1 + f(month,
  replicate = manzana, model = "rw1",
  cyclic = TRUE, constr = FALSE,
  scale.model = FALSE, hyper = precision.prior, diagonal = 1
) +
  f(manzana,
    model = "bym2", replicate = year_index,
    graph = nb_manzanas_inla, adjust.for.con.comp = TRUE,
    scale.model = TRUE, hyper = precision.prior, diagonal = 1
  )

## 3.5.1 Adding fixed effects
formula_manzanas_spat <- update(formula_manzanas.0, ~ . + density + kids +
  women + higher_ed + gas + sewage + garbage)

formula_manzanas_temp <- update(formula_manzanas.0, ~ . + basis_mean_temp +
  basis_ndvi + basis_tot_precip +
  basis_wet_days + basis_days_ov_32)

formula_manzanas_spte <- update(formula_manzanas.0, ~ . + density + kids +
  women + higher_ed + gas + sewage + garbage + basis_mean_temp + basis_ndvi +
  basis_tot_precip + basis_wet_days +
  basis_days_ov_32)

### 3.6 Model implementation
fit_inla <- function(data, formula) {
  model <- inla(
    formula = formula,
    data = data,
    family = "nbinomial",
    offset = log(population),
    control.inla = list(strategy = "adaptive", int.strategy = "ccd"),
    control.compute = list(
      dic = TRUE, config = TRUE,
      cpo = TRUE, return.marginals = TRUE,
      waic = TRUE
    ),
    control.fixed = list(
      correlation.matrix = TRUE,
      prec.intercept = 1, prec = 1
    ),
    control.predictor = list(link = 1, compute = TRUE),
    verbose = TRUE,
    debug = TRUE
  )
  return(model)
}

manzanas_spat <- fit_inla(manzanas_inla, formula_manzanas_spat)
manzanas_temp <- fit_inla(manzanas_inla, formula_manzanas_temp)
manzanas_spte <- fit_inla(manzanas_inla, formula_manzanas_spte)

## 3.7 Model Selection ##

manzanas_compared_df <- data.frame(
  model = c("manzanas_spat", "manzanas_temp", "manzanas_spte"),
  DIC = c(0, 0, 0),
  WAIC = c(0, 0, 0)
)
for (i in c("manzanas_spat", "manzanas_temp", "manzanas_spte")) {
  temp <- get(i)
  manzanas_compared_df[manzanas_compared_df$model == i, "DIC"] <- temp$dic$dic
  manzanas_compared_df[manzanas_compared_df$model == i, "WAIC"] <- temp$waic$waic
}
manzanas_compared_df$model <- c("Spatial", "Temporal", "Spatiotemporal")

# Best model is the one containing all predictors (remove others)
rm(temp)
rm(manzanas_spat)
rm(manzanas_temp)

### 4. Results report ###

## 4.1 Temporal variagbles (Lagged) ##
coeficients <- manzanas_spte$summary.fixed$mean
covariance <- manzanas_spte$misc$lincomb.derived.covariance.matrix

# Mean temperature
col_mean_temp <- grep("basis_mean_temp", manzanas_spte$names.fixed)
predicted_mean_temp <- crosspred(basis_mean_temp,
  coef = coeficients[col_mean_temp],
  vcov = covariance[col_mean_temp, col_mean_temp],
  model.link = "log", bylag = 0.25, cen = round(mean(manzanas$mean_temp), 0)
)
png("figures/inla_lag/mean_temp_manzanas.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
plot_mean_temp <- filled.contour(seq(0, nlag, 0.25),
  predicted_mean_temp$predvar,
  t(predicted_mean_temp$matRRfit),
  family = windowsFont(family = "Times New Roman"),
  plot.axes = {
    axis(1,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 1, 2, 3, 4, 5, 6), at = c(0, 1, 2, 3, 4, 5, 6)
    )
    axis(2,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(20, 25, 30, 35, 40), at = c(20, 25, 30, 35, 40)
    )
  }, key.axes = axis(4,
    cex.axis = 2.5,
    family = windowsFont(family = "Times New Roman")
  ),
  plot.title = {
    title(
      main = "Manzanas", ylab = "",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Total precipitation
col_tot_precip <- grep("basis_tot_precip", manzanas_spte$names.fixed)
predicted_tot_precip <- crosspred(basis_tot_precip,
  coef = coeficients[col_tot_precip],
  vcov = covariance[col_tot_precip, col_tot_precip],
  model.link = "log", bylag = 0.25, cen = round(mean(manzanas$tot_precip), 0)
)
png("figures/inla_lag/tot_precip_manzanas.png",
  width = 2200, height = 1700, res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
plot_tot_precip <- plot_tot_precip <- filled.contour(seq(0, nlag, 0.25),
  predicted_tot_precip$predvar,
  t(predicted_tot_precip$matRRfit),
  plot.axes = {
    axis(1,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 1, 2, 3, 4, 5, 6), at = c(0, 1, 2, 3, 4, 5, 6)
    )
    axis(2,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 2, 4, 6, 8, 10), at = c(0, 2, 4, 6, 8, 10)
    )
  }, key.axes = axis(4,
    cex.axis = 2.5,
    family = windowsFont(family = "Times New Roman")
  ),
  plot.title = {
    title(
      main = "", ylab = "",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# NDVI
col_ndvi <- grep("basis_ndvi", manzanas_spte$names.fixed)
predicted_ndvi <- crosspred(basis_ndvi,
  coef = coeficients[col_ndvi], vcov = covariance[col_ndvi, col_ndvi],
  model.link = "log", bylag = 0.25, cen = round(mean(manzanas$ndvi), 0)
)
png("figures/inla_lag/ndvi_manzanas.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_ndvi$predvar,
  t(predicted_ndvi$matRRfit),
  family = windowsFont(family = "Times New Roman"),
  plot.axes = {
    axis(1,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 1, 2, 3, 4, 5, 6), at = c(0, 1, 2, 3, 4, 5, 6)
    )
    axis(2,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), at = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    )
  }, key.axes = axis(4,
    cex.axis = 2.5,
    family = windowsFont(family = "Times New Roman")
  ),
  plot.title = {
    title(
      main = "", ylab = "",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Wet days
col_wet_days <- grep("basis_wet_days", manzanas_spte$names.fixed)
predicted_wet_days <- crosspred(basis_wet_days,
  coef = coeficients[col_wet_days],
  vcov = covariance[col_wet_days, col_wet_days],
  model.link = "log", bylag = 0.25, cen = round(mean(manzanas$wet_days), 0)
)
png("figures/inla_lag/wet_days_manzanas.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_wet_days$predvar,
  t(predicted_wet_days$matRRfit),
  plot.axes = {
    axis(1,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 1, 2, 3, 4, 5, 6), at = c(0, 1, 2, 3, 4, 5, 6)
    )
    axis(2,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 5, 10, 15, 20), at = c(0, 5, 10, 15, 20)
    )
  }, key.axes = axis(4,
    cex.axis = 2.5,
    family = windowsFont(family = "Times New Roman")
  ),
  plot.title = {
    title(
      main = "", ylab = "",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Days Over 32
col_days_ov_32 <- grep("basis_days_ov_32", manzanas_spte$names.fixed)
predicted_days_ov_32 <- crosspred(basis_days_ov_32,
  coef = coeficients[col_days_ov_32],
  vcov = covariance[col_days_ov_32, col_days_ov_32],
  model.link = "log", bylag = 0.25, cen = round(mean(manzanas$days_ov_32), 0)
)
png("figures/inla_lag/days_ov_32_manzanas.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_days_ov_32$predvar,
  t(predicted_days_ov_32$matRRfit),
  plot.axes = {
    axis(1,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 1, 2, 3, 4, 5, 6), at = c(0, 1, 2, 3, 4, 5, 6)
    )
    axis(2,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 5, 10, 15, 20, 25, 30), at = c(0, 5, 10, 15, 20, 25, 30)
    )
  }, key.axes = axis(4,
    cex.axis = 2.5,
    family = windowsFont(family = "Times New Roman")
  ),
  plot.title = {
    title(
      main = "", ylab = "",
      xlab = "Lag [Months]", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

## 4.2 Spatial variables (non-lagged) ##
df_results <- manzanas_spte$summary.fixed[1:8, 3:5]
rownames(df_results) <- c(
  "intercept", "density", "kids", "women", "higher_ed", "gas", "sewage",
  "garbage"
)
df_results$rr <- ifelse(sign(df_results$"0.5quant") > 0,
  exp(df_results$"0.5quant"),
  1 - exp(df_results$"0.5quant")
)
df_results <- round(df_results, 3)
df_results[1, "rr"] <- NaN
manzanas_spatial <- df_results

### 5. Save results ###
manzanas_results <- list(
  spatial = manzanas_spatial,
  compared = manzanas_compared_df,
  summary.fitted = manzanas_spte$summary.fitted.values
)
saveRDS(manzanas_results, "data/workspace/manzanas_results.rds")
save.image(file = "data/workspace/manzanas_inla.RData")
