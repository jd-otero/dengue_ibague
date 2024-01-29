## ---------- INLA ----------##

### 1. Set up ###

# Libraries
rm(list = ls())
packages <- c(
  "sf", "dplyr", "tmap", "tsModel", "INLA", "dlnm", "lubridate", "Hmisc",
  "dichromat", "RColorBrewer", "png"
)
lapply(packages, require, character.only = TRUE)
inla.setOption(inla.mode = "classic")

# Import
secciones <- st_read("data/secciones.shp") %>% st_drop_geometry()
secciones <- secciones[order(
  secciones$year,
  secciones$month,
  secciones$seccion
), ]

# Removing geometries without population (therefore no demographic data)
secciones <- secciones[which(secciones$population != 0), ]

# Renaming secciones (consecutive)
n_secciones <- nrow(secciones) / (12 * 6)
new_secciones <- rep(1:n_secciones, 12 * 6)
secciones$seccion <- new_secciones

### 2. Lagged variables (environmental) ###
# Number of lags
nlag <- 6

# Lag knots
lagknot <- equalknots(0:nlag, 2)

# Lag weather variables (up to 6 months)
lag_mean_temp <- tsModel::Lag(secciones$mean_temp,
  group = secciones$seccion,
  k = 0:nlag
)
lag_mean_temp <- lag_mean_temp[as.Date(paste(secciones$year,
  secciones$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_tot_precip <- tsModel::Lag(secciones$tot_precip,
  group = secciones$seccion,
  k = 0:nlag
)
lag_tot_precip <- lag_tot_precip[as.Date(paste(secciones$year,
  secciones$month,
  "1",
  sep = "-"
)) > "2013-6-30", ] # Over june 2013


lag_ndvi <- tsModel::Lag(secciones$ndvi,
  group = secciones$seccion,
  k = 0:nlag
)
lag_ndvi <- lag_ndvi[as.Date(paste(secciones$year,
  secciones$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_wet_days <- tsModel::Lag(secciones$wet_days,
  group = secciones$seccion,
  k = 0:nlag
)
lag_wet_days <- lag_wet_days[as.Date(paste(secciones$year,
  secciones$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_days_ov_32 <- tsModel::Lag(secciones$days_ov_32,
  group = secciones$seccion,
  k = 0:nlag
)
lag_days_ov_32 <- lag_days_ov_32[as.Date(paste(secciones$year,
  secciones$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

# Redefining time window for the rest of the data
secciones <- secciones[as.Date(paste(secciones$year,
  secciones$month, "1",
  sep = "-"
)) >
  "2013-6-30", ] # Over jun-2013

## 2.1. Crossbasis matrix ##
basis_mean_temp <- crossbasis(lag_mean_temp,
  argvar = list(fun = "ns", knots = equalknots(secciones$mean_temp, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_mean_temp) <- paste0(
  "basis_mean_temp.",
  colnames(basis_mean_temp)
)

basis_tot_precip <- crossbasis(lag_tot_precip,
  argvar = list(fun = "ns", knots = equalknots(secciones$tot_precip, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_tot_precip) <- paste0(
  "basis_tot_precip.",
  colnames(basis_tot_precip)
)

basis_ndvi <- crossbasis(lag_ndvi,
  argvar = list(fun = "ns", knots = equalknots(secciones$ndvi, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_ndvi) <- paste0(
  "basis_ndvi.",
  colnames(basis_ndvi)
)

basis_wet_days <- crossbasis(lag_wet_days,
  argvar = list(fun = "ns", knots = equalknots(secciones$wet_days, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_wet_days) <- paste0(
  "basis_wet_days.",
  colnames(basis_wet_days)
)
basis_days_ov_32 <- crossbasis(lag_days_ov_32,
  argvar = list(fun = "ns", knots = equalknots(secciones$days_ov_32, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_days_ov_32) <- paste0(
  "basis_days_ov_32.",
  colnames(basis_days_ov_32)
)

### 3. INLA ###

## 3.1 General data ##
secciones_inla <- secciones %>%
  summarise(
    seccion = seccion,
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

secciones_inla <- secciones_inla[order(
  secciones_inla$year_index,
  secciones_inla$month,
  secciones_inla$seccion
), ]


### 3.2 Prior definition ##
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

## 3.4 Spatial data import ##

nb_secciones_inla <- inla.read.graph("maps/secciones.adj") # INLA geometries

## 3.5 Formulas ##
# Base formula demographic (random effects)
formula_secciones.0 <- cases ~ 1 + f(month,
  replicate = seccion, model = "rw1",
  cyclic = TRUE, constr = FALSE,
  scale.model = FALSE, hyper = precision.prior, diagonal = 1
) +
  f(seccion,
    model = "bym2", replicate = year_index,
    graph = nb_secciones_inla, adjust.for.con.comp = TRUE,
    scale.model = TRUE, hyper = precision.prior, diagonal = 1
  )

## 3.5.1 Adding fixed effects
formula_secciones_spat <- update(formula_secciones.0, ~ . + density + kids +
  women + higher_ed + gas + sewage + garbage)

formula_secciones_temp <- update(formula_secciones.0, ~ . + basis_mean_temp +
  basis_ndvi + basis_tot_precip +
  basis_wet_days + basis_days_ov_32)

formula_secciones_spte <- update(formula_secciones.0, ~ . + density + kids +
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

secciones_spat <- fit_inla(secciones_inla, formula_secciones_spat)
secciones_temp <- fit_inla(secciones_inla, formula_secciones_temp)
secciones_spte <- fit_inla(secciones_inla, formula_secciones_spte)

## 3.7 Model Selection ##

secciones_compared_df <- data.frame(
  model = c("secciones_spat", "secciones_temp", "secciones_spte"),
  DIC = c(0, 0, 0),
  WAIC = c(0, 0, 0)
)
for (i in c("secciones_spat", "secciones_temp", "secciones_spte")) {
  temp <- get(i)
  secciones_compared_df[secciones_compared_df$model == i, "DIC"] <- temp$dic$dic
  secciones_compared_df[secciones_compared_df$model == i, "WAIC"] <- temp$waic$waic
}
secciones_compared_df$model <- c("Spatial", "Temporal", "Spatiotemporal")

# Best model is the one containing all predictors (remove others)
rm(temp)
rm(secciones_spat)
rm(secciones_temp)

### 4. Results report ###

## 4.1 Temporal variagbles (Lagged) ##
coeficients <- secciones_spte$summary.fixed$mean
covariance <- secciones_spte$misc$lincomb.derived.covariance.matrix

# Mean temperature
col_mean_temp <- grep("basis_mean_temp", secciones_spte$names.fixed)
predicted_mean_temp <- crosspred(basis_mean_temp,
  coef = coeficients[col_mean_temp],
  vcov = covariance[col_mean_temp, col_mean_temp],
  model.link = "log", bylag = 0.25, cen = round(mean(secciones$mean_temp), 0)
)
png("figures/inla_lag/mean_temp_secciones.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
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
      main = "Secciones", ylab = "",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Total precipitation
col_tot_precip <- grep("basis_tot_precip", secciones_spte$names.fixed)
predicted_tot_precip <- crosspred(basis_tot_precip,
  coef = coeficients[col_tot_precip],
  vcov = covariance[col_tot_precip, col_tot_precip],
  model.link = "log", bylag = 0.25, cen = round(mean(secciones$tot_precip), 0)
)
png("figures/inla_lag/tot_precip_secciones.png",
  width = 2200, height = 1700, res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_tot_precip$predvar,
  t(predicted_tot_precip$matRRfit),
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

# NDVI
col_ndvi <- grep("basis_ndvi", secciones_spte$names.fixed)
predicted_ndvi <- crosspred(basis_ndvi,
  coef = coeficients[col_ndvi], vcov = covariance[col_ndvi, col_ndvi],
  model.link = "log", bylag = 0.25, cen = round(mean(secciones$ndvi), 0)
)
png("figures/inla_lag/ndvi_secciones.png",
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
col_wet_days <- grep("basis_wet_days", secciones_spte$names.fixed)
predicted_wet_days <- crosspred(basis_wet_days,
  coef = coeficients[col_wet_days],
  vcov = covariance[col_wet_days, col_wet_days],
  model.link = "log", bylag = 0.25, cen = round(mean(secciones$wet_days), 0)
)
png("figures/inla_lag/wet_days_secciones.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_wet_days$predvar,
  t(predicted_wet_days$matRRfit),
  family = windowsFont(family = "Times New Roman"),
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
col_days_ov_32 <- grep("basis_days_ov_32", secciones_spte$names.fixed)
predicted_days_ov_32 <- crosspred(basis_days_ov_32,
  coef = coeficients[col_days_ov_32],
  vcov = covariance[col_days_ov_32, col_days_ov_32],
  model.link = "log", bylag = 0.25, cen = round(mean(secciones$days_ov_32), 0)
)
png("figures/inla_lag/days_ov_32_secciones.png",
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
df_results <- secciones_spte$summary.fixed[1:8, 3:5]
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
secciones_spatial <- df_results

### 5. Save results ###
secciones_results <- list(
  spatial = secciones_spatial,
  compared = secciones_compared_df,
  summary.fitted = secciones_spte$summary.fitted.values
)
saveRDS(secciones_results, "data/workspace/secciones_results.rds")
save.image(file = "data/workspace/secciones_inla.RData")
