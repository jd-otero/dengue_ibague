## ---------- INLA ----------##

### 1. Set up ###

# Libraries
rm(list = ls())
packages <- c(
  "sf", "dplyr", "tsModel", "INLA", "dlnm", "lubridate", "Hmisc",
  "dichromat", "RColorBrewer", "png"
)
lapply(packages, require, character.only = TRUE)
inla.setOption(inla.mode = "classic")

# Import
sectores <- st_read("data/sectores.shp") %>% st_drop_geometry()
sectores <- sectores[order(
  sectores$year,
  sectores$month,
  sectores$sector
), ]

# Removing geometries without population (therefore no demographic data)
sectores <- sectores[which(sectores$population != 0), ]

# Renaming sectores (consecutive)
n_sectores <- nrow(sectores) / (12 * 6)
new_sectores <- rep(1:n_sectores, 12 * 6)
sectores$sector <- new_sectores

### 2. Lagged variables (environmental) ###
# Number of lags
nlag <- 6

# Lag knots
lagknot <- equalknots(0:nlag, 2)

# Lag weather variables (up to 6 months)
lag_mean_temp <- tsModel::Lag(sectores$mean_temp,
  group = sectores$sector,
  k = 0:nlag
)
lag_mean_temp <- lag_mean_temp[as.Date(paste(sectores$year,
  sectores$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_tot_precip <- tsModel::Lag(sectores$tot_precip,
  group = sectores$sector,
  k = 0:nlag
)
lag_tot_precip <- lag_tot_precip[as.Date(paste(sectores$year,
  sectores$month,
  "1",
  sep = "-"
)) > "2013-6-30", ] # Over june 2013


lag_ndvi <- tsModel::Lag(sectores$ndvi,
  group = sectores$sector,
  k = 0:nlag
)
lag_ndvi <- lag_ndvi[as.Date(paste(sectores$year,
  sectores$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_wet_days <- tsModel::Lag(sectores$wet_days,
  group = sectores$sector,
  k = 0:nlag
)
lag_wet_days <- lag_wet_days[as.Date(paste(sectores$year,
  sectores$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_days_ov_32 <- tsModel::Lag(sectores$days_ov_32,
  group = sectores$sector,
  k = 0:nlag
)
lag_days_ov_32 <- lag_days_ov_32[as.Date(paste(sectores$year,
  sectores$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

# Redefining time window for the rest of the data
sectores <- sectores[as.Date(paste(sectores$year,
  sectores$month, "1",
  sep = "-"
)) >
  "2013-6-30", ] # Over jun-2013

## 2.1. Crossbasis matrix ##
basis_mean_temp <- crossbasis(lag_mean_temp,
  argvar = list(fun = "ns", knots = equalknots(sectores$mean_temp, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_mean_temp) <- paste0(
  "basis_mean_temp.",
  colnames(basis_mean_temp)
)

basis_tot_precip <- crossbasis(lag_tot_precip,
  argvar = list(fun = "ns", knots = equalknots(sectores$tot_precip, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_tot_precip) <- paste0(
  "basis_tot_precip.",
  colnames(basis_tot_precip)
)

basis_ndvi <- crossbasis(lag_ndvi,
  argvar = list(fun = "ns", knots = equalknots(sectores$ndvi, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_ndvi) <- paste0(
  "basis_ndvi.",
  colnames(basis_ndvi)
)

basis_wet_days <- crossbasis(lag_wet_days,
  argvar = list(fun = "ns", knots = equalknots(sectores$wet_days, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_wet_days) <- paste0(
  "basis_wet_days.",
  colnames(basis_wet_days)
)
basis_days_ov_32 <- crossbasis(lag_days_ov_32,
  argvar = list(fun = "ns", knots = equalknots(sectores$days_ov_32, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_days_ov_32) <- paste0(
  "basis_days_ov_32.",
  colnames(basis_days_ov_32)
)

### 3. INLA ###

## 3.1 General data ##
sectores_inla <- sectores %>%
  summarise(
    sector = sector,
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

sectores_inla <- sectores_inla[order(
  sectores_inla$year_index,
  sectores_inla$month,
  sectores_inla$sector
), ]


### 3.2 Prior definition ##
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

## 3.4 Spatial data import ##

nb_sectores_inla <- inla.read.graph("maps/sectores.adj") # INLA geometries

## 3.5 Formulas ##
# Base formula demographic (random effects)
formula_sectores.0 <- cases ~ 1 + f(month,
  replicate = sector, model = "rw1",
  cyclic = TRUE, constr = FALSE,
  scale.model = FALSE, hyper = precision.prior, diagonal = 1
) +
  f(sector,
    model = "bym2", replicate = year_index,
    graph = nb_sectores_inla, adjust.for.con.comp = TRUE,
    scale.model = TRUE, hyper = precision.prior, diagonal = 1
  )

## 3.5.1 Adding fixed effects
formula_sectores_spat <- update(formula_sectores.0, ~ . + density + kids +
  women + higher_ed + gas + sewage + garbage)

formula_sectores_temp <- update(formula_sectores.0, ~ . + basis_mean_temp +
  basis_ndvi + basis_tot_precip +
  basis_wet_days + basis_days_ov_32)

formula_sectores_spte <- update(formula_sectores.0, ~ . + density + kids +
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

sectores_spat <- fit_inla(sectores_inla, formula_sectores_spat)
sectores_temp <- fit_inla(sectores_inla, formula_sectores_temp)
sectores_spte <- fit_inla(sectores_inla, formula_sectores_spte)

## 3.7 Model Selection ##

sectores_compared_df <- data.frame(
  model = c("sectores_spat", "sectores_temp", "sectores_spte"),
  DIC = c(0, 0, 0),
  WAIC = c(0, 0, 0)
)
for (i in c("sectores_spat", "sectores_temp", "sectores_spte")) {
  temp <- get(i)
  sectores_compared_df[sectores_compared_df$model == i, "DIC"] <- temp$dic$dic
  sectores_compared_df[sectores_compared_df$model == i, "WAIC"] <- temp$waic$waic
}
sectores_compared_df$model <- c("Spatial", "Temporal", "Spatiotemporal")

# Best model is the one containing all predictors (remove others)
rm(temp)
rm(sectores_spat)
rm(sectores_temp)

### 4. Results report ###

## 4.1 Temporal variagbles (Lagged) ##
coeficients <- sectores_spte$summary.fixed$mean
covariance <- sectores_spte$misc$lincomb.derived.covariance.matrix

# Mean temperature
col_mean_temp <- grep("basis_mean_temp", sectores_spte$names.fixed)
predicted_mean_temp <- crosspred(basis_mean_temp,
  coef = coeficients[col_mean_temp],
  vcov = covariance[col_mean_temp, col_mean_temp],
  model.link = "log", bylag = 0.25, cen = round(mean(sectores$mean_temp), 0)
)
png("figures/inla_lag/mean_temp_sectores.png",
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
      main = "Sectores", ylab = "",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Total precipitation
col_tot_precip <- grep("basis_tot_precip", sectores_spte$names.fixed)
predicted_tot_precip <- crosspred(basis_tot_precip,
  coef = coeficients[col_tot_precip],
  vcov = covariance[col_tot_precip, col_tot_precip],
  model.link = "log", bylag = 0.25, cen = round(mean(sectores$tot_precip), 0)
)
png("figures/inla_lag/tot_precip_sectores.png",
  width = 2200, height = 1700, res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_tot_precip$predvar,
  t(predicted_tot_precip$matRRfit),
  family = windowsFont(family = "Times New Roman"),
  plot.axes = {
    axis(1,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 1, 2, 3, 4, 5, 6), at = c(0, 1, 2, 3, 4, 5, 6)
    )
    axis(2,
      cex.axis = 2.5, family = windowsFont(family = "Times New Roman"),
      labels = c(0, 20, 40, 60, 80, 100), at = c(0, 20, 40, 60, 80, 100)
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
col_ndvi <- grep("basis_ndvi", sectores_spte$names.fixed)
predicted_ndvi <- crosspred(basis_ndvi,
  coef = coeficients[col_ndvi], vcov = covariance[col_ndvi, col_ndvi],
  model.link = "log", bylag = 0.25, cen = round(mean(sectores$ndvi), 0)
)
png("figures/inla_lag/ndvi_sectores.png",
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
col_wet_days <- grep("basis_wet_days", sectores_spte$names.fixed)
predicted_wet_days <- crosspred(basis_wet_days,
  coef = coeficients[col_wet_days],
  vcov = covariance[col_wet_days, col_wet_days],
  model.link = "log", bylag = 0.25, cen = round(mean(sectores$wet_days), 0)
)
png("figures/inla_lag/wet_days_sectores.png",
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
col_days_ov_32 <- grep("basis_days_ov_32", sectores_spte$names.fixed)
predicted_days_ov_32 <- crosspred(basis_days_ov_32,
  coef = coeficients[col_days_ov_32],
  vcov = covariance[col_days_ov_32, col_days_ov_32],
  model.link = "log", bylag = 0.25, cen = round(mean(sectores$days_ov_32), 0)
)
png("figures/inla_lag/days_ov_32_sectores.png",
  width = 2200, height = 1700,
  res = 300
)
par(mar = c(6, 6.7, 3, 0), mgp = c(4.5, 1.5, 0))
filled.contour(seq(0, nlag, 0.25),
  predicted_days_ov_32$predvar,
  t(predicted_days_ov_32$matRRfit),
  family = windowsFont(family = "Times New Roman"),
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
df_results <- sectores_spte$summary.fixed[1:8, 3:5]
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
sectores_spatial <- df_results

### 5. Save results ###
sectores_results <- list(
  spatial = sectores_spatial,
  compared = sectores_compared_df,
  summary.fitted = sectores_spte$summary.fitted.values
)
saveRDS(sectores_results, "data/workspace/sectores_results.rds")
save.image(file = "data/workspace/sectores_inla.RData")
