## ---------- INLA ----------##

### 1. Set up ###

# Libraries
rm(list = ls())
packages <- c(
  "sf", "dplyr", "tsModel", "INLA", "dlnm", "lubridate", "Hmisc",
  "dichromat", "RColorBrewer", "png", "extrafont"
)
lapply(packages, require, character.only = TRUE)
inla.setOption(inla.mode = "classic")

# Import
comunas <- st_read("data/comunas.shp") %>% st_drop_geometry()
comunas <- comunas[order(
  comunas$year,
  comunas$month,
  comunas$comuna
), ]

# Removing geometries without population (therefore no demographic data)
comunas <- comunas[which(comunas$population != 0), ]

# Renaming comunas (consecutive)
n_comunas <- nrow(comunas) / (12 * 6)
new_comunas <- rep(1:n_comunas, 12 * 6)
comunas$comuna <- new_comunas

### 2. Lagged variables (environmental) ###
# Number of lags
nlag <- 6

# Lag knots
lagknot <- equalknots(0:nlag, 2)

# Lag weather variables (up to 6 months)
lag_mean_temp <- tsModel::Lag(comunas$mean_temp,
  group = comunas$comuna,
  k = 0:nlag
)
lag_mean_temp <- lag_mean_temp[as.Date(paste(comunas$year,
  comunas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_tot_precip <- tsModel::Lag(comunas$tot_precip,
  group = comunas$comuna,
  k = 0:nlag
)
lag_tot_precip <- lag_tot_precip[as.Date(paste(comunas$year,
  comunas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ] # Over june 2013


lag_ndvi <- tsModel::Lag(comunas$ndvi,
  group = comunas$comuna,
  k = 0:nlag
)
lag_ndvi <- lag_ndvi[as.Date(paste(comunas$year,
  comunas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_wet_days <- tsModel::Lag(comunas$wet_days,
  group = comunas$comuna,
  k = 0:nlag
)
lag_wet_days <- lag_wet_days[as.Date(paste(comunas$year,
  comunas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

lag_days_ov_32 <- tsModel::Lag(comunas$days_ov_32,
  group = comunas$comuna,
  k = 0:nlag
)
lag_days_ov_32 <- lag_days_ov_32[as.Date(paste(comunas$year,
  comunas$month,
  "1",
  sep = "-"
)) > "2013-6-30", ]

# Redefining time window for the rest of the data
comunas <- comunas[as.Date(paste(comunas$year,
  comunas$month, "1",
  sep = "-"
)) >
  "2013-6-30", ] # Over jun-2013

## 2.1. Crossbasis matrix ##
basis_mean_temp <- crossbasis(lag_mean_temp,
  argvar = list(fun = "ns", knots = equalknots(comunas$mean_temp, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_mean_temp) <- paste0(
  "basis_mean_temp.",
  colnames(basis_mean_temp)
)

basis_tot_precip <- crossbasis(lag_tot_precip,
  argvar = list(fun = "ns", knots = equalknots(comunas$tot_precip, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_tot_precip) <- paste0(
  "basis_tot_precip.",
  colnames(basis_tot_precip)
)

basis_ndvi <- crossbasis(lag_ndvi,
  argvar = list(fun = "ns", knots = equalknots(comunas$ndvi, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_ndvi) <- paste0(
  "basis_ndvi.",
  colnames(basis_ndvi)
)

basis_wet_days <- crossbasis(lag_wet_days,
  argvar = list(fun = "ns", knots = equalknots(comunas$wet_days, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_wet_days) <- paste0(
  "basis_wet_days.",
  colnames(basis_wet_days)
)
basis_days_ov_32 <- crossbasis(lag_days_ov_32,
  argvar = list(fun = "ns", knots = equalknots(comunas$days_ov_32, 2)),
  arglag = list(fun = "ns", knots = lagknot)
)
colnames(basis_days_ov_32) <- paste0(
  "basis_days_ov_32.",
  colnames(basis_days_ov_32)
)

### 3. INLA ###

## 3.1 General data ##
comunas_inla <- comunas %>%
  summarise(
    comuna = comuna,
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

comunas_inla <- comunas_inla[order(
  comunas_inla$year_index,
  comunas_inla$month,
  comunas_inla$comuna
), ]


### 3.2 Prior definition ##
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

## 3.4 Spatial data import ##

nb_comunas_inla <- inla.read.graph("maps/comunas.adj") # INLA geometries

## 3.5 Formulas ##
# Base formula demographic (random effects)
formula_comunas.0 <- cases ~ 1 + f(month,
  replicate = comuna, model = "rw1",
  cyclic = TRUE, constr = FALSE,
  scale.model = FALSE, hyper = precision.prior, diagonal = 1
) +
  f(comuna,
    model = "bym2", replicate = year_index,
    graph = nb_comunas_inla, adjust.for.con.comp = TRUE,
    scale.model = TRUE, hyper = precision.prior, diagonal = 1
  )

## 3.5.1 Adding fixed effects
formula_comunas_spat <- update(formula_comunas.0, ~ . + density + kids +
  women + higher_ed + gas + sewage + garbage)

formula_comunas_temp <- update(formula_comunas.0, ~ . + basis_mean_temp +
  basis_ndvi + basis_tot_precip +
  basis_wet_days + basis_days_ov_32)

formula_comunas_spte <- update(formula_comunas.0, ~ . + density + kids +
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

comunas_spat <- fit_inla(comunas_inla, formula_comunas_spat)
comunas_temp <- fit_inla(comunas_inla, formula_comunas_temp)
comunas_spte <- fit_inla(comunas_inla, formula_comunas_spte)

## 3.7 Model Selection ##

comunas_compared_df <- data.frame(
  model = c("comunas_spat", "comunas_temp", "comunas_spte"),
  DIC = c(0, 0, 0),
  WAIC = c(0, 0, 0)
)
for (i in c("comunas_spat", "comunas_temp", "comunas_spte")) {
  temp <- get(i)
  comunas_compared_df[comunas_compared_df$model == i, "DIC"] <- temp$dic$dic
  comunas_compared_df[comunas_compared_df$model == i, "WAIC"] <- temp$waic$waic
}
comunas_compared_df$model <- c("Spatial", "Temporal", "Spatiotemporal")

# Best model is the one containing all predictors (remove others)
rm(temp)
rm(comunas_spat)
rm(comunas_temp)

### 4. Results report ###

## 4.1 Temporal variagbles (Lagged) ##
coeficients <- comunas_spte$summary.fixed$mean
covariance <- comunas_spte$misc$lincomb.derived.covariance.matrix

# Mean temperature
col_mean_temp <- grep("basis_mean_temp", comunas_spte$names.fixed)
predicted_mean_temp <- crosspred(basis_mean_temp,
  coef = coeficients[col_mean_temp],
  vcov = covariance[col_mean_temp, col_mean_temp],
  model.link = "log", bylag = 0.25, cen = round(mean(comunas$mean_temp), 0)
)
png("figures/inla_lag/mean_temp_comunas.png",
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
      main = "Comunas", ylab = "Mean Temp. [°C]",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Total precipitation
col_tot_precip <- grep("basis_tot_precip", comunas_spte$names.fixed)
predicted_tot_precip <- crosspred(basis_tot_precip,
  coef = coeficients[col_tot_precip],
  vcov = covariance[col_tot_precip, col_tot_precip],
  model.link = "log", bylag = 0.25, cen = round(mean(comunas$tot_precip), 0)
)
png("figures/inla_lag/tot_precip_comunas.png",
  width = 2200, height = 1700,
  res = 300
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
      labels = c(0, 20, 40, 60, 80, 100, 120),
      at = c(0, 20, 40, 60, 80, 100, 120)
    )
  }, key.axes = axis(4,
    cex.axis = 2.5, family =
      windowsFont(family = "Times New Roman")
  ),
  plot.title = {
    title(
      main = "", ylab = "Total Prec. [mm/month]",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# NDVI
col_ndvi <- grep("basis_ndvi", comunas_spte$names.fixed)
predicted_ndvi <- crosspred(basis_ndvi,
  coef = coeficients[col_ndvi], vcov = covariance[col_ndvi, col_ndvi],
  model.link = "log", bylag = 0.25, cen = round(mean(comunas$ndvi), 0)
)
png("figures/inla_lag/ndvi_comunas.png", width = 2200, height = 1700, res = 300)
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
      main = "", ylab = "NDVI",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Wet days
col_wet_days <- grep("basis_wet_days", comunas_spte$names.fixed)
predicted_wet_days <- crosspred(basis_wet_days,
  coef = coeficients[col_wet_days],
  vcov = covariance[col_wet_days, col_wet_days],
  model.link = "log", bylag = 0.25, cen = round(mean(comunas$wet_days), 0)
)
png("figures/inla_lag/wet_days_comunas.png",
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
      main = "", ylab = "Wet Days",
      xlab = "", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

# Days Over 32
col_days_ov_32 <- grep("basis_days_ov_32", comunas_spte$names.fixed)
predicted_days_ov_32 <- crosspred(basis_days_ov_32,
  coef = coeficients[col_days_ov_32],
  vcov = covariance[col_days_ov_32, col_days_ov_32],
  model.link = "log", bylag = 0.25, cen = round(mean(comunas$days_ov_32), 0)
)
png("figures/inla_lag/days_ov_32_comunas.png",
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
      main = "", ylab = "Days over 32°C",
      xlab = "Lag [Months]", cex.lab = 3, cex.main = 3, font.main = 1,
      family = windowsFont(family = "Times New Roman")
    )
  }
)
dev.off()

## 4.2 Spatial variables (non-lagged) ##
df_results <- comunas_spte$summary.fixed[1:8, 3:5]
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
comunas_spatial <- df_results

### 5. Save results ###
comunas_results <- list(
  spatial = comunas_spatial,
  compared = comunas_compared_df,
  summary.fitted = comunas_spte$summary.fitted.values
)
saveRDS(comunas_results, "data/workspace/comunas_results.rds")
save.image(file = "data/workspace/comunas_inla.RData")
