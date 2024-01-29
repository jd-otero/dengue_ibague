## ---------- Geographical Weighted Regression ---------- ##

### 1. Set up ###

# Libraries
library(GWmodel)
library(sf)
library(dplyr)
library(sp)

# Import data
comunas <- st_read("data/comunas.shp") %>% st_drop_geometry()
sectores <- st_read("data/sectores.shp") %>% st_drop_geometry()
secciones <- st_read("data/secciones.shp") %>% st_drop_geometry()
manzanas <- st_read("data/manzanas.shp") %>% st_drop_geometry()

### 2. Data preparation ###

# Transform data (group all)
group_data <- function(df) {
  df_base <- df[!duplicated(df[, 1]), ] %>%
    subset(select = -c(
      mean_temp, days_ov_32, tot_precip,
      wet_days, dry_days, ndvi, elevation, year,
      month, cases
    )) %>%
    group_by_at(1) %>%
    summarise(
      latitude = latitude,
      longitude = longitude,
      density = density,
      population = population,
      houses = houses,
      homes = homes / houses,
      electric = electric / houses,
      aqueduct = aqueduct / houses,
      sewage = sewage / houses,
      gas = gas / houses,
      garbage = garbage / houses,
      internet = internet / houses,
      strata = strata,
      higher_ed = higher_ed / population,
      middle_ed = middle_ed / population,
      no_ed = no_ed / population,
      kids = kids / population,
      adults = adults / population,
      men = men / population,
      women = women / population
    ) %>%
    replace(is.na(.), 0)
  df_cases <- df %>%
    group_by_at(1) %>%
    summarise(cases = sum(cases)) %>%
    mutate(cases = cases / df_base$population) %>%
    replace(is.na(.), 0)
  df_grouped <- merge(x = df_base, y = df_cases)
  return(df_grouped)
}

# Get groupped and normalized data
comunas_g <- group_data(comunas)
sectores_g <- group_data(sectores)
secciones_g <- group_data(secciones)
manzanas_g <- group_data(manzanas)

### 3. Geographically Wwighted Regression ###

# Change sf to sp and fit gwr
fit_gwr <- function(df, formula) {
  coords <- df[, c("latitude", "longitude")] # coordinates
  data <- df %>% select(-c(latitude, longitude))
  crs <- CRS("+init=epsg:4686")

  df_sp <- SpatialPointsDataFrame(
    coords = coords,
    data = data,
    proj4string = crs
  )

  gwr <- GWmodel::gwr.basic(
    formula = formula,
    data = df_sp,
    bw = 5,
    kernel = "gaussian",
    longlat = T,
    adaptive = F
  )
  return(gwr)
}

# Formulas
formula_soc <- as.formula(cases ~ density + homes + electric + aqueduct +
  sewage + gas + garbage + internet + strata)

formula_dem <- as.formula(cases ~ higher_ed + no_ed + kids + women)

# Fitting gwr
gwr_comunas_dem <- fit_gwr(comunas_g, formula_dem)
gwr_sectores_dem <- fit_gwr(sectores_g, formula_dem)
gwr_secciones_dem <- fit_gwr(secciones_g, formula_dem)
gwr_manzanas_dem <- fit_gwr(manzanas_g, formula_dem)

gwr_comunas_soc <- fit_gwr(comunas_g, formula_soc)
gwr_sectores_soc <- fit_gwr(sectores_g, formula_soc)
gwr_secciones_soc <- fit_gwr(secciones_g, formula_soc)
gwr_manzanas_soc <- fit_gwr(manzanas_g, formula_soc)

### 4. Significant variables ###
sign_var <- c()

extract_sign_var <- function(summary) {
  vars <- names(summary)[
    which(summary < 0.05)
  ]
  return(vars)
}
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_comunas_soc$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_comunas_dem$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_sectores_soc$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_sectores_dem$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_secciones_soc$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_secciones_dem$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_manzanas_soc$lm
)$coefficients[, 4]))
sign_var <- append(sign_var, extract_sign_var(summary(
  gwr_manzanas_dem$lm
)$coefficients[, 4]))

sign_var <- sign_var[!duplicated(sign_var)]

# Extract significance of each variable and level
significance <- function(gwr) {
  names <- names(gwr$lm$coefficients)[-1]
  df_significance <- data.frame(matrix(rep(0, length(names) * nrow(gwr$SDF)),
    ncol = length(names)
  ))
  names(df_significance) <- names
  for (column in names) {
    t <- abs(gwr$SDF[, paste0(column)][[1]]) - 2 *
      gwr$SDF[paste(column, "SE", sep = "_")][[1]]
    sig <- t > 0
    df_significance[, column] <- sig
  }
  return(df_significance)
}

sig_comunas_dem <- significance(gwr_comunas_dem)
sig_comunas_soc <- significance(gwr_comunas_soc)
significance_comunas <- cbind(sig_comunas_dem, sig_comunas_soc)

sig_sectores_dem <- significance(gwr_sectores_dem)
sig_sectores_soc <- significance(gwr_sectores_soc)
significance_sectores <- cbind(sig_sectores_dem, sig_sectores_soc)

sig_secciones_dem <- significance(gwr_secciones_dem)
sig_secciones_soc <- significance(gwr_secciones_soc)
significance_secciones <- cbind(sig_secciones_dem, sig_secciones_soc)

sig_manzanas_dem <- significance(gwr_manzanas_dem)
sig_manzanas_soc <- significance(gwr_manzanas_soc)
significance_manzanas <- cbind(sig_manzanas_dem, sig_manzanas_soc)

# Define variables
significance_global <- data.frame(
  matrix(rep(0, ncol(significance_comunas) * 4),
    ncol = 4
  ),
  row.names = names(significance_comunas)
)
names(significance_global) <- c("comunas", "sectores", "secciones", "manzanas")
significance_global$comunas <- colSums(significance_comunas) /
  nrow(comunas_g)
significance_global$sectores <- colSums(significance_sectores) /
  nrow(sectores_g)
significance_global$secciones <- colSums(significance_secciones) /
  nrow(secciones_g)
significance_global$manzanas <- colSums(significance_manzanas) /
  nrow(manzanas_g)
significance_global$total <- apply(significance_global,
  function(i) ifelse(any(i > 0.5), T, F),
  MARGIN = 1
)

# Group significance
significance <- list(
  comunas = cbind(comuna = comunas_g$comuna, significance_comunas),
  sectores = cbind(sector = sectores_g$sector, significance_sectores),
  secciones = cbind(seccion = secciones_g$seccion, significance_secciones),
  manzanas = cbind(manzana = manzanas_g$manzana, significance_manzanas),
  global = significance_global
)

### 5. Save workspace ###
saveRDS(significance, "data/workspace/gwr_significance.rds")
save.image("data/workspace/gwr.RData")
