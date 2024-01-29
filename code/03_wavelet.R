#### ---------- Wavelet Coherenec ---------- ####

### 1. Set up ###

# Libraries
library(sf)
library(dplyr)
library(biwavelet)
library(cowplot)
library(extrafont)
library(grid)
library(gridExtra)
library(png)

# Import data
comunas <- st_read("data/comunas.shp") %>% st_drop_geometry()

### 2. Data preparation ###

# Group time data
data_temp <- comunas %>%
  group_by(year, month) %>%
  summarise(
    mean_temp = mean(mean_temp),
    tot_precip = sum(tot_precip),
    ndvi = mean(ndvi),
    wet_days = mean(wet_days),
    days_ov_32 = mean(days_ov_32),
    cases = sum(cases)
  ) %>%
  mutate(date = as.Date(paste(year,
    month,
    "1",
    sep = "-"
  )))

### 3. Wavelet Coherence ###

data_temp <- data_temp[order(data_temp$date), ]
data_temp$month <- seq(1, nrow(data_temp))

# Cases
cases_wc <- as.matrix(data_temp[, c("month", "cases")])

# Mean temperature
mean_temp_wc <- as.matrix(data_temp[, c("month", "mean_temp")])
wtc_temp <- wtc(cases_wc, mean_temp_wc, nrands = 100, max.scale = 30)

png("figures/wavelets/mean_temp.png", width = 650, height = 550, res = 300)
tryCatch(
  {
    par(mar = c(2.2, 2, 1.3, 3.2), mgp = c(1.2, 0.4, 0))
    plot(wtc_temp,
      plot.phase = TRUE, lty.coi = 0.25, col.coi = "grey", lwd.coi = 0.8,
      lwd.sig = 0.8, arrow.lwd = 0.01, arrow.len = 0.1, ylab = "Period",
      xlab = "Month", plot.cb = TRUE, main = "Mean Temperature", font.main = 1,
      legend.loc = c(0.72, 0.78, 0.25, 0.85),
      family = windowsFont(family = "Times New Roman"),
      cex.main = 0.81,
      cex.lab = 0.81, cex.axis = 0.81,
    )
  },
  error = function(e) {
    # NA
  }
)
dev.off()

# Total Precipitation
tot_precip_wc <- as.matrix(data_temp[, c("month", "tot_precip")])
wtc_precip <- wtc(cases_wc, tot_precip_wc, nrands = 100, max.scale = 30)

png("figures/wavelets/tot_precip.png", width = 650, height = 550, res = 300)
tryCatch(
  {
    par(mar = c(2.2, 2, 1.3, 3.2), mgp = c(1.2, 0.4, 0))
    plot(wtc_precip,
      plot.phase = TRUE, lty.coi = 0.25, col.coi = "grey", lwd.coi = 0.8,
      lwd.sig = 0.8, arrow.lwd = 0.01, arrow.len = 0.1, ylab = "Period",
      xlab = "Month", plot.cb = TRUE, main = "Total Precipitation",
      font.main = 1, legend.loc = c(0.72, 0.78, 0.25, 0.85),
      family = windowsFont(family = "Times New Roman"), cex.main = 0.81,
      cex.lab = 0.81, cex.axis = 0.81
    )
  },
  error = function(e) {
    # NA
  }
)
dev.off()

# NDVI
ndvi_wc <- as.matrix(data_temp[, c("month", "ndvi")])
wtc_ndvi <- wtc(cases_wc, ndvi_wc, nrands = 100, max.scale = 30)

png("figures/wavelets/ndvi.png", width = 650, height = 550, res = 300)
tryCatch(
  {
    par(mar = c(2.2, 2, 1.3, 3.2), mgp = c(1.2, 0.4, 0))
    plot(wtc_ndvi,
      plot.phase = TRUE, lty.coi = 0.25, col.coi = "grey", lwd.coi = 0.8,
      lwd.sig = 0.8, arrow.lwd = 0.01, arrow.len = 0.1, ylab = "Period",
      xlab = "Month", plot.cb = TRUE, main = "NDVI", font.main = 1,
      legend.loc = c(0.72, 0.78, 0.25, 0.85),
      family = windowsFont(family = "Times New Roman"), cex.main = 0.81,
      cex.lab = 0.81, cex.axis = 0.81
    )
  },
  error = function(e) {
    # NA
  }
)
dev.off()

# Wet Days
wet_days_wc <- as.matrix(data_temp[, c("month", "wet_days")])
wtc_wet_days <- wtc(cases_wc, wet_days_wc, nrands = 100, max.scale = 30)

png("figures/wavelets/wet_days.png", width = 650, height = 550, res = 300)
tryCatch(
  {
    par(mar = c(2.2, 2, 1.3, 3.2), mgp = c(1.2, 0.4, 0))
    plot(wtc_wet_days,
      plot.phase = TRUE, lty.coi = 0.25, col.coi = "grey", lwd.coi = 0.8,
      lwd.sig = 0.8, arrow.lwd = 0.01, arrow.len = 0.1, ylab = "Period",
      xlab = "Month", plot.cb = TRUE, main = "Wet Days", font.main = 1,
      legend.loc = c(0.72, 0.78, 0.25, 0.85),
      family = windowsFont(family = "Times New Roman"), cex.main = 0.81,
      cex.lab = 0.81, cex.axis = 0.81
    )
  },
  error = function(e) {
    # NA
  }
)
dev.off()

# Days ov 32
days_ov_32_wc <- as.matrix(data_temp[, c("month", "days_ov_32")])
wtc_days_ov_32 <- wtc(cases_wc, days_ov_32_wc, nrands = 100, max.scale = 30)

png("figures/wavelets/days_ov_32.png", width = 650, height = 550, res = 300)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
tryCatch(
  {
    par(mar = c(2.2, 2, 1.3, 3.2), mgp = c(1.2, 0.4, 0))
    plot(wtc_days_ov_32,
      plot.phase = TRUE, lty.coi = 0.25, col.coi = "grey", lwd.coi = 0.8,
      lwd.sig = 0.8, arrow.lwd = 0.01, arrow.len = 0.1, ylab = "Period",
      xlab = "Month", plot.cb = TRUE, main = "Days Over 32 C", font.main = 1,
      legend.loc = c(0.72, 0.78, 0.25, 0.85),
      family = windowsFont(family = "Times New Roman"), cex.main = 0.81,
      cex.lab = 0.81, cex.axis = 0.81
    )
  },
  error = function(e) {
    # NA
  }
)
dev.off()

## 4. Plot ###
plot_mean_temp <- readPNG("figures/wavelets/mean_temp.png")
plot_tot_precip <- readPNG("figures/wavelets/tot_precip.png")
plot_ndvi <- readPNG("figures/wavelets/ndvi.png")
plot_wet_days <- readPNG("figures/wavelets/wet_days.png")
plot_days_ov_32 <- readPNG("figures/wavelets/days_ov_32.png")

# As list
rl <- list(
  plot_mean_temp, plot_tot_precip, plot_ndvi, plot_wet_days,
  plot_days_ov_32
)
gl <- lapply(rl, rasterGrob)
par(mfrow = c(1, 1))
mat <- rbind(
  c(1, 1, 2, 2, 3, 3),
  c(NA, 4, 4, 5, 5, NA)
)
png("figures/wavelets.png", width = 1950, height = 1100, res = 300)
grid.arrange(grobs = gl, layout_matrix = mat)
dev.off()
