## ---------- Plot non lagged variables ----------##

### 1. Set up ###

# Libraries
library(ggplot2)
library(dplyr)
library(extrafont)
library(png)
library(sf)
library(lubridate)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(forcats)

# Data from INLA models
comunas_results <- readRDS("data/workspace/comunas_results.rds")
sectores_results <- readRDS("data/workspace/sectores_results.rds")
secciones_results <- readRDS("data/workspace/secciones_results.rds")
manzanas_results <- readRDS("data/workspace/manzanas_results.rds")

### 2. Spatial effects ###

# Results data
results <- data.frame(var = rep(c(
  "intercept", "density", "sewage", "gas",
  "garbage", "higher_ed", "kids", "women"
), 4))
results$level <- c(
  rep("manzanas", 8), rep("secciones", 8), rep("sectores", 8),
  rep("comunas", 8)
)
results$median <- NaN
results$upper <- NaN
results$lower <- NaN
results$rr <- NaN

# Extract needed information
for (i in 1:nrow(results)) {
  if (results[i, 2] == "manzanas") {
    if (results[i, 1] %in% rownames(manzanas_results$spatial)) {
      row_manzanas <- which(rownames(manzanas_results$spatial) == results[
        i,
        1
      ])
      results[i, "median"] <- manzanas_results$spatial[
        row_manzanas,
        "0.5quant"
      ]
      results[i, "upper"] <- manzanas_results$spatial[
        row_manzanas,
        "0.975quant"
      ]
      results[i, "lower"] <- manzanas_results$spatial[
        row_manzanas,
        "0.025quant"
      ]
      results[i, "rr"] <- manzanas_results$spatial[
        row_manzanas,
        "rr"
      ]
    }
  }
  if (results[i, 2] == "secciones") {
    if (results[i, 1] %in% rownames(secciones_results$spatial)) {
      row_secciones <- which(rownames(secciones_results$spatial) == results[
        i,
        1
      ])
      results[i, "median"] <- secciones_results$spatial[
        row_secciones,
        "0.5quant"
      ]
      results[i, "upper"] <- secciones_results$spatial[
        row_secciones,
        "0.975quant"
      ]
      results[i, "lower"] <- secciones_results$spatial[
        row_secciones,
        "0.025quant"
      ]
      results[i, "rr"] <- secciones_results$spatial[
        row_secciones,
        "rr"
      ]
    }
  }
  if (results[i, 2] == "sectores") {
    if (results[i, 1] %in% rownames(sectores_results$spatial)) {
      row_sectores <- which(rownames(sectores_results$spatial) == results[
        i,
        1
      ])
      results[i, "median"] <- sectores_results$spatial[
        row_sectores,
        "0.5quant"
      ]
      results[i, "upper"] <- sectores_results$spatial[
        row_sectores,
        "0.975quant"
      ]
      results[i, "lower"] <- sectores_results$spatial[
        row_sectores,
        "0.025quant"
      ]
      results[i, "rr"] <- sectores_results$spatial[row_sectores, "rr"]
    }
  }
  if (results[i, 2] == "comunas") {
    if (results[i, 1] %in% rownames(comunas_results$spatial)) {
      row_comunas <- which(rownames(comunas_results$spatial) == results[
        i,
        1
      ])
      results[i, "median"] <- comunas_results$spatial[
        row_comunas,
        "0.5quant"
      ]
      results[i, "upper"] <- comunas_results$spatial[
        row_comunas,
        "0.975quant"
      ]
      results[i, "lower"] <- comunas_results$spatial[
        row_comunas,
        "0.025quant"
      ]
      results[i, "rr"] <- comunas_results$spatial[
        row_comunas,
        "rr"
      ]
    }
  }
}

results$level <- factor(results$level, levels = c(
  "comunas", "sectores",
  "secciones", "manzanas"
))
results <- results %>% mutate(var = fct_relevel(var, c(
  "intercept", "density", "sewage", "gas",
  "garbage", "higher_ed", "kids", "women"
)))
results_no_inter <- results[which(results$var != "intercept"), ]

# Plot
coeff_plot <- ggplot(data = results_no_inter) +
  geom_pointrange(aes(
    y = median, x = var, ymin = lower, ymax = upper,
    color = level,
  ), size = 0.1, position = position_dodge(0.3)) +
  scale_x_discrete(labels = c(
    "Density", "Sewage", "Gas", "Garbage",
    "Higher Ed.", "Kids", "Women"
  )) +
  scale_color_manual(
    values = c(
      "comunas" = "#8A9045", "sectores" = "#155F83",
      "secciones" = "#B24745", "manzanas" = "#DF8F44"
    ),
    labels = c("Comunas", "Sectores", "Secciones", "Manzanas")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.x = element_text(
      size = 10, hjust = 0.5, vjust = -1,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.text.x = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    axis.text.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    legend.text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman"), color = "black"
    )
  ) +
  ylab("Coefficient") +
  xlab("Variable")
ggsave2("figures/coefficients.jpg", coeff_plot, width = 6.5, height = 3.7)

### 3. Temporal effects ###
# Extracted from previous individual plots
plt_names <- c()
i <- 0
rl <- list()
for (var in c("mean_temp", "tot_precip", "ndvi", "wet_days", "days_ov_32")) {
  for (level in c("comunas", "sectores", "secciones", "manzanas")) {
    i <- i + 1
    plt_names <- c(plt_names, paste(var, level, sep = "_"))
    temp_plot <- readPNG(paste0(
      "figures/inla_lag/",
      paste(var, level, sep = "_"),
      ".png"
    ))
    rl[[i]] <- temp_plot
  }
}
gl <- lapply(rl, rasterGrob)
png("figures/lag.png", width = 1950, height = 1880, res = 300)
grid.arrange(grobs = gl, ncol = 4)
dev.off()

### 4. Comparison ###

# Comunas
comunas <- st_read("data/comunas.shp") %>% st_drop_geometry()
comunas <- comunas %>% filter(population > 0)
comunas <- comunas[order(
  comunas$year,
  comunas$month,
  comunas$comuna
), ]
comunas <- comunas %>%
  mutate(date = ymd(paste(year, month, "1", sep = "-"))) %>%
  filter(date > "2013-06-01") %>%
  select(date, cases)
comunas <- cbind(
  comunas, comunas_results$summary.fitted$mean,
  comunas_results$summary.fitted$sd
)
names(comunas) <- c("date", "cases", "estimated", "sd")
comunas <- comunas %>%
  group_by(date) %>%
  summarise(
    cases = sum(cases),
    estimated = sum(estimated),
    sd = sqrt(sum(sd^2))
  )

# Sectores
sectores <- st_read("data/sectores.shp") %>% st_drop_geometry()
sectores <- sectores %>% filter(population > 0)
sectores <- sectores[order(
  sectores$year,
  sectores$month,
  sectores$sector
), ]
sectores <- sectores %>%
  mutate(date = ymd(paste(year, month, "1", sep = "-"))) %>%
  filter(date > "2013-06-01") %>%
  select(date, cases)
sectores <- cbind(
  sectores, sectores_results$summary.fitted$mean,
  sectores_results$summary.fitted$sd
)
names(sectores) <- c("date", "cases", "estimated", "sd")
sectores <- sectores %>%
  group_by(date) %>%
  summarise(
    cases = sum(cases),
    estimated = sum(estimated),
    sd = sqrt(sum(sd^2))
  )

# Secciones
secciones <- st_read("data/secciones.shp") %>% st_drop_geometry()
secciones <- secciones %>% filter(population > 0)
secciones <- secciones[order(
  secciones$year,
  secciones$month,
  secciones$seccion
), ]
secciones <- secciones %>%
  mutate(date = ymd(paste(year, month, "1", sep = "-"))) %>%
  filter(date > "2013-06-01") %>%
  select(date, cases)
secciones <- cbind(
  secciones, secciones_results$summary.fitted$mean,
  secciones_results$summary.fitted$sd
)
names(secciones) <- c("date", "cases", "estimated", "sd")
secciones <- secciones %>%
  group_by(date) %>%
  summarise(
    cases = sum(cases),
    estimated = sum(estimated),
    sd = sqrt(sum(sd^2))
  )

# Manzanas
manzanas <- st_read("data/manzanas.shp") %>% st_drop_geometry()
manzanas <- manzanas %>% filter(population > 0)
manzanas <- manzanas[order(
  manzanas$year,
  manzanas$month,
  manzanas$manzana
), ]
manzanas <- manzanas %>%
  mutate(date = ymd(paste(year, month, "1", sep = "-"))) %>%
  filter(date > "2013-06-01") %>%
  select(date, cases)
manzanas <- cbind(
  manzanas, manzanas_results$summary.fitted$mean,
  manzanas_results$summary.fitted$sd
)
names(manzanas) <- c("date", "cases", "estimated", "sd")
manzanas <- manzanas %>%
  group_by(date) %>%
  summarise(
    cases = sum(cases),
    estimated = sum(estimated),
    sd = sqrt(sum(sd^2))
  )

all <- cbind(
  comunas, sectores$estimated, sectores$sd,
  secciones$estimated, secciones$sd, manzanas$estimated, manzanas$sd
)
names(all) <- c(
  "date", "cases", "comunas", "comunas_sd", "sectores",
  "sectores_sd", "secciones", "secciones_sd", "manzanas",
  "manzanas_sd"
)


fit_comunas <- ggplot(data = all) +
  geom_ribbon(aes(
    x = date, ymin = comunas - 2 * comunas_sd,
    ymax = comunas + 2 * comunas_sd
  ), fill = "#FBD2D0") +
  geom_line(aes(x = date, y = comunas, colour = "Fitted"), lwd = 0.5) +
  geom_line(aes(x = date, y = cases, colour = "Observed"), lwd = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.text.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    axis.text.x = element_blank(),
    plot.title = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    legend.text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    )
  ) +
  scale_x_date(breaks = scales::breaks_pretty(12)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_color_manual(name = "", values = c(
    "Observed" = "#1B2D2A",
    "Fitted" = "#FB4B4E"
  )) +
  ylab("Cases") +
  ggtitle("Comunas")

# Plot sectores
fit_sectores <- ggplot(data = all) +
  geom_ribbon(aes(
    x = date, ymin = sectores - 2 * sectores_sd,
    ymax = sectores + 2 * sectores_sd
  ), fill = "#FBD2D0") +
  geom_line(aes(x = date, y = sectores, colour = "Fitted"), lwd = 0.5) +
  geom_line(aes(x = date, y = cases, colour = "Observed"), lwd = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    )
  ) +
  scale_x_date(breaks = scales::breaks_pretty(12)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_color_manual(name = "", values = c(
    "Observed" = "#1B2D2A",
    "Fitted" = "#FB4B4E"
  )) +
  ggtitle("Sectores")

# Plot secciones
fit_secciones <- ggplot(data = all) +
  geom_ribbon(aes(
    x = date, ymin = secciones - 2 * secciones_sd,
    ymax = secciones + 2 * secciones_sd
  ), fill = "#FBD2D0") +
  geom_line(aes(x = date, y = secciones, colour = "Fitted"), lwd = 0.5) +
  geom_line(aes(x = date, y = cases, colour = "Observed"), lwd = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.x = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.text.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    axis.text.x = element_text(
      size = 10, hjust = 1, angle = 45,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    plot.title = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    )
  ) +
  scale_x_date(breaks = scales::breaks_pretty(12)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_color_manual(name = "", values = c(
    "Observed" = "#1B2D2A",
    "Fitted" = "#FB4B4E"
  )) +
  ylab("Cases") +
  xlab("Date") +
  ggtitle("Secciones")

# Plot manzanas
fit_manzanas <- ggplot(data = all) +
  geom_ribbon(aes(
    x = date, ymin = manzanas - 2 * manzanas_sd,
    ymax = manzanas + 2 * manzanas_sd
  ), fill = "#FBD2D0") +
  geom_line(aes(x = date, y = manzanas, colour = "Fitted"), lwd = 0.4) +
  geom_line(aes(x = date, y = cases, colour = "Observed"), lwd = 0.4) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    text = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.x = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      size = 10, hjust = 1, angle = 45,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    axis.text.y = element_blank(),
    plot.title = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman")
    )
  ) +
  scale_x_date(breaks = scales::breaks_pretty(12)) +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_color_manual(name = "", values = c(
    "Observed" = "#1B2D2A",
    "Fitted" = "#FB4B4E"
  )) +
  ylab("Cases") +
  xlab("Date") +
  ggtitle("Manzanas")

# Plot grid
legend <- ggpubr::get_legend(fit_comunas)
png("figures/comparison.png", width = 1950, height = 1500, res = 300)
ggarrange(
  ggarrange(fit_comunas + theme(legend.position = "none"),
    fit_sectores,
    fit_secciones,
    fit_manzanas,
    widths = c(1.1, 1), heights = c(1, 1.28)
  ),
  legend,
  nrow = 2, heights = c(1, 0.06)
)
dev.off()

### 5. RMSE ###
rmse_comunas <- sqrt((1 / nrow(all)) * sum((all$cases - all$comunas)^2))
rmse_sectores <- sqrt((1 / nrow(all)) * sum((all$cases - all$sectores)^2))
rmse_secciones <- sqrt((1 / nrow(all)) * sum((all$cases - all$secciones)^2))
rmse_manzanas <- sqrt((1 / nrow(all)) * sum((all$cases - all$manzanas)^2))
