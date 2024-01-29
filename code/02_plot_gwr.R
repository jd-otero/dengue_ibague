## ---------- Plot GWR results ----------##

### 1. Set up ###

# Libraries
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)
library(extrafont)
library(stringr)

significance <- readRDS("data/workspace/gwr_significance.rds")

# Import data
comunas_sf <- st_read("data/comunas.shp") %>% select(
  comuna,
  latitude,
  longitude,
  geometry
)
sectores_sf <- st_read("data/sectores.shp") %>% select(
  sector,
  latitude,
  longitude,
  geometry
)
secciones_sf <- st_read("data/secciones.shp") %>% select(
  seccion,
  latitude,
  longitude,
  geometry
)
manzanas_sf <- st_read("data/manzanas.shp") %>% select(
  manzana,
  latitude,
  longitude,
  geometry
)

comunas_sf <- comunas_sf[!duplicated(comunas_sf$comuna), ]
sectores_sf <- sectores_sf[!duplicated(sectores_sf$sector), ]
secciones_sf <- secciones_sf[!duplicated(secciones_sf$seccion), ]
manzanas_sf <- manzanas_sf[!duplicated(manzanas_sf$manzana), ]

comunas <- merge(significance$comunas, comunas_sf) %>% st_as_sf()
sectores <- merge(significance$sectores, sectores_sf) %>% st_as_sf()
secciones <- merge(significance$secciones, secciones_sf) %>% st_as_sf()
manzanas <- merge(significance$manzanas, manzanas_sf) %>% st_as_sf()

### 2. Plot ###

# Selected covariates
vars <- c(
  "density", "sewage", "gas", "garbage", "higher_ed",
  "kids", "women"
)

var_names <- list(
  density = expression("Density"),
  sewage = expression("Sewage"),
  gas = expression("Gas"),
  garbage = expression("Garbage"),
  higher_ed = expression("Higher Ed."),
  kids = expression("Kids"),
  women = expression("Women")
)

# Define colors
colors <- list(
  density = c("TRUE" = "#49997c", "FALSE" = "#e6e6e6"),
  sewage = c("TRUE" = "#027ab0", "FALSE" = "#e6e6e6"),
  gas = c("TRUE" = "#ae3918", "FALSE" = "#e6e6e6"),
  garbage = c("TRUE" = "#4f5b67", "FALSE" = "#e6e6e6"),
  higher_ed = c("TRUE" = "#d19c2f", "FALSE" = "#e6e6e6"),
  kids = c("TRUE" = "#6b4a74", "FALSE" = "#e6e6e6"),
  women = c("TRUE" = "#556b2f", "FALSE" = "#e6e6e6")
)

# List for all plots
plot_list <- list()

for (var in vars) {
  for (level in c("comunas", "sectores", "secciones", "manzanas")) {
    data <- get(level)
    color <- colors[var][[1]]
    var_plot <- ggplot(data = data) +
      geom_sf(mapping = aes(fill = .data[[var]]), color = NA) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.2
        )
      ) +
      scale_fill_manual(
        values = color,
        labels = c(TRUE, FALSE)
      )
    if (var == "density") {
      var_plot <- var_plot + ggtitle(str_to_title(level)) +
        theme(plot.title = element_text(
          size = 10, hjust = 0.5,
          family = windowsFont(family = "Times New Roman")
        ))
    }
    if (level == "comunas") {
      var_plot <- var_plot + ylab(var_names[[var]]) +
        theme(axis.title.y = element_text(
          size = 10, vjust = 0.5, hjust = 0.5,
          family = windowsFont(family = "Times New Roman")
        ))
    }
    plot_list[[paste("p", var, level, sep = "_")]] <- var_plot
  }
}

# Plot grid
plot_grid <- plot_grid(
  plotlist = plot_list, ncol = 4,
  rel_heights = c(1.24, 1, 1, 1, 1, 1, 1),
  rel_widths = c(1.08, 1, 1, 1)
)
ggsave2("figures/gwr.jpg", plot_grid, width = 6.5, height = 8)
