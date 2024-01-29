## ---------- Plot Dengue Cases ----------##

### 1. Set up ###

# Libraries
library(lubridate)
library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(extrafont)
library(dichromat)

# Import and transform data
ibague_shp <- st_read("data/raw/ibague_shp/ibague.shp")
ibague_shp <- st_transform(ibague_shp, 4686)

cases <- read_excel("data/raw/CasosGeorreferenciados.xlsx",
  sheet = "Registros Depurados"
)
cases <- cases %>%
  mutate(
    cx = as.numeric(CX),
    cy = as.numeric(CY),
    date = dmy(`FECHA DE INICIO DE SINTOM`),
    year = year(date),
    month = month(date),
    case = 1
  ) %>%
  filter(
    year >= 2013 & year <= 2018,
    comuna != "R"
  ) %>%
  select(date, year, month, cx, cy, case)
cases <- cases[!is.na(cases$cx) & !is.na(cases$cy), ] # Delete cases without coordinates
cases <- st_as_sf(cases, coords = c("cx", "cy"), remove = FALSE)
cases <- st_set_crs(cases, 4686)

filtered_cases <- st_intersects(cases, ibague_shp)
cases <- cases[which(lengths(filtered_cases) != 0), ]

### 2. Line plot ###
cases_grouped <- cases %>%
  mutate(yearmon = as.Date(paste(year, month, "1", sep = "-"))) %>%
  group_by(yearmon) %>%
  summarise(cases = sum(case)) %>%
  st_drop_geometry()

line_plot <- ggplot(data = cases_grouped) +
  geom_line(mapping = aes(x = yearmon, y = cases), colour = "#ae3918") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(
      size = 10, hjust = 0.5, vjust = -1,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.title.y = element_text(
      size = 10, hjust = 0.5, vjust = 2,
      family = windowsFont(family = "Times New Roman")
    ),
    axis.text.y = element_text(
      size = 10, hjust = 0.5,
      family = windowsFont(family = "Times New Roman"), color = "black"
    ),
    axis.text.x = element_text(
      size = 10, hjust = 1, angle = 45,
      family = windowsFont(family = "Times New Roman"), color = "black"
    )
  ) +
  scale_x_date(breaks = scales::breaks_pretty(24)) +
  xlab("Date") +
  ylab("Cases")


### 3. Heatmaps ###

pal <- colorRampPalette(c("#FFFFC8", "#F8C964", "#F17D00", "#CF2C00", "#7D0025"))
for (y in 2013:2018) {
  count_year <- cases %>% filter(year == y)
  assign(paste("count", y, sep = "_"), ggplot() +
    geom_sf(data = ibague_shp, fill = "white", color = "black") +
    stat_density2d(
      data = count_year,
      aes(x = cx, y = cy, fill = ..level..), h = 0.008, alpha = 0.4,
      geom = "polygon"
    ) +
    scale_fill_gradientn(
      colours = pal(30),
      breaks = c(0, 1800), labels = c("Low", "High"),
      limits = c(0, 1800)
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(
        size = 10, hjust = 0.5,
        family = windowsFont(family = "Times New Roman")
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 10, hjust = 0.5,
        family = windowsFont(family = "Times New Roman")
      )
    ) +
    ggtitle(paste(y)))
}

heatmaps <- ggarrange(count_2013, count_2014, count_2015,
  count_2016, count_2017, count_2018,
  common.legend = TRUE,
  legend = "bottom",
  nrow = 2,
  ncol = 3
)

### 4. Grid plot ###
plot_grid <- plot_grid(
  plotlist = list(line_plot, heatmaps), ncol = 1,
  labels = "AUTO", rel_heights = c(0.4, 0.6)
)
ggsave2("figures/Fig1.jpg", plot_grid, width = 6.5, height = 5.7)
