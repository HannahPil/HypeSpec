#!/usr/bin/env Rscript
# =============================================================================
# Figure: PC1 vs PC2 of plot-level reflectance, coloured by scan date
# In : output/pca_scores.csv, output/pca_summary.csv, output/bz_plot_meta.csv
# Out: output/pca_pc1_pc2.png
# Run: Rscript scripts/03_plot_pca.R   (from the repo root)
# =============================================================================
suppressMessages(library(ggplot2))

sc <- read.csv("output/pca_scores.csv")
sm <- read.csv("output/pca_summary.csv")
meta <- read.csv("output/bz_plot_meta.csv")
sc$date <- meta$date[match(sc$PlotID, meta$PlotID)]
days <- unique(sc$date[order(as.Date(sc$date, "%m.%d.%y"))])
sc$date <- factor(sc$date, levels = days)
# ordinal single-hue ramp (blue 250 -> 700), one step per scan date
ramp <- c("#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b")
pct <- function(k) sprintf("PC%d (%.1f%% of variance)", k, 100 * sm$var_explained[k])

p <- ggplot(sc, aes(PC1, PC2, colour = date)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_colour_manual(values = setNames(ramp[seq_along(days)], days), name = "Scan date (2025)") +
  labs(
    x = pct(1), y = pct(2),
    title = "Principal components 1 and 2 of plot-level leaf reflectance, CLY25-B5 block 2",
    subtitle = sprintf(paste0(
      "PCA (centred, unscaled) of splice-corrected reflectance, 350–2500 nm, mean of 5 scans per plot; ",
      "n = %d plots, leaf 2 above the ear"), nrow(sc)),
    caption = paste0(
      "PC1 loads mainly on shortwave-infrared bands (peak 1893 nm); ",
      "PC2 loads mainly on near-infrared and red-edge bands (peak 715 nm).")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey35", size = 9),
    plot.caption = element_text(colour = "grey35", hjust = 0, size = 9),
    plot.title.position = "plot", plot.caption.position = "plot"
  )
ggsave("output/pca_pc1_pc2.png", p, width = 9, height = 6.5, dpi = 200)
