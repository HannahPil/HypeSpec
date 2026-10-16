#!/usr/bin/env Rscript
# =============================================================================
# Figure: mean leaf spectrum with the FieldSpec4 detector ranges and the main
# leaf absorption features marked (for presentation)
# In : output/bz_plot_reflectance.csv.gz, output/outlier_group_plots.csv
# Out: output/spectrum_regions.png
# Run: Rscript scripts/09_spectrum_regions_figure.R   (from the repo root)
# =============================================================================
suppressMessages(library(ggplot2))

spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
spec <- spec[!spec$PlotID %in% read.csv("output/outlier_group_plots.csv")$PlotID, ]
X <- as.matrix(spec[, grep("^X", names(spec))])
d <- data.frame(wavelength = 350:2500, R = colMeans(X), sd = apply(X, 2, sd))

det <- data.frame(
  xmin = c(350, 1000, 1800), xmax = c(1000, 1800, 2500),
  label = c("VNIR detector\n350–1000 nm", "SWIR1 detector\n1001–1800 nm", "SWIR2 detector\n1801–2500 nm"),
  fill = c("#cde2fb", "#b7d3f6", "#9ec5f4")
)
feat <- data.frame(
  x = c(450, 550, 670, 800, 880, 1200, 1450, 1940, 2100),
  label = c("Blue absorption\n(chlorophyll,\ncarotenoids)", "Green\npeak", "Red\nabsorption\n(chlorophyll)",
            "Red edge", "NIR plateau\n(leaf structure)", "Water\n~1200", "Water\n~1450",
            "Water\n~1940", "Dry matter\n(protein,\ncellulose, lignin)"),
  y = c(0.20, 0.27, 0.025, 0.33, 0.53, 0.52, 0.27, 0.17, 0.30)
)

p <- ggplot(d, aes(wavelength, R)) +
  geom_rect(data = det, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = I(fill)),
            inherit.aes = FALSE, alpha = 0.45) +
  geom_vline(xintercept = c(1000, 1800), colour = "#256abf", linetype = "dashed", linewidth = 0.4) +
  geom_text(data = det, aes(x = (xmin + xmax) / 2, y = 0.66, label = label), inherit.aes = FALSE,
            size = 3.3, fontface = "bold", colour = "#0d366b", lineheight = 0.9) +
  geom_ribbon(aes(ymin = R - sd, ymax = R + sd), fill = "grey40", alpha = 0.2) +
  geom_line(linewidth = 0.8, colour = "#222222") +
  geom_text(data = feat, aes(x = x, y = y, label = label), size = 2.7, colour = "#333333", lineheight = 0.85) +
  scale_x_continuous(breaks = seq(500, 2500, 250), expand = c(0.01, 0)) +
  coord_cartesian(ylim = c(-0.01, 0.72)) +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = "Mean maize leaf reflectance spectrum with spectroradiometer detector ranges",
       subtitle = sprintf(paste0(
         "ASD FieldSpec4 with contact probe, 350–2500 nm at 1 nm (2,151 bands); mean ± 1 SD of %d plot-level spectra ",
         "\n(leaf 2 above the ear, CLY25-B5 block 2, senesced plots excluded); dashed lines = detector joins (splice-corrected)"), nrow(X)),
       caption = "Feature positions after Curran (1989), Ustin et al. (2009) and Ceccato et al. (2001).") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey35", hjust = 0, size = 8.5), plot.caption.position = "plot")
ggsave("output/spectrum_regions.png", p, width = 11, height = 5.8, dpi = 220)
