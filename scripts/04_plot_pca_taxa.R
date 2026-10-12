#!/usr/bin/env Rscript
# =============================================================================
# Figure: PC1 vs PC2 of plot-level reflectance, coloured by taxon
# Chal, Dura, Mesa and Nobo are grouped as Zea mays ssp. mexicana. B73
# (recurrent parent) in near-black, Purple checks in violet, the other groups
# take categorical hues in fixed order.
# In : output/pca_scores.csv, output/pca_summary.csv, output/bz_plot_meta.csv
# Out: output/pca_pc1_pc2_taxa.png
# Run: Rscript scripts/04_plot_pca_taxa.R   (from the repo root)
# =============================================================================
suppressMessages(library(ggplot2))
set.seed(1) # fixed random draw order

MEXICANA <- c("Chal", "Dura", "Mesa", "Nobo")
TAG <- Sys.getenv("TAG", "") # input/output suffix, see 02_pca.R
tg <- function(p) sub("\\.(csv|png)$", paste0(TAG, ".\\1"), p)
sc <-read.csv(tg("output/pca_scores.csv"))
sm <- read.csv(tg("output/pca_summary.csv"))
meta <- read.csv("output/bz_plot_meta.csv")
sc$group <- meta$taxa_mex[match(sc$PlotID, meta$PlotID)]
lev <- c("B73", "mexicana", "Bals", "Zdip", "Zlux", "Hueh", "Purple")
stopifnot(all(sc$group %in% lev))
n <- table(factor(sc$group, levels = lev))
cols <- setNames(c("#222222", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"), lev)
labs_n <- sprintf("%s (n = %d)", ifelse(lev == "mexicana",
  sprintf("mexicana: %s", paste(MEXICANA, collapse = ", ")), lev), n)
sc$group <- factor(sc$group, levels = lev)
pct <- function(k) sprintf("PC%d (%.1f%% of variance)", k, 100 * sm$var_explained[k])

p <- ggplot(sc[sample(nrow(sc)), ], aes(PC1, PC2, colour = group)) + # random draw order
  geom_point(size = 1.7, alpha = 0.85) +
  scale_colour_manual(values = cols, labels = labs_n, name = "Taxon") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    x = pct(1), y = pct(2),
    title = "Principal components 1 and 2 of plot-level leaf reflectance by taxon, CLY25-B5 block 2",
    subtitle = sprintf(paste0(
      "PCA (centred, unscaled) of splice-corrected reflectance, 350–2500 nm, mean of 5 scans per plot; ",
      "n = %d plots%s"), nrow(sc),
      if (nzchar(TAG)) "; senesced (high-PC1) plots excluded" else "")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey35", size = 9),
    plot.title.position = "plot"
  )
ggsave(tg("output/pca_pc1_pc2_taxa.png"), p, width = 9.5, height = 6.5, dpi = 200)
