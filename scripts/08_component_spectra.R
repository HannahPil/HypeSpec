#!/usr/bin/env Rscript
# =============================================================================
# Spectral signature of selected components (for QTL interpretation)
# -----------------------------------------------------------------------------
# PCs have loadings; autoencoder LVs do not. So every component is described
# the same way: Pearson r between its plot scores and reflectance at each
# wavelength, across the plots it was fit on. For PCs the loadings are also
# summarised (share of squared loading per 50-nm band).
#
# Env: TAG (default _green), COMPONENTS (default HSGPC4,HSGRLV5,HSGPC6,HSGPC2)
# In : output/bz_plot_reflectance.csv.gz, output/pca_scores<TAG>.csv,
#      output/pca_loadings<TAG>.csv, output/pca_summary<TAG>.csv,
#      output/ae_10lv<TAG>/LVs_Reflectance_Data_bz_10lv<TAG>.csv,
#      output/tross_pretrained_lvs.csv, output/pca_loadings.csv (all-plot PC8)
# Out: output/component_spectra<TAG>.csv (wavelength x component r),
#      output/component_spectra<TAG>.png, printed summary
# Run: Rscript scripts/08_component_spectra.R   (from the repo root)
# =============================================================================
suppressMessages(library(ggplot2))
TAG <- Sys.getenv("TAG", "_green")
COMP <- strsplit(Sys.getenv("COMPONENTS", "HSGPC4,HSGRLV5,HSGPC6,HSGPC2"), ",")[[1]]

spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
pcs <- read.csv(sprintf("output/pca_scores%s.csv", TAG))
L <- read.csv(sprintf("output/pca_loadings%s.csv", TAG))
sm <- read.csv(sprintf("output/pca_summary%s.csv", TAG))
rlv <- read.csv(sprintf("output/ae_10lv%s/LVs_Reflectance_Data_bz_10lv%s.csv", TAG, TAG))
tlv <- read.csv("output/tross_pretrained_lvs.csv")
spec <- spec[match(pcs$PlotID, spec$PlotID), ]
X <- as.matrix(spec[, grep("^X", names(spec))])
wav <- as.integer(sub("^X", "", colnames(X)))

score <- function(cname) { # HSG?PCk / HSG?RLVk / HSG?TLVk -> plot scores
  k <- sub("^HSG?(PC|RLV|TLV)", "", cname)
  src <- switch(sub("^HSG?(PC|RLV|TLV).*", "\\1", cname),
                PC = pcs, RLV = rlv, TLV = tlv)
  col <- if (grepl("PC", cname)) paste0("PC", k) else paste0("LV", k)
  src[[col]][match(pcs$PlotID, src$PlotID)]
}
S <- sapply(COMP, score)
R <- cor(S, X) # components x wavelengths

bands <- cut(wav, seq(349, 2549, 50), dig.lab = 4)
for (cn in COMP) {
  cat(sprintf("\n== %s\n", cn))
  if (grepl("PC", cn)) {
    k <- paste0("PC", sub("^HSG?PC", "", cn))
    cat(sprintf("  variance explained %.3f%% (rank %s of 10)\n",
                100 * sm$var_explained[sm$PC == k], sub("PC", "", k)))
    s <- tapply(L[[k]]^2, bands, sum); m <- tapply(L[[k]], bands, mean)
    o <- order(-s)[1:5]
    cat("  loadings, top 50-nm bands (share of squared loading, mean sign):\n")
    cat(sprintf("    %-14s %.3f  %+.3f\n", names(s)[o], s[o], m[o]), sep = "")
  }
  r <- R[cn, ]
  o <- order(-abs(r))
  cat(sprintf("  score-reflectance r: max |r| %.2f at %d nm; |r| >= 0.5 at %d of 2151 bands\n",
              max(abs(r)), wav[o[1]], sum(abs(r) >= 0.5)))
  rb <- tapply(r, bands, function(v) v[which.max(abs(v))])
  ob <- order(-abs(rb))[1:5]
  cat("  strongest bands (peak r within band):\n")
  cat(sprintf("    %-14s %+.2f\n", names(rb)[ob], rb[ob]), sep = "")
}

cat("\nSimilarity of score-reflectance profiles (Pearson r across 2151 wavelengths):\n")
print(round(cor(t(R)), 2))
cat("\nScore correlations across plots:\n")
print(round(cor(S), 2))

# all-plot PC8 vs this set's PC6 loadings
if ("HSGPC6" %in% COMP) {
  L0 <- read.csv("output/pca_loadings.csv")
  cat(sprintf("\nLoading correlation, all-plot PC8 vs %s PC6: %.2f\n", TAG, cor(L0$PC8, L$PC6)))
}

out <- data.frame(wavelength = wav, t(R))
write.csv(out, sprintf("output/component_spectra%s.csv", TAG), row.names = FALSE)
long <- data.frame(wavelength = rep(wav, length(COMP)), component = rep(COMP, each = length(wav)), r = c(t(R)))
long$component <- factor(long$component, levels = COMP)
p <- ggplot(long, aes(wavelength, r)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(colour = "#2a78d6", linewidth = 0.6) +
  facet_wrap(~component, ncol = 1) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Wavelength (nm)", y = "Pearson r (component score vs reflectance)",
       title = "Correlation of component scores with reflectance by wavelength, CLY25-B5 block 2",
       subtitle = sprintf("Plot-level splice-corrected reflectance, 350–2500 nm; n = %d non-senesced plots", nrow(X))) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 8.5))
ggsave(sprintf("output/component_spectra%s.png", TAG), p, width = 8, height = 8, dpi = 200)
