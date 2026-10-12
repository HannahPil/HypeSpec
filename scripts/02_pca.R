#!/usr/bin/env Rscript
# =============================================================================
# PCA baseline (Tross et al. 2024): 10 PCs of the plot-level reflectance
# -----------------------------------------------------------------------------
# Tross used scikit-learn PCA with defaults: bands centred, not scaled. prcomp
# with center = TRUE, scale. = FALSE is the same decomposition (centred SVD);
# only the sign of each PC is arbitrary.
#
# Also reports, per PC, Spearman rho (Tross's annotation statistic) with plot
# mean Predicted_N, scan day and field Range. No modelling of genotype here;
# heritability comes later with SpATS.
#
# In : output/bz_plot_reflectance.csv.gz, output/bz_plot_meta.csv
# Out: output/pca_scores.csv    (PlotID, PC1..PC10)
#      output/pca_loadings.csv  (wavelength, PC1..PC10)
#      output/pca_summary.csv   (per PC: variance explained, Spearman rhos)
# Run: Rscript scripts/02_pca.R   (from the repo root)
# =============================================================================

NPC <- 10
spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
meta <- read.csv("output/bz_plot_meta.csv")
X <- as.matrix(spec[, grep("^X", names(spec))])
wav <- as.integer(sub("^X", "", colnames(X)))

pca <- prcomp(X, center = TRUE, scale. = FALSE, rank. = NPC)
ve <- pca$sdev^2 / sum(pca$sdev^2)

scores <- data.frame(PlotID = spec$PlotID, pca$x)
write.csv(scores, "output/pca_scores.csv", row.names = FALSE)
write.csv(data.frame(wavelength = wav, pca$rotation), "output/pca_loadings.csv", row.names = FALSE)

m <- meta[match(spec$PlotID, meta$PlotID), ]
day <- as.numeric(as.Date(m$date, "%m.%d.%y"))
rho <- function(v) apply(pca$x, 2, function(pc) cor(pc, v, method = "spearman"))
summ <- data.frame(
  PC = colnames(pca$x), var_explained = ve[1:NPC], cum_var = cumsum(ve)[1:NPC],
  rho_predN = rho(m$predN_mean), rho_scan_day = rho(day), rho_range = rho(m$Range)
)
write.csv(summ, "output/pca_summary.csv", row.names = FALSE)
print(format(summ, digits = 3), row.names = FALSE)
