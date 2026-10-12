#!/usr/bin/env Rscript
# =============================================================================
# Generalized heritability per component (SpATS, genotype random)
# -----------------------------------------------------------------------------
# Tross used repeatability sigma2_G / (sigma2_G + sigma2_e / 2) from two reps.
# BZ has one plot per line, so genetic and plot-level variance are separated
# with the replicated checks instead: B73 (125 plots) and the Purple check
# (2 plots) are genotype levels like any line, and their repeated plots inform
# the residual. Field trend, which here is also scanning drift (scan day vs
# Range r = 0.98), is absorbed by the spatial surface.
#
# Model: the SpATS recipe of ZeaL_QTL/scripts/zeal_predn_pheno.R
# (PSANOVA(Row, Range) surface, Row and Range factors random) with
# genotype.as.random = TRUE; H2 from SpATS::getHeritability (generalized
# heritability, Oakey et al. 2006).
#
# Traits: 10 PCs, 10 Tross-encoder LVs, 10 retrained-autoencoder LVs, and
# plot-mean Predicted_N as a reference.
#
# In : output/pca_scores.csv, output/tross_pretrained_lvs.csv,
#      output/ae_10lv/LVs_Reflectance_Data_bz_10lv.csv, output/bz_plot_meta.csv
# Out: output/heritability.csv  (method, component, H2, n_plots, n_genotypes)
# Run: Rscript scripts/09_heritability.R   (from the repo root)
# =============================================================================
suppressMessages(library(SpATS))

meta <- read.csv("output/bz_plot_meta.csv")
sets <- list(
  PCA = read.csv("output/pca_scores.csv"),
  `Tross encoder` = read.csv("output/tross_pretrained_lvs.csv"),
  `Retrained autoencoder` = read.csv("output/ae_10lv/LVs_Reflectance_Data_bz_10lv.csv"),
  Reference = data.frame(PlotID = meta$PlotID, Predicted_N = meta$predN_mean)
)

fit_h2 <- function(y, plot_id) {
  m <- meta[match(plot_id, meta$PlotID), ]
  d <- data.frame(Genotype = factor(m$Genotype), Row = as.numeric(m$Row),
                  Range = as.numeric(m$Range), y = y)
  d$Rf <- factor(d$Range)
  d$Cf <- factor(d$Row)
  fit <- SpATS(
    response = "y", genotype = "Genotype", genotype.as.random = TRUE,
    spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf,
    data = d, control = list(monitoring = 0, maxit = 50)
  )
  c(H2 = unname(getHeritability(fit)), n_plots = nrow(d), n_genotypes = nlevels(d$Genotype))
}

rows <- list()
for (s in names(sets)) {
  d <- sets[[s]]
  for (k in setdiff(names(d), c("PlotID", "trt"))) {
    r <- fit_h2(d[[k]], d$PlotID)
    rows[[length(rows) + 1]] <- data.frame(method = s, component = k, H2 = round(r[["H2"]], 3),
                                           n_plots = r[["n_plots"]], n_genotypes = r[["n_genotypes"]])
    cat(sprintf("%-22s %-12s H2 = %.3f\n", s, k, r[["H2"]]))
  }
}
write.csv(do.call(rbind, rows), "output/heritability.csv", row.names = FALSE)
