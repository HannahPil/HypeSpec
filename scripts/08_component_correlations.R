#!/usr/bin/env Rscript
# =============================================================================
# Spearman correlations of every component (PCA, Tross encoder, retrained
# autoencoder) with plot-mean Predicted_N, scan day and field Range
# -----------------------------------------------------------------------------
# Same statistic as Tross's LV annotation (Spearman rho). Latent variables have
# no loadings or variance-explained, so those columns from the PCA table are
# not repeated here.
#
# In : output/pca_scores.csv, output/tross_pretrained_lvs.csv,
#      output/ae_10lv/LVs_Reflectance_Data_bz_10lv.csv, output/bz_plot_meta.csv
# Out: output/component_correlations.csv
# Run: Rscript scripts/08_component_correlations.R   (from the repo root)
# =============================================================================

source("scripts/component_sets.R") # meta, sets, tg(); TAG / EXCLUDE env
rho_p <- function(x, v) {
  t <- suppressWarnings(cor.test(x, v, method = "spearman", exact = FALSE))
  c(unname(t$estimate), t$p.value)
}

rows <- list()
for (s in names(sets)) {
  d <- sets[[s]]
  m <- meta[match(d$PlotID, meta$PlotID), ]
  stopifnot(!anyNA(m$PlotID))
  day <- as.numeric(as.Date(m$date, "%m.%d.%y"))
  for (k in grep("^(PC|LV)[0-9]+$", names(d), value = TRUE)) {
    n <- rho_p(d[[k]], m$predN_mean)
    dd <- rho_p(d[[k]], day)
    r <- rho_p(d[[k]], m$Range)
    rows[[length(rows) + 1]] <- data.frame(
      method = s, component = k, sd = signif(sd(d[[k]]), 3),
      spearman_rho_predicted_N = round(n[1], 3), p_predicted_N = signif(n[2], 2),
      spearman_rho_scan_day = round(dd[1], 3), p_scan_day = signif(dd[2], 2),
      spearman_rho_range = round(r[1], 3), p_range = signif(r[2], 2)
    )
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, tg("output/component_correlations.csv"), row.names = FALSE)
print(tab[, c("method", "component", "spearman_rho_predicted_N", "spearman_rho_scan_day",
              "spearman_rho_range")], row.names = FALSE)
