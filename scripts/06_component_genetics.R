#!/usr/bin/env Rscript
# =============================================================================
# Per component (10 PCs, 10 Tross-encoder LVs, 10 retrained-autoencoder LVs):
#   1. Spearman rho with plot-mean Predicted_N, scan day and Range
#   2. generalized heritability (SpATS, genotype random)
#   3. genotype BLUEs (SpATS, genotype fixed) for R/qtl in ZeaL_QTL
# -----------------------------------------------------------------------------
# Heritability: Tross used repeatability sigma2_G / (sigma2_G + sigma2_e / 2)
# from two reps. BZ has one plot per line, so genetic and plot-level variance
# are separated with the replicated checks instead (B73, 125 plots; Purple
# check, 2 plots). Field trend, which here is also scanning drift (scan day vs
# Range r = 0.98), is absorbed by the spatial surface. H2 from
# SpATS::getHeritability (generalized heritability, Oakey et al. 2006).
# Plot-mean Predicted_N is included as a reference trait for H2.
#
# BLUE trait names: HS{PC,TLV,RLV}<k> (all plots) or HSG{PC,TLV,RLV}<k> (tagged
# run, so staged traits never collide). Files follow zeal_predn_pheno.R:
# pheno_<trait>_cly_blue.csv with Genotype, <TRAIT>_cly25, <TRAIT>_CLY_mean.
#
# Env: TAG / EXCLUDE (lib.R)
# In : output/pca_scores<TAG>.csv, output/tross_pretrained_lvs.csv,
#      output/ae_10lv<TAG>/LVs_Reflectance_Data_bz_10lv<TAG>.csv,
#      output/bz_plot_meta.csv
# Out: output/component_correlations<TAG>.csv, output/heritability<TAG>.csv,
#      output/blues<TAG>/pheno_<trait>_cly_blue.csv (30 files)
# Run: Rscript scripts/06_component_genetics.R   (from the repo root, ~3 min)
# =============================================================================
source("scripts/lib.R")
cs <- load_component_sets()
meta <- cs$meta
sets <- cs$sets
comps <- function(d) grep("^(PC|LV)[0-9]+$", names(d), value = TRUE)

# --- 1. correlations ---------------------------------------------------------
cor_rows <- list()
for (s in names(sets)) {
  d <- sets[[s]]
  m <- meta[match(d$PlotID, meta$PlotID), ]
  day <- as.numeric(as.Date(m$date, "%m.%d.%y"))
  for (k in comps(d)) {
    n <- rho_p(d[[k]], m$predN_mean)
    dd <- rho_p(d[[k]], day)
    r <- rho_p(d[[k]], m$Range)
    cor_rows[[length(cor_rows) + 1]] <- data.frame(
      method = s, component = k, sd = signif(sd(d[[k]]), 3),
      spearman_rho_predicted_N = round(n[[1]], 3), p_predicted_N = signif(n[[2]], 2),
      spearman_rho_scan_day = round(dd[[1]], 3), p_scan_day = signif(dd[[2]], 2),
      spearman_rho_range = round(r[[1]], 3), p_range = signif(r[[2]], 2)
    )
  }
}
write.csv(do.call(rbind, cor_rows), tg("output/component_correlations.csv"), row.names = FALSE)

# --- 2. heritability ---------------------------------------------------------
h2_sets <- c(sets, list(Reference = data.frame(PlotID = meta$PlotID, Predicted_N = meta$predN_mean)))
h2_rows <- list()
for (s in names(h2_sets)) {
  d <- h2_sets[[s]]
  m <- meta[match(d$PlotID, meta$PlotID), ]
  for (k in setdiff(names(d), c("PlotID", "trt"))) {
    h <- spats_h2(m, d[[k]])
    h2_rows[[length(h2_rows) + 1]] <- data.frame(method = s, component = k, H2 = round(h, 3),
                                                 n_plots = nrow(m), n_genotypes = length(unique(m$Genotype)))
    cat(sprintf("%-22s %-12s H2 = %.3f\n", s, k, h))
  }
}
write.csv(do.call(rbind, h2_rows), tg("output/heritability.csv"), row.names = FALSE)

# --- 3. BLUEs ----------------------------------------------------------------
prefix <- setNames(paste0("HS", if (nzchar(TAG)) "G" else "", c("PC", "TLV", "RLV")), names(sets))
for (s in names(sets)) {
  d <- sets[[s]]
  m <- meta[match(d$PlotID, meta$PlotID), ]
  for (k in comps(d)) {
    trait <- paste0(prefix[[s]], sub("^(PC|LV)", "", k)) # e.g. HSPC1
    n <- write_blues(m, d[[k]], trait, paste0("output/blues", TAG))
    cat(sprintf("%-8s %d genotype BLUEs\n", trait, n))
  }
}
