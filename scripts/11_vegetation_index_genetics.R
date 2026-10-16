#!/usr/bin/env Rscript
# =============================================================================
# Vegetation indices: generalized H2 and genotype BLUEs for R/qtl
# -----------------------------------------------------------------------------
# Part A, NDVI / NDRE from the plot-mean reflectance (first pass):
#   NDVI = (R800 - R670) / (R800 + R670), NDRE = (R790 - R720) / (R790 + R720)
#   Tokens HSNDVI, HSNDRE (all plots), HSGNDVI, HSGNDRE (non-senesced).
# Part B, all 28 Table-1 indices from 10_vegetation_indices.py:
#   Tokens HSV<ACR> (all plots) and HSGV<ACR> (non-senesced), <ACR> the acronym
#   upper-cased with non-alphanumerics removed ("WI/NDVI" -> WINDVI). The V
#   keeps them apart from Part A, whose NDVI uses R670 rather than R680.
# SpATS model and BLUE file layout as in lib.R / 06_component_genetics.R.
#
# In : output/bz_plot_reflectance.csv.gz, output/bz_plot_meta.csv,
#      output/outlier_group_plots.csv, output/vegetation_indices_table1.csv
# Out: output/vegetation_indices.csv (PlotID, NDVI, NDRE, senesced)
#      output/blues_indices/ (4 files), output/blues_vi/ (56 files + trait_key.csv)
# Run: Rscript scripts/11_vegetation_index_genetics.R   (from the repo root, ~3 min)
# =============================================================================
source("scripts/lib.R")
meta <- read.csv("output/bz_plot_meta.csv")
sen <- read.csv("output/outlier_group_plots.csv")$PlotID
rows_for <- function(ids) meta[match(ids, meta$PlotID), ]

# --- Part A: NDVI / NDRE -----------------------------------------------------
spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
nd <- function(a, b) (spec[[paste0("X", a)]] - spec[[paste0("X", b)]]) /
  (spec[[paste0("X", a)]] + spec[[paste0("X", b)]])
vi <- data.frame(PlotID = spec$PlotID, NDVI = nd(800, 670), NDRE = nd(790, 720),
                 senesced = spec$PlotID %in% sen)
write.csv(vi, "output/vegetation_indices.csv", row.names = FALSE)
for (set in c("all plots", "non-senesced")) {
  v <- if (set == "all plots") vi else vi[!vi$senesced, ]
  m <- rows_for(v$PlotID)
  cat(sprintf("\n== %s (n = %d)\n", set, nrow(v)))
  for (k in c("NDVI", "NDRE")) {
    q <- quantile(v[[k]], c(0, .25, .5, .75, 1))
    cat(sprintf("  %s  median %.3f  IQR %.3f-%.3f  range %.3f-%.3f | rho with Predicted_N %+.2f | H2 %.2f\n",
                k, q[3], q[2], q[4], q[1], q[5], cor(v[[k]], m$predN_mean, method = "spearman"),
                spats_h2(m, v[[k]])))
    write_blues(m, v[[k]], paste0(if (set == "all plots") "HS" else "HSG", k), "output/blues_indices")
  }
}

# --- Part B: Table-1 indices -------------------------------------------------
t1 <- read.csv("output/vegetation_indices_table1.csv", check.names = FALSE)
t1$senesced <- t1$senesced %in% c("True", "TRUE")
idx <- names(t1)[(which(names(t1) == "Predicted_N_mean") + 1):ncol(t1)]
key <- list()
for (set in c("all", "nonsenesced")) {
  d0 <- if (set == "all") t1 else t1[!t1$senesced, ]
  for (k in idx) {
    token <- paste0(if (set == "all") "HSV" else "HSGV", toupper(gsub("[^A-Za-z0-9]", "", k)))
    n <- write_blues(d0, d0[[k]], token, "output/blues_vi")
    key[[length(key) + 1]] <- data.frame(token = paste0(token, "_CLY"), index = k, plots = set,
                                         n_plots = nrow(d0), n_genotypes = n,
                                         H2 = round(spats_h2(d0, d0[[k]]), 2))
  }
}
key <- do.call(rbind, key)
write.csv(key, "output/blues_vi/trait_key.csv", row.names = FALSE)
cat(sprintf("\nTable-1 indices: %d BLUE files in output/blues_vi/\n", nrow(key)))
