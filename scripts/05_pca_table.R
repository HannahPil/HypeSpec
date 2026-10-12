#!/usr/bin/env Rscript
# =============================================================================
# Summary table of the 10 PCs for presentation
# -----------------------------------------------------------------------------
# Per PC: variance explained; share of the squared loadings in each detector /
# spectral region (each PC's squared loadings sum to 1); wavelength of the
# largest absolute loading and the leaf absorption feature that wavelength
# falls in (lookup below, after Curran 1989, Ustin et al. 2009, Ceccato et al.
# 2001); Spearman rho and p with plot-mean Predicted_N, scan day and Range.
# The feature column is a wavelength lookup, not a measured trait.
#
# In : output/pca_loadings.csv, output/pca_scores.csv, output/pca_summary.csv,
#      output/bz_plot_meta.csv
# Out: output/pca_component_table.csv
# Run: Rscript scripts/05_pca_table.R   (from the repo root)
# =============================================================================

TAG <- Sys.getenv("TAG", "") # input/output suffix, see 02_pca.R
tg <- function(p) sub("\\.(csv|png)$", paste0(TAG, ".\\1"), p)
L <- read.csv(tg("output/pca_loadings.csv"))
sc <- read.csv(tg("output/pca_scores.csv"))
sm <- read.csv(tg("output/pca_summary.csv"))
meta <- read.csv("output/bz_plot_meta.csv")
m <- meta[match(sc$PlotID, meta$PlotID), ]
day <- as.numeric(as.Date(m$date, "%m.%d.%y"))

# leaf reflectance features by wavelength (nm), [lo, hi)
feat <- data.frame(
  lo = c(350, 400, 500, 600, 690, 750, 950, 990, 1180, 1220, 1300, 1500, 1800, 2000),
  hi = c(400, 500, 600, 690, 750, 950, 990, 1180, 1220, 1300, 1500, 1800, 2000, 2501),
  label = c(
    "Ultraviolet, lower edge of instrument range (low signal)",
    "Blue absorption (chlorophyll a/b, carotenoids)",
    "Green reflectance peak (chlorophyll)",
    "Red absorption (chlorophyll)",
    "Red edge (chlorophyll, nitrogen)",
    "Near-infrared plateau (leaf structure)",
    "Weak water absorption (~970 nm)",
    "Near-infrared plateau (leaf structure)",
    "Weak water absorption (~1200 nm)",
    "Near-infrared plateau (leaf structure)",
    "Water absorption (~1450 nm)",
    "Shortwave-infrared dry matter (protein, cellulose, lignin)",
    "Water absorption (~1940 nm)",
    "Shortwave-infrared dry matter (protein, cellulose, lignin)"
  )
)
feature_at <- function(w) feat$label[w >= feat$lo & w < feat$hi]
reg <- cut(L$wavelength, c(349, 700, 1300, 1800, 2500),
           c("VIS_350_700", "NIR_701_1300", "SWIR1_1301_1800", "SWIR2_1801_2500"))
rho_p <- function(pc, v) {
  t <- suppressWarnings(cor.test(pc, v, method = "spearman", exact = FALSE))
  c(rho = unname(t$estimate), p = t$p.value)
}

rows <- lapply(paste0("PC", 1:10), function(k) {
  share <- tapply(L[[k]]^2, reg, sum)
  peak <- L$wavelength[which.max(abs(L[[k]]))]
  n <- rho_p(sc[[k]], m$predN_mean)
  d <- rho_p(sc[[k]], day)
  r <- rho_p(sc[[k]], m$Range)
  data.frame(
    component = k,
    variance_explained_pct = round(100 * sm$var_explained[sm$PC == k], 3),
    cumulative_variance_pct = round(100 * sm$cum_var[sm$PC == k], 3),
    loading_share_VIS_350_700 = round(share[["VIS_350_700"]], 2),
    loading_share_NIR_701_1300 = round(share[["NIR_701_1300"]], 2),
    loading_share_SWIR1_1301_1800 = round(share[["SWIR1_1301_1800"]], 2),
    loading_share_SWIR2_1801_2500 = round(share[["SWIR2_1801_2500"]], 2),
    peak_loading_nm = peak,
    feature_at_peak = feature_at(peak),
    spearman_rho_predicted_N = round(n[["rho"]], 3),
    p_predicted_N = signif(n[["p"]], 2),
    spearman_rho_scan_day = round(d[["rho"]], 3),
    p_scan_day = signif(d[["p"]], 2),
    spearman_rho_range = round(r[["rho"]], 3),
    p_range = signif(r[["p"]], 2)
  )
})
tab <- do.call(rbind, rows)
write.csv(tab, tg("output/pca_component_table.csv"), row.names = FALSE)
print(tab[, c("component", "variance_explained_pct", "peak_loading_nm", "feature_at_peak",
              "spearman_rho_predicted_N", "spearman_rho_scan_day")], row.names = FALSE)
