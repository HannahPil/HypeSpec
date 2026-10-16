#!/usr/bin/env Rscript
# =============================================================================
# PCA (Tross et al. 2024 baseline): 10 PCs of the plot-level reflectance,
# component table, and PC1 vs PC2 figures
# -----------------------------------------------------------------------------
# Tross used scikit-learn PCA with defaults: bands centred, not scaled. prcomp
# with center = TRUE, scale. = FALSE is the same decomposition (centred SVD);
# only the sign of each PC is arbitrary.
#
# Component table, per PC: variance explained; share of the squared loadings
# in each spectral region (each PC's squared loadings sum to 1); wavelength of
# the largest absolute loading and the leaf feature it falls in (lookup in
# lib.R, not a measured trait); Spearman rho and p with plot-mean Predicted_N,
# scan day and Range.
#
# Env: TAG / EXCLUDE (lib.R), e.g. EXCLUDE=output/outlier_group_plots.csv
#      TAG=_green for the non-senesced plots.
# In : output/bz_plot_reflectance.csv.gz, output/bz_plot_meta.csv
# Out: output/pca_scores<TAG>.csv, pca_loadings<TAG>.csv, pca_summary<TAG>.csv,
#      pca_component_table<TAG>.csv, pca_pc1_pc2_taxa<TAG>.png,
#      pca_pc1_pc2.png (scan-date colouring; untagged run only)
# Run: Rscript scripts/02_pca.R   (from the repo root)
# =============================================================================
source("scripts/lib.R")
suppressMessages(library(ggplot2))
NPC <- 10

spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
meta <- read.csv("output/bz_plot_meta.csv")
spec <- spec[!spec$PlotID %in% excluded(), ]
if (nzchar(EXCLUDE)) cat(sprintf("EXCLUDE: %d plots remain\n", nrow(spec)))
X <- as.matrix(spec[, grep("^X", names(spec))])
wav <- as.integer(sub("^X", "", colnames(X)))

# --- PCA ---------------------------------------------------------------------
pca <- prcomp(X, center = TRUE, scale. = FALSE, rank. = NPC)
ve <- pca$sdev^2 / sum(pca$sdev^2)
sc <- data.frame(PlotID = spec$PlotID, pca$x)
L <- data.frame(wavelength = wav, pca$rotation)
write.csv(sc, tg("output/pca_scores.csv"), row.names = FALSE)
write.csv(L, tg("output/pca_loadings.csv"), row.names = FALSE)

m <- meta[match(spec$PlotID, meta$PlotID), ]
day <- as.numeric(as.Date(m$date, "%m.%d.%y"))
rho <- function(v) apply(pca$x, 2, function(pc) cor(pc, v, method = "spearman"))
sm <- data.frame(
  PC = colnames(pca$x), var_explained = ve[1:NPC], cum_var = cumsum(ve)[1:NPC],
  rho_predN = rho(m$predN_mean), rho_scan_day = rho(day), rho_range = rho(m$Range)
)
write.csv(sm, tg("output/pca_summary.csv"), row.names = FALSE)

# --- component table ---------------------------------------------------------
reg <- cut(L$wavelength, c(349, 700, 1300, 1800, 2500),
           c("VIS_350_700", "NIR_701_1300", "SWIR1_1301_1800", "SWIR2_1801_2500"))
tab <- do.call(rbind, lapply(paste0("PC", 1:NPC), function(k) {
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
}))
write.csv(tab, tg("output/pca_component_table.csv"), row.names = FALSE)
print(tab[, c("component", "variance_explained_pct", "peak_loading_nm", "feature_at_peak",
              "spearman_rho_predicted_N", "spearman_rho_scan_day")], row.names = FALSE)

# --- figures -----------------------------------------------------------------
pct <- function(k) sprintf("PC%d (%.1f%% of variance)", k, 100 * sm$var_explained[k])
sub_txt <- sprintf(paste0(
  "PCA (centred, unscaled) of splice-corrected reflectance, 350–2500 nm, mean of 5 scans per plot; ",
  "n = %d plots%s"), nrow(sc), if (nzchar(TAG)) "; senesced (high-PC1) plots excluded" else "")
theme_pca <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption = element_text(colour = "grey35", hjust = 0, size = 9),
        plot.caption.position = "plot")

# by taxon (Chal, Dura, Mesa, Nobo grouped as mexicana via taxa_mex)
set.seed(1) # fixed random draw order
MEXICANA <- c("Chal", "Dura", "Mesa", "Nobo")
sc$group <- m$taxa_mex
lev <- c("B73", "mexicana", "Bals", "Zdip", "Zlux", "Hueh", "Purple")
stopifnot(all(sc$group %in% lev))
n <- table(factor(sc$group, levels = lev))
cols <- setNames(c("#222222", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"), lev)
labs_n <- sprintf("%s (n = %d)", ifelse(lev == "mexicana",
  sprintf("mexicana: %s", paste(MEXICANA, collapse = ", ")), lev), n)
sc$group <- factor(sc$group, levels = lev)
p <- ggplot(sc[sample(nrow(sc)), ], aes(PC1, PC2, colour = group)) +
  geom_point(size = 1.7, alpha = 0.85) +
  scale_colour_manual(values = cols, labels = labs_n, name = "Taxon") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = pct(1), y = pct(2), subtitle = sub_txt,
       title = "Principal components 1 and 2 of plot-level leaf reflectance by taxon, CLY25-B5 block 2") +
  theme_pca
ggsave(tg("output/pca_pc1_pc2_taxa.png"), p, width = 9.5, height = 6.5, dpi = 200)

# by scan date (all-plot run only)
if (!nzchar(TAG)) {
  sc$date <- m$date
  days <- unique(sc$date[order(as.Date(sc$date, "%m.%d.%y"))])
  sc$date <- factor(sc$date, levels = days)
  ramp <- c("#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b")
  p <- ggplot(sc, aes(PC1, PC2, colour = date)) +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_colour_manual(values = setNames(ramp[seq_along(days)], days), name = "Scan date (2025)") +
    labs(x = pct(1), y = pct(2), subtitle = paste0(sub_txt, ", leaf 2 above the ear"),
         title = "Principal components 1 and 2 of plot-level leaf reflectance, CLY25-B5 block 2",
         caption = sprintf("PC1 peak loading at %d nm (%s); PC2 peak loading at %d nm (%s).",
                           tab$peak_loading_nm[1], tolower(tab$feature_at_peak[1]),
                           tab$peak_loading_nm[2], tolower(tab$feature_at_peak[2]))) +
    theme_pca
  ggsave("output/pca_pc1_pc2.png", p, width = 9, height = 6.5, dpi = 200)
}
