#!/usr/bin/env Rscript
# =============================================================================
# Narrowband vegetation indices per plot: NDVI and NDRE
# -----------------------------------------------------------------------------
#   NDVI = (R800 - R670) / (R800 + R670)   red vs near-infrared
#   NDRE = (R790 - R720) / (R790 + R720)   red edge vs near-infrared
# Single 1-nm bands from the splice-corrected plot means. Reports, for all
# plots and for the non-senesced plots: distribution, Spearman rho with
# plot-mean Predicted_N, and SpATS generalized H2 (same model as
# scripts/09_heritability.R).
#
# In : output/bz_plot_reflectance.csv.gz, output/bz_plot_meta.csv,
#      output/outlier_group_plots.csv
# Out: output/vegetation_indices.csv (PlotID, NDVI, NDRE, senesced)
# Run: Rscript scripts/15_vegetation_indices.R   (from the repo root)
# =============================================================================
suppressMessages(library(SpATS))

spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
meta <- read.csv("output/bz_plot_meta.csv")
sen <- read.csv("output/outlier_group_plots.csv")$PlotID
nd <- function(a, b) (spec[[paste0("X", a)]] - spec[[paste0("X", b)]]) /
  (spec[[paste0("X", a)]] + spec[[paste0("X", b)]])
vi <- data.frame(PlotID = spec$PlotID, NDVI = nd(800, 670), NDRE = nd(790, 720),
                 senesced = spec$PlotID %in% sen)
write.csv(vi, "output/vegetation_indices.csv", row.names = FALSE)

h2 <- function(y, ids) {
  m <- meta[match(ids, meta$PlotID), ]
  d <- data.frame(Genotype = factor(m$Genotype), Row = as.numeric(m$Row),
                  Range = as.numeric(m$Range), y = y)
  d$Rf <- factor(d$Range); d$Cf <- factor(d$Row)
  fit <- SpATS("y", genotype = "Genotype", genotype.as.random = TRUE,
               spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf, data = d,
               control = list(monitoring = 0, maxit = 50))
  unname(getHeritability(fit))
}

for (set in c("all plots", "non-senesced")) {
  v <- if (set == "all plots") vi else vi[!vi$senesced, ]
  pn <- meta$predN_mean[match(v$PlotID, meta$PlotID)]
  cat(sprintf("\n== %s (n = %d)\n", set, nrow(v)))
  for (k in c("NDVI", "NDRE")) {
    q <- quantile(v[[k]], c(0, .25, .5, .75, 1))
    cat(sprintf("  %s  median %.3f  IQR %.3f-%.3f  range %.3f-%.3f | rho with Predicted_N %+.2f | H2 %.2f\n",
                k, q[3], q[2], q[4], q[1], q[5], cor(v[[k]], pn, method = "spearman"),
                h2(v[[k]], v$PlotID)))
  }
}
cat(sprintf("\nSenesced vs other plots, median NDVI %.3f vs %.3f, NDRE %.3f vs %.3f\n",
            median(vi$NDVI[vi$senesced]), median(vi$NDVI[!vi$senesced]),
            median(vi$NDRE[vi$senesced]), median(vi$NDRE[!vi$senesced])))

# --- genotype BLUEs for R/qtl (same model and file layout as 10_blues.R) -----
# HSNDVI / HSNDRE: all plots; HSGNDVI / HSGNDRE: non-senesced plots.
dir.create("output/blues_indices", showWarnings = FALSE)
blue <- function(y, ids, trait) {
  m <- meta[match(ids, meta$PlotID), ]
  d <- data.frame(Genotype = factor(m$Genotype), Row = as.numeric(m$Row),
                  Range = as.numeric(m$Range), y = y)
  d$Rf <- factor(d$Range); d$Cf <- factor(d$Row)
  fit <- SpATS("y", genotype = "Genotype", genotype.as.random = FALSE,
               spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf, data = d,
               control = list(monitoring = 0, maxit = 50))
  pr <- predict(fit, which = "Genotype")
  vc <- grep("predicted", names(pr), ignore.case = TRUE, value = TRUE)[1]
  pr <- pr[!is.na(pr[[vc]]), ]
  out <- data.frame(Genotype = as.character(pr$Genotype), v = pr[[vc]])
  out <- out[order(out$Genotype), ]
  names(out)[2] <- paste0(trait, "_cly25")
  out[[paste0(trait, "_CLY_mean")]] <- out[[2]]
  f <- sprintf("output/blues_indices/pheno_%s_cly_blue.csv", tolower(trait))
  write.csv(out, f, row.names = FALSE)
  cat(sprintf("%-8s %d genotype BLUEs -> %s\n", trait, nrow(out), f))
}
g <- vi[!vi$senesced, ]
for (k in c("NDVI", "NDRE")) {
  blue(vi[[k]], vi$PlotID, paste0("HS", k))
  blue(g[[k]], g$PlotID, paste0("HSG", k))
}
