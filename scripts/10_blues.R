#!/usr/bin/env Rscript
# =============================================================================
# Genotype BLUEs per component (SpATS, genotype fixed) for R/qtl in ZeaL_QTL
# -----------------------------------------------------------------------------
# Same model as ZeaL_QTL/scripts/zeal_predn_pheno.R: genotype FIXED,
# PSANOVA(Row, Range) surface, Row and Range factors random, one fit over the
# block-2 grid. Output files follow that script's format so they can be staged
# into ZeaL_QTL/data/zeal/ unchanged:
#   pheno_<trait>_blue.csv with columns Genotype, <TRAIT>_cly25, <TRAIT>_CLY_mean
#
# Trait names (all CLY25 block 2):
#   HSPC1..10  PCA components
#   HSTLV1..10 Tross et al. 2024 trained encoder latent variables
#   HSRLV1..10 autoencoder retrained on ZEAL, latent variables
#
# In : output/pca_scores.csv, output/tross_pretrained_lvs.csv,
#      output/ae_10lv/LVs_Reflectance_Data_bz_10lv.csv, output/bz_plot_meta.csv
# Out: output/blues/pheno_<trait>_cly_blue.csv  (30 files)
# Run: Rscript scripts/10_blues.R   (from the repo root)
# Next (in ZeaL_QTL, after staging): TRAITS=HSPC1_CLY,... Rscript scripts/zeal_rqtl_scan.R
# =============================================================================
suppressMessages(library(SpATS))

meta <- read.csv("output/bz_plot_meta.csv")
sets <- list(
  HSPC = read.csv("output/pca_scores.csv"),
  HSTLV = read.csv("output/tross_pretrained_lvs.csv"),
  HSRLV = read.csv("output/ae_10lv/LVs_Reflectance_Data_bz_10lv.csv")
)
dir.create("output/blues", showWarnings = FALSE)

for (prefix in names(sets)) {
  d <- sets[[prefix]]
  m <- meta[match(d$PlotID, meta$PlotID), ]
  for (k in grep("^(PC|LV)[0-9]+$", names(d), value = TRUE)) {
    trait <- paste0(prefix, sub("^(PC|LV)", "", k)) # e.g. HSPC1
    fd <- data.frame(Genotype = factor(m$Genotype), Row = as.numeric(m$Row),
                     Range = as.numeric(m$Range), y = d[[k]])
    fd$Rf <- factor(fd$Range)
    fd$Cf <- factor(fd$Row)
    fit <- SpATS(
      response = "y", genotype = "Genotype", genotype.as.random = FALSE,
      spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf,
      data = fd, control = list(monitoring = 0, maxit = 50)
    )
    pr <- predict(fit, which = "Genotype")
    vc <- grep("predicted", names(pr), ignore.case = TRUE, value = TRUE)[1]
    pr <- pr[!is.na(pr[[vc]]), ]
    out <- data.frame(Genotype = as.character(pr$Genotype), v = pr[[vc]])
    out <- out[order(out$Genotype), ]
    names(out)[2] <- paste0(trait, "_cly25")
    out[[paste0(trait, "_CLY_mean")]] <- out[[2]]
    f <- sprintf("output/blues/pheno_%s_cly_blue.csv", tolower(trait))
    write.csv(out, f, row.names = FALSE)
    cat(sprintf("%-8s %d genotype BLUEs -> %s\n", trait, nrow(out), basename(f)))
  }
}
