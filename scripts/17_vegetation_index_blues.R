#!/usr/bin/env Rscript
# =============================================================================
# Genotype BLUEs (and generalized H2) for all 28 Table-1 vegetation indices,
# all plots and non-senesced plots, in ZeaL_QTL phenotype format
# -----------------------------------------------------------------------------
# Same SpATS model as 10_blues.R (genotype fixed for BLUEs, random for H2;
# PSANOVA(Row, Range), Row and Range factors random).
# Trait tokens: HSV<ACR> (all 921 plots) and HSGV<ACR> (799 non-senesced),
# where <ACR> is the acronym upper-cased with non-alphanumerics removed
# (e.g. "WI/NDVI" -> WINDVI, "CIred-edge" -> CIREDEDGE, "Red-edge" -> REDEDGE).
# The V keeps them apart from HSNDVI / HSNDRE of 15_vegetation_indices.R,
# whose NDVI uses R670 rather than Table 1's R680.
#
# In : output/vegetation_indices_table1.csv, output/bz_plot_meta.csv
# Out: output/blues_vi/pheno_<token>_cly_blue.csv  (56 files)
#      output/blues_vi/trait_key.csv  (token, index, plot set, n, H2)
# Run: Rscript scripts/17_vegetation_index_blues.R   (from the repo root)
# =============================================================================
suppressMessages(library(SpATS))

vi <- read.csv("output/vegetation_indices_table1.csv", check.names = FALSE)
vi$senesced <- vi$senesced %in% c("True", "TRUE")
idx <- names(vi)[(which(names(vi) == "Predicted_N_mean") + 1):ncol(vi)]
dir.create("output/blues_vi", showWarnings = FALSE)

fit <- function(d0, y, random) {
  d <- data.frame(Genotype = factor(d0$Genotype), Row = as.numeric(d0$Row),
                  Range = as.numeric(d0$Range), y = y)
  d$Rf <- factor(d$Range); d$Cf <- factor(d$Row)
  SpATS("y", genotype = "Genotype", genotype.as.random = random,
        spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf, data = d,
        control = list(monitoring = 0, maxit = 50))
}

key <- list()
for (set in c("all", "nonsenesced")) {
  d0 <- if (set == "all") vi else vi[!vi$senesced, ]
  for (k in idx) {
    token <- paste0(if (set == "all") "HSV" else "HSGV", toupper(gsub("[^A-Za-z0-9]", "", k)))
    pr <- predict(fit(d0, d0[[k]], FALSE), which = "Genotype")
    vc <- grep("predicted", names(pr), ignore.case = TRUE, value = TRUE)[1]
    pr <- pr[!is.na(pr[[vc]]), ]
    out <- data.frame(Genotype = as.character(pr$Genotype), v = pr[[vc]])
    out <- out[order(out$Genotype), ]
    names(out)[2] <- paste0(token, "_cly25")
    out[[paste0(token, "_CLY_mean")]] <- out[[2]]
    write.csv(out, sprintf("output/blues_vi/pheno_%s_cly_blue.csv", tolower(token)), row.names = FALSE)
    key[[length(key) + 1]] <- data.frame(
      token = paste0(token, "_CLY"), index = k, plots = set, n_plots = nrow(d0),
      n_genotypes = nrow(out), H2 = round(unname(getHeritability(fit(d0, d0[[k]], TRUE))), 2))
  }
}
key <- do.call(rbind, key)
write.csv(key, "output/blues_vi/trait_key.csv", row.names = FALSE)
print(key, row.names = FALSE)
