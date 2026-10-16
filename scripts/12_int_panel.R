#!/usr/bin/env Rscript
# =============================================================================
# Pre-specified trait panel, rank-based inverse-normal transformed BLUEs
# -----------------------------------------------------------------------------
# Panel fixed before computing (2026-10-16), all on the 799 non-senesced plots:
#   NDRE, PSRI, NDWI_2                (one index each: red edge, senescence, water)
#   PC1-PC4                           (top four PCA components by variance)
#   PC6                               (confirmatory test of the chr1 locus
#                                      S1_252335612, not discovery)
# Transform: Blom rank-INT, qnorm((rank - 3/8) / (n + 1/4)), ties averaged,
# applied to the SpATS genotype BLUEs. No plots or genotypes are removed.
#
# In : output/blues_indices/, output/blues_vi/, output/blues_green/
# Out: output/blues_int_panel/pheno_<token>_cly_blue.csv  (8 files, tokens
#      <original>INT, e.g. HSGNDREINT), output/blues_int_panel/panel_key.csv
# Run: Rscript scripts/12_int_panel.R   (from the repo root)
# =============================================================================
panel <- data.frame(
  trait = c("NDRE", "PSRI", "NDWI_2", "PC1", "PC2", "PC3", "PC4", "PC6"),
  role = c("red edge / chlorophyll", "senescence", "water status",
           rep("unsupervised (top 4 PCs by variance)", 4), "confirmatory: chr1 S1_252335612"),
  source = c("blues_indices/pheno_hsgndre", "blues_vi/pheno_hsgvpsri", "blues_vi/pheno_hsgvndwi2",
             sprintf("blues_green/pheno_hsgpc%d", c(1:4, 6)))
)
skew <- function(v) mean((v - mean(v))^3) / sd(v)^3
int <- function(v) qnorm((rank(v, ties.method = "average") - 3 / 8) / (length(v) + 1 / 4))
dir.create("output/blues_int_panel", showWarnings = FALSE)

key <- list()
for (i in seq_len(nrow(panel))) {
  x <- read.csv(sprintf("output/%s_cly_blue.csv", panel$source[i]))
  orig <- sub("_cly25$", "", names(x)[2])
  token <- paste0(orig, "INT")
  y <- int(x[[2]])
  out <- data.frame(Genotype = x$Genotype, y, y)
  names(out)[2:3] <- c(paste0(token, "_cly25"), paste0(token, "_CLY_mean"))
  write.csv(out, sprintf("output/blues_int_panel/pheno_%s_cly_blue.csv", tolower(token)), row.names = FALSE)
  key[[i]] <- data.frame(token = paste0(token, "_CLY"), trait = panel$trait[i], role = panel$role[i],
                         from = paste0(orig, "_CLY"), n_genotypes = nrow(x),
                         skew_before = round(skew(x[[2]]), 2), skew_after = round(skew(y), 2),
                         spearman_with_original = round(cor(x[[2]], y, method = "spearman"), 3))
}
key <- do.call(rbind, key)
write.csv(key, "output/blues_int_panel/panel_key.csv", row.names = FALSE)
print(key, row.names = FALSE)
