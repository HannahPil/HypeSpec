#!/usr/bin/env Rscript
# =============================================================================
# All plots (921) vs non-senesced plots (799, "_green"): heritability and
# Predicted_N correlation per method, side by side
# -----------------------------------------------------------------------------
# Components are refit in the _green set (PCA, retrained autoencoder), so
# component k in one set is not component k in the other; comparisons are per
# method, not per component. The Tross encoder is fixed, so its LVs are the
# same functions in both sets.
#
# In : output/heritability{,_green}.csv, output/component_correlations{,_green}.csv
# Out: output/compare_sets_by_method.csv, output/compare_sets_by_component.csv
# Run: Rscript scripts/12_compare_sets.R   (from the repo root)
# =============================================================================

load_set <- function(tag, label) {
  h <- read.csv(sprintf("output/heritability%s.csv", tag))
  r <- read.csv(sprintf("output/component_correlations%s.csv", tag))
  x <- merge(h, r[, c("method", "component", "spearman_rho_predicted_N", "spearman_rho_scan_day")],
             by = c("method", "component"), all.x = TRUE)
  x$set <- label
  x
}
x <- rbind(load_set("", "all_921"), load_set("_green", "nonsenesced_799"))
x <- x[x$method != "Reference", ]

by_method <- do.call(rbind, lapply(split(x, list(x$method, x$set), drop = TRUE), function(d) {
  b <- which.max(abs(d$spearman_rho_predicted_N))
  data.frame(method = d$method[1], set = d$set[1], n_plots = d$n_plots[1],
             H2_mean = round(mean(d$H2), 2), H2_min = min(d$H2), H2_max = max(d$H2),
             n_H2_ge_0.5 = sum(d$H2 >= 0.5),
             best_predN_component = d$component[b],
             best_predN_rho = d$spearman_rho_predicted_N[b],
             best_predN_H2 = d$H2[b],
             median_abs_rho_scan_day = round(median(abs(d$spearman_rho_scan_day)), 2))
}))
by_method <- by_method[order(match(by_method$method, c("PCA", "Tross encoder", "Retrained autoencoder")),
                             by_method$set), ]
ref <- rbind(read.csv("output/heritability.csv"), read.csv("output/heritability_green.csv"))
ref <- ref[ref$method == "Reference", ]
write.csv(by_method, "output/compare_sets_by_method.csv", row.names = FALSE)
x <- x[order(x$set, match(x$method, c("PCA", "Tross encoder", "Retrained autoencoder")),
             as.numeric(gsub("[^0-9]", "", x$component))), ]
write.csv(x, "output/compare_sets_by_component.csv", row.names = FALSE)
print(by_method, row.names = FALSE)
cat(sprintf("\nPredicted_N H2: all %.2f, non-senesced %.2f\n", ref$H2[1], ref$H2[2]))
