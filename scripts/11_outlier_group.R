#!/usr/bin/env Rscript
# =============================================================================
# Characterise the high-PC1 group (PC1 > 5): spectra, field position,
# composition, Predicted_N / M-distance, and whether membership clusters by
# accession beyond chance.
# -----------------------------------------------------------------------------
# Accession test: observed variance of the per-accession membership fraction
# vs 10,000 permutations of membership across lines (B73 and Purple checks
# excluded). Clustering by accession would point to a genetic basis; a random
# scatter would not rule one out (lines within an accession differ).
#
# In : output/bz_plot_reflectance.csv.gz, output/pca_scores.csv,
#      output/bz_plot_meta.csv
# Out: output/outlier_group_spectra.png, output/outlier_group_fieldmap.png,
#      output/outlier_group_summary.txt, output/outlier_group_plots.csv
# Run: Rscript scripts/11_outlier_group.R   (from the repo root)
# =============================================================================
suppressMessages(library(ggplot2))
set.seed(1)
CUT <- 5

spec <- read.csv("output/bz_plot_reflectance.csv.gz", check.names = FALSE)
pc <- read.csv("output/pca_scores.csv")
meta <- read.csv("output/bz_plot_meta.csv")
meta$PC1 <- pc$PC1[match(meta$PlotID, pc$PlotID)]
meta$group <- ifelse(meta$PC1 > CUT, "High-PC1 group", "Other plots")
X <- as.matrix(spec[, grep("^X", names(spec))])
wav <- as.integer(sub("^X", "", colnames(X)))
g <- meta$group[match(spec$PlotID, meta$PlotID)]

sink("output/outlier_group_summary.txt", split = TRUE)
cat(sprintf("High-PC1 group: PC1 > %g, %d of %d plots\n\n", CUT, sum(meta$PC1 > CUT), nrow(meta)))

# --- spectra -----------------------------------------------------------------
mean_by <- sapply(c("High-PC1 group", "Other plots"), function(k) colMeans(X[g == k, ]))
sd_by <- sapply(c("High-PC1 group", "Other plots"), function(k) apply(X[g == k, ], 2, sd))
at <- c(550, 680, 800, 1450, 1650, 1940, 2200)
cat("Mean reflectance at selected wavelengths (group / other / difference):\n")
print(data.frame(nm = at, group = round(mean_by[match(at, wav), 1], 4),
                 other = round(mean_by[match(at, wav), 2], 4),
                 diff = round(mean_by[match(at, wav), 1] - mean_by[match(at, wav), 2], 4)),
      row.names = FALSE)
sd_ <- data.frame(wavelength = rep(wav, 2), group = rep(colnames(mean_by), each = length(wav)),
                  mean = c(mean_by), sd = c(sd_by))
ps <- ggplot(sd_, aes(wavelength, mean, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("High-PC1 group" = "#eb6834", "Other plots" = "#2a78d6"), name = NULL) +
  scale_fill_manual(values = c("High-PC1 group" = "#eb6834", "Other plots" = "#2a78d6"), name = NULL) +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = "Mean leaf reflectance of high-PC1 plots and other plots, CLY25-B5 block 2",
       subtitle = sprintf("Splice-corrected plot means (5 scans per plot); high-PC1 group = PC1 > %g (n = %d), other plots n = %d; band = ±1 SD",
                          CUT, sum(g == "High-PC1 group"), sum(g == "Other plots"))) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave("output/outlier_group_spectra.png", ps, width = 9, height = 5.5, dpi = 200)

# --- field position ----------------------------------------------------------
pf <- ggplot(meta, aes(Range, Row, fill = group)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("High-PC1 group" = "#eb6834", "Other plots" = "#d9d9d9"), name = NULL) +
  scale_y_reverse() + coord_equal() +
  labs(title = "Field position of high-PC1 plots, CLY25-B5 block 2",
       subtitle = sprintf("One tile per scanned plot (n = %d); high-PC1 group = PC1 > %g; scanning proceeded in Range order, 2025-07-18 to 07-26",
                          nrow(meta), CUT)) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.grid = element_blank(),
        legend.position = "top", plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave("output/outlier_group_fieldmap.png", pf, width = 10, height = 6, dpi = 200)

# --- composition -------------------------------------------------------------
tab <- function(v, grp = meta$group) {
  t <- table(v, grp)
  data.frame(level = rownames(t), in_group = t[, 1], total = rowSums(t),
             pct_in_group = round(100 * t[, 1] / rowSums(t), 1), row.names = NULL)
}
cat("\nBy taxon (taxa_mex):\n"); print(tab(meta$taxa_mex), row.names = FALSE)
cat("\nBy scan date:\n"); print(tab(meta$date), row.names = FALSE)
cat("\nNotes in group:\n"); print(meta[meta$group == "High-PC1 group" & meta$notes != "", c("PlotID", "notes")], row.names = FALSE)

cat("\nPredicted_N and M-distance (median [IQR]):\n")
for (v in c("predN_mean", "Mdist_mean", "Mdist_max")) {
  q <- tapply(meta[[v]], meta$group, quantile, c(.25, .5, .75))
  cat(sprintf("  %-11s group %.2f [%.2f-%.2f]   other %.2f [%.2f-%.2f]   Wilcoxon p = %.2g\n", v,
              q[[1]][2], q[[1]][1], q[[1]][3], q[[2]][2], q[[2]][1], q[[2]][3],
              wilcox.test(meta[[v]] ~ meta$group)$p.value))
}

# --- clustering by accession -------------------------------------------------
lines <- meta[!meta$Genotype %in% c("B73", "Purple check"), ]
memb <- lines$group == "High-PC1 group"
stat <- function(m) var(tapply(m, lines$Accession, mean))
obs <- stat(memb)
null <- replicate(10000, stat(sample(memb)))
cat(sprintf("\nAccession clustering (lines only, %d lines, %d accessions): observed var of per-accession fraction %.4f, null median %.4f, permutation p = %.4f\n",
            nrow(lines), length(unique(lines$Accession)), obs, median(null), mean(null >= obs)))
acc <- tab(lines$Accession, lines$group)
acc <- acc[acc$total >= 5, ]
cat("Accessions (>= 5 lines) with the highest fraction in the group:\n")
print(head(acc[order(-acc$pct_in_group), ], 10), row.names = FALSE)
sink()

write.csv(meta[meta$group == "High-PC1 group",
               c("PlotID", "date", "Row", "Range", "Genotype", "Accession", "Taxa", "notes",
                 "predN_mean", "Mdist_mean", "Mdist_max", "PC1")],
          "output/outlier_group_plots.csv", row.names = FALSE)
