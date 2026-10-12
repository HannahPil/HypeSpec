# =============================================================================
# Shared loader for scripts 08-10: plot metadata and the three component sets
# -----------------------------------------------------------------------------
# Env:
#   TAG      suffix of the PCA / retrained-autoencoder run to read and of every
#            output written (e.g. _green). Default "" = all 921 plots.
#   EXCLUDE  CSV with a PlotID column; those plots are dropped from metadata and
#            every set (e.g. output/outlier_group_plots.csv for _green).
# The Tross encoder is fixed, so its LVs are the same per plot in every run and
# are only subset.
# =============================================================================
TAG <- Sys.getenv("TAG", "")
EXCLUDE <- Sys.getenv("EXCLUDE", "")
tg <- function(p) sub("\\.(csv|png)$", paste0(TAG, ".\\1"), p)

meta <- read.csv("output/bz_plot_meta.csv")
if (nzchar(EXCLUDE)) meta <- meta[!meta$PlotID %in% read.csv(EXCLUDE)$PlotID, ]

keep_plots <- function(d) d[d$PlotID %in% meta$PlotID, ]
sets <- lapply(list(
  PCA = read.csv(tg("output/pca_scores.csv")),
  `Tross encoder` = read.csv("output/tross_pretrained_lvs.csv"),
  `Retrained autoencoder` = read.csv(sprintf("output/ae_10lv%s/LVs_Reflectance_Data_bz_10lv%s.csv", TAG, TAG))
), keep_plots)
stopifnot(all(vapply(sets, nrow, 1L) == nrow(meta)))
cat(sprintf("component sets%s: %d plots\n", if (nzchar(TAG)) paste0(" (", TAG, ")") else "", nrow(meta)))
