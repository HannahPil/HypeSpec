# =============================================================================
# Shared helpers for the R scripts (source("scripts/lib.R") from the repo root)
# -----------------------------------------------------------------------------
# Env:
#   TAG      suffix for PCA / retrained-autoencoder inputs and every output
#            (e.g. _green). Default "" = all 921 plots.
#   EXCLUDE  CSV with a PlotID column; those plots are dropped
#            (e.g. output/outlier_group_plots.csv for _green).
# =============================================================================
suppressMessages(library(SpATS))

TAG <- Sys.getenv("TAG", "")
EXCLUDE <- Sys.getenv("EXCLUDE", "")
tg <- function(p) sub("\\.(csv|png)$", paste0(TAG, ".\\1"), p)
excluded <- function() if (nzchar(EXCLUDE)) read.csv(EXCLUDE)$PlotID else integer(0)

# --- leaf reflectance features by wavelength (nm), [lo, hi) ------------------
# After Curran 1989, Ustin et al. 2009, Ceccato et al. 2001.
FEATURES <- data.frame(
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
feature_at <- function(w) FEATURES$label[w >= FEATURES$lo & w < FEATURES$hi]

rho_p <- function(x, v) {
  t <- suppressWarnings(cor.test(x, v, method = "spearman", exact = FALSE))
  c(rho = unname(t$estimate), p = t$p.value)
}

# --- SpATS: the recipe of ZeaL_QTL/scripts/zeal_predn_pheno.R ----------------
# d0: one row per plot with Genotype, Row, Range; y: trait in the same order.
spats_fit <- function(d0, y, random) {
  d <- data.frame(Genotype = factor(d0$Genotype), Row = as.numeric(d0$Row),
                  Range = as.numeric(d0$Range), y = y)
  d$Rf <- factor(d$Range)
  d$Cf <- factor(d$Row)
  SpATS(response = "y", genotype = "Genotype", genotype.as.random = random,
        spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf,
        data = d, control = list(monitoring = 0, maxit = 50))
}
# generalized heritability (genotype random)
spats_h2 <- function(d0, y) unname(getHeritability(spats_fit(d0, y, TRUE)))
# genotype BLUEs (genotype fixed) written in ZeaL_QTL phenotype format:
# <dir>/pheno_<trait>_cly_blue.csv with Genotype, <TRAIT>_cly25, <TRAIT>_CLY_mean
write_blues <- function(d0, y, trait, dir) {
  pr <- predict(spats_fit(d0, y, FALSE), which = "Genotype")
  vc <- grep("predicted", names(pr), ignore.case = TRUE, value = TRUE)[1]
  pr <- pr[!is.na(pr[[vc]]), ]
  out <- data.frame(Genotype = as.character(pr$Genotype), v = pr[[vc]])
  out <- out[order(out$Genotype), ]
  names(out)[2] <- paste0(trait, "_cly25")
  out[[paste0(trait, "_CLY_mean")]] <- out[[2]]
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(out, sprintf("%s/pheno_%s_cly_blue.csv", dir, tolower(trait)), row.names = FALSE)
  nrow(out)
}

# --- component sets: plot metadata + PCA, Tross encoder, retrained AE --------
# The Tross encoder is fixed, so its LVs are the same per plot in every run.
load_component_sets <- function() {
  meta <- read.csv("output/bz_plot_meta.csv")
  meta <- meta[!meta$PlotID %in% excluded(), ]
  keep <- function(d) d[d$PlotID %in% meta$PlotID, ]
  sets <- lapply(list(
    PCA = read.csv(tg("output/pca_scores.csv")),
    `Tross encoder` = read.csv("output/tross_pretrained_lvs.csv"),
    `Retrained autoencoder` = read.csv(sprintf("output/ae_10lv%s/LVs_Reflectance_Data_bz_10lv%s.csv", TAG, TAG))
  ), keep)
  stopifnot(all(vapply(sets, nrow, 1L) == nrow(meta)))
  cat(sprintf("component sets%s: %d plots\n", if (nzchar(TAG)) paste0(" (", TAG, ")") else "", nrow(meta)))
  list(meta = meta, sets = sets)
}
