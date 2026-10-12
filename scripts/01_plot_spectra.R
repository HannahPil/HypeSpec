#!/usr/bin/env Rscript
# =============================================================================
# BZ block 2 leaf spectra -> one splice-corrected reflectance spectrum per plot
# -----------------------------------------------------------------------------
# Builds the autoencoder/PCA input in the format the Tross et al. 2024 code
# expects (scripts/tross2024): one row per plot, PlotID, trt, X350..X2500.
#
# PIPELINE
#   1. scan list from the BZ-remap sheet of the collection workbook; plots
#      scanned twice (985, 1441: measurements 1-10) keep the latter set (6-10)
#   2. read each scan's .spc (Indico Pro "Asd to Spc" output, log(1/R),
#      350-2500 nm, 1 nm) and back-transform: R = 10^-y. Cross-checked against
#      DN / white reference from the .asd files: identical.
#   3. splice correction at the detector joins (1000 nm VNIR/SWIR1, 1800 nm
#      SWIR1/SWIR2, per the .asd header) with prospectr::spliceCorrection.
#      Tross's published spectra show no steps at these joins and ours do; the
#      paper does not say how they were removed.
#   4. plain mean of the scans in each plot (Tross / Wijewardane et al. 2023)
#   5. QC as in Tross: flag plots with any reflectance < 0 or > 1
#
# Predicted_N and M_Distance are joined at PLOT level only: the predicted-N
# file is sorted by plot, not scan order, so scans cannot be matched one to one.
# They are carried as ground truth / a "looks wrong" flag, never as a filter.
# Field notes are carried too, for filtering experiments later; none applied.
#
# In : data/CLY25-B5-octexp-spectrum-collection.xlsx (sheets BZ-remap, BZ-all)
#      data/B5_block2_predictedN.csv
#      $SPC_DIR (default ../2025 spectra/BZ) — the per-scan .spc files
# Out: output/bz_plot_reflectance.csv.gz  (PlotID, trt = scan date, X350..X2500)
#      output/bz_plot_meta.csv             (one row per plot, see bottom)
# Run: Rscript scripts/01_plot_spectra.R   (from the repo root)
# =============================================================================

WB <- "data/CLY25-B5-octexp-spectrum-collection.xlsx"
PREDN <- "data/B5_block2_predictedN.csv"
SPC_DIR <- Sys.getenv("SPC_DIR", "../2025 spectra/BZ")
OUT_SPEC <- "output/bz_plot_reflectance.csv.gz"
OUT_META <- "output/bz_plot_meta.csv"
WAV <- 350:2500
SPLICE <- c(1000, 1800)
dir.create("output", showWarnings = FALSE)

# --- 1. scan list ------------------------------------------------------------
remap <- as.data.frame(readxl::read_excel(WB, sheet = "BZ-remap"))
notes <- as.data.frame(readxl::read_excel(WB, sheet = "BZ-all"))
n_per_plot <- table(remap$plot)
double <- as.numeric(names(n_per_plot[n_per_plot == 10]))
stopifnot(all(n_per_plot %in% c(5, 10)))
scans <- remap[!(remap$plot %in% double & remap$measurement <= 5), ]
stopifnot(all(table(scans$plot) == 5))
scans$spc <- file.path(SPC_DIR, sub("\\.asd$", ".spc", scans$spectra_name.asd))
cat(sprintf("%d scans in remap; %d kept (plots scanned twice: %s, latter set kept)\n",
            nrow(remap), nrow(scans), paste(double, collapse = ", ")))

# --- 2. read .spc (GRAMS SPC new format, single evenly spaced spectrum) -----
read_spc <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  stopifnot(b[2] == as.raw(0x4B), b[4] == as.raw(0x80)) # new format, IEEE float
  npts <- readBin(b[5:8], "integer", 1, size = 4)
  x <- readBin(b[9:24], "double", 2, size = 8)
  stopifnot(npts == length(WAV), x[1] == min(WAV), x[2] == max(WAV))
  readBin(b[545:(544 + 4 * npts)], "double", npts, size = 4) # 512 hdr + 32 subhdr
}
logir <- t(vapply(scans$spc, read_spc, numeric(length(WAV)), USE.NAMES = FALSE))
R <- 10^(-logir)
colnames(R) <- paste0("X", WAV)

# --- 3. splice correction ----------------------------------------------------
step_ratio <- function(M) { # median |step| at each join / median |step| nearby
  sapply(SPLICE, function(s) {
    j <- which(WAV == s)
    nb <- abs(M[, (j - 9):(j + 11)] - M[, (j - 10):(j + 10)])[, -11] # w -> w+1, minus the join
    median(abs(M[, j + 1] - M[, j])) / median(nb)
  })
}
before <- step_ratio(R)
R <- prospectr::spliceCorrection(R, WAV, splice = SPLICE)
after <- step_ratio(R)
cat(sprintf("join step / neighbour step at %s nm: before %s, after %s\n",
            paste(SPLICE, collapse = "/"), paste(round(before, 1), collapse = "/"),
            paste(round(after, 1), collapse = "/")))

# --- 4. plot means -----------------------------------------------------------
plots <- unique(scans$plot) # field collection order
P <- t(vapply(plots, function(p) colMeans(R[scans$plot == p, , drop = FALSE]),
              numeric(length(WAV))))
date <- tapply(scans$date, scans$plot, unique)[as.character(plots)]
stopifnot(lengths(date) == 1)

# --- 5. QC -------------------------------------------------------------------
qc_pass <- apply(P, 1, function(r) all(r >= 0 & r <= 1))
cat(sprintf("%d plots; %d fail the 0-1 reflectance check\n", length(plots), sum(!qc_pass)))

# --- outputs -----------------------------------------------------------------
spec <- data.frame(PlotID = plots, trt = unlist(date), P, check.names = FALSE)
con <- gzfile(OUT_SPEC, "w")
write.csv(spec[qc_pass, ], con, row.names = FALSE)
close(con)

pn <- read.csv(PREDN)
agg <- function(fun, col) tapply(pn[[col]], pn$Plot, fun)[as.character(plots)]
first <- function(col) agg(function(v) v[1], col)
note <- notes$notes[match(plots, notes$plot)]
meta <- data.frame(
  PlotID = plots, date = unlist(date), Row = first("Row"), Range = first("Range"),
  Genotype = first("Genotype"), Accession = first("Accession"), Taxa = first("Taxa"),
  notes = ifelse(is.na(note), "", note), n_scans = 5L,
  predN_mean = agg(mean, "Predicted_N"), Mdist_mean = agg(mean, "M_Distance"),
  Mdist_max = agg(max, "M_Distance"), qc_pass = qc_pass
)
write.csv(meta, OUT_META, row.names = FALSE)
cat(sprintf("wrote %s (%d plots) and %s\n", OUT_SPEC, sum(qc_pass), OUT_META))
