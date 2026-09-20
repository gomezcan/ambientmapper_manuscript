###################################################################################################
## 5_1_build_acr_matrix.R   --   Step 5.1: consensus-ACR x cell matrix (SEACells input)
##
## Builds an ACR x cell sparse matrix of raw Tn5 insertion counts against a FROZEN consensus peak
## set. Mirrors 4_0_build_gene_body_matrix.R exactly -- same Tn5 BED, same cellID convention, same
## sparseMatrix accumulation -- ONLY the feature set changes (consensus ACRs instead of gene windows).
##
## WHY THIS EXISTS: the Socrates embedding uses FIXED 500 bp TILES (1_1:461 -> peaks=F) and the
## annotation uses gene body +/-500 bp, so the called MACS2 peaks have never been used as features
## anywhere -- they feed FRiP/QC only.
##
## NO GenomicRanges DEPENDENCY. The consensus BED is pre-merged, hence its intervals are DISJOINT,
## so each 1 bp Tn5 site falls in at most ONE peak. That turns the overlap into a per-chromosome
## interval lookup (findInterval on sorted starts), which is both exact and faster than
## findOverlaps -- and runs on the bare local R (Matrix + data.table only).
##
## Usage:
##   Rscript 5_1_build_acr_matrix.R <out_prefix> <peaks.bed> <tn5.bed.gz> [meta.txt]
##     out_prefix : matrix written to <out_prefix>.acrs.sparse.rds
##     peaks.bed  : BED4 consensus set from SM2v2_plate/_consensus_peaks/ (MERGED = disjoint)
##     tn5.bed.gz : Tn5 BED (chr, start, end, cellID, strand)
##     meta.txt   : OPTIONAL clustered metadata; if given, columns are restricted to its cellIDs
##                  (v6 QC-passing cells). STRONGLY RECOMMENDED -- the BED carries every barcode,
##                  including tens of thousands of non-cells.
##
## Env:
##   MIN_CELLS_FRAC : feature kept if detected in >= this FRACTION of cells   (default 0.005)
##   MIN_CELLS_ABS  : absolute floor on the same filter                       (default 10)
##       kept  <=>  n_cells_detected >= max(MIN_CELLS_ABS, MIN_CELLS_FRAC * ncol)
##
##   The filter is FRACTIONAL by design. Socrates' cleanData(min.c=50) is an ABSOLUTE cell count.
##   MEASURED on the consensus ACR matrices: at t=50 At retains 157 / 39,923 ACRs (0.4%)
##   while maize retains 42,706 / 284,132 (15.0%) -- a 272x difference in surviving feature count at
##   the SAME threshold. Not merely "stricter on the minor genome": near-annihilating for At
##   (the min.c=50 failure mode). A retention sweep is printed at several thresholds so
##   the choice is VISIBLE, not silent: this is a bandwidth-like parameter and must be declared and
##   swept, not tuned to taste.
##
##   FEATURE_WHITELIST : file of ACR names (one per line), typically <Pre>.acrs.features.txt
##     THE FEATURE SET MUST BE FROZEN ACROSS STAGES. Applying the detection filter per object makes
##     the feature set stage-dependent -- measured: At kept 7,105 (Pre) / 6,917 (wd) / 8,684 (nd),
##     a 25% swing, purely because nd cells are ~65% deeper so each ACR clears the threshold in more
##     cells. That would make the meta-cell distance metric a function of stage: exactly the confound
##     the frozen peak set exists to prevent. So: build Pre WITHOUT a whitelist, then pass its
##     emitted .acrs.features.txt when building wd and nd. Same rule as the frozen cluster config
##     (plan decision B) and the Pre-primary consensus peak set.
###################################################################################################

suppressMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: Rscript 5_1_build_acr_matrix.R <out_prefix> <peaks.bed> <tn5.bed.gz> [meta.txt]")

out_prefix <- args[1]
peaks_f    <- args[2]
tn5_f      <- args[3]
meta_f     <- if (length(args) >= 4) args[4] else NA_character_

MIN_FRAC <- as.numeric(Sys.getenv("MIN_CELLS_FRAC", "0.005"))
MIN_ABS  <- as.integer(Sys.getenv("MIN_CELLS_ABS",  "10"))
WHITELIST <- Sys.getenv("FEATURE_WHITELIST", "")

message(" - out_prefix = ", out_prefix)
message(" - peaks      = ", basename(peaks_f))
message(" - tn5 bed    = ", basename(tn5_f))
message(" - metadata   = ", if (is.na(meta_f)) "(none -- ALL barcodes kept)" else basename(meta_f))
message(" - filter     = n_cells >= max(", MIN_ABS, ", ", MIN_FRAC, " * ncol)")

## (1) consensus peaks ------------------------------------------------------------
pk <- fread(peaks_f, header = FALSE, sep = "\t", select = 1:4,
            col.names = c("chrom", "start", "end", "name"))
pk[, chrom := as.character(chrom)]
setorder(pk, chrom, start)
pk[, fid := .I]                                   # row index in the output matrix
n_feat <- nrow(pk)
message("   consensus ACRs = ", n_feat, " (", round(sum(as.numeric(pk$end - pk$start)) / 1e6, 1), " Mb)")

## Disjointness is the assumption that licenses the findInterval shortcut. Verify, do not assume.
bad <- pk[, sum(start[-1] < end[-.N]), by = chrom]$V1
if (sum(bad) > 0)
  stop("peaks.bed is NOT disjoint (", sum(bad), " overlapping pairs) -- re-merge before use")

## split peaks by chromosome for the lookup
pk_by <- split(pk[, .(start, end, fid)], pk$chrom)

## (2) optional cell whitelist -----------------------------------------------------
keep_cells <- NULL
if (!is.na(meta_f)) {
  md <- read.table(meta_f, header = TRUE, sep = "\t", row.names = 1,
                   check.names = FALSE, comment.char = "", quote = "")
  keep_cells <- rownames(md)
  message("   whitelist  = ", length(keep_cells), " QC-passing cells")
}

## (3) Tn5 BED -> (feature, cell) pairs ---------------------------------------------
## Only chr/start/cell are needed: a Tn5 site is 1 bp, so `start` (0-based) IS the position.
message(" - reading Tn5 BED ...")
bed <- fread(cmd = paste("gzip -dc", shQuote(tn5_f)), header = FALSE, sep = "\t",
             select = c(1, 2, 4), col.names = c("chrom", "pos", "cell"))
bed[, chrom := as.character(chrom)]
message("   insertions = ", format(nrow(bed), big.mark = ","))

if (!is.null(keep_cells)) {
  bed <- bed[cell %chin% keep_cells]
  message("   after cell whitelist = ", format(nrow(bed), big.mark = ","))
  if (!nrow(bed)) stop("no insertions left after whitelisting -- cellID convention mismatch?")
}

## per-chromosome interval lookup on DISJOINT peaks:
##   j = last peak whose start <= pos ; in-peak  <=>  j > 0 AND pos < end[j]
bed[, fid := NA_integer_]
for (ch in intersect(names(pk_by), unique(bed$chrom))) {
  idx <- which(bed$chrom == ch)
  if (!length(idx)) next
  P <- pk_by[[ch]]
  j <- findInterval(bed$pos[idx], P$start)
  ok <- j > 0L
  ok[ok] <- bed$pos[idx][ok] < P$end[j[ok]]
  bed$fid[idx[ok]] <- P$fid[j[ok]]
}
hits <- bed[!is.na(fid)]
message("   in-peak insertions = ", format(nrow(hits), big.mark = ","),
        "  (FRiP-equivalent ", round(100 * nrow(hits) / nrow(bed), 1), "%)")
rm(bed); gc()

## (4) sparse matrix ----------------------------------------------------------------
cf <- factor(hits$cell)
## sparseMatrix() sums x over duplicated (i,j) pairs -> per-(ACR, cell) insertion counts
m <- sparseMatrix(i = hits$fid, j = as.integer(cf), x = 1,
                  dims = c(n_feat, nlevels(cf)),
                  dimnames = list(pk$name, levels(cf)))
rm(hits, cf); gc()
message(" - raw matrix: ", nrow(m), " ACRs x ", ncol(m), " cells (nnz = ",
        format(length(m@x), big.mark = ","), ")")

## (5) feature filter, with the retention sweep printed ------------------------------
det <- Matrix::rowSums(m > 0)                      # n cells in which the ACR is detected
thr <- max(MIN_ABS, ceiling(MIN_FRAC * ncol(m)))
message(" - detection-per-ACR quantiles: ",
        paste(sprintf("%s=%d", c("min", "25%", "50%", "75%", "90%", "max"),
                      as.integer(quantile(det, c(0, .25, .5, .75, .9, 1)))), collapse = "  "))
message(" - retention sweep (n ACRs surviving at threshold t):")
for (t in sort(unique(c(1L, 5L, 10L, 20L, 50L, 100L, as.integer(thr)))))
  message(sprintf("     t = %4d cells (%5.2f%% of cells) -> %7d ACRs (%5.1f%%)%s",
                  t, 100 * t / ncol(m), sum(det >= t), 100 * mean(det >= t),
                  if (t == thr) "   <-- APPLIED" else ""))

## A whitelist OVERRIDES the computed threshold: the feature space must be identical across stages
## (see header). Pre emits the list; wd/nd consume it.
if (nzchar(WHITELIST)) {
  wl <- readLines(WHITELIST)
  miss <- setdiff(wl, rownames(m))
  if (length(miss))
    stop("whitelist has ", length(miss), " ACRs absent from this matrix -- peak sets differ?")
  m <- m[wl, , drop = FALSE]
  message(" - FEATURE_WHITELIST applied: ", nrow(m), " ACRs from ", basename(WHITELIST),
          "  (computed threshold ", thr, " would have kept ", sum(det >= thr), " -- OVERRIDDEN)")
} else {
  m <- m[det >= thr, , drop = FALSE]
  message(" - threshold applied: ", nrow(m), " ACRs (no whitelist)")
}
message(" - final matrix: ", nrow(m), " ACRs x ", ncol(m), " cells (nnz = ",
        format(length(m@x), big.mark = ","), ")")

## (6) save --------------------------------------------------------------------------
saveRDS(m, file = paste0(out_prefix, ".acrs.sparse.rds"))
message(" - wrote ", out_prefix, ".acrs.sparse.rds")
writeLines(rownames(m), paste0(out_prefix, ".acrs.features.txt"))
message(" - wrote ", out_prefix, ".acrs.features.txt  (pass as FEATURE_WHITELIST to the other stages)")

stats <- data.frame(
  out_prefix = out_prefix, peaks = basename(peaks_f), tn5 = basename(tn5_f),
  consensus_acrs = n_feat, cells = ncol(m), acrs_kept = nrow(m),
  min_cells_frac = MIN_FRAC, min_cells_abs = MIN_ABS, threshold_applied = thr,
  nnz = length(m@x), mean_acrs_per_cell = round(mean(Matrix::colSums(m > 0)), 1),
  mean_counts_per_cell = round(mean(Matrix::colSums(m)), 1))
write.table(stats, paste0(out_prefix, ".acrs.stats.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message(" - wrote ", out_prefix, ".acrs.stats.tsv")
