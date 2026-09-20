###################################################################################################
## 5_0_make_consensus_peaks.R  --  Step 5.0: frozen consensus ACR set + blind-spot report
##
## Merges each stage's whole-object bulk MACS2 calls and emits the FROZEN feature space that
## Step 5 (meta-cells) and Step 6 (DAR) both consume.
##
## THE RULE: ONE peak set, frozen across all three stages. Stage-specific peak sets confound a
## "lost" peak with a peak-calling threshold shift -- the same logic as the frozen cluster config
## (plan decision B).
##
## PRIMARY = Pre. Pre defines the measurement space; wd and nd are measured IN it. Letting the
## treatments shape the feature space they are then evaluated in is the failure mode this project
## keeps catching. Second reason: at ~888 fragments/cell (At), adding near-empty features degrades
## the distance metric more than the marginal signal helps.
## The union is emitted too, for the reverse-direction SUPPLEMENTARY arm (start from wd/nd, go back
## to Pre) -- the sensitivity check on this very decision.
##
## NOT pseudobulk-per-cluster. DAR would normally want that, but Step 5 feeds the BOUNDARY TEST,
## which asks whether cluster boundaries are supported. A feature space built from cluster
## pseudobulks would make that test circular by construction. Bulk-only is REQUIRED, not convenient.
##
## No GenomicRanges dependency (bare local R): interval merge + overlap are done directly.
##
## Usage:
##   Rscript 5_0_make_consensus_peaks.R <outdir> <genome> <pre.narrowPeak> <wd.narrowPeak> <nd.narrowPeak>
###################################################################################################

suppressMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  stop("Usage: Rscript 5_0_make_consensus_peaks.R <outdir> <genome> <pre.np> <wd.np> <nd.np>")
outdir <- args[1]; genome <- args[2]
files  <- c(Pre = args[3], wd = args[4], nd = args[5])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_np <- function(f) {
  d <- fread(f, header = FALSE, sep = "\t", select = 1:3,
             col.names = c("chrom", "start", "end"))
  d[, chrom := as.character(chrom)][order(chrom, start)]
}

## merge overlapping / book-ended intervals within one set, per chromosome
merge_iv <- function(d) {
  setorder(d, chrom, start)
  out <- vector("list", length(unique(d$chrom))); k <- 0L
  for (ch in unique(d$chrom)) {
    x <- d[chrom == ch]; s <- x$start; e <- x$end
    ks <- s[1]; ke <- e[1]; S <- integer(0); E <- integer(0)
    if (length(s) > 1) for (i in 2:length(s)) {
      if (s[i] <= ke) ke <- max(ke, e[i]) else { S <- c(S, ks); E <- c(E, ke); ks <- s[i]; ke <- e[i] }
    }
    k <- k + 1L
    out[[k]] <- data.table(chrom = ch, start = c(S, ks), end = c(E, ke))
  }
  rbindlist(out)[order(chrom, start)]
}

## does each interval of A overlap ANY interval of B?
## per chrom: j = last B-start <= A-end ; overlap <=> cummax(B-end)[j] >= A-start
has_ov <- function(A, B) {
  res <- logical(nrow(A))
  for (ch in unique(A$chrom)) {
    ia <- which(A$chrom == ch); ib <- which(B$chrom == ch)
    if (!length(ib)) next
    bs <- B$start[ib]; be <- B$end[ib]
    o <- order(bs); bs <- bs[o]; M <- cummax(be[o])
    j <- findInterval(A$end[ia], bs)
    ok <- j > 0L
    res[ia][ok] <- M[j[ok]] >= A$start[ia][ok]
  }
  res
}

bpsum <- function(d) sum(as.numeric(d$end - d$start))
wr <- function(d, f) {
  d[, name := sprintf("%s_%d_%d", chrom, start, end)]
  fwrite(d[, .(chrom, start, end, name)], f, sep = "\t", col.names = FALSE, quote = FALSE)
  message(sprintf("   wrote %-34s %7d intervals  %6.1f Mb", basename(f), nrow(d), bpsum(d) / 1e6))
}

message("== ", genome)
P <- lapply(files, function(f) merge_iv(read_np(f)))
for (n in names(P))
  message(sprintf("   %-4s %7d merged intervals  %6.1f Mb", n, nrow(P[[n]]), bpsum(P[[n]]) / 1e6))

## blind-spot report: what does a Pre-only feature space miss?
rep <- list()
for (n in c("wd", "nd")) {
  ov  <- has_ov(P[[n]], P$Pre)
  ov2 <- has_ov(P$Pre, P[[n]])
  uq  <- P[[n]][!ov]
  w   <- if (nrow(uq)) uq$end - uq$start else NA_integer_
  message(sprintf("   %s: %d/%d (%.2f%%) %s-unique | %.2f Mb | width med %s max %s | reverse %d/%d (%.2f%%) Pre-unique",
                  n, sum(!ov), nrow(P[[n]]), 100 * mean(!ov), n, bpsum(uq) / 1e6,
                  ifelse(all(is.na(w)), "-", median(w)), ifelse(all(is.na(w)), "-", max(w)),
                  sum(!ov2), nrow(P$Pre), 100 * mean(!ov2)))
  rep[[n]] <- data.table(genome = genome, stage = n,
                         n_stage = nrow(P[[n]]), n_unique_vs_Pre = sum(!ov),
                         pct_unique_vs_Pre = round(100 * mean(!ov), 3),
                         mb_unique = round(bpsum(uq) / 1e6, 3),
                         width_median = ifelse(all(is.na(w)), NA, median(w)),
                         width_max = ifelse(all(is.na(w)), NA, max(w)),
                         n_Pre_unique_vs_stage = sum(!ov2),
                         pct_Pre_unique_vs_stage = round(100 * mean(!ov2), 3))
}

U <- merge_iv(rbindlist(P))
covPre <- has_ov(U, P$Pre)
message(sprintf("   UNION = %d intervals, %.1f Mb | Pre covers %.2f%% -> BLIND SPOT %.2f%% (%d)",
                nrow(U), bpsum(U) / 1e6, 100 * mean(covPre), 100 * mean(!covPre), sum(!covPre)))

wr(copy(P$Pre), file.path(outdir, sprintf("%s.consensus_Pre.bed", genome)))
wr(copy(U),     file.path(outdir, sprintf("%s.union_PreWdNd.bed", genome)))

rr <- rbindlist(rep)
rr[, `:=`(union_intervals = nrow(U),
          pct_union_covered_by_Pre = round(100 * mean(covPre), 3),
          blind_spot_pct = round(100 * mean(!covPre), 3))]
fwrite(rr, file.path(outdir, sprintf("%s.peak_overlap_report.tsv", genome)), sep = "\t")
message("   wrote ", genome, ".peak_overlap_report.tsv")
