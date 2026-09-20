#!/usr/bin/env Rscript
# =============================================================================
# figS7_fripnorm.R  -  peak-count-normalised FRiP for Fig S7 panel D (and the top-K control):
#   plate-split design, per genome, three stages (PreClean / WD / ND). Recomputes per-cell FRiP
#   from reads x peaks (data.table foverlaps), then rarefies FRiP as a function of #peaks.
# Inputs  <scratch>/{TAIR10,B73v5}_cells.txt + {TAIR10,B73v5}_{pre,wd,nd}_reads.tsv  (figS7_prep_beds.sh)
#         data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step0_qc/<base>_macs2_temp/<base>_peaks_combined_peaks.narrowPeak
# Output  figures/supplementary/figS7/S7_frip_{rarefaction,full_points,percell}[.top].tsv  (caches read by figS7.R)
# Run     Rscript analysis/supplementary/figS7_fripnorm.R <scratch_dir>                    (random draws, ~12 min)
#         PEAK_SELECT=top Rscript analysis/supplementary/figS7_fripnorm.R <scratch_dir>    (top-K control, ~4 min)
# =============================================================================
#
# WHY THIS EXISTS, AND WHY IT IS EXPECTED TO BE A NEAR-NO-OP HERE.
#   In the COMBINED-GENOME co-projection (main Fig 5 QC panels) the At peak set
#   collapses 186,774 -> 65,379 (-65%) on cleaning, so naive FRiP is dominated
#   by peak-set territory and no peak-matching rule is neutral (Fig 5F).
#   On this PLATE-SPLIT design the peak sets are nearly stable across stages:
#       At ->TAIR10  39,923 -> wd 39,381 (-1.4%) / nd 39,382 (-1.4%)
#       B73->B73v5  284,132 -> wd 283,554 (-0.2%) / nd 268,826 (-5.4%)
#   so the correction should be small. Showing that it IS small is the result:
#   it bounds the peak-count confound to the co-projection. The one arm where
#   it could bite is B73/nd.
#
# wd and nd are INDEPENDENT treatments of the same raw input, NOT a chain.
# Cells are the set SHARED by all three stages, so every comparison is exactly
#   paired. The stage cell sets are NOT nested (cleaning both drops and GAINS
#   cells), so an unpaired comparison would confound QC change with cell-set
#   change -- on At/nd that flips the apparent direction of pTSS. See figS7.R.
# The scratch reads are filtered against the PRE-GATE metadata (v1). v6 keeps
#   only qc_check == 1 cells, which would make the QC comparison circular -- see
#   the header of figS7.R. Re-run the prep with MDVER=v1 if in any doubt.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
SC <- if (length(args) >= 1) args[1] else
  stop("usage: figS7_fripnorm.R <scratch_dir>  (run figS7_prep_beds.sh first)")

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SOC    <- file.path(DATA, "socrates")
QC     <- file.path(SOC, "SM2v2_plate", "step0_qc")
OUTDIR <- "figures/supplementary/figS7"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

R_DRAWS <- 12        # peak-subsampling repeats per K
MAXCELL <- 10000     # cap cells per genome (B73 has ~18.3k shared); disclosed on the panel

# HOW PEAKS ARE SUBSAMPLED IS A REAL METHODOLOGICAL CHOICE, NOT A DETAIL.
#   "random": draw K of N uniformly. This is what the earlier combined-genome QC
#     panel does. It answers "what if we threw away N-K peaks at random?" -- but
#     that is NOT the counterfactual we want. A stage that had genuinely called
#     only K peaks would have MACS2's K STRONGEST, which capture more reads than a
#     random K. So random subsampling UNDERSTATES the larger-peak-set stage and
#     biases the comparison toward whichever stage has the FEWEST peaks (it is the
#     anchor and is never subsampled at all -- it draws K of K).
#   "top": take the K highest-scoring peaks (narrowPeak column 5). This is the
#     honest "if this stage had called only K peaks" counterfactual.
#   Set PEAK_SELECT=top in the environment to run the control. The bias matters
#   most where the peak sets differ most: here <=5.4%, but in the combined-genome
#   panels At is 186,774 vs 65,379 (-65%), where the sign flips with the rule.
PEAK_SELECT <- Sys.getenv("PEAK_SELECT", "random")
if (!PEAK_SELECT %in% c("random", "top")) stop("PEAK_SELECT must be 'random' or 'top'")
SUFFIX <- if (PEAK_SELECT == "random") "" else paste0(".", PEAK_SELECT)
cat(sprintf("[config] peak subsampling = %s%s\n", PEAK_SELECT,
            if (PEAK_SELECT == "top") "  (control run; writes *.top.tsv)" else ""))

STAGES <- c(pre = "SM2", wd = "Clean.SM2v2wd", nd = "Clean.SM2v2")
GEN <- list(
  TAIR10 = list(suf = "At_TAIR10",  label = "Arabidopsis (At) -> TAIR10"),
  B73v5  = list(suf = "B73_B73v5",  label = "Maize (B73) -> B73v5")
)

peak_file <- function(g, stg) {
  base <- sprintf("%s_%s", STAGES[[stg]], GEN[[g]]$suf)
  file.path(QC, paste0(base, "_macs2_temp"), paste0(base, "_peaks_combined_peaks.narrowPeak"))
}

# --- one (genome, stage): read x peak overlaps ---------------------------------
# `reads` is keyed in place (no copy) to keep peak memory down on the B73 arm;
# `rid` is assigned BEFORE the key so it survives the reorder, and totals are
# counted by barcode, which is order-independent.
build_stage <- function(readfile, peakfile, keep_bc) {
  reads <- fread(readfile, col.names = c("chrom", "start", "end", "bc"))
  reads <- reads[bc %in% keep_bc]
  reads[, rid := .I]
  total <- reads[, .(total = .N), by = bc]
  pk <- fread(peakfile, select = c(1, 2, 3, 5),
              col.names = c("chrom", "start", "end", "score"))
  pk[, pid := .I]
  score_by_pid <- pk$score      # pid == .I here, so this vector is pid-indexed;
  setkey(pk, chrom, start, end) # capture it BEFORE setkey reorders the table
  setkey(reads, chrom, start, end)
  ov <- foverlaps(reads, pk, type = "any", which = TRUE, nomatch = NULL)
  # keep ALL overlaps: a read over >1 peak must stay in-peak if ANY sampled peak
  # covers it. Using only the first peak biases the large-peak-set arm down.
  hits <- unique(data.table(rid = reads$rid[ov$xid], bc = reads$bc[ov$xid], pid = ov$yid),
                 by = c("rid", "pid"))
  list(hits = hits, total = total, npeak = nrow(pk), bp = pk[, sum(end - start)],
       score = score_by_pid)
}

per_cell <- function(st, pids_in) {
  acr <- unique(st$hits[pids_in[pid], .(rid, bc)])[, .(acrs = .N), by = bc]
  f <- merge(st$total, acr, by = "bc", all.x = TRUE)
  f[is.na(acrs), acrs := 0][, FRiP := acrs / total]
  f[]
}
draw_K <- function(st, K) {
  s <- logical(st$npeak)
  if (PEAK_SELECT == "top") s[order(st$score, decreasing = TRUE)[seq_len(K)]] <- TRUE
  else                      s[sample.int(st$npeak, K)] <- TRUE
  s
}
# "top" is deterministic, so repeating it just recomputes the same answer
N_DRAWS <- function() if (PEAK_SELECT == "top") 1L else R_DRAWS
med_at_K <- function(st, K, R = N_DRAWS())
  vapply(seq_len(R), function(i) median(per_cell(st, draw_K(st, K))$FRiP), numeric(1))
percell_at_K <- function(st, K, R = N_DRAWS())
  rbindlist(lapply(seq_len(R), function(i) per_cell(st, draw_K(st, K))[, .(bc, FRiP)]
              ))[, .(frip_norm = mean(FRiP)), by = bc]

rare_all <- list(); pts_all <- list(); pc_all <- list()

for (g in names(GEN)) {
  cells <- fread(file.path(SC, sprintf("%s_cells.txt", g)), header = FALSE)$V1
  n_shared <- length(cells)
  if (n_shared > MAXCELL) cells <- sort(sample(cells, MAXCELL))
  message(sprintf("[%s] %d shared cells (%d used)", g, n_shared, length(cells)))

  # peak counts first: K_common is the min across the three stages
  npk <- vapply(names(STAGES), function(s) {
    f <- peak_file(g, s); if (!file.exists(f)) stop("missing peaks: ", f)
    as.integer(nrow(fread(f, select = 1, col.names = "chrom")))
  }, integer(1))
  Kc <- min(npk)
  cat(sprintf("  peaks: pre %s | wd %s | nd %s  -> common K = %s\n",
              format(npk["pre"], big.mark = ","), format(npk["wd"], big.mark = ","),
              format(npk["nd"], big.mark = ","), format(Kc, big.mark = ",")))

  # K grid for the curve, shared across stages so the three lines are comparable
  grid <- sort(unique(pmin(c(2000, 5000, 10000, 20000, 40000, 80000, 150000, 250000, Kc), Kc)))

  for (s in names(STAGES)) {
    st <- build_stage(file.path(SC, sprintf("%s_%s_reads.tsv", g, s)), peak_file(g, s), cells)
    naive <- per_cell(st, rep(TRUE, st$npeak))          # own full peak set
    pcn   <- percell_at_K(st, Kc)                        # normalized to common K
    rare_all[[paste(g, s)]] <- rbindlist(lapply(grid, function(K) {
      m <- med_at_K(st, K); data.table(genome = g, stage = s, K = K, med = mean(m), sd = sd(m))
    }))
    pts_all[[paste(g, s)]] <- data.table(
      genome = g, stage = s, npeak = st$npeak, bp = st$bp,
      med_naive = median(naive$FRiP), Kcommon = Kc)
    pc_all[[paste(g, s)]] <- merge(naive[, .(bc, total, frip_naive = FRiP)], pcn, by = "bc")[
      , `:=`(genome = g, stage = s)]
    cat(sprintf("    %-3s  %s pk / %.1f Mbp | naive FRiP %.4f | norm@K %.4f | %d cells\n",
                s, format(st$npeak, big.mark = ","), st$bp / 1e6,
                median(naive$FRiP), median(pcn$frip_norm), nrow(naive)))
    rm(st); invisible(gc(FALSE))
  }
}

rare <- rbindlist(rare_all); pts <- rbindlist(pts_all); pc <- rbindlist(pc_all)
fwrite(rare, file.path(OUTDIR, sprintf("S7_frip_rarefaction%s.tsv", SUFFIX)), sep = "\t")
fwrite(pts,  file.path(OUTDIR, sprintf("S7_frip_full_points%s.tsv", SUFFIX)), sep = "\t")
fwrite(pc,   file.path(OUTDIR, sprintf("S7_frip_percell%s.tsv",     SUFFIX)), sep = "\t")

cat("\n### FRiP at the common peak count (each stage subsampled to K_common)\n")
print(merge(rare, unique(pts[, .(genome, Kcommon)]), by = "genome")[K == Kcommon][order(genome, stage)])
cat("\n[done] wrote S7_frip_{rarefaction,full_points,percell}", SUFFIX, ".tsv to ", OUTDIR, "/\n", sep = "")
