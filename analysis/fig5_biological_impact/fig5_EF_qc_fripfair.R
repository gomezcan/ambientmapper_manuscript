#!/usr/bin/env Rscript
# =============================================================================
# fig5_EF_qc_fripfair.R  -  Fig 5 panel F: peak-matched per-cell FRiP under two matching rules
#   (random K peaks / top K peaks), PreClean vs PostClean per species, on the combined-genome
#   co-projection (SM2 -> Clean.SM2v2) fixed cell set.
# Inputs  <scratch>/{at,b73}_{pre,post}_reads.tsv   (written by fig5_EF_prefilter_beds.sh)
#         data/processed/scifiATAC_B73_Arabidopsis/socrates/_data/_PeakFiles/{SM2,Clean.SM2v2}_{At,B73}/*_peaks.narrowPeak
# Output  figures/main/fig5/Fig5_P1_QC_FRiPfair.{pdf,png} + 7 Fig5_P1_QC_FRiPfair_*.tsv
# Run     bash analysis/fig5_biological_impact/fig5_EF_prefilter_beds.sh <scratch>   (heavy, once)
#         Rscript analysis/fig5_biological_impact/fig5_EF_qc_fripfair.R <scratch>    (from the repo root)
# =============================================================================
#
# WHY THIS EXISTS
#   The earlier panel (fig5_EF_qc.R -> Fig5_P1_QC_FRiPnorm) rarefies to the COMMON peak count
#   K = min(N_pre, N_post) using a RANDOM draw of K peaks. Two things are wrong with that as a
#   matched comparison:
#   1. THE ANCHOR IS NEVER SUBSAMPLED. On At, Post has exactly 65,379 peaks and K = 65,379, so
#      Post draws K of K (its full set) while Pre throws away 65% of its 186,774 peaks. A
#      "+0.162 at equal peak count" compares Pre's RANDOM 65,379 peaks against Post's ACTUAL
#      BEST 65,379.
#   2. LOWERING K DOES NOT FIX IT. Measured on the existing curve, Post/Pre is ~1.6x at EVERY K
#      (1.42 at 5k, 1.66 at 10k, 1.59 at 50k, 1.61 at 65,379), because K is a different SHARE
#      of each peak set (10k = 15% of Post's set but 5% of Pre's). A random draw keeps favouring
#      the smaller set at all K. Maize is the proof: its peak sets are near-equal (325,754 vs
#      320,580) and the Post-Pre gap is 0.0003-0.008 at every K.
#
# THE RESULT: THERE IS NO UNBIASED CHOICE, AND THE ANSWER FLIPS.
#   `top` is NOT a fix. It carries the OPPOSITE bias:
#     random  Pre draws 65,379 of its 186,774 peaks (a 35% sample, so it keeps ~35% of its
#             in-peak reads) while Post draws 65,379 of 65,379 and keeps 100%.
#                                => FAVOURS THE SMALLER PEAK SET (Post).
#     top     Pre picks its BEST 65,379 out of 186,774 candidates (its top 35%) while Post must
#             use all 65,379 including its weakest.
#                                => FAVOURS THE LARGER PEAK SET (Pre).
#   The two rules therefore BRACKET the truth -- and on At they bracket ZERO:
#     At, K = 65,379:   random  +0.161      top  -0.068
#   So the sign of the At FRiP effect is NOT IDENTIFIABLE from these data. A "+0.162 gain at
#   equal peak count" is one end of that bracket, not a measurement.
#   The `common` control does not rescue it either, because neither peak set is neutral: each
#   was CALLED ON the reads of one stage, so it favours that stage. Scored on Post's peaks At
#   reads give +0.259; on Pre's peaks, -0.110.
#   Maize is the clean negative control throughout: its two peak sets are near-equal
#   (325,754 vs 320,580), so every rule agrees at ~0.00 and there is nothing to identify.
#   WHAT WOULD FIX IT: an EXTERNAL reference peak set, called on neither stage (a published ACR
#   atlas, or peaks from an independent dataset). That is exactly what makes pTSS confound-free
#   -- a fixed annotation -- which is why the pTSS panel (E), not this one, carries the
#   data-quality claim.
#
# THREE COMPARISONS, none neutral, reported together so the spread is visible:
#   random  K of N, uniform         -- the earlier method; biased toward Post
#   top     K highest-scoring       -- biased toward Pre
#   common  both stages on ONE peak set (Pre's, then Post's) -- biased toward whichever stage
#           authored the peaks
#
# NOT a circularity risk (unlike the Fig S7 step0_qc trap): the `fixedSet` cell list is frozen
#   with thresholds derived from PRE (mu/sd_pTSS_pre) and applied identically to both stages,
#   and the loader filters on `present`, not on `qc_check`.
# Significance, if ever added, is ACROSS CELLS -- never across peak draws. The draws are
#   resamples of one dataset; n would be a number we chose.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
SC <- if (length(args) >= 1) args[1] else
  stop("usage: fig5_EF_qc_fripfair.R <scratch_dir>   (run fig5_EF_prefilter_beds.sh first)")

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SOC    <- file.path(DATA, "socrates")
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ONE stem for every output, so the figure and its data sort together in OUTDIR:
# Part-1 outputs share the stem `Fig5_P1*`.
STEM    <- "Fig5_P1_QC_FRiPfair"
K_FAIR  <- 10000                                   # headline matched peak count
K_GRID  <- c(2000, 5000, 10000, 20000, 40000, 65379)
R_DRAWS <- 12                                      # random draws (top-K is deterministic)
SP_LEV  <- c("Maize (B73)", "Arabidopsis (At)")    # maize first, as in the assembled panel
ST_LEV  <- c("Pre", "Post")
ST_COL  <- c(Pre = "#FF83FA", Post = "#43CD80")    # Fig 3 stage colours

PK <- function(...) file.path(SOC, "_data/_PeakFiles", ...)
SPP <- list(
  B73 = list(prefix = "b73", label = "Maize (B73)",
             pre = PK("SM2_B73/SM2_B73_peaks.narrowPeak"),
             post = PK("Clean.SM2v2_B73/Clean.SM2v2_B73_peaks.narrowPeak")),
  At  = list(prefix = "at",  label = "Arabidopsis (At)",
             pre = PK("SM2_At/SM2_At_peaks.narrowPeak"),
             post = PK("Clean.SM2v2_At/Clean.SM2v2_At_peaks.narrowPeak"))
)

# --- reads x peaks -------------------------------------------------------------
# `score` is captured BEFORE setkey reorders the table, so it stays pid-indexed.
build <- function(prefix, stg, peakfile) {
  rf <- file.path(SC, sprintf("%s_%s_reads.tsv", prefix, stg))
  if (!file.exists(rf)) stop("missing prefiltered reads: ", rf,
                             "\n  run: bash analysis/fig5_biological_impact/fig5_EF_prefilter_beds.sh ", SC)
  reads <- fread(rf, col.names = c("chrom", "start", "end", "bc"))
  reads[, rid := .I]
  total <- reads[, .(total = .N), by = bc]
  pk <- fread(peakfile, select = c(1, 2, 3, 5),
              col.names = c("chrom", "start", "end", "score"))
  pk[, pid := .I]
  score_by_pid <- pk$score
  setkey(pk, chrom, start, end); setkey(reads, chrom, start, end)
  ov <- foverlaps(reads, pk, type = "any", which = TRUE, nomatch = NULL)
  # keep ALL overlaps: a read over >1 peak stays in-peak if ANY sampled peak
  # covers it. Using only the first peak biases the larger-peak-set arm down.
  hits <- unique(data.table(rid = reads$rid[ov$xid], bc = reads$bc[ov$xid], pid = ov$yid),
                 by = c("rid", "pid"))
  list(hits = hits, total = total, npeak = nrow(pk),
       bp = pk[, sum(end - start)], score = score_by_pid)
}

per_cell <- function(st, keep) {
  acr <- unique(st$hits[keep[pid], .(rid, bc)])[, .(acrs = .N), by = bc]
  f <- merge(st$total, acr, by = "bc", all.x = TRUE)
  f[is.na(acrs), acrs := 0][, FRiP := acrs / total]
  f[]
}
sel_random <- function(st, K) { s <- logical(st$npeak); s[sample.int(st$npeak, K)] <- TRUE; s }
sel_top    <- function(st, K) { s <- logical(st$npeak)
                                s[order(st$score, decreasing = TRUE)[seq_len(K)]] <- TRUE; s }

# per-cell FRiP at K: mean over draws for random, single deterministic pass for top
pc_at_K <- function(st, K, how) {
  if (how == "top") return(per_cell(st, sel_top(st, K))[, .(bc, FRiP)])
  rbindlist(lapply(seq_len(R_DRAWS), function(i) per_cell(st, sel_random(st, K))[, .(bc, FRiP)]
    ))[, .(FRiP = mean(FRiP)), by = bc]
}

# =============================================================================
# COMPUTE
# =============================================================================
# ONE SPECIES AT A TIME: the B73 arm carries ~35M reads per stage, so holding all
# four read x peak tables at once is what exhausts memory. Each species' four
# results are collected, then its tables are dropped before the next.
stg_lab <- c(pre = "Pre", post = "Post")
inv_l <- pc_l <- ksw_l <- com_l <- list()

# The compute is ~20 min (31M B73 reads x 325k peaks, ~170 per-cell passes), so
# the raw tables are cached and reused for plot/text iteration.
# QC_RECOMPUTE=1 forces a fresh run.
out <- function(suffix, ext = "tsv") file.path(OUTDIR, sprintf("%s_%s.%s", STEM, suffix, ext))
F_INV <- out("inventory"); F_PC  <- out("percell")
F_KSW <- out("ksweep_raw"); F_COM <- out("common_raw")
RECOMPUTE <- nzchar(Sys.getenv("QC_RECOMPUTE")) ||
             !all(file.exists(c(F_INV, F_PC, F_KSW, F_COM)))

if (RECOMPUTE) {
for (sp in names(SPP)) {
  message("[", sp, "] building read x peak overlaps ...")
  S   <- list(pre  = build(SPP[[sp]]$prefix, "pre",  SPP[[sp]]$pre),
              post = build(SPP[[sp]]$prefix, "post", SPP[[sp]]$post))
  lab <- SPP[[sp]]$label

  inv_l[[sp]] <- rbindlist(lapply(c("pre", "post"), function(stg) {
    st <- S[[stg]]
    data.table(species = lab, stage = stg_lab[[stg]], npeak = st$npeak,
               Mbp = round(st$bp / 1e6, 1),
               own_FRiP = round(median(per_cell(st, rep(TRUE, st$npeak))$FRiP), 4),
               cells = nrow(st$total))
  }))

  pc_l[[sp]] <- rbindlist(lapply(c("random", "top"), function(how)
    rbindlist(lapply(c("pre", "post"), function(stg)
      pc_at_K(S[[stg]], K_FAIR, how)[, `:=`(species = lab, stage = stg_lab[[stg]],
                                            selection = how, K = K_FAIR)]))))

  Kmax <- min(S$pre$npeak, S$post$npeak)
  ksw_l[[sp]] <- rbindlist(lapply(c("random", "top"), function(how)
    rbindlist(lapply(K_GRID[K_GRID <= Kmax], function(K) {
      m <- vapply(c("pre", "post"), function(stg)
        median(pc_at_K(S[[stg]], K, how)$FRiP), numeric(1))
      data.table(species = lab, selection = how, K = K, Pre = m[["pre"]], Post = m[["post"]])
    }))))

  com_l[[sp]] <- rbindlist(lapply(c("pre", "post"), function(ref)
    rbindlist(lapply(c("pre", "post"), function(stg) {
      # a stage on its OWN peaks is already built; only the cross terms are new
      st <- if (stg == ref) S[[stg]] else build(SPP[[sp]]$prefix, stg, SPP[[sp]][[ref]])
      r <- data.table(species = lab, peak_set = paste(stg_lab[[ref]], "peaks"),
                      stage = stg_lab[[stg]], npeak = st$npeak,
                      med = median(per_cell(st, rep(TRUE, st$npeak))$FRiP))
      if (stg != ref) { rm(st); invisible(gc(FALSE)) }
      r
    }))))

  rm(S); invisible(gc(FALSE))
}
  inv <- rbindlist(inv_l); percell <- rbindlist(pc_l)
  ksw <- rbindlist(ksw_l); common <- rbindlist(com_l)
  # saved RAW (pre-factor) so the cache branch can apply levels identically
  fwrite(inv, F_INV); fwrite(percell, F_PC); fwrite(ksw, F_KSW); fwrite(common, F_COM)
} else {
  message("[cache] reusing saved tables -- QC_RECOMPUTE=1 to force a fresh ~20 min compute")
  inv <- fread(F_INV); percell <- fread(F_PC); ksw <- fread(F_KSW); common <- fread(F_COM)
}

cat("=== peak sets (combined-genome) ===\n")
print(inv)
cat("\n  K_FAIR =", format(K_FAIR, big.mark = ","),
    "-- well below every peak set above, so NEITHER stage is scored at its ceiling.\n")

# --- headline: per-cell FRiP at K_FAIR, both selection rules -------------------
percell[, species   := factor(species, levels = SP_LEV)]
percell[, stage     := factor(stage,   levels = ST_LEV)]
percell[, selection := factor(selection, levels = c("random", "top"),
                              labels = c("random K peaks", "top K peaks"))]

summ <- percell[, .(cells = .N, med = median(FRiP)), by = .(species, selection, stage)]
w <- dcast(summ, species + selection ~ stage, value.var = "med")
w[, `:=`(diff = round(Post - Pre, 4), ratio = round(Post / Pre, 3))]
cat(sprintf("\n=== HEADLINE: median per-cell FRiP at K = %s ===\n", format(K_FAIR, big.mark = ",")))
print(w)
cat("\n  NOTE: NEITHER RULE IS NEUTRAL, AND THEY BIAS IN OPPOSITE DIRECTIONS:\n",
    "      random -> favours the SMALLER peak set (K is a larger share of it)\n",
    "      top    -> favours the LARGER peak set (more candidates to pick the best from)\n",
    "    So they BRACKET the truth. Where they straddle zero, the sign of the effect is\n",
    "    NOT identifiable from these data -- do not quote either end as the answer.\n",
    "    Maize, whose peak sets are near-equal, agrees under both rules: nothing to identify.\n",
    sep = "")

# --- K sweep, both rules -------------------------------------------------------
ksw[, `:=`(diff = round(Post - Pre, 4), ratio = round(Post / Pre, 3))]
cat("\n=== K sweep: is the effect stable in K? ===\n"); print(ksw[order(species, selection, K)])

# --- common-peak-set control: score BOTH stages on the SAME peaks --------------
# The most decisive test. If Post's reads are genuinely more concentrated, Post
# beats Pre even when the peak set is held fixed -- no subsampling involved.
cw <- dcast(common, species + peak_set + npeak ~ stage, value.var = "med")
cw[, `:=`(diff = round(Post - Pre, 4), ratio = round(Post / Pre, 3))]
cat("\n=== COMMON PEAK SET control (no subsampling; peak set held fixed) ===\n")
print(cw[order(species, peak_set)])
cat("\n  NOTE: THIS CONTROL IS ALSO NOT NEUTRAL -- and it shows why most clearly.\n",
    "    Each peak set was CALLED ON the reads of one stage, so it favours that stage:\n",
    "    on Post's peaks At reads give a large POSITIVE diff, on Pre's peaks a NEGATIVE\n",
    "    one. Same cells, same reads, opposite conclusions from the choice of peaks alone.\n",
    "    A neutral test needs an EXTERNAL peak set called on neither stage. Until then the\n",
    "    per-cell quality claim belongs to pTSS, which uses a fixed TSS annotation.\n",
    "    Maize moves by <=0.004 either way -- the negative control behaves.\n", sep = "")

# =============================================================================
# PANEL F -- violin/box of per-cell FRiP, in the style of the TSS panel
# =============================================================================
med_lab <- summ[, .(species, selection, stage, med)]
pFair <- ggplot(percell, aes(stage, FRiP, fill = stage)) +
  geom_violin(trim = FALSE, alpha = 0.75, width = 0.9, linewidth = 0.2, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.55, fill = "white", linewidth = 0.3) +
  geom_text(data = med_lab, aes(y = med, label = sprintf("%.3f", med)), vjust = -0.7,
            size = 2.5, fontface = "bold", colour = "grey15") +
  facet_grid(species ~ selection, scales = "free_y") +
  scale_fill_manual(values = ST_COL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7.6, colour = "grey35", lineheight = 1.25)) +
  labs(title = sprintf("Peak-matched FRiP at K = %s: the two matching rules disagree in SIGN",
                       format(K_FAIR, big.mark = ",")),
       subtitle = paste0(
         "Both stages cut to the SAME number of peaks, K well below either peak set so neither is\n",
         "scored at its ceiling. LEFT: K drawn at random (the rule the earlier panel uses) -- favours\n",
         "the SMALLER peak set, because K is a larger share of it. RIGHT: each stage's K STRONGEST\n",
         "peaks -- favours the LARGER set, because it has more candidates to pick the best from.\n",
         "Neither is neutral; they BRACKET the truth. On Arabidopsis they straddle zero (Post is\n",
         "higher on the left, lower on the right), so the sign of the effect is NOT identifiable.\n",
         "Maize, whose two peak sets are near-equal, agrees under both rules -- the negative control.\n",
         "y scales are per-panel; compare Pre vs Post WITHIN a panel, never across."),
       x = NULL, y = sprintf("Per-cell FRiP at K = %s", format(K_FAIR, big.mark = ",")))

written <- character()
for (ext in c("pdf", "png")) {
  f <- file.path(OUTDIR, paste0(STEM, ".", ext))
  # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
  if (ext == "pdf") ggsave(f, pFair, width = 6.4, height = 5.4, bg = "white")
  else              ggsave(f, pFair, width = 6.4, height = 5.4, dpi = 300, bg = "white")
  written <- c(written, f)
}
fwrite(w,  out("summary"), sep = "\t")
fwrite(ksw, out("ksweep"), sep = "\t")
fwrite(cw, out("common"),  sep = "\t")
written <- c(written, F_INV, F_PC, F_KSW, F_COM,
             out("summary"), out("ksweep"), out("common"))

cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz <- if (file.exists(f)) file.size(f) else NA_integer_
  good <- !is.na(sz) && sz > (if (grepl("\\.tsv$", f)) 50 else 1000)
  ok <- ok && good
  cat(sprintf("  %-4s %-34s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
cat("\n[done] ", STEM, ".{pdf,png} + 7 TSVs -> ", OUTDIR, "/\n", sep = "")
