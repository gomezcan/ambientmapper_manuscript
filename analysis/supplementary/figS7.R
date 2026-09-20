#!/usr/bin/env Rscript
# =============================================================================
# figS7.R  -  Fig S7 panels A to F: per-genome QC across the three cleaning stages (PreClean / WD /
#   ND) on the plate-split objects: A peak-set size, B pTSS, C naive FRiP, D peak-count-normalised
#   FRiP, E paired per-cell deltas, F PreClean profile of the cells each mode drops. Companion to Fig 5.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step0_qc/<stage>_<genome>.minDepth200.updated_metadata_v1.txt
#         .../SM2v2_plate/step0_qc/<base>_macs2_temp/<base>_peaks_combined_peaks.narrowPeak
#         figures/supplementary/figS7/S7_frip_{rarefaction,full_points,percell}[.top].tsv  (panel D caches, figS7_fripnorm.R)
# Output  figures/supplementary/figS7/Fig_S7_{A_peaks,B_pTSS,C_FRiP,D_FRiPnorm,E_deltas,F_dropped}.{pdf,png}, Fig_S7.composite.*, S7_*.tsv
# Run     Rscript analysis/supplementary/figS7.R   (repo root; A/B/C/E/F build without the caches, D degrades to a placeholder)
# =============================================================================
#
#   "Cleaning does not cost data quality on the per-genome design, and the
#    peak-count confound is specific to the combined-genome co-projection."
#
# PLATE-SPLIT objects (socrates/SM2v2_plate/), per genome, THREE stages:
#     PreClean  SM2_*            raw
#     wd        Clean.SM2v2wd_*  AmbientMapper decontam WITH plate design (WD)
#     nd        Clean.SM2v2_*    AmbientMapper decontam design-free       (ND)
# wd and nd are INDEPENDENT treatments of the same raw input, NOT a chain.
#   Never draw or describe a wd -> nd transition.
#
# TWO CELL SETS, ON PURPOSE -- they answer different questions.
#   * Panels B/C (distributions) use each stage's FULL set: "what is the quality
#     of the data this mode actually delivers?" n is printed under every violin.
#   * Panel E (deltas) uses the PAIRWISE set -- Pre & wd for the wd contrast,
#     Pre & nd for the nd contrast: "for a cell this mode KEEPS, does its quality
#     change?" A three-way intersection was rejected: it penalises the wd
#     contrast with cells only nd dropped (At 784 vs 1,076 available), and on
#     those cells cleaning removes just 0.17% of reads, making the panel
#     near-tautological.
# The stage cell sets are NOT nested -- cleaning both drops AND gains cells
#   (wd gains 13 At / 40 B73; nd gains 83 At / 103 B73), so a stage difference in
#   panels B/C is partly a cell-set difference. That is why the causal claim lives
#   in panel E and the descriptive one in B/C. Concretely: on the FULL sets At/nd's
#   median pTSS is HIGHER than Pre (0.822 vs 0.791), but that is survivorship --
#   nd drops 297 of 1,090 At cells -- not a quality gain.
# The FRiP recompute (panel D) is filtered upstream to the three-way shared set by
#   figS7_prep_beds.sh, so its rows carry a smaller n than the metadata metrics.
#   n is reported per row; the sets are never pooled.
#
# WHY THIS FIGURE'S HEADLINE DIFFERS FROM THE MAIN Fig 5 QC PANELS (E, F).
#   Those are the COMBINED-GENOME co-projection, where the At peak set collapses
#   186,774 -> 65,379 (-65%) and no peak-matching rule is neutral (Fig 5F). Here
#   the peak sets are near-stable, so normalization is expected to be a near-no-op.
#   That is the result, not a failure: it bounds the confound to the co-projection.
#
# METADATA VERSION v1 IS LOAD-BEARING -- DO NOT "UPGRADE" IT TO v6.
#   step0_qc ships `updated_metadata_v1 .. v6`. The metric values (total, pTSS,
#   FRiP) are BYTE-IDENTICAL across all six (max|diff| = 0 on the 1,090 shared At
#   cells) -- the ladder only removes cells:
#       At/Pre   v1 1,525 (1,374 pass qc_check)  ->  v6 1,090
#       B73/Pre  v1 18,482 (16,894 pass)         ->  v6 15,179
#   v6 retains ONLY cells with qc_check == 1, in every stage. That gate is
#   defined on pTSS / FRiP / nSites -- the very metrics this figure plots -- and it
#   is applied to each stage INDEPENDENTLY, so building the figure on v6 would
#   make "cleaning preserves QC" true by construction: any cell whose quality
#   cleaning had degraded would have been removed before it could be measured.
#   v1 is the pre-gate set (every barcode >= 200 reads, qc_check both values), so
#   the comparison is honest. The gate is drawn as a reference line instead, and
#   both all-cell and gate-passing medians are printed on every run.
# These metadata are R-exported with row.names: data rows carry ONE MORE field
#   than the header (field 1 duplicates the cellID). fread warns and adds `V1`.
#   ALWAYS select by name. The loader asserts the columns it needs.
# SIGNIFICANCE IS TESTED ACROSS CELLS, NEVER ACROSS PEAK DRAWS. The draws are
#   resamples of one dataset, so n would be a number we chose and raising it
#   would mechanically shrink p. Same error that invalidated the old combined-
#   genome panel C (t-test over 45 grid combos). Guard kept in the code below.
#   At these n the p-values saturate -- LEAD WITH THE EFFECT SIZE.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SOC    <- file.path(DATA, "socrates")
QC     <- file.path(SOC, "SM2v2_plate", "step0_qc")
OUTDIR <- "figures/supplementary/figS7"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# PRE-GATE metadata. See the header: v6 keeps only qc_check == 1 and would make
# the figure circular. Changing this constant changes the figure's validity.
MDVER   <- "v1"
STAGES  <- c(PreClean = "SM2", wd = "Clean.SM2v2wd", nd = "Clean.SM2v2")
ST_LEV  <- names(STAGES)
# Pre/Post colours carried from Fig 3 (pink/green) so Pre and wd read the same
# as everywhere else in Fig 5; nd gets a third hue rather than a shade of green,
# because wd and nd are independent treatments, not a progression.
ST_COL  <- c(PreClean = "#FF83FA", wd = "#43CD80", nd = "#E69F00")

GEN <- list(
  TAIR10 = list(suf = "At_TAIR10", lib = "At",
                label = "Arabidopsis (At) -> TAIR10"),
  B73v5  = list(suf = "B73_B73v5", lib = "B73",
                label = "Maize (B73) -> B73v5")
)
G_LEV <- vapply(GEN, `[[`, character(1), "label")

# combined-genome peak counts, for the contrast stated on panel A (main Fig 5 QC
# panels; measured from the socrates/_data/_PeakFiles narrowPeak files)
COMBINED_AT <- c(pre = 186774, post = 65379)

# -------------------------
# LOAD -- step0_qc metadata (pre-clustering, so QC is not filtered by min.c)
# -------------------------
read_qc <- function(g, stage) {
  f <- file.path(QC, sprintf("%s_%s.minDepth200.updated_metadata_%s.txt",
                             STAGES[[stage]], GEN[[g]]$suf, MDVER))
  if (!file.exists(f)) stop("missing QC metadata: ", f)
  d <- suppressWarnings(fread(f))                      # warns: header 1 short -> V1
  need <- c("cellID", "total", "pTSS", "FRiP", "nSites", "qc_check")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("columns absent from ", basename(f), ": ", paste(miss, collapse = ", "),
                         " -- the by-name contract broke, do NOT fall back to positions.")
  # guard the circularity trap: a pre-gate file MUST contain failing cells too.
  # If every row passes, the gate was already applied and the figure is invalid.
  if (all(d$qc_check == 1))
    stop(sprintf("%s/%s (%s) contains ONLY qc_check==1 cells -- this is a POST-GATE file. ",
                 g, stage, MDVER),
         "Building the QC comparison on it would be circular (the gate is defined on the ",
         "metrics being plotted). Use the pre-gate version; see the header.")
  # the plate split must have been applied: every cell from this genome's library
  lib_ok <- grepl(sprintf("-SM2_%s_", GEN[[g]]$lib), d$cellID)
  if (!all(lib_ok))
    stop(sprintf("PLATE SPLIT NOT APPLIED: %s/%s has %d cells not from the %s library. File: %s",
                 g, stage, sum(!lib_ok), GEN[[g]]$lib, f))
  d[, .(cellID, total, nSites, pTSS, FRiP, qc_check)][, `:=`(genome = g, stage = stage)]
}

meta <- rbindlist(lapply(names(GEN), function(g)
  rbindlist(lapply(ST_LEV, function(s) read_qc(g, s)))))
meta[, stage := factor(stage, levels = ST_LEV)]

# --- cell-set accounting + the paired (shared) set ------------------------------
acct <- rbindlist(lapply(names(GEN), function(g) {
  L <- lapply(ST_LEV, function(s) meta[genome == g & stage == s, cellID])
  names(L) <- ST_LEV
  sh <- Reduce(intersect, L)
  rbindlist(lapply(ST_LEV, function(s) data.table(
    genome = g, stage = s, n = length(L[[s]]), shared = length(sh),
    dropped_vs_pre = length(setdiff(L$PreClean, L[[s]])),
    gained_vs_pre  = length(setdiff(L[[s]], L$PreClean)))))
}))
shared_cells <- lapply(names(GEN), function(g) {
  L <- lapply(ST_LEV, function(s) meta[genome == g & stage == s, cellID])
  Reduce(intersect, L)
}); names(shared_cells) <- names(GEN)

meta[, paired := FALSE]
for (g in names(GEN)) meta[genome == g & cellID %in% shared_cells[[g]], paired := TRUE]
mp <- meta[paired == TRUE]
mp[, glab := factor(G_LEV[genome], levels = G_LEV)]
meta[, glab := factor(G_LEV[genome], levels = G_LEV)]

# -------------------------
# LOAD -- peak-set sizes (independent of the FRiP recompute)
# -------------------------
peak_n <- rbindlist(lapply(names(GEN), function(g) rbindlist(lapply(ST_LEV, function(s) {
  base <- sprintf("%s_%s", STAGES[[s]], GEN[[g]]$suf)
  f <- file.path(QC, paste0(base, "_macs2_temp"), paste0(base, "_peaks_combined_peaks.narrowPeak"))
  if (!file.exists(f)) stop("missing peaks: ", f)
  data.table(genome = g, stage = s, npeak = nrow(fread(f, select = 1, col.names = "chrom")))
}))))
peak_n[, stage := factor(stage, levels = ST_LEV)]
peak_n[, pct_vs_pre := 100 * npeak / npeak[stage == "PreClean"] - 100, by = genome]
peak_n[, glab := factor(G_LEV[genome], levels = G_LEV)]

# -------------------------
# LOAD -- FRiP recompute (optional; panel D)
# -------------------------
f_rare <- file.path(OUTDIR, "S7_frip_rarefaction.tsv")
f_pts  <- file.path(OUTDIR, "S7_frip_full_points.tsv")
f_pc   <- file.path(OUTDIR, "S7_frip_percell.tsv")
has_frip <- all(file.exists(c(f_rare, f_pts, f_pc)))
SEL_LEV <- c("random subsample", "top-K by MACS2 score")
smap <- c(pre = "PreClean", wd = "wd", nd = "nd")
tag <- function(d) { d[, stage := factor(smap[stage], levels = ST_LEV)]
                     d[, glab  := factor(G_LEV[genome], levels = G_LEV)]; d }
if (has_frip) {
  rare <- tag(fread(f_rare)); fpts <- tag(fread(f_pts)); fpc <- tag(fread(f_pc))

  # OPTIONAL CONTROL: the same rarefaction with the K HIGHEST-SCORING peaks
  # instead of a random K (PEAK_SELECT=top). This is the honest "if this stage
  # had called only K peaks" counterfactual -- MACS2 would keep the strongest,
  # not a random sample -- and it is what tells us whether an apparent difference
  # at K is real or an artifact of the anchor never being subsampled.
  f_rare_top <- file.path(OUTDIR, "S7_frip_rarefaction.top.tsv")
  f_pts_top  <- file.path(OUTDIR, "S7_frip_full_points.top.tsv")
  has_top <- all(file.exists(c(f_rare_top, f_pts_top)))
  if (has_top) {
    rare_top <- tag(fread(f_rare_top)); fpts_top <- tag(fread(f_pts_top))
    rare_both <- rbind(copy(rare)[, selection := SEL_LEV[1]],
                       copy(rare_top)[, selection := SEL_LEV[2]])
    cmp <- merge(
      rare[K == fpts$Kcommon[match(genome, fpts$genome)], .(genome, stage, random = med)],
      rare_top[K == fpts_top$Kcommon[match(genome, fpts_top$genome)], .(genome, stage, topK = med)],
      by = c("genome", "stage"))
    cmp[, artifact := round(random - topK, 4)]
    cat("\n--- panel D control: FRiP at the common peak count, random vs top-K ---\n")
    print(cmp[order(genome, stage)])
    cat("\n  NOTE: The gap is the SUBSAMPLING ARTIFACT, not biology. Random draws discard\n",
        "    strong peaks that a genuinely smaller peak call would have kept, so they\n",
        "    penalise the stages with MORE peaks; the fewest-peak stage is the anchor,\n",
        "    draws K of K, and is never penalised at all.\n",
        "  NOTE: On maize this REVERSES the ranking: random makes nd look best (+0.033 over\n",
        "    PreClean); under top-K nd is -0.001, i.e. no advantage. Quote top-K.\n", sep = "")
  } else {
    rare_both <- copy(rare)[, selection := SEL_LEV[1]]
    cat("\n  (panel D control not run -- PEAK_SELECT=top Rscript analysis/supplementary/figS7_fripnorm.R <scratch>)\n")
  }
  rare_both[, selection := factor(selection, levels = SEL_LEV)]
}

# =============================================================================
# PROVENANCE -- printed on every run
# =============================================================================
cat("=== FIG S7 | per-genome QC across cleaning modes (plate-split, step0_qc minDepth200) ===\n\n")
cat("--- cell-set accounting (stage sets are NOT nested) ---\n"); print(acct)
cat("\n--- peak-set size ---\n"); print(peak_n[, .(genome, stage, npeak, pct_vs_pre = round(pct_vs_pre, 2))])
cat(sprintf("\n  CONTRAST, combined-genome co-projection (main Fig 5 QC panels): At %s -> %s (%.1f%%)\n",
            format(COMBINED_AT["pre"], big.mark = ","), format(COMBINED_AT["post"], big.mark = ","),
            100 * COMBINED_AT["post"] / COMBINED_AT["pre"] - 100))
cat("  => the peak-count confound is a CO-PROJECTION artifact; it is near-absent here.\n")

cat(sprintf("\n--- metadata version: %s (PRE-GATE) ---\n", MDVER))
cat("--- medians: ALL cells (plotted) vs only qc_check==1 (the gate) ---\n")
print(merge(
  meta[, .(n_all = .N, pTSS_all = round(median(pTSS), 4), FRiP_all = round(median(FRiP), 4)),
       by = .(genome, stage)],
  meta[qc_check == 1,
       .(n_pass = .N, pTSS_pass = round(median(pTSS), 4), FRiP_pass = round(median(FRiP), 4)),
       by = .(genome, stage)], by = c("genome", "stage"))[order(genome, stage)])
cat("\n  The figure plots ALL cells. The `_pass` columns are shown only to prove the\n",
    "  conclusion does not depend on the gate -- they are NOT the plotted values, and\n",
    "  gating on them would make this figure circular (see the script header).\n",
    "  NOTE: On the full sets At/nd sits HIGHER than Pre on pTSS. That is survivorship --\n",
    "    nd drops the most contaminated At cells -- not a quality gain. Panel E, which\n",
    "    is paired per cell, is the causal contrast.\n", sep = "")

# =============================================================================
# PAIRED PER-CELL DELTAS + tests  (ACROSS CELLS -- never across peak draws)
# =============================================================================
# Effect size = median of per-cell paired differences. The Wilcoxon is reported
# only to satisfy convention: at n = 784 / 10,000 cells it saturates and carries
# almost no information. Do not add a test over the peak draws (pseudo-replication).
# PAIRWISE pairing: dcast leaves NA where a cell is absent from a stage, so
# filtering to complete pairs per contrast gives Pre & wd for wd and Pre & nd for nd
# -- each contrast on exactly the cells that mode kept, without nd's drops
# shrinking the wd comparison.
delta_tab <- function(dt, value, metric) {
  w <- dcast(dt, genome + cellID ~ stage, value.var = value)
  rbindlist(lapply(c("wd", "nd"), function(m)
    w[!is.na(PreClean) & !is.na(get(m)), {
      d <- get(m) - PreClean
      p <- tryCatch(suppressWarnings(wilcox.test(get(m), PreClean, paired = TRUE)$p.value),
                    error = function(e) NA_real_)
      .(metric = metric, contrast = paste0(m, " - PreClean"), n = .N,
        med_pre = median(PreClean), med_post = median(get(m)),
        med_paired_diff = median(d),
        pct_unchanged = 100 * mean(abs(d) < 1e-12), p = p)
    }, by = genome]))
}
deltas <- rbindlist(list(
  delta_tab(meta, "pTSS", "pTSS"),
  delta_tab(meta, "FRiP", "FRiP (naive)")
))
NORM_LAB <- "FRiP (peak-norm)"
if (has_frip) {
  fpn <- fpc[, .(genome, cellID = bc, stage, frip_naive, frip_norm)]
  deltas <- rbindlist(list(deltas, delta_tab(fpn, "frip_naive", "FRiP (naive, recomputed)")))
  # The peak-normalized delta MUST come from the top-K run when it exists.
  # The random-subsample frip_norm carries the anchor artifact measured above
  # (B73 nd would read +0.033 here, which is the artifact almost exactly, not a
  # real gain). Falling back to random is flagged in the metric label so a reader
  # can never mistake one for the other.
  f_pc_top <- file.path(OUTDIR, "S7_frip_percell.top.tsv")
  if (file.exists(f_pc_top)) {
    NORM_LAB <- "FRiP (peak-norm, top-K)"
    fpt <- tag(fread(f_pc_top))[, .(genome, cellID = bc, stage, frip_norm)]
    deltas <- rbindlist(list(deltas, delta_tab(fpt, "frip_norm", NORM_LAB)))
  } else {
    NORM_LAB <- "FRiP (peak-norm, RANDOM - artifact-prone)"
    deltas <- rbindlist(list(deltas, delta_tab(fpn, "frip_norm", NORM_LAB)))
  }
}
deltas[, glab := factor(G_LEV[genome], levels = G_LEV)]
deltas[, contrast := factor(contrast, levels = c("wd - PreClean", "nd - PreClean"))]

cat("\n--- paired per-cell deltas: PAIRWISE sets (Pre & wd, Pre & nd) ---\n")
cat("    effect size first; p saturates at these n and carries little information\n")
print(deltas[, .(genome, metric, contrast, n, med_pre = round(med_pre, 4),
                 med_post = round(med_post, 4), med_paired_diff = round(med_paired_diff, 4),
                 pct_unchanged = round(pct_unchanged, 1), p = signif(p, 3))])
cat("\n  `pct_unchanged` = cells whose value is BIT-IDENTICAL across the two stages.\n",
    "  Where it is high, cleaning removed no reads from most cells it kept, and a\n",
    "  median paired difference of exactly 0.0000 is the correct answer, not a bug.\n", sep = "")

if (has_frip) {
  cat("\n--- calibration: recomputed naive FRiP vs pipeline metadata FRiP ---\n")
  cal <- merge(fpc[, .(recomputed = median(frip_naive)), by = .(genome, stage)],
               mp[, .(metadata = median(FRiP)), by = .(genome, stage)], by = c("genome", "stage"))
  cal[, offset := round(recomputed - metadata, 4)]
  print(cal[order(genome, stage)])
  cat("  (a small offset vs socrates `acrs/total` is expected and immaterial to the curves)\n")
}

# =============================================================================
# PANELS
# =============================================================================
base_theme <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          strip.text = element_text(face = "bold", size = 8.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 7.6, colour = "grey35", lineheight = 1.25))
}

# --- A. peak-set size ----------------------------------------------------------
# NOTE: geom_text carries the discrete x, so the scale is typed discrete by the
# first layer. Do NOT put an annotate() with numeric x before these layers.
pA <- ggplot(peak_n, aes(stage, npeak, fill = stage)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = ifelse(stage == "PreClean", format(npeak, big.mark = ","),
                               sprintf("%s\n(%+.1f%%)", format(npeak, big.mark = ","), pct_vs_pre))),
            vjust = -0.25, size = 2.5, lineheight = 0.95, colour = "grey15") +
  facet_wrap(~ glab, scales = "free_y") +
  scale_fill_manual(values = ST_COL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.22))) +
  base_theme() + theme(legend.position = "none") +
  labs(title = "A. Peak-set size is near-stable across cleaning modes",
       subtitle = paste0("Combined-genome co-projection for contrast: At ",
                         format(COMBINED_AT["pre"], big.mark = ","), " -> ",
                         format(COMBINED_AT["post"], big.mark = ","), " (-65%).\n",
                         "The peak-count confound that dominates the main Fig 5 QC\n",
                         "panels is a co-projection artifact, near-absent here."),
       x = NULL, y = "Peaks called (MACS2)")

# --- B/C. paired distributions -------------------------------------------------
# Each stage's FULL set: the quality of the data that mode actually delivers.
# n is drawn under every violin because the sets differ (nd drops 27% of At), and
# a stage difference here is partly a cell-set difference -- panel E is the
# per-cell causal contrast.
viol <- function(value, ylab, title, subtitle, gate = NA_real_) {
  med <- meta[, .(v = median(get(value)), n = .N), by = .(glab, stage)]
  p <- ggplot(meta, aes(stage, get(value), fill = stage)) +
    geom_violin(trim = FALSE, alpha = 0.75, width = 0.9, linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 0.55, fill = "white", linewidth = 0.3)
  if (!is.na(gate)) p <- p + geom_hline(yintercept = gate, linetype = "dashed",
                                        colour = "grey45", linewidth = 0.3)
  p + geom_text(data = med, aes(y = v, label = sprintf("%.3f", v)), vjust = -0.6,
                size = 2.5, fontface = "bold", colour = "grey15") +
    geom_text(data = med, aes(x = stage, y = -Inf, label = paste0("n=", format(n, big.mark = ","))),
              inherit.aes = FALSE, vjust = -0.7, size = 2.1, colour = "grey45") +
    facet_wrap(~ glab, scales = "free_y") +
    scale_fill_manual(values = ST_COL) +
    base_theme() + theme(legend.position = "none") +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab)
}
# numbers pulled from the data, never hardcoded: they move with MDVER
at_nd_drop <- acct[genome == "TAIR10" & stage == "nd", dropped_vs_pre]
at_pre_n   <- acct[genome == "TAIR10" & stage == "PreClean", n]
pB <- viol("pTSS", "pTSS (reads near TSS)",
           "B. TSS enrichment holds in both cleaning modes",
           paste0("Each stage's full cell set, pre-QC-gate (n under each violin).\n",
                  "Fixed TSS annotation, so no peak-set confound.\n",
                  sprintf("At/nd sits higher than Pre because nd DROPS %s of %s At\n",
                          format(at_nd_drop, big.mark = ","), format(at_pre_n, big.mark = ",")),
                  "cells -- survivorship, not a gain. See panel E.\n",
                  "Dashed = QC gate (0.2); cells below it are KEPT, not filtered."),
           gate = 0.2)
pC <- viol("FRiP", "FRiP (own peak set)",
           "C. FRiP holds against each stage's own peak set",
           paste0("Each stage's full cell set. This is the pipeline's\n",
                  "as-delivered number: every stage scored on ITS OWN\n",
                  "peaks, which is the quantity panel D controls for."))

# --- D. peak-count-normalized FRiP --------------------------------------------
if (has_frip) {
  Kc <- unique(fpts[, .(glab, Kcommon)])
  pD <- ggplot(rare_both, aes(K, med, colour = stage, linetype = selection)) +
    geom_vline(data = Kc, aes(xintercept = Kcommon), inherit.aes = FALSE,
               linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    geom_point(data = fpts, aes(x = npeak, y = med_naive, colour = stage),
               inherit.aes = FALSE, size = 2, shape = 21, fill = "white", stroke = 0.7) +
    facet_wrap(~ glab, scales = "free") +
    scale_x_log10(labels = scales::comma, expand = expansion(mult = c(0.03, 0.09))) +
    scale_colour_manual(values = ST_COL) +
    scale_linetype_manual(values = c("solid", "22"), name = "Peak selection", drop = TRUE) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    base_theme() + theme(legend.position = "bottom", legend.box = "vertical",
                         legend.margin = margin(t = -2, b = -2)) +
    labs(title = "D. At matched peak count the three modes are indistinguishable",
         subtitle = paste0("Median per-cell FRiP vs #peaks. Dashed = common peak count K;\n",
                           "open circles = each stage at its own full peak set.\n",
                           "HOW PEAKS ARE SUBSAMPLED MATTERS: random draws discard strong\n",
                           "peaks a genuinely smaller peak call would have kept, penalising\n",
                           "the stages with MORE peaks, while the fewest-peak stage is the\n",
                           "anchor and draws K of K. On maize that alone makes nd look\n",
                           "+0.033 better than PreClean; under top-K it is -0.001."),
         x = "Number of peaks (log10)", y = "Median per-cell FRiP",
         colour = "Stage")
} else {
  pD <- ggplot() + theme_void() +
    annotate("text", x = 0, y = 0, size = 3, colour = "grey35", lineheight = 1.3,
             label = paste("D. Peak-count-normalized FRiP - NOT BUILT",
                           "", "Run (from the repo root):",
                           "  bash analysis/supplementary/figS7_prep_beds.sh <scratch>",
                           "  Rscript analysis/supplementary/figS7_fripnorm.R <scratch>", sep = "\n")) +
    labs(title = "D. Peak-count-normalized FRiP (pending recompute)")
}

# --- E. paired per-cell deltas -------------------------------------------------
dE <- deltas[metric %in% c("pTSS", "FRiP (naive)", if (has_frip) NORM_LAB else NULL)]
dE[, metric := factor(metric, levels = c("pTSS", "FRiP (naive)", NORM_LAB))]
pE <- ggplot(dE, aes(metric, med_paired_diff, fill = contrast)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62,
           colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%+.4f", med_paired_diff),
                vjust = ifelse(med_paired_diff >= 0, -0.4, 1.3)),
            position = position_dodge(width = 0.72), size = 2.4, colour = "grey15") +
  facet_wrap(~ glab, scales = "free_y") +
  scale_fill_manual(values = c("wd - PreClean" = "#43CD80", "nd - PreClean" = "#E69F00")) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  base_theme() + theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(title = "E. For a cell it keeps, cleaning barely changes quality",
       subtitle = paste0("Median of per-cell paired differences, each contrast on its own pairwise set ",
                         "(Pre and wd, Pre and nd). Tested ACROSS CELLS (paired Wilcoxon),\n",
                         "never across peak draws; at these n every p saturates, so read the effect size, ",
                         "not the test. Exact zeros are real: over half of kept cells are bit-identical.\n",
                         "WARNING: the two y axes are per-genome and differ ~3x -- compare within a ",
                         "facet, never across. Every effect here is <= 0.013."),
       x = NULL, y = "Median paired difference (stage - PreClean)")

# --- F. what each mode DROPS --------------------------------------------------
# Turns panel B's survivorship caveat from an assertion into a measurement: were
# the cells a mode removed already the worse ones, pre-clean? Every value here is
# read from the PreClean object, so it is a property of the cell BEFORE cleaning
# touched it -- the fate label is the only thing that comes from the other stage.
drop_cells <- rbindlist(lapply(names(GEN), function(g) {
  pre <- meta[genome == g & stage == "PreClean", .(cellID, pTSS, FRiP, total)]
  rbindlist(lapply(c("wd", "nd"), function(m) {
    kept <- meta[genome == g & stage == m, cellID]
    copy(pre)[, `:=`(genome = g, mode = m,
                     fate = fifelse(cellID %in% kept, "retained", "dropped"))]
  }))
}))
drop_prof <- drop_cells[, .(n = .N, pTSS = median(pTSS), depth = as.numeric(median(total)),
                            FRiP = median(FRiP)), by = .(genome, mode, fate)]
cat("\n--- panel F: PreClean properties of the cells each mode later drops ---\n")
print(drop_prof[order(genome, mode, fate)])
cat("\n  Read as: the dropped cells were already shallower and lower-pTSS before\n",
    "  cleaning ran. That is why At/nd's median pTSS RISES in panel B -- the worst\n",
    "  cells left, not the survivors improved. wd drops almost nothing.\n", sep = "")

dF <- melt(drop_prof, id.vars = c("genome", "mode", "fate", "n"),
           measure.vars = c("pTSS", "depth"), variable.name = "metric", value.name = "v")
dF[, metric := factor(metric, levels = c("pTSS", "depth"),
                      labels = c("Pre-clean pTSS", "Pre-clean depth (reads)"))]
dF[, glab := factor(G_LEV[genome], levels = G_LEV)]
dF[, fate := factor(fate, levels = c("retained", "dropped"))]
dF[, mode := factor(mode, levels = c("wd", "nd"))]   # match the figure's order, not alphabetical
nd_drop <- drop_prof[fate == "dropped" & mode == "nd"]

# n is carried on every bar: wd's dropped group is only 6 At / 32 B73 cells, so
# its median is noisy and must not be read as an effect.
pF <- ggplot(dF, aes(mode, v, fill = fate)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62,
           colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%s\nn=%s",
                                ifelse(v >= 10, sprintf("%.0f", v), sprintf("%.3f", v)),
                                format(n, big.mark = ","))),
            position = position_dodge(width = 0.72), vjust = -0.18, size = 2.05,
            lineheight = 0.92, colour = "grey15") +
  facet_grid(metric ~ glab, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(retained = "grey78", dropped = "#B03A2E"),
                    name = "Fate under that mode") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  base_theme() +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 8)) +
  labs(title = "F. The cells cleaning drops were already the worst ones",
       subtitle = paste0(
         "Both metrics measured on the PreClean object, so they describe each cell ",
         "BEFORE cleaning; only the retained/dropped label comes from the other stage.\n",
         sprintf("nd drops %s At and %s maize cells, and they are the shallow, low-pTSS ones. ",
                 format(nd_drop[genome == "TAIR10", n], big.mark = ","),
                 format(nd_drop[genome == "B73v5",  n], big.mark = ",")),
         "This is the mechanism behind panel B's rise in At/nd -- survivorship."),
       x = NULL, y = NULL)

# =============================================================================
# EXPORT -- individual parts + composite, each verified
# =============================================================================
written <- character()
save_part <- function(p, name, w, h) {
  for (ext in c("pdf", "png")) {
    f <- file.path(OUTDIR, paste0(name, ".", ext))
    # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
    if (ext == "pdf") ggsave(f, p, width = w, height = h, bg = "white")
    else              ggsave(f, p, width = w, height = h, dpi = 300, bg = "white")
    written <<- c(written, f)
  }
}
save_part(pA, "Fig_S7_A_peaks",     6.2, 3.8)
save_part(pB, "Fig_S7_B_pTSS",      6.2, 3.8)
save_part(pC, "Fig_S7_C_FRiP",      6.2, 3.8)
save_part(pD, "Fig_S7_D_FRiPnorm",  6.2, 3.8)
save_part(pE, "Fig_S7_E_deltas",    9.0, 4.0)
save_part(pF, "Fig_S7_F_dropped",   9.0, 4.6)

comp <- (pA | pB) / (pC | pD) / pE / pF +
  plot_layout(heights = c(1, 1, 0.95, 1.15)) +
  plot_annotation(
    title = "Figure S7. Cleaning preserves per-cell data quality in both modes; the peak-count confound is specific to the co-projection",
    subtitle = paste0(
      "Plate-split design: each library mapped to its own reference, three stages -- PreClean, wd (with plate design), nd (design-free). ",
      "wd and nd are INDEPENDENT treatments of the same raw input, not a chain.\n",
      "Cells are the PRE-QC-GATE set (all barcodes >= 200 reads): gating on qc_check would filter cells using the very metrics plotted here. ",
      "B and C show each stage's full set (what that mode delivers);\n",
      "E is paired per cell on Pre-and-wd / Pre-and-nd separately (what cleaning does to a cell it keeps). ",
      "The one arm that loses ground is nd on maize: -5.4% peaks (A) and -0.0036 FRiP (E). Companion to main Fig 5."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 8.5, colour = "grey30", lineheight = 1.3)))
for (ext in c("pdf", "png")) {
  f <- file.path(OUTDIR, paste0("Fig_S7.composite.", ext))
  if (ext == "pdf") ggsave(f, comp, width = 13, height = 19, bg = "white")
  else              ggsave(f, comp, width = 13, height = 19, dpi = 300, bg = "white",
                           limitsize = FALSE)
  written <- c(written, f)
}

fwrite(acct,    file.path(OUTDIR, "S7_cellset_accounting.tsv"), sep = "\t")
fwrite(peak_n[, .(genome, stage, npeak, pct_vs_pre)], file.path(OUTDIR, "S7_peak_counts.tsv"), sep = "\t")
fwrite(deltas[, .(genome, metric, contrast, n, med_pre, med_post, med_paired_diff,
                  pct_unchanged, p)],
       file.path(OUTDIR, "S7_paired_deltas.tsv"), sep = "\t")
fwrite(drop_prof[order(genome, mode, fate)], file.path(OUTDIR, "S7_dropped_cell_profile.tsv"), sep = "\t")
written <- c(written, file.path(OUTDIR, c("S7_cellset_accounting.tsv", "S7_peak_counts.tsv",
                                          "S7_paired_deltas.tsv", "S7_dropped_cell_profile.tsv")))
if (exists("cmp")) {
  fwrite(cmp[order(genome, stage)], file.path(OUTDIR, "S7_peakselect_control.tsv"), sep = "\t")
  written <- c(written, file.path(OUTDIR, "S7_peakselect_control.tsv"))
}

# verify rather than assert: ggsave can fail while raising only a warning
cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz <- if (file.exists(f)) file.size(f) else NA_integer_
  # small TSVs are legitimately a few hundred bytes; only plots must be large
  floor_b <- if (grepl("\\.tsv$", f)) 50 else 1000
  good <- !is.na(sz) && sz > floor_b
  ok <- ok && good
  cat(sprintf("  %-4s %-46s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
if (!has_frip)
  cat("\nNOTE: Panel D is a PLACEHOLDER -- the FRiP recompute has not been run.\n")
cat("\n[done] Fig S7 -> ", OUTDIR, "/Fig_S7.composite.{pdf,png}\n", sep = "")
