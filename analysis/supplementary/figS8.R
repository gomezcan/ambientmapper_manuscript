#!/usr/bin/env Rscript
# =============================================================================
# figS8.R  -  Fig S8 panels A to C: meta-cell reproducibility across the 100 SEACells seeds on the
#   plate-split objects. A seed-pair ARI/AMI violins (descriptive), B stage contrasts with run-level
#   bootstrap CIs (inferential), C per-meta-cell annotation support (descriptive). Companion to Fig 5.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/consensus/
#           {At,maize}.ari_ami_ci.tsv, {At,maize}.seedpair_ami.tsv, <prefix>.consensus.seedpairs.tsv
#         .../SM2v2_plate/step5_metacell/rZ_annotation/consensus/<prefix>.seed_label_sweep.metacell.tsv, {At,maize}.seed_label_sweep.summary.tsv
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/supplementary/figS8/Fig_S8.{pdf,png}, Fig_S8_{violins,diffs,support}.{pdf,png}, Fig_S8_{summary,ci,support_summary}.tsv
# Run     Rscript analysis/supplementary/figS8.R      (from the repo root, ~30 s)
# =============================================================================
#
# WHAT IS SHOWN
#   A. Violins of the seed-pair agreement values: for every pair of the 100
#      SEACells runs (C(100,2) = 4,950 pairs per genome x stage), the ARI and
#      AMI between the two runs' meta-cell partitions. DESCRIPTIVE ONLY.
#   B. Stage contrasts (wd-Pre, nd-Pre, wd-nd) with 95% percentile CIs from
#      the upstream run-level bootstrap -- the INFERENTIAL panel, and the
#      numbers a manuscript sentence should quote.
#   C. Per-meta-cell ANNOTATION support: for each consensus meta-cell, the
#      fraction of the 100 seeds whose own annotation gives it the same type
#      as the consensus call. DESCRIPTIVE ONLY. This is the former Part-3
#      support panel, rebuilt here from the upstream sweep files with the same
#      asserts.
#
# THE STATISTICAL UNIT (this project's recurring trap -- do not undo):
#   Seed PAIRS are not independent: each run appears in 99 pairs. Testing
#   across the 4,950 pairs would be pseudo-replication (same failure mode as
#   the retired grid-combo t-test and the meta-cell-unit tests). Therefore:
#     - NO test is computed across pairs anywhere in this script.
#     - All inference (points + CIs) comes from the UPSTREAM bootstrap that
#       resamples RUNS (B = 10,000; consensus/{At,maize}.ari_ami_ci.tsv),
#       the correct replication unit.
#     - Panel C's unit is the CONSENSUS META-CELL and it is descriptive for
#       the same reason: meta-cells share one partition and one set of 100
#       seeds, so they are not independent replicates of a stage. No test
#       across meta-cells; the stage median is UNWEIGHTED (never weight a
#       meta-cell stage summary by cell count -- point size shows size).
#   With n = 100 runs, tiny effects are significant -- LEAD WITH EFFECT SIZE.
#     At wd-Pre ARI +0.25 is large; maize -0.035 excludes 0 but is ~7x smaller
#     and negative: the minor/major asymmetry, stated as such on the panel.
#   REPRODUCIBILITY / AGREEMENT, never "accuracy" -- partitions are compared
#     between runs, never to a ground truth.
#   wd and nd are INDEPENDENT treatments of the same raw input, not a chain.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
OUTDIR <- "figures/supplementary/figS8"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source("analysis/_helpers/fig5_part2_helpers.R")   # CFG, STAGES, PLATE

CONS <- file.path(PLATE, "step5_metacell", "consensus")
RZC  <- file.path(PLATE, "step5_metacell", "rZ_annotation", "consensus")
N_SEED  <- 100L
N_PAIRS <- N_SEED * (N_SEED - 1L) / 2L          # 4,950

GKEY  <- c(B73v5 = "maize", TAIR10 = "At")      # helper key -> upstream key
G_LEV <- vapply(CFG, `[[`, character(1), "label")
S_LEV <- names(STAGES)                          # PreClean, wd, nd
STAGE_UP <- c(PreClean = "Pre", wd = "wd", nd = "nd")
CMP_LEV  <- c("wd-Pre", "nd-Pre", "wd-nd")
ST_COL   <- c(PreClean = "#FF83FA", wd = "#43CD80", nd = "#E69F00")
# Display labels only: the manuscript nomenclature is WD/ND. Data values stay
# wd/nd everywhere -- they join against upstream files.
STAGE_LAB <- c(PreClean = "PreClean", wd = "WD", nd = "ND")
CMP_LAB   <- c("wd-Pre" = "WD-Pre", "nd-Pre" = "ND-Pre", "wd-nd" = "WD-ND")

prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)

# =============================================================================
# LOAD + CONTRACT ASSERTIONS
# =============================================================================
# --- upstream bootstrap estimates and CIs (the inferential layer) -----------
ci <- rbindlist(lapply(unname(GKEY), function(gu) {
  f <- file.path(CONS, sprintf("%s.ari_ami_ci.tsv", gu))
  if (!file.exists(f)) stop("missing upstream CI table: ", f)
  fread(f)
}))
if (!all(ci$B == 10000L) ||
    !all(ci[kind == "stage", n_pairs] == N_PAIRS))   # diff rows carry n_pairs = NA
  stop("ari_ami_ci.tsv: unexpected B or n_pairs -- re-check the upstream run")
if (!setequal(unique(ci$kind), c("stage", "diff")))
  stop("ari_ami_ci.tsv: expected kinds {stage, diff}, got: ",
       paste(unique(ci$kind), collapse = ", "))
if (!setequal(ci[kind == "stage", unique(name)], unname(STAGE_UP)) ||
    !setequal(ci[kind == "diff",  unique(name)], CMP_LEV))
  stop("ari_ami_ci.tsv: unexpected stage/comparison names")
if (nrow(ci) != 2L * 2L * (3L + 3L))            # genomes x metrics x (stages+diffs)
  stop("ari_ami_ci.tsv: expected 24 rows total, got ", nrow(ci))

# --- seed-pair ARI (per arm) + AMI (per genome file): the descriptive layer --
pairs <- rbindlist(c(
  lapply(names(CFG), function(g) rbindlist(lapply(S_LEV, function(s) {
    f <- file.path(CONS, sprintf("%s.consensus.seedpairs.tsv", prefix_for(g, s)))
    if (!file.exists(f)) stop("missing seed-pair ARI: ", f)
    d <- fread(f)
    if (nrow(d) != N_PAIRS) stop(f, ": expected ", N_PAIRS, " pairs, got ", nrow(d))
    d[, .(genome = GKEY[[g]], stage = s, metric = "ARI", value = ari)]
  }))),
  lapply(unname(GKEY), function(gu) {
    f <- file.path(CONS, sprintf("%s.seedpair_ami.tsv", gu))
    if (!file.exists(f)) stop("missing seed-pair AMI: ", f)
    d <- fread(f)
    if (nrow(d) != 3L * N_PAIRS) stop(f, ": expected ", 3L * N_PAIRS, " rows")
    if (!setequal(unique(d$stage), unname(STAGE_UP)))
      stop(f, ": unexpected stages: ", paste(unique(d$stage), collapse = ", "))
    d[, .(genome = gu, stage = names(STAGE_UP)[match(stage, STAGE_UP)],
          metric = "AMI", value = ami)]
  })))

# --- reconcile the two layers: mean over pairs must equal the CI estimate ----
rec <- merge(pairs[, .(pair_mean = mean(value)), by = .(genome, stage, metric)],
             ci[kind == "stage",
                .(genome, stage = names(STAGE_UP)[match(name, STAGE_UP)],
                  metric, estimate)],
             by = c("genome", "stage", "metric"))
if (rec[abs(pair_mean - estimate) > 5e-3, .N] > 0) {
  print(rec[abs(pair_mean - estimate) > 5e-3])
  stop("seed-pair means do not reconcile with the upstream CI estimates")
}

# --- per-meta-cell label support over the 100 seeds (panel C) ----------------
sup <- rbindlist(lapply(names(CFG), function(g) rbindlist(lapply(S_LEV, function(s) {
  f <- file.path(RZC, sprintf("%s.seed_label_sweep.metacell.tsv", prefix_for(g, s)))
  if (!file.exists(f)) stop("missing seed-label sweep: ", f)
  d <- fread(f)
  need <- c("stage", "metacell", "n_cells", "label_support", "labels_agree")
  if (length(setdiff(need, names(d))) > 0)
    stop(f, ": missing columns: ", paste(setdiff(need, names(d)), collapse = ", "))
  # upstream stage naming is Pre/wd/nd, the figures' is PreClean/wd/nd --
  # assert the file really is the requested arm instead of joining raw strings
  if (uniqueN(d$stage) != 1L || d$stage[1] != STAGE_UP[[s]])
    stop(f, ": stage column is not ", STAGE_UP[[s]])
  if (d[!is.finite(label_support) | label_support < 0 | label_support > 1, .N] > 0)
    stop(f, ": label_support outside [0, 1]")
  if (d[n_cells <= 0L, .N] > 0) stop(f, ": non-positive meta-cell size")
  d[, .(genome = GKEY[[g]], stage = s, metacell, n_cells, label_support,
        labels_agree)]
}))))
# k is pinned per genome: all three stages annotate a same-size consensus partition
if (sup[, .N, by = .(genome, stage)][, uniqueN(N), by = genome][, any(V1 != 1L)])
  stop("meta-cell count differs across stages within a genome -- k should be pinned")

# reconcile against the upstream sweep summary before plotting anything
sup_med <- sup[, .(n_mc = .N, median_support = median(label_support),
                   frac_agree = mean(labels_agree)), by = .(genome, stage)]
for (gu in unname(GKEY)) {
  f <- file.path(RZC, sprintf("%s.seed_label_sweep.summary.tsv", gu))
  if (!file.exists(f)) stop("missing upstream sweep summary: ", f)
  up <- fread(f)
  if (!all(up$n_seeds == N_SEED)) stop(f, ": sweep is not over ", N_SEED, " seeds")
  chk <- merge(sup_med[genome == gu][, stage_up := STAGE_UP[stage]],
               up[, .(stage_up = stage, up_med = median_mc_label_support,
                      up_agree = frac_mc_label_matches_modal)],
               by = "stage_up")
  if (nrow(chk) != length(S_LEV))
    stop(gu, ": upstream sweep summary does not cover all three stages")
  if (chk[, max(abs(median_support - up_med))] > 1e-3 ||
      chk[, max(abs(frac_agree - up_agree))] > 1e-3)
    stop(gu, ": computed support disagrees with the upstream sweep summary")
}

lab_g <- function(d) {
  d[, glab := factor(fifelse(genome == "maize", G_LEV[["B73v5"]], G_LEV[["TAIR10"]]),
                     levels = unname(G_LEV))]
  d
}
pairs <- lab_g(pairs)[, stage := factor(stage, levels = S_LEV)]
pairs[, metric := factor(metric, levels = c("ARI", "AMI"))]
ci_st <- lab_g(ci[kind == "stage"])[
  , stage := factor(names(STAGE_UP)[match(name, STAGE_UP)], levels = S_LEV)]
ci_st[, metric := factor(metric, levels = c("ARI", "AMI"))]
ci_df <- lab_g(ci[kind == "diff"])[, cmp := factor(name, levels = CMP_LEV)]
ci_df[, metric := factor(metric, levels = c("ARI", "AMI"))]
sup     <- lab_g(sup)[,     stage := factor(stage, levels = S_LEV)]
sup_med <- lab_g(sup_med)[, stage := factor(stage, levels = S_LEV)]

# =============================================================================
# REPORT
# =============================================================================
cat("=== FIG S8 | run-level meta-cell reproducibility (", N_SEED, " seeds) ===\n\n",
    sep = "")
cat("source: ", CONS, "\n",
    "reconciliation seed-pair mean vs bootstrap estimate: max |diff| = ",
    format(rec[, max(abs(pair_mean - estimate))], digits = 2), "  (gate 5e-3, PASS)\n\n",
    sep = "")
cat("--- stage estimates [95% CI], bootstrap over RUNS (B = 10,000) ---\n")
print(dcast(ci_st, genome + name ~ metric, value.var = c("estimate", "ci_lo", "ci_hi"))[
  order(genome, factor(name, levels = unname(STAGE_UP)))])
cat("\n--- stage contrasts [95% CI] -- the quotable numbers ---\n")
print(ci_df[order(genome, metric, cmp),
            .(genome, metric, cmp, estimate, ci_lo, ci_hi,
              excludes_0 = (ci_lo > 0 | ci_hi < 0))])
cat("\n  NOTE: Unit = RUN (upstream bootstrap). No test is computed across seed pairs\n",
    "  here -- pairs are non-independent (each run sits in 99 pairs) and n_pairs is\n",
    "  a design choice; that would be the retired pseudo-replication error.\n",
    "  NOTE: Lead with EFFECT SIZE. All At contrasts exclude 0 (wd-Pre ARI +0.25, and\n",
    "  wd beats nd by +0.07). maize wd-Pre and nd-Pre exclude 0 but are ~7x smaller\n",
    "  and negative (-0.035; the minor/major asymmetry, consistent with S7 FRiP and\n",
    "  Part-3 agreement), while maize wd-nd does NOT exclude 0 -- the two cleaning\n",
    "  modes are indistinguishable on maize.\n",
    sep = "")

cat("\n--- panel C: label support by stage (unit = META-CELL, descriptive only) ---\n")
print(sup_med[order(genome, stage),
              .(genome, stage, n_mc, median_support = round(median_support, 3),
                pct_agree_modal = round(100 * frac_agree, 1))])
cat("  NOTE: No test across meta-cells (one partition, one seed set). The label-level\n",
    "  readout mirrors the partition level: At wd median support 0.62 -> 0.83 with\n",
    "  93% of calls modal, maize flat-to-down (nd worst, 0.61).\n", sep = "")

# =============================================================================
# PANEL A -- seed-pair distributions + run-bootstrap mean and CI
# =============================================================================
pA <- ggplot(pairs, aes(stage, value)) +
  geom_violin(aes(fill = stage), colour = "grey30", linewidth = 0.25,
              alpha = 0.75) +
  geom_linerange(data = ci_st, aes(x = stage, ymin = ci_lo, ymax = ci_hi),
                 inherit.aes = FALSE, linewidth = 0.9, colour = "grey10") +
  geom_point(data = ci_st, aes(x = stage, y = estimate),
             inherit.aes = FALSE, shape = 23, size = 1.9,
             fill = "white", colour = "grey10", stroke = 0.7) +
  facet_grid(metric ~ glab, scales = "free_y") +
  scale_fill_manual(values = ST_COL) +
  scale_x_discrete(labels = STAGE_LAB) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7.4, colour = "grey35", lineheight = 1.3)) +
  labs(title = "A. Agreement between independent SEACells runs",
       subtitle = paste0(
         "One value per PAIR of runs (n = 4,950 = C(100,2) per violin): ARI / AMI between the two runs' meta-cell partitions of the same object.\n",
         "White diamond + bar = mean and 95% CI from the run-level bootstrap (B = 10,000, resampling RUNS - the correct unit; ",
         "pairs are NOT independent, each run sits in\n99 pairs, so no test is computed across pairs). ",
         "REPRODUCIBILITY, not accuracy: runs are compared to each other, never to a ground truth. ",
         "WD and ND are independent\ntreatments of the same raw input, not a chain."),
       x = NULL, y = "agreement between runs (seed-pair value)")

# =============================================================================
# PANEL B -- stage contrasts with run-bootstrap CIs (the inferential panel)
# =============================================================================
pB <- ggplot(ci_df, aes(cmp, estimate)) +
  geom_hline(yintercept = 0, linetype = "22", colour = "grey50", linewidth = 0.35) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.14,
                linewidth = 0.55, colour = "grey10") +
  geom_point(shape = 23, size = 2.1, fill = "white", colour = "grey10", stroke = 0.8) +
  geom_text(aes(label = sprintf("%+.3f", estimate)),
            hjust = -0.28, vjust = 0.5, size = 2.35, colour = "grey25") +
  facet_grid(metric ~ glab, scales = "free_y") +
  scale_x_discrete(labels = CMP_LAB, expand = expansion(add = c(0.45, 0.75))) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7.4, colour = "grey35", lineheight = 1.3)) +
  labs(title = "B. Stage contrasts, bootstrap over runs",
       subtitle = paste0(
         "Difference of stage means with 95% percentile CIs from the run-level bootstrap (B = 10,000). ",
         "y is shared between the genomes within each metric row, so the\nAt vs maize magnitude difference is visible. ",
         "Read the EFFECT SIZE - at n = 100 runs even tiny effects reach significance. All At contrasts exclude 0:\n",
         "cleaning raises At run-to-run agreement (WD-Pre ARI +0.25) and WD beats ND (+0.07). ",
         "maize WD-Pre and ND-Pre exclude 0 but are ~7x smaller and negative\n",
         "(-0.035) - the minor/major asymmetry seen independently in Fig S7 (per-cell FRiP) and Part 3 - ",
         "while maize WD-ND does NOT exclude 0: the two cleaning modes\nare indistinguishable on maize."),
       x = NULL, y = "difference between stages (with 95% CI)")

# =============================================================================
# PANEL C -- per-meta-cell annotation support (former Part-3 support panel)
# =============================================================================
# One point = ONE CONSENSUS META-CELL -- descriptive only (see header). The
#   At violins are kernel densities over 14 points: read the points, not the
#   shape; the violin is drawn for visual parity with maize, nothing more.
k_at <- sup_med[genome == "At", unique(n_mc)]
pC <- ggplot(sup, aes(stage, label_support)) +
  geom_violin(aes(fill = stage), colour = "grey35", linewidth = 0.3,
              width = 0.85, trim = TRUE, alpha = 0.5) +
  geom_jitter(aes(size = n_cells),
              position = position_jitter(width = 0.13, height = 0, seed = 1),
              colour = "grey20", alpha = 0.5, stroke = 0) +
  geom_errorbar(data = sup_med, aes(y = median_support, ymin = median_support,
                                    ymax = median_support),
                width = 0.62, colour = "grey10", linewidth = 0.6) +
  geom_text(data = sup_med, aes(y = median_support,
                                label = sprintf("%.2f", median_support)),
            vjust = -0.7, hjust = -0.15, size = 2.9, fontface = "bold",
            colour = "grey10") +
  geom_text(data = sup_med, aes(y = -0.045, label = sprintf("n=%d", n_mc)),
            size = 2.3, colour = "grey40") +
  facet_wrap(~ glab, nrow = 1) +
  scale_fill_manual(values = ST_COL, guide = "none") +
  scale_x_discrete(labels = STAGE_LAB) +
  scale_size_continuous(range = c(0.35, 2.6), name = "cells in\nmeta-cell") +
  scale_y_continuous(limits = c(-0.07, 1.06), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8.5),
        legend.position = "right", legend.key.size = unit(0.32, "cm"),
        legend.title = element_text(size = 6.4),
        legend.text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7.4, colour = "grey35", lineheight = 1.3)) +
  labs(title = "C. Reproducibility of the meta-cell annotation across the 100 seeds",
       subtitle = paste0(
         "One point = one consensus meta-cell; y = the fraction of the ", N_SEED,
         " SEACells seeds whose own annotation gives that meta-cell the same type as the consensus does.\n",
         "Black bar = the stage median (UNWEIGHTED - a meta-cell stage summary is never weighted by cell count; point size shows the size instead). ",
         "Values reconcile against\nthe upstream sweep summary to <1e-3. ",
         "DESCRIPTIVE ONLY: no test is run across meta-cells - they share one partition and one set of seeds, so they are not independent\n",
         "replicates of a stage; the replication unit for a stage-level claim is the RUN (panel B's bootstrap). ",
         sprintf("At densities are over %d points: read the points, not the shape.", k_at)),
       x = NULL, y = "label support (fraction of the 100 seeds)")

# =============================================================================
# EXPORT + VERIFY
# =============================================================================
written <- character()
save_pair <- function(plot, stem, w, h) {
  for (ext in c("pdf", "png")) {
    f <- file.path(OUTDIR, sprintf("%s.%s", stem, ext))
    # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
    if (ext == "pdf") ggsave(f, plot, width = w, height = h, bg = "white")
    else              ggsave(f, plot, width = w, height = h, dpi = 300, bg = "white")
    written <<- c(written, f)
  }
}
save_pair(pA, "Fig_S8_violins", 8.6, 5.4)
save_pair(pB, "Fig_S8_diffs",   8.6, 4.4)
save_pair(pC, "Fig_S8_support", 8.6, 4.6)
comp <- pA / pB / pC + plot_layout(heights = c(1.35, 1, 1.1))
save_pair(comp, "Fig_S8", 8.6, 14.2)

summ <- pairs[, .(n_pairs = .N, median = round(median(value), 4),
                  q25 = round(quantile(value, .25), 4),
                  q75 = round(quantile(value, .75), 4),
                  mean = round(mean(value), 4)),
              by = .(genome, stage, metric)][order(genome, metric, stage)]
f <- file.path(OUTDIR, "Fig_S8_summary.tsv"); fwrite(summ, f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, "Fig_S8_ci.tsv")
fwrite(ci[order(genome, kind, metric, name),
          .(genome, kind, name, metric, estimate, ci_lo, ci_hi, B, n_pairs)],
       f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, "Fig_S8_support_summary.tsv")
fwrite(sup_med[order(genome, stage),
               .(genome, stage, n_mc, median_support = round(median_support, 4),
                 frac_agree_modal = round(frac_agree, 4))],
       f, sep = "\t")
written <- c(written, f)

cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz <- if (file.exists(f)) file.size(f) else NA_integer_
  good <- !is.na(sz) && sz > (if (grepl("\\.tsv$", f)) 50 else 1000)
  ok <- ok && good
  cat(sprintf("  %-4s %-28s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
cat("\n[done] Fig_S8.{pdf,png} (A/B/C) + panels + 3 TSVs -> ", OUTDIR, "/\n", sep = "")
