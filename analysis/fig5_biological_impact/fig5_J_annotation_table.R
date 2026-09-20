#!/usr/bin/env Rscript
# =============================================================================
# fig5_J_annotation_table.R  -  the rZ annotation of the CONSENSUS meta-cells with each call's support
#   across the 100 SEACells seeds. Its load-bearing product is the per-meta-cell TABLE
#   Fig5_P3B_annotation_metacells.tsv (the seed-support gate of fig5_K_examples.R, panel K).
#   Its three figures (rZ heatmap "B", support violins "C", marker-gene rZ matrix "D") are NOT
#   manuscript panels: the support panel was moved to Fig S8C (analysis/supplementary/figS8.R).
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/
#           rZ_annotation/consensus/<prefix>.{permetacell_type_rZ.geom,seed_label_sweep.metacell,metacell_table,permetacell_marker_rZ.geom}.tsv
#           rZ_annotation/consensus/{maize,At}.seed_label_sweep.summary.tsv, consensus/<prefix>.consensus.metacells.tsv
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P3B_annotation.{pdf,png} + _metacells.tsv + _types.tsv, Fig5_P3C_support.* + _summary.tsv, Fig5_P3D_marker_rZ.* + _genes.tsv
# Run     Rscript analysis/fig5_biological_impact/fig5_J_annotation_table.R   (repo root, ~20 s; run BEFORE fig5_K_examples.R)
# =============================================================================
#
# WHAT rZ IS (upstream 4_3p): a RECIPROCAL z-score on the raw `.plate.perkb` gene matrix
# (accessibility normalised per kb of gene length -- NOT the 4_3f-smoothed matrix, because
# pooling cells into a meta-cell is already smoothing). Zi = z across meta-cells, Zj = z across
# genes; `geom` (the primary metric) is their AND-like geometric combination, `euclid` the
# other. Per meta-cell it yields one score per cell type, and the argmax is the call --
# structurally the same table as the Leiden `cluster_annotation.tsv`.
#
# PANEL FORM: per arm, meta-cells on x GROUPED BY THEIR CALL (then by descending rZ within a
# group), types on y in a fixed order, so the calls trace a staircase down the diagonal and
# the full rZ vector behind each call stays visible -- that vector is the evidence for how
# specific the call is. An orange segment marks each call; a bar track above gives its seed
# support.
#
# THE ONE THING NOT TO OVER-READ -- THESE ARE ARGMAX CALLS. Every meta-cell gets a type because
#   argmax always returns something (the Leiden cluster annotation has exactly the same
#   property). The expression-matched gene null (upstream 4_3za, bgnull_consensus) finds no
#   significant meta-cell x type pair at q < 0.05 in five of the six objects (0 of 252 in every
#   At object, 0 of 5,720 in maize PreClean and WD) and 2 of 5,720 in maize ND (q = 0.048); its
#   fixture gate passes (planted set called, expression-matched and high-expression decoys
#   rejected) and the p floor clears the BH rank-1 bar, so the near-zero is readable, not
#   underpowered.
#   => Use these labels as DESCRIPTIVE annotation, weighted by support. Never write "this
#      meta-cell IS type T" as a tested claim.
#
# OTHER TRAPS:
#   Zi CAP. Zi is a z across meta-cells, so |Zi| <= (N-1)/sqrt(N): 3.47 on At (only 14
#     meta-cells) vs 16.85 on maize. On At that cap makes `euclid` Zj-dominated -- `geom` is
#     primary. METRIC is a constant here; do not switch it to euclid for At without re-reading
#     the upstream 4_3p warning.
#   DEPTH RESIDUAL. Meta-cells are scaled to a common total, which removes the SCALE difference
#     exactly but NOT the PRECISION difference (the ACR matrix was rarefied upstream, not this
#     gene matrix). A deeper stage yields intrinsically peakier profiles and a mildly inflated
#     Zi, so small cross-stage rZ differences are not biology.
#   rZ is NOT comparable across genomes -- different marker panels (275 At / 221 maize markers)
#     and different N. The fill scale is per genome row.
#   Partitions differ per stage: meta-cell IDs, orderings and x positions are per panel. Never
#     compare cMC-7 across stages.
#   wd and nd are INDEPENDENT treatments of the same raw input, NOT a chain.
#   Support is IN-SAMPLE in the same sense as F: the consensus was built from the same 100
#     seeds the support is measured against. The non-circular stability check remains the
#     split-half number stamped on the consensus-matrix panel (I).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source("analysis/_helpers/fig5_part2_helpers.R")   # SOC, PLATE, CFG, STAGES

STEM   <- "Fig5_P3B_annotation"    # panel B: the rZ heatmap
STEM_C <- "Fig5_P3C_support"       # panel C: the support distribution (now Fig S8C)
STEM_D <- "Fig5_P3D_marker_rZ"     # panel D: the marker-gene matrix
MC     <- file.path(PLATE, "step5_metacell")
RZDIR  <- file.path(MC, "rZ_annotation", "consensus")
CONS   <- file.path(MC, "consensus")
METRIC <- "geom"      # primary (AND). See the Zi-cap note above.
N_SEED <- 100L

G_LEV    <- vapply(CFG, `[[`, character(1), "label")
S_LEV    <- names(STAGES)
STAGE_UP <- c(PreClean = "Pre", wd = "wd", nd = "nd")
FIXED    <- c("metacell", "seed", "cluster", "n_cells", "top_type", "top_rZ", "margin")
TOL      <- 1e-4      # upstream values are rounded to 4 dp

# rZ fill: light -> dark (mako). Support: the SAME white->maroon ramp the consensus-matrix
# panel uses, because it is the same quantity -- a fraction of the 100 seeds.
RZ_PAL  <- viridisLite::viridis(256, option = "G", direction = -1)
SUP_PAL <- colorRampPalette(c("#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33",
                              "#7F0000", "#250000"))(256)
CALLCOL <- "#FF7F00"   # the call marker: orange reads against mako at every level
SEPCOL  <- "grey88"    # separator between call groups

prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)

# =============================================================================
# LOAD + VALIDATE one arm
# =============================================================================
read_arm <- function(g, s) {
  pf <- prefix_for(g, s)
  f_rz  <- file.path(RZDIR, sprintf("%s.permetacell_type_rZ.%s.tsv", pf, METRIC))
  f_sw  <- file.path(RZDIR, sprintf("%s.seed_label_sweep.metacell.tsv", pf))
  f_tab <- file.path(RZDIR, sprintf("%s.metacell_table.tsv", pf))
  f_mc  <- file.path(CONS,  sprintf("%s.consensus.metacells.tsv", pf))
  for (f in c(f_rz, f_sw, f_tab, f_mc))
    if (!file.exists(f)) stop("missing upstream input: ", f)

  rz <- fread(f_rz); sw <- fread(f_sw); tb <- fread(f_tab); mcs <- fread(f_mc)

  # ---- schema ---------------------------------------------------------------
  if (length(setdiff(FIXED, names(rz))) > 0)
    stop(pf, ": rZ table is missing ", paste(setdiff(FIXED, names(rz)), collapse = ", "))
  tcol <- setdiff(names(rz), FIXED)
  if (length(tcol) < 2L) stop(pf, ": no type columns in the rZ table")

  # ---- the call really is the argmax of the row (and margin the runner-up gap)
  M  <- as.matrix(rz[, ..tcol])
  o2 <- t(apply(M, 1L, function(v) order(v, decreasing = TRUE)[1:2]))
  rz[, `:=`(chk_top  = tcol[o2[, 1]],
            chk_rZ   = M[cbind(seq_len(.N), o2[, 1])],
            sec_type = tcol[o2[, 2]],
            sec_rZ   = M[cbind(seq_len(.N), o2[, 2])])]
  if (rz[chk_top != top_type, .N] > 0)
    stop(pf, ": top_type is not the argmax of the type columns")
  if (rz[abs(chk_rZ - top_rZ) > TOL, .N] > 0)
    stop(pf, ": top_rZ disagrees with the argmax value")
  if (rz[abs((top_rZ - sec_rZ) - margin) > TOL, .N] > 0)
    stop(pf, ": margin is not top_rZ - second_rZ")

  # ---- the label sweep annotates the SAME consensus meta-cells --------------
  if (!identical(sort(sw$metacell), sort(rz$metacell)))
    stop(pf, ": seed_label_sweep and rZ table cover different meta-cells")
  if (uniqueN(sw$stage) != 1L || sw$stage[1] != STAGE_UP[[s]])
    stop(pf, ": seed_label_sweep stage is not ", STAGE_UP[[s]])
  d <- merge(rz, sw[, .(metacell, sw_n = n_cells, consensus_label, label_support,
                        modal_seed_label, modal_seed_freq, labels_agree)],
             by = "metacell")
  if (nrow(d) != nrow(rz)) stop(pf, ": meta-cell join lost rows")
  if (d[consensus_label != top_type, .N] > 0)
    stop(pf, ": the sweep's consensus_label disagrees with the rZ argmax call")
  if (d[sw_n != n_cells, .N] > 0)  stop(pf, ": n_cells disagree between rZ and sweep")
  if (d[!is.finite(label_support) | label_support < 0 | label_support > 1, .N] > 0)
    stop(pf, ": label_support outside [0, 1]")

  # ---- cross-check against the consensus itself (same partition, same sizes) -
  d <- merge(d, tb[, .(metacell, cluster_tb = cluster, purity, pooled_perkb_mass)],
             by = "metacell")
  d[, cmc := as.integer(sub("^cMC-", "", metacell))]
  if (anyNA(d$cmc)) stop(pf, ": cannot parse a cMC-<n> meta-cell id")
  d <- merge(d, mcs[, .(cmc = consensus_mc, cons_n = n_cells, conf_meanF)], by = "cmc")
  if (nrow(d) != nrow(rz))
    stop(pf, ": meta-cells do not match the published consensus partition")
  if (d[cons_n != n_cells, .N] > 0)
    stop(pf, ": meta-cell sizes disagree with consensus.metacells.tsv")

  d[, `:=`(genome = g, stage = s)]
  list(d = d[], types = tcol)
}

# =============================================================================
# LOAD ALL SIX ARMS
# =============================================================================
arms <- list(); types_by_g <- list()
for (g in names(CFG)) for (s in S_LEV) {
  a <- read_arm(g, s)
  arms[[prefix_for(g, s)]] <- a$d
  # the type panel must be identical across the stages of a genome, else the
  # y axis is not shared and the three panels are not comparable
  if (is.null(types_by_g[[g]])) types_by_g[[g]] <- sort(a$types)
  else if (!identical(types_by_g[[g]], sort(a$types)))
    stop(g, ": the marker type panel differs between stages")
}

# =============================================================================
# ORDER + LONG FORM
# =============================================================================
long <- list(); mcout <- list(); report <- list()

for (g in names(CFG)) for (s in S_LEV) {
  pf <- prefix_for(g, s); d <- arms[[pf]]; tl <- types_by_g[[g]]
  d[, ypos := match(top_type, tl)]
  setorder(d, ypos, -top_rZ)
  d[, xpos := seq_len(.N)]
  if (anyNA(d$ypos)) stop(pf, ": a call is not in the type panel")

  L <- melt(d[, c("xpos", "metacell", ..tl)],
            id.vars = c("xpos", "metacell"),
            variable.name = "type", value.name = "rZ", variable.factor = FALSE)
  L[, `:=`(ty = match(type, tl), genome = g, stage = s)]
  if (anyNA(L$ty)) stop(pf, ": melted a column that is not a type")
  long[[pf]] <- L

  mcout[[pf]] <- d[, .(genome, stage, metacell, xpos, n_cells, cluster = cluster_tb,
                       purity, conf_meanF, pooled_perkb_mass,
                       top_type, top_rZ, second_type = sec_type, second_rZ = sec_rZ,
                       margin, label_support, modal_seed_label, modal_seed_freq,
                       labels_agree)]

  report[[pf]] <- data.table(
    genome = g, stage = s, metacells = nrow(d), cells = d[, sum(n_cells)],
    types_called = uniqueN(d$top_type), types_avail = length(tl),
    rZ_med = round(d[, median(top_rZ)], 3), rZ_max = round(d[, max(top_rZ)], 3),
    margin_med = round(d[, median(margin)], 3),
    support_med = round(d[, median(label_support)], 3),
    pct_agree_modal = round(100 * d[, mean(labels_agree)], 1))
}

mcall  <- rbindlist(mcout)
lall   <- rbindlist(long)

# Fill scale, per genome. The raw range is set by a handful of very high cells
# (maize median rZ 0.84 vs max 4.28), which flattens everything else to one pale
# tone and hides exactly what the panel is for -- how specific each call is
# against the rest of its column. So the scale is CLIPPED at the 99th percentile
# and values above it are squished to the top colour: a display choice on the
# COLOUR only, never on the data (the TSVs carry the full values, and the true
# max is printed on each colourbar).
RZ_CLIP <- 0.99
rZ_rng  <- lall[, .(lo = min(rZ), hi = quantile(rZ, RZ_CLIP), true_hi = max(rZ)),
                by = genome]

# per genome x stage x type: how much of the object each call accounts for
types_out <- mcall[, .(n_metacells = .N, n_cells = sum(n_cells),
                       mean_rZ = round(mean(top_rZ), 4),
                       mean_margin = round(mean(margin), 4),
                       mean_support = round(mean(label_support), 4)),
                   by = .(genome, stage, type = top_type)][order(genome, stage, -n_cells)]
types_out[, pct_cells := round(100 * n_cells / sum(n_cells), 2), by = .(genome, stage)]

# =============================================================================
# PANELS
# =============================================================================
mk_panel <- function(g, s, ylab_on) {
  pf <- prefix_for(g, s)
  d <- mcall[genome == g & stage == s]; L <- long[[pf]]
  tl <- types_by_g[[g]]; n <- nrow(d)
  rng <- unlist(rZ_rng[genome == g, .(lo, hi)])

  # call groups: contiguous by construction (x was ordered by the call)
  grp <- d[, .(x0 = min(xpos), x1 = max(xpos), ypos = match(top_type[1], tl)),
           by = top_type]
  if (grp[, sum(x1 - x0 + 1L)] != n) stop(pf, ": call groups are not contiguous")
  seps <- grp$x1[-nrow(grp)] + 0.5

  # --- support track (position, not colour: it is a quantity worth reading) ---
  # The arm title + stamp live on THIS plot, not on a plot_annotation() of the
  # inner patchwork: patchwork DROPS a nested patchwork's annotations when it is
  # composed into a parent, so an inner plot_annotation silently loses the
  # per-panel numbers.
  medsup <- d[, median(label_support)]
  sub <- paste0(
    sprintf("%s meta-cells, %s cells | %d of %d types called",
            format(n, big.mark = ","), format(d[, sum(n_cells)], big.mark = ","),
            uniqueN(d$top_type), length(tl)),
    "\n",
    sprintf("median support %.2f | %.0f%% of calls match the modal seed label",
            medsup, 100 * d[, mean(labels_agree)]))
  top <- ggplot(d, aes(x = xpos, y = label_support)) +
    geom_col(width = 1, fill = "grey45") +
    geom_hline(yintercept = medsup, colour = CALLCOL, linewidth = 0.35) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", ".5", "1"), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    labs(y = "support", title = sprintf("%s - %s", G_LEV[[g]], s), subtitle = sub) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_text(size = 5.2, colour = "grey30"),
          axis.text.y = element_text(size = 4.6, colour = "grey45"),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5),
          plot.subtitle = element_text(size = 5.4, colour = "grey35", hjust = 0.5,
                                       lineheight = 1.25,
                                       margin = margin(b = 3)),
          plot.margin = margin(2, 3, 0, 3))
  if (!ylab_on) top <- top + theme(axis.title.y = element_blank(),
                                   axis.text.y = element_blank())

  # --- the rZ heatmap --------------------------------------------------------
  hm <- ggplot(L, aes(x = xpos, y = ty, fill = rZ)) +
    geom_raster() +
    geom_vline(xintercept = seps, colour = SEPCOL, linewidth = 0.15) +
    geom_segment(data = grp, inherit.aes = FALSE,
                 aes(x = x0 - 0.5, xend = x1 + 0.5, y = ypos, yend = ypos),
                 colour = CALLCOL, linewidth = 0.55) +
    scale_fill_gradientn(colours = RZ_PAL, limits = rng, guide = "none",
                         oob = scales::squish) +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(tl), labels = sub("^[A-Za-z]+:", "", tl),
                       expand = c(0, 0)) +
    labs(x = sprintf("%s consensus meta-cells, grouped by call", n)) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size = 5.2, colour = "grey30"),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 4.4, colour = "grey25"),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3),
          plot.margin = margin(0, 3, 2, 3))
  if (!ylab_on) hm <- hm + theme(axis.text.y = element_blank())

  (top / hm) + plot_layout(heights = c(0.20, 1))
}

# hand-built colourbar per genome row (rZ is not comparable across genomes)
mk_cbar <- function(g) {
  rng   <- unlist(rZ_rng[genome == g, .(lo, hi)])
  thi   <- rZ_rng[genome == g, true_hi]
  labs4 <- sprintf("%.2f", seq(rng[1], rng[2], length.out = 4))
  labs4[4] <- paste0(">=", labs4[4])   # squished: the top colour is an open bin
  ggplot() +
    annotation_raster(as.raster(matrix(rev(RZ_PAL), ncol = 1)),
                      xmin = 0, xmax = 1, ymin = 0, ymax = 1, interpolate = TRUE) +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = NA, colour = "grey25", linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1.35, y = c(0, 1/3, 2/3, 1),
             yend = c(0, 1/3, 2/3, 1), linewidth = 0.3, colour = "grey25") +
    annotate("text", x = 1.6, y = c(0, 1/3, 2/3, 1), label = labs4,
             hjust = 0, size = 2.1, colour = "grey20") +
    coord_cartesian(xlim = c(-0.2, 5.2), ylim = c(-0.03, 1.03), expand = FALSE) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 5.6, colour = "grey20", hjust = 0),
          plot.margin = margin(30, 2, 34, 2)) +
    labs(title = sprintf("rZ (%s)\n%s\ncolour clipped at\np%g; max %.2f", METRIC,
                         if (g == "TAIR10") "Arabidopsis" else "maize",
                         100 * RZ_CLIP, thi))
}

P <- list()
for (g in names(CFG)) for (s in S_LEV)
  P[[prefix_for(g, s)]] <- mk_panel(g, s, ylab_on = (s == "PreClean"))

# =============================================================================
# REPORT
# =============================================================================
cat("=== FIG 5 PART 3 | rZ annotation of the CONSENSUS meta-cells (table for panel K) ===\n\n")
cat("rZ source   : ", RZDIR, "  (metric ", METRIC, ")\n", sep = "")
cat("support     : seed_label_sweep, ", N_SEED, " seeds\n\n", sep = "")
print(rbindlist(report)[order(genome, factor(stage, levels = S_LEV))])
cat("\n-- calls accounting for the most cells, per arm --\n")
print(types_out[, head(.SD, 3), by = .(genome, stage)][
  , .(genome, stage, type, n_metacells, n_cells, pct_cells, mean_rZ, mean_support)])
cat("\n  ARGMAX CALLS. The expression-matched gene null finds no significant meta-cell x\n",
    "  type pair at q<0.05 in five of the six objects (0/252 in every At object, 0/5720 in\n",
    "  maize PreClean and WD) and 2/5720 in maize ND; the fixture gate passes and the p floor\n",
    "  clears the BH rank-1 bar, so the near-zero is readable. Treat these labels as\n",
    "  descriptive annotation weighted by support, never as tested identity.\n",
    "  rZ is not comparable across genomes (different marker panels, different N) and\n",
    "  carries a depth residual, so small cross-stage differences are not biology.\n",
    sep = "")

# =============================================================================
# ASSEMBLY: maize row | At row, one colourbar per row
# =============================================================================
k <- function(g, s) prefix_for(g, s)
comp <- (P[[k("B73v5", "PreClean")]]  | P[[k("B73v5", "wd")]]  | P[[k("B73v5", "nd")]]  | mk_cbar("B73v5")) /
        (P[[k("TAIR10", "PreClean")]] | P[[k("TAIR10", "wd")]] | P[[k("TAIR10", "nd")]] | mk_cbar("TAIR10")) +
  plot_layout(widths = c(1, 1, 1, 0.17)) +
  plot_annotation(
    title = "B. Annotating the consensus meta-cells - reciprocal z-score (rZ) on per-kb gene accessibility",
    subtitle = paste0(
      "The Leiden-cluster annotation strategy applied to the CONSENSUS meta-cells - the partition built from the ", N_SEED,
      " SEACells seeds via the co-assignment matrix F of panel A.\n",
      "HEATMAP: rZ of every cell type (y) for every meta-cell (x). The ORANGE segment marks each meta-cell's call; meta-cells are grouped by that call and ordered by\n",
      "descending rZ within a group, so the calls trace a staircase and the rest of each column shows how specific the call was. BAR TRACK: that call's SUPPORT =\n",
      "the fraction of the ", N_SEED, " seeds whose own annotation gives the meta-cell the same type (orange line = arm median). Colour is clipped at the 99th percentile.\n",
      "*** ARGMAX CALLS. Every meta-cell gets a type, exactly as the cluster annotation does. The expression-matched gene null finds no significant meta-cell x type\n",
      "pair at q<0.05 in five of the six objects (0/252 in every At object, 0/5720 in maize PreClean and WD) and 2/5720 in maize ND - fixture gate passes, p floor\n",
      "clears the BH bar. The labels are DESCRIPTIVE annotation to be weighted by support, never a tested identity claim for any pair. ***\n",
      "rZ is NOT comparable across genomes (275 At vs 221 maize markers; 14 vs 286 units) - the fill scale is per row. On At the Zi cap is (N-1)/sqrt(N) = 3.47, which is\n",
      "why `geom` is primary. A depth residual remains, so small cross-stage rZ differences are not biology. Partitions and x positions are PER PANEL; wd and nd are\n",
      "independent treatments, not a chain. Support is in-sample against the seeds that built the consensus - the non-circular check is panel A's split-half number."),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 6.3, colour = "grey35",
                                               lineheight = 1.28)))

# =============================================================================
# PANEL C: the support distribution, all three stages, both genomes (now Fig S8C)
# =============================================================================
# STATISTICAL UNIT. One point = ONE CONSENSUS META-CELL, which is the unit
# `label_support` is defined at (n = 14 At / 286 maize per stage). This panel is
# DESCRIPTIVE ONLY -- there is NO test across meta-cells anywhere here, for the
# same reason S8 runs no test across seed pairs: meta-cells within a stage share
# one partition and one set of 100 seeds, so they are not independent replicates
# of the stage. A stage-level claim's replication unit is the RUN. (Same class
# of error as the retired 45-grid-combo t-test and the meta-cell-unit tests.)
#
# The median is UNWEIGHTED. Never weight a meta-cell stage summary by n_cells --
# point SIZE shows the size so the reader can see the weighting question without
# the statistic silently answering it.
#
# At violins are kernel densities over 14 points. Read the points, not the
# shape. The violin is drawn for visual parity with maize, nothing more.
set.seed(1)   # deterministic jitter

# --- reconcile against the upstream stage summary before plotting anything ----
GEN_KEY <- c(B73v5 = "maize", TAIR10 = "At")
sup_med <- mcall[, .(n_mc = .N, median_support = median(label_support),
                     frac_agree = mean(labels_agree)), by = .(genome, stage)]
for (g in names(CFG)) {
  f <- file.path(RZDIR, sprintf("%s.seed_label_sweep.summary.tsv", GEN_KEY[[g]]))
  if (!file.exists(f)) stop("missing upstream sweep summary: ", f)
  up <- fread(f)
  if (!all(up$n_seeds == N_SEED)) stop(f, ": sweep is not over ", N_SEED, " seeds")
  # upstream labels the stages Pre / wd / nd; this figure uses the pipeline's
  # PreClean / wd / nd. Join on the mapped name, never on the raw string.
  chk <- merge(sup_med[genome == g][, stage_up := STAGE_UP[stage]],
               up[, .(stage_up = stage, up_med = median_mc_label_support,
                      up_agree = frac_mc_label_matches_modal)],
               by = "stage_up")
  if (nrow(chk) != length(S_LEV))
    stop(g, ": upstream sweep summary does not cover all three stages")
  d1 <- chk[, max(abs(median_support - up_med))]
  d2 <- chk[, max(abs(frac_agree - up_agree))]
  if (d1 > 1e-3 || d2 > 1e-3)
    stop(sprintf("%s: computed support disagrees with upstream sweep summary (med %.2e, agree %.2e)",
                 g, d1, d2))
  message(sprintf("  [reconcile] %-6s vs %s: max|diff| median %.1e, agree %.1e",
                  g, basename(f), d1, d2))
}

ST_COL <- c(PreClean = "#FF83FA", wd = "#43CD80", nd = "#E69F00")   # as S7 / S8
vd <- copy(mcall)[, `:=`(stage = factor(stage, levels = S_LEV),
                         gl    = factor(G_LEV[genome], levels = unname(G_LEV)))]
vm <- sup_med[, `:=`(stage = factor(stage, levels = S_LEV),
                     gl    = factor(G_LEV[genome], levels = unname(G_LEV)))]

panelC <- ggplot(vd, aes(x = stage, y = label_support)) +
  geom_violin(aes(fill = stage), colour = "grey35", linewidth = 0.3,
              width = 0.85, trim = TRUE, alpha = 0.5) +
  geom_jitter(aes(size = n_cells), width = 0.13, height = 0,
              colour = "grey20", alpha = 0.5, stroke = 0) +
  geom_errorbar(data = vm, aes(y = median_support, ymin = median_support,
                               ymax = median_support),
                width = 0.62, colour = "grey10", linewidth = 0.6) +
  geom_text(data = vm, aes(y = median_support, label = sprintf("%.2f", median_support)),
            vjust = -0.7, hjust = -0.15, size = 2.9, fontface = "bold", colour = "grey10") +
  geom_text(data = vm, aes(y = -0.045, label = sprintf("n=%d", n_mc)),
            size = 2.3, colour = "grey40") +
  facet_wrap(~ gl, nrow = 1) +
  scale_fill_manual(values = ST_COL, guide = "none") +
  scale_size_continuous(range = c(0.35, 2.6), name = "cells in\nmeta-cell") +
  scale_y_continuous(limits = c(-0.07, 1.06), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "label support (fraction of the 100 seeds)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold", size = 8.4),
        axis.title.y = element_text(size = 7.6),
        legend.position = "right", legend.key.size = unit(0.32, "cm"),
        legend.title = element_text(size = 6.4), legend.text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 6.5, colour = "grey35", lineheight = 1.3)) +
  labs(title = "C. Reproducibility of the meta-cell annotation across the 100 seeds",
       subtitle = paste0(
         "One point = one consensus meta-cell; y = the fraction of the ", N_SEED,
         " SEACells seeds whose own annotation gives that meta-cell the same type as the consensus does.\n",
         "Black bar = the stage median (UNWEIGHTED - a meta-cell stage summary is never weighted by cell count; point size shows the size instead). ",
         "Values reconcile\nagainst the upstream sweep summary to <1e-3. ",
         "*** DESCRIPTIVE ONLY: no test is run across meta-cells. They share one partition and one set of seeds, so they are not\n",
         "independent replicates of a stage - the replication unit for a stage-level claim is the RUN (that is what S8's bootstrap uses). *** ",
         "At densities are over 14 points:\nread the points, not the shape. wd and nd are independent treatments of the same raw input, not a chain."))

cat("\n-- panel C: label support by stage (unit = meta-cell, DESCRIPTIVE) --\n")
print(sup_med[order(genome, factor(stage, levels = S_LEV))][
  , .(genome, stage, n_mc, median_support = round(median_support, 3),
      pct_agree_modal = round(100 * frac_agree, 1))])

# =============================================================================
# PANEL D: the standardised gene accessibility BEHIND the calls
# =============================================================================
# Panel B shows one aggregated score per (meta-cell, type). This is the layer
# under it: the same reciprocal z, per MARKER GENE -- `permetacell_marker_rZ`,
# a genes x meta-cells matrix. x reuses panel B's ORDERING EXACTLY (meta-cells
# grouped by call), so the two panels align column for column; y is the marker
# panel grouped by the type each gene marks, in panel B's type order.
#
# It is standardised accessibility, NOT raw. The per-kb gene matrix itself is a
#   large .rds and is never loaded here (fig5_J_access_cache.R does that). What
#   is plotted is the gene-level rZ that panel B aggregates over.
# MARKER PANEL != SCORED TYPES. Some types carry too few markers to be scored:
#   maize has 221 markers over 39 type labels but only 20 types are scored
#   (194 genes); At has 275 markers over 20 labels, 18 scored. The unscored
#   genes are KEPT and parked in a separate trailing block -- dropping them
#   silently would overstate how completely the panel is covered -- but they can
#   never form a diagonal block because no meta-cell can be called their type.
# GENE SET IS NOT IDENTICAL ACROSS STAGES: At/nd is missing AT2G08610
#   (At:phloem). y must be the same genes in all three panels of a row or they
#   are not comparable, so the row is INTERSECTED to the shared set and the
#   drop is reported. Never plot a per-stage gene set on a shared axis.
# Gene order is symmetric across stages by construction (score = the mean over
#   the three stages of the gene's max rZ), so no stage sets the layout.
BOXCOL_D <- "#2166AC"

read_markers <- function(g, s) {
  f <- file.path(RZDIR, sprintf("%s.permetacell_marker_rZ.%s.tsv",
                                prefix_for(g, s), METRIC))
  if (!file.exists(f)) stop("missing upstream marker matrix: ", f)
  m <- fread(f)
  idc <- c("geneID", "name", "type_label")
  if (length(setdiff(idc, names(m))) > 0)
    stop(basename(f), ": marker matrix is missing an id column")
  mc <- setdiff(names(m), idc)
  exp_mc <- mcall[genome == g & stage == s, metacell]
  if (!setequal(mc, exp_mc))
    stop(basename(f), ": marker matrix meta-cells differ from the rZ table")
  if (anyDuplicated(m$geneID) > 0) stop(basename(f), ": duplicated geneID")
  # melt warns when some meta-cell columns are integer and others double; that
  # coercion is lossless, but a non-numeric column would silently become NA.
  L <- suppressWarnings(
    melt(m, id.vars = idc, variable.name = "metacell", value.name = "rZ",
         variable.factor = FALSE))
  if (!is.numeric(L$rZ) || anyNA(L$rZ))
    stop(basename(f), ": marker matrix has non-numeric or missing rZ values")
  L[, `:=`(genome = g, stage = s)][]
}

mk_long <- rbindlist(lapply(names(CFG), function(g)
  rbindlist(lapply(S_LEV, function(s) read_markers(g, s)))))

# --- intersect each genome's gene set across its three stages ------------------
gene_keep <- list(); lost_note <- character()
for (g in names(CFG)) {
  sets   <- lapply(S_LEV, function(s) mk_long[genome == g & stage == s, unique(geneID)])
  shared <- Reduce(intersect, sets)
  lost   <- setdiff(Reduce(union, sets), shared)
  gene_keep[[g]] <- shared
  if (length(lost) > 0) {
    lab <- unique(mk_long[genome == g & geneID %in% lost, paste0(geneID, " (", type_label, ")")])
    lost_note <- c(lost_note, sprintf("%s is missing %s",
                                      G_LEV[[g]], paste(lab, collapse = ", ")))
    message(sprintf("  [panel D] %s: %d gene(s) absent from at least one stage, dropped from the shared y axis: %s",
                    g, length(lost), paste(lab, collapse = ", ")))
  }
}
keep_dt <- rbindlist(lapply(names(CFG), function(g)
  data.table(genome = g, geneID = gene_keep[[g]])))
mk_long <- merge(mk_long, keep_dt, by = c("genome", "geneID"))

# --- gene ordering: type group (panel B order, unscored last), then strength ---
gord <- list(); ggrp <- list()
for (g in names(CFG)) {
  tl <- types_by_g[[g]]
  gi <- unique(mk_long[genome == g, .(geneID, name, type_label)])
  # symmetric strength: mean over stages of the gene's max rZ across meta-cells
  st <- mk_long[genome == g, .(mx = max(rZ)), by = .(geneID, stage)][
    , .(strength = mean(mx)), by = geneID]
  gi <- merge(gi, st, by = "geneID")
  gi[, scored := type_label %in% tl]
  gi[, tord := fifelse(scored, match(type_label, tl), length(tl) + 1L)]
  setorder(gi, tord, type_label, -strength)
  gi[, ypos := seq_len(.N)]
  gord[[g]] <- gi
  ggrp[[g]] <- gi[, .(y0 = min(ypos), y1 = max(ypos), n = .N,
                      scored = scored[1]),
                  by = .(grp = fifelse(scored, type_label, "(types with too few markers to score)"))]
}

mk_panelD <- function(g, s, ylab_on) {
  gi <- gord[[g]]; grp <- ggrp[[g]]; tl <- types_by_g[[g]]
  L  <- merge(mk_long[genome == g & stage == s], gi[, .(geneID, ypos)], by = "geneID")
  d  <- mcall[genome == g & stage == s]
  L  <- merge(L, d[, .(metacell, xpos)], by = "metacell")
  if (nrow(L) != nrow(gi) * nrow(d)) stop(prefix_for(g, s), ": marker matrix is not complete")
  n  <- nrow(d); ng <- nrow(gi)
  rng <- unlist(rZ_rngD[genome == g, .(lo, hi)])

  # call groups on x (same rule as panel B) and their matching gene block on y
  cg <- d[, .(x0 = min(xpos), x1 = max(xpos)), by = top_type]
  diagb <- merge(cg, grp[scored == TRUE], by.x = "top_type", by.y = "grp")

  ggplot(L, aes(x = xpos, y = ypos, fill = rZ)) +
    geom_raster() +
    geom_hline(yintercept = grp$y1[-nrow(grp)] + 0.5, colour = "grey80", linewidth = 0.12) +
    geom_vline(xintercept = cg[order(x1), x1][-nrow(cg)] + 0.5,
               colour = "grey80", linewidth = 0.12) +
    geom_rect(data = diagb, inherit.aes = FALSE,
              aes(xmin = x0 - 0.5, xmax = x1 + 0.5, ymin = y0 - 0.5, ymax = y1 + 0.5),
              fill = NA, colour = BOXCOL_D, linewidth = 0.3) +
    scale_fill_gradientn(colours = RZ_PAL, limits = rng, guide = "none",
                         oob = scales::squish) +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, ng + 0.5), expand = c(0, 0),
                       breaks = grp[, (y0 + y1) / 2],
                       labels = sub("^[A-Za-z]+:", "", grp$grp)) +
    labs(x = sprintf("%s meta-cells, grouped by call (panel B order)", n),
         y = NULL, title = sprintf("%s - %s", G_LEV[[g]], s),
         subtitle = sprintf("%d marker genes shared by all 3 stages | %d in a scored type",
                            ng, gi[, sum(scored)])) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size = 5.2, colour = "grey30"),
          axis.text.x = element_blank(),
          axis.text.y = if (ylab_on) element_text(size = 4.2, colour = "grey25")
                        else element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5),
          plot.subtitle = element_text(size = 5.2, colour = "grey35", hjust = 0.5),
          plot.margin = margin(2, 4, 2, 4))
}

rZ_rngD <- mk_long[, .(lo = min(rZ), hi = quantile(rZ, RZ_CLIP), true_hi = max(rZ)),
                   by = genome]
PD <- list()
for (g in names(CFG)) for (s in S_LEV)
  PD[[prefix_for(g, s)]] <- mk_panelD(g, s, ylab_on = (s == "PreClean"))

mk_cbarD <- function(g) {
  rng <- unlist(rZ_rngD[genome == g, .(lo, hi)]); thi <- rZ_rngD[genome == g, true_hi]
  l4 <- sprintf("%.2f", seq(rng[1], rng[2], length.out = 4)); l4[4] <- paste0(">=", l4[4])
  ggplot() +
    annotation_raster(as.raster(matrix(rev(RZ_PAL), ncol = 1)),
                      xmin = 0, xmax = 1, ymin = 0, ymax = 1, interpolate = TRUE) +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = NA,
             colour = "grey25", linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1.35, y = c(0, 1/3, 2/3, 1),
             yend = c(0, 1/3, 2/3, 1), linewidth = 0.3, colour = "grey25") +
    annotate("text", x = 1.6, y = c(0, 1/3, 2/3, 1), label = l4, hjust = 0,
             size = 2.1, colour = "grey20") +
    coord_cartesian(xlim = c(-0.2, 5.2), ylim = c(-0.03, 1.03), expand = FALSE) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 5.6, colour = "grey20", hjust = 0),
          plot.margin = margin(30, 2, 34, 2)) +
    labs(title = sprintf("gene rZ (%s)\n%s\nclipped at p%g\nmax %.2f", METRIC,
                         if (g == "TAIR10") "Arabidopsis" else "maize",
                         100 * RZ_CLIP, thi))
}

compD <- (PD[[k("B73v5", "PreClean")]]  | PD[[k("B73v5", "wd")]]  | PD[[k("B73v5", "nd")]]  | mk_cbarD("B73v5")) /
         (PD[[k("TAIR10", "PreClean")]] | PD[[k("TAIR10", "wd")]] | PD[[k("TAIR10", "nd")]] | mk_cbarD("TAIR10")) +
  plot_layout(widths = c(1, 1, 1, 0.17)) +
  plot_annotation(
    title = "D. The standardised gene accessibility behind the calls - marker-gene rZ per consensus meta-cell",
    subtitle = paste0(
      "The layer under panel B: the SAME reciprocal z, but per MARKER GENE rather than aggregated to a type score. x is panel B's ordering exactly (meta-cells\n",
      "grouped by their call), so the two panels align column for column; y is the marker panel grouped by the type each gene marks, in panel B's type order.\n",
      "BLUE BOXES mark each call group against its own markers - that is where a correct call should put its signal, so read the panel as a block diagonal.\n",
      "Standardised, NOT raw accessibility: the per-kb gene matrix is never loaded here; this is the gene-level rZ that panel B aggregates over. Colour is clipped\n",
      "at the 99th percentile, per genome row - rZ is not comparable across genomes.\n",
      "*** The marker panel is WIDER than the scored type set: maize carries 221 markers over 39 type labels but only 20 types have enough markers to be scored\n",
      "(194 genes), At 275 markers over 20 labels with 18 scored. The unscored genes are kept, parked in the trailing block, and can never form a diagonal box\n",
      "because no meta-cell can be called their type. ",
      if (length(lost_note))
        paste0("The marker set is not identical across stages (", paste(lost_note, collapse = "; "),
               "), so each row is\nintersected to the genes shared by all three of its stages")
      else "Each row is intersected to the genes shared by all three of its stages",
      " - a per-stage gene set on a shared y axis would not be comparable.\nGene order within a group is symmetric across stages by\n",
      "construction (mean over stages of the gene's max rZ), so no stage sets the layout. ***"),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 6.3, colour = "grey35",
                                               lineheight = 1.28)))

cat("\n-- panel D: marker-gene matrix --\n")
print(rbindlist(lapply(names(CFG), function(g)
  data.table(genome = g, genes_shared = length(gene_keep[[g]]),
             genes_scored_type = gord[[g]][, sum(scored)],
             type_labels = uniqueN(gord[[g]]$type_label),
             scored_types = length(types_by_g[[g]]),
             metacells = mcall[genome == g & stage == "wd", .N]))))

# =============================================================================
# EXPORT + VERIFY
# =============================================================================
written <- character()
for (ext in c("pdf", "png")) {
  # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
  f <- file.path(OUTDIR, sprintf("%s.%s", STEM, ext))
  if (ext == "pdf") ggsave(f, comp, width = 13.2, height = 10.4, bg = "white")
  else              ggsave(f, comp, width = 13.2, height = 10.4, dpi = 300, bg = "white")
  written <- c(written, f)
  f <- file.path(OUTDIR, sprintf("%s.%s", STEM_C, ext))
  if (ext == "pdf") ggsave(f, panelC, width = 8.6, height = 4.6, bg = "white")
  else              ggsave(f, panelC, width = 8.6, height = 4.6, dpi = 300, bg = "white")
  written <- c(written, f)
  f <- file.path(OUTDIR, sprintf("%s.%s", STEM_D, ext))
  if (ext == "pdf") ggsave(f, compD, width = 13.2, height = 11.0, bg = "white")
  else              ggsave(f, compD, width = 13.2, height = 11.0, dpi = 300, bg = "white")
  written <- c(written, f)
}
f <- file.path(OUTDIR, paste0(STEM, "_metacells.tsv")); fwrite(mcall, f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, paste0(STEM, "_types.tsv"));     fwrite(types_out, f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, paste0(STEM_C, "_summary.tsv"))
fwrite(sup_med[order(genome, factor(stage, levels = S_LEV))], f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, paste0(STEM_D, "_genes.tsv"))
fwrite(rbindlist(lapply(names(CFG), function(g) copy(gord[[g]])[, genome := g]))[
  , .(genome, geneID, name, type_label, scored, strength = round(strength, 4), ypos)],
  f, sep = "\t")
written <- c(written, f)

cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz   <- if (file.exists(f)) file.size(f) else NA_integer_
  good <- !is.na(sz) && sz > (if (grepl("\\.tsv$", f)) 50 else 1000)
  ok   <- ok && good
  cat(sprintf("  %-4s %-42s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
cat("\n[done] ", STEM, " (B) + ", STEM_C, " (C) + ", STEM_D,
    " (D), each .{pdf,png}, + 4 TSVs -> ", OUTDIR, "/\n", sep = "")
