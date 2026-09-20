#!/usr/bin/env Rscript
# =============================================================================
# fig5_J_typeaccess.R  -  Fig 5 panel J: for each arm (2 genomes x PreClean / WD / ND) the type-level
#   reciprocal z-score (rZ) of the consensus meta-cells beside the pooled marker ACCESSIBILITY
#   aggregated to the same (cell type x meta-cell) grid, aligned row for row and column for column.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/rZ_annotation/consensus/
#           <prefix>.permetacell_{type,marker}_rZ.geom.tsv
#         figures/main/fig5/Fig5_P3D_marker_access_{pooled,genes}.tsv   (cache from fig5_J_access_cache.R)
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P3BD_rZ_access.{pdf,png} + _typeaccess.tsv + _agreement.tsv
# Run     Rscript analysis/fig5_biological_impact/fig5_J_typeaccess.R   (repo root, ~25 s; run fig5_J_access_cache.R first)
# Part 3 internal lettering kept in titles and cross-references: "B | D" = this pair = manuscript J.
# =============================================================================
#
# THE AGGREGATION (the one real choice in this figure):
#   value(type, mc) = MEAN over the type's marker genes of the per-gene min-max-scaled (0-1)
#   pooled accessibility. Each gene is scaled BEFORE averaging because within a type the
#   markers differ in absolute per-kb level by orders of magnitude -- a mean of raw values
#   would be the profile of the type's single most accessible gene. Scaling first weights
#   every marker equally (the standard dot-plot / module-score convention). The mean of the
#   raw pooled values is exported alongside so nothing is hidden.
#
# DISPLAY RESCALE: those means never reach 1 (a type's markers do not all peak in the same
#   meta-cell), so for COLOUR each TYPE ROW is additionally min-max rescaled 0-1 across its
#   three stage panels. The colourbar is then a true 0-1 and every row uses the full range.
#   The stamps and TSV stats stay on the PRE-rescale means -- after a per-type rescale a
#   cross-type argmax is meaningless by construction (every type hits 1 somewhere).
#
# Pipeline for the accessibility half (identical to the gene-level panel):
#   perkb genes x cells -> pooled (sum) into the CONSENSUS meta-cells -> each meta-cell scaled
#   to a common total over ALL genes (depth out) -> per gene min-max 0-1 across the genome's
#   THREE stages together -> mean over each scored type's shared markers.
# It READS THE CACHE written by fig5_J_access_cache.R (no perkb load here) and reconciles the
# per-gene scale ranges against that script's genes TSV.
#
# TRAPS (the panel-specific ones first):
#   COLOUR IS RELATIVE THREE TIMES OVER on the access half: per-gene min-max, mean over
#     markers, then the per-type-row 0-1 display rescale. It shows WHERE a type's markers are
#     most active, never HOW accessible -- and cross-ROW colour comparison is meaningless by
#     construction. Raw pooled and pre-rescale means are in the TSV.
#   The orange marks are the rZ CALL drawn on BOTH halves -- on the access half they are NOT
#     the accessibility argmax. That is the point: the reader checks whether the called row is
#     also the accessible row. The agreement rate is stamped per panel and exported.
#   ARGMAX CALLS, still. The expression-matched gene null finds no significant meta-cell x
#     type pair at q < 0.05 in five of the six objects (0 of 252 in every At object, 0 of 5,720
#     in maize PreClean and WD) and 2 of 5,720 in maize ND, so the labels are descriptive
#     annotation weighted by support, never tested identity.
#   Scored types only: the y axis must match panel B, so markers of unscored types (no row in
#     B) are NOT here -- they stay in the gene-level panel.
#   Gene sets intersected across a genome's stages (At/nd lacks AT2G08610).
#   Partitions differ per stage -- a column is not the same meta-cell across panels. wd and nd
#     are independent treatments, not a chain. rZ is not comparable across genomes; every
#     scale here is per genome row.
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

STEM    <- "Fig5_P3BD_rZ_access"
MCD     <- file.path(PLATE, "step5_metacell")
RZDIR   <- file.path(MCD, "rZ_annotation", "consensus")
METRIC  <- "geom"
RZ_CLIP <- 0.99
CACHE   <- file.path(OUTDIR, "Fig5_P3D_marker_access_pooled.tsv")
GENES_T <- file.path(OUTDIR, "Fig5_P3D_marker_access_genes.tsv")

G_LEV <- vapply(CFG, `[[`, character(1), "label")
S_LEV <- names(STAGES)
FIXED <- c("metacell", "seed", "cluster", "n_cells", "top_type", "top_rZ", "margin")
TOL   <- 1e-4

RZ_PAL  <- viridisLite::viridis(256, option = "G", direction = -1)  # as panel B
AC_PAL  <- colorRampPalette(c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6",
                              "#08519C", "#08214A"))(256)           # as panel D
CALLCOL <- "#FF7F00"
SEPCOL  <- "grey85"

prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)

if (!file.exists(CACHE))
  stop("pooled-accessibility cache missing -- run fig5_J_access_cache.R first: ", CACHE)

# =============================================================================
# 1. type-level rZ (panel B's data + ordering rule) and the marker id map
# =============================================================================
calls <- list(); rzl <- list(); types_by_g <- list(); mk_ids <- list()

for (g in names(CFG)) for (s in S_LEV) {
  pf <- prefix_for(g, s)
  f  <- file.path(RZDIR, sprintf("%s.permetacell_type_rZ.%s.tsv", pf, METRIC))
  if (!file.exists(f)) stop("missing type rZ table: ", f)
  rz <- fread(f)
  tcol <- setdiff(names(rz), FIXED)
  if (length(tcol) < 2L) stop(pf, ": no type columns")
  if (is.null(types_by_g[[g]])) types_by_g[[g]] <- sort(tcol)
  else if (!identical(types_by_g[[g]], sort(tcol)))
    stop(g, ": scored type set differs between stages")
  M <- as.matrix(rz[, ..tcol])
  am <- tcol[max.col(M, ties.method = "first")]
  if (any(am != rz$top_type & abs(M[cbind(seq_len(nrow(M)), match(rz$top_type, tcol))] -
                                  apply(M, 1, max)) > TOL))
    stop(pf, ": top_type is not the argmax of the type columns")

  tl <- types_by_g[[g]]
  d  <- rz[, .(metacell, n_cells, top_type, top_rZ)]
  d[, ypos_t := match(top_type, tl)]
  if (anyNA(d$ypos_t)) stop(pf, ": a call is outside the scored type set")
  setorder(d, ypos_t, -top_rZ)                      # panel B's ordering rule
  d[, `:=`(xpos = seq_len(.N), genome = g, stage = s)]
  calls[[pf]] <- d

  L <- melt(rz[, c("metacell", ..tcol)], id.vars = "metacell",
            variable.name = "type", value.name = "rZ", variable.factor = FALSE)
  rzl[[pf]] <- merge(L, d[, .(metacell, xpos)], by = "metacell")[
    , `:=`(genome = g, stage = s)]

  fm <- file.path(RZDIR, sprintf("%s.permetacell_marker_rZ.%s.tsv", pf, METRIC))
  if (!file.exists(fm)) stop("missing marker table: ", fm)
  mk_ids[[pf]] <- fread(fm, select = c("geneID", "name", "type_label"))[
    , `:=`(genome = g, stage = s)]
}
calls_all <- rbindlist(calls)
rzl_all   <- rbindlist(rzl)

# marker genes shared by a genome's three stages, with their type label
gmap <- list()
for (g in names(CFG)) {
  sets   <- lapply(S_LEV, function(s) mk_ids[[prefix_for(g, s)]]$geneID)
  shared <- Reduce(intersect, sets)
  mk <- unique(rbindlist(mk_ids)[genome == g & geneID %in% shared,
                                 .(geneID, name, type_label)])
  if (anyDuplicated(mk$geneID) > 0)
    stop(g, ": a shared marker carries two type labels")
  gmap[[g]] <- mk[, genome := g]
}
gmap <- rbindlist(gmap)

# =============================================================================
# 2. pooled accessibility from the cache -> per-gene 0-1 -> mean per type
# =============================================================================
pooled <- fread(CACHE)
if (pooled[!is.finite(access), .N] > 0) stop("cache has non-finite access values")
for (g in names(CFG)) {
  if (!setequal(pooled[genome == g, unique(geneID)], gmap[genome == g, geneID]))
    stop("cache gene set does not match the shared marker set for ", g,
         " -- rerun fig5_J_access_cache.R")
  for (s in S_LEV)
    if (!setequal(pooled[genome == g & stage == s, unique(metacell)],
                  calls[[prefix_for(g, s)]]$metacell))
      stop("cache meta-cells do not match the consensus for ", prefix_for(g, s),
           " -- rerun fig5_J_access_cache.R")
}

rngs <- pooled[, .(lo = min(access), hi = max(access)), by = .(genome, geneID)]
rngs[, flat := hi <= lo]
if (file.exists(GENES_T)) {          # reconcile vs the gene-level panel's export
  gt <- fread(GENES_T)[, .(genome, geneID, scale_lo, scale_hi)]
  ck <- merge(rngs, gt, by = c("genome", "geneID"))
  if (nrow(ck) != nrow(rngs)) stop("genes TSV does not cover the cached gene set")
  dmax <- ck[, max(abs(lo - scale_lo), abs(hi - scale_hi))]
  if (dmax > TOL) stop(sprintf("per-gene scale ranges disagree with %s (max|diff| %.2e)",
                               basename(GENES_T), dmax))
  message(sprintf("  [reconcile] per-gene scale ranges vs %s: max|diff| = %.1e",
                  basename(GENES_T), dmax))
} else message("  [reconcile] ", basename(GENES_T), " absent -- skipping the range cross-check")

pooled <- merge(pooled, rngs, by = c("genome", "geneID"))
pooled[, scaled := fifelse(flat, 0, (access - lo) / (hi - lo))]

# aggregate to (type x meta-cell); SCORED types only (the y axis must match B)
tacc <- merge(pooled, gmap[, .(genome, geneID, type_label)],
              by = c("genome", "geneID"))
tacc <- tacc[mapply(function(gg, tt) tt %in% types_by_g[[gg]], genome, type_label)]
tacc <- tacc[, .(mean_scaled = mean(scaled), mean_pooled = mean(access),
                 n_markers = .N),
             by = .(genome, stage, metacell, type = type_label)]

nmk <- unique(tacc[, .(genome, type, n_markers)])
for (g in names(CFG)) {
  missing_t <- setdiff(types_by_g[[g]], nmk[genome == g, type])
  if (length(missing_t) > 0)
    stop(g, ": scored type(s) with no shared marker gene: ",
         paste(missing_t, collapse = ", "))
  for (s in S_LEV) {
    n_exp <- length(types_by_g[[g]]) * nrow(calls[[prefix_for(g, s)]])
    if (tacc[genome == g & stage == s, .N] != n_exp)
      stop(prefix_for(g, s), ": type-access grid is not complete")
  }
}

tacc <- merge(tacc, rzl_all, by = c("genome", "stage", "metacell", "type"))
tacc <- merge(tacc, calls_all[, .(genome, stage, metacell, top_type)],
              by = c("genome", "stage", "metacell"))
tacc[, is_call := type == top_type]
if (tacc[is.na(rZ) | is.na(xpos), .N] > 0) stop("type-access join lost rZ or xpos")

# display rescale: each TYPE ROW 0-1 across the genome's three stages together
# (the same joint-across-stages rule as the per-gene scaling, one level up)
trng <- tacc[, .(t_lo = min(mean_scaled), t_hi = max(mean_scaled)),
             by = .(genome, type)]
trng[, t_flat := t_hi <= t_lo]
if (trng[, sum(t_flat)] > 0)
  message(sprintf("  [scale] %d type row(s) flat across all stages -> set to 0",
                  trng[, sum(t_flat)]))
tacc <- merge(tacc, trng, by = c("genome", "type"))
tacc[, scaled01 := fifelse(t_flat, 0, (mean_scaled - t_lo) / (t_hi - t_lo))]

# =============================================================================
# 3. agreement between the rZ call and the accessibility argmax  (descriptive)
# =============================================================================
acc_top <- tacc[, .SD[which.max(mean_scaled)],
                by = .(genome, stage, metacell)][, .(genome, stage, metacell,
                                                     acc_top = type)]
agree <- merge(unique(tacc[, .(genome, stage, metacell, top_type)]), acc_top,
               by = c("genome", "stage", "metacell"))
stats <- tacc[, .(spearman = cor(rZ, mean_scaled, method = "spearman"),
                  n_pairs = .N), by = .(genome, stage)]
stats <- merge(stats,
               agree[, .(argmax_agree = mean(acc_top == top_type), n_mc = .N),
                     by = .(genome, stage)], by = c("genome", "stage"))

# =============================================================================
# 4. PANELS -- for each arm a pair: rZ (mako) | mean marker access (blue)
# =============================================================================
rz_rng <- rzl_all[, .(lo = min(rZ), hi = quantile(rZ, RZ_CLIP), true_hi = max(rZ)),
                  by = genome]
# access half needs no range table: the per-type display rescale makes its
# scale a fixed 0-1 by construction (no clip, no squish)

mk_half <- function(g, s, kind, ylab_on) {
  pf <- prefix_for(g, s)
  d  <- calls[[pf]]; tl <- types_by_g[[g]]; n <- nrow(d)
  L  <- if (kind == "rz") copy(rzl[[pf]])[, val := rZ]
        else copy(tacc[genome == g & stage == s])[, val := scaled01]
  L[, ty := match(type, tl)]
  if (anyNA(L$ty)) stop(pf, ": a plotted type is outside the panel")
  rng <- if (kind == "rz") unlist(rz_rng[genome == g, .(lo, hi)]) else c(0, 1)
  pal <- if (kind == "rz") RZ_PAL else AC_PAL

  grp <- d[, .(x0 = min(xpos), x1 = max(xpos), ypos = match(top_type[1], tl)),
           by = top_type]
  if (grp[, sum(x1 - x0 + 1L)] != n) stop(pf, ": call groups not contiguous")
  seps <- grp[order(x1), x1][-nrow(grp)] + 0.5
  cnt  <- nmk[genome == g][match(tl, type), n_markers]
  ylab <- sprintf("%s (%d)", sub("^[A-Za-z]+:", "", tl), cnt)
  st   <- stats[genome == g & stage == s]

  ggplot(L, aes(x = xpos, y = ty, fill = val)) +
    geom_raster() +
    geom_vline(xintercept = seps, colour = SEPCOL, linewidth = 0.15) +
    geom_segment(data = grp, inherit.aes = FALSE,
                 aes(x = x0 - 0.5, xend = x1 + 0.5, y = ypos, yend = ypos),
                 colour = CALLCOL, linewidth = 0.5) +
    scale_fill_gradientn(colours = pal, limits = rng, guide = "none",
                         oob = scales::squish) +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, length(tl) + 0.5), expand = c(0, 0),
                       breaks = seq_along(tl), labels = ylab) +
    labs(x = sprintf("%d meta-cells", n), y = NULL,
         title = sprintf("%s | %s", s, if (kind == "rz") "rZ" else "marker access"),
         subtitle = if (kind == "rz") sprintf("the call (orange) = rZ argmax")
                    else sprintf("call = access argmax: %.0f%% | rho %.2f",
                                 100 * st$argmax_agree, st$spearman)) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size = 5.0, colour = "grey30"),
          axis.text.x = element_blank(),
          axis.text.y = if (ylab_on) element_text(size = 4.6, colour = "grey25")
                        else element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 7.2, hjust = 0.5),
          plot.subtitle = element_text(size = 4.9, colour = "grey35", hjust = 0.5),
          plot.margin = if (kind == "rz") margin(2, 1, 2, 2) else margin(2, 9, 2, 1))
}

mk_cbar <- function(g, kind) {
  pal <- if (kind == "rz") RZ_PAL else AC_PAL
  gnm <- if (g == "TAIR10") "Arabidopsis" else "maize"
  if (kind == "rz") {
    r   <- rz_rng[genome == g]
    l4  <- sprintf("%.2f", seq(r$lo, r$hi, length.out = 4)); l4[4] <- paste0(">=", l4[4])
    ttl <- sprintf("rZ (%s)\n%s\nclip p%g\nmax %.2f", METRIC, gnm, 100 * RZ_CLIP, r$true_hi)
  } else {
    l4  <- c("0", "0.33", "0.67", "1")
    ttl <- sprintf("marker access\n%s\nper-type 0-1\n(display rescale)", gnm)
  }
  ggplot() +
    annotation_raster(as.raster(matrix(rev(pal), ncol = 1)),
                      xmin = 0, xmax = 1, ymin = 0, ymax = 1, interpolate = TRUE) +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = NA,
             colour = "grey25", linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1.35, y = c(0, 1/3, 2/3, 1),
             yend = c(0, 1/3, 2/3, 1), linewidth = 0.3, colour = "grey25") +
    annotate("text", x = 1.6, y = c(0, 1/3, 2/3, 1), label = l4, hjust = 0,
             size = 2.0, colour = "grey20") +
    coord_cartesian(xlim = c(-0.2, 5.6), ylim = c(-0.03, 1.03), expand = FALSE) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 5.2, colour = "grey20", hjust = 0),
          plot.margin = margin(14, 1, 10, 1)) +
    labs(title = ttl)
}

row_of <- function(g) {
  ps <- list()
  for (s in S_LEV) {
    ps[[paste0(s, "_rz")]] <- mk_half(g, s, "rz",
                                      ylab_on = (s == S_LEV[1]))
    ps[[paste0(s, "_ac")]] <- mk_half(g, s, "access", ylab_on = FALSE)
  }
  # the two colourbars STACK in one narrow column -- side by side their titles
  # overflow into each other at this panel width
  (ps[[1]] | ps[[2]] | ps[[3]] | ps[[4]] | ps[[5]] | ps[[6]] |
     (mk_cbar(g, "rz") / mk_cbar(g, "access"))) +
    plot_layout(widths = c(1, 1, 1, 1, 1, 1, 0.16))
}

comp <- row_of("B73v5") / row_of("TAIR10") +
  plot_annotation(
    title = "Panel B beside its evidence - type-level rZ | pooled marker accessibility, per consensus meta-cell",
    subtitle = paste0(
      "Each stage is a PAIR on the same (cell type x meta-cell) grid: LEFT the reciprocal z-score panel B plots; RIGHT the pooled per-kb accessibility of the same\n",
      "markers, aggregated to the type level. Rows and columns align exactly (x = panel B's ordering; y = the scored types, with the marker count per type in the\n",
      "label). ORANGE = the rZ call, drawn on BOTH halves - on the access half it is NOT the argmax, so the reader can check whether the called row is also the\n",
      "accessible row; the per-panel stamp gives that agreement rate and the Spearman correlation between the two matrices.\n",
      "Access half: cells pooled (sum) into their consensus meta-cell; each meta-cell scaled to a common total over ALL genes (depth out); each gene min-max\n",
      "scaled 0-1 across the genome's three stages together; the MEAN over the type's shared markers (genes scaled before averaging - raw per-kb levels differ\n",
      "by orders of magnitude within a type, so a raw mean would be the profile of the single most accessible marker); finally each TYPE ROW is min-max rescaled\n",
      "0-1 across its three stage panels, so the colourbar is a true 0-1 and every row uses the full range. Colour is therefore RELATIVE THREE TIMES OVER: it says\n",
      "WHERE a type's markers are most active, never how accessible, and cross-ROW colour comparison is meaningless by construction (every type hits 1 somewhere).\n",
      "The stamps (argmax agreement, Spearman) are computed on the PRE-rescale means - a cross-type argmax after per-type rescaling would be meaningless for the\n",
      "same reason. Raw pooled and pre-rescale means are in the TSV. The rZ half is clipped at p99 per genome row.\n",
      "ARGMAX calls throughout (the expression-matched gene null finds no significant meta-cell x type pair in five of the six objects and 2 of 5,720 in maize ND):\n",
      "descriptive annotation, weighted by seed support. Markers of unscored types have no row in panel B and are shown only in the gene-level panel. Partitions\n",
      "differ per stage - a column is never the same meta-cell across panels. wd and nd are independent treatments of the same raw input, not a chain. rZ is not\n",
      "comparable across genomes; every scale here is per genome row."),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 6.3, colour = "grey35",
                                               lineheight = 1.28)))

# =============================================================================
# REPORT
# =============================================================================
cat("=== FIG 5 PANEL J (Part 3 B | D) | type-level rZ vs pooled marker accessibility ===\n\n")
cat("access cache : ", CACHE, "\n", sep = "")
cat("aggregation  : mean of per-gene min-max-scaled pooled accessibility over each\n",
    "               scored type's shared markers (equal gene weighting)\n",
    "display      : each type row rescaled 0-1 across its three stages for colour;\n",
    "               all stats below are on the PRE-rescale means\n\n", sep = "")
# stats already carries n_mc (from the agreement merge); assert it matches the
# call tables rather than merging a duplicate column in.
chk_n <- merge(stats, calls_all[, .(n_calls = .N), by = .(genome, stage)],
               by = c("genome", "stage"))
if (chk_n[n_mc != n_calls, .N] > 0) stop("agreement n_mc disagrees with the call tables")
rp <- stats[order(genome, factor(stage, levels = S_LEV))]
print(rp[, .(genome, stage, n_mc, types = n_pairs / n_mc,
             spearman = round(spearman, 3),
             pct_call_is_access_argmax = round(100 * argmax_agree, 1))])
cat("\n  DESCRIPTIVE: no test anywhere (meta-cells share one partition and one seed\n",
    "  set). The agreement rate asks: in what fraction of meta-cells is the rZ call\n",
    "  also the most-accessible type? rZ standardises each gene across meta-cells\n",
    "  before combining, the access half min-max scales and averages - they need not\n",
    "  agree, and where they diverge the divergence is the information.\n", sep = "")

# =============================================================================
# EXPORT + VERIFY
# =============================================================================
written <- character()
for (ext in c("pdf", "png")) {
  f <- file.path(OUTDIR, sprintf("%s.%s", STEM, ext))
  # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
  if (ext == "pdf") ggsave(f, comp, width = 16.0, height = 9.6, bg = "white")
  else              ggsave(f, comp, width = 16.0, height = 9.6, dpi = 300, bg = "white")
  written <- c(written, f)
}
f <- file.path(OUTDIR, paste0(STEM, "_typeaccess.tsv"))
fwrite(tacc[order(genome, factor(stage, levels = S_LEV), xpos, type),
            .(genome, stage, metacell, xpos, type, n_markers,
              mean_scaled = round(mean_scaled, 5), mean_pooled = round(mean_pooled, 5),
              display_scaled01 = round(scaled01, 5), rZ, is_call)], f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, paste0(STEM, "_agreement.tsv"))
fwrite(rp[, .(genome, stage, n_mc, spearman = round(spearman, 4),
              argmax_agree = round(argmax_agree, 4))], f, sep = "\t")
written <- c(written, f)

cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz   <- if (file.exists(f)) file.size(f) else NA_integer_
  good <- !is.na(sz) && sz > (if (grepl("\\.tsv$", f)) 50 else 1000)
  ok   <- ok && good
  cat(sprintf("  %-4s %-44s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
cat("\n[done] ", STEM, ".{pdf,png} + 2 TSVs -> ", OUTDIR, "/\n", sep = "")
