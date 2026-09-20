#!/usr/bin/env Rscript
# =============================================================================
# fig5_J_access_cache.R  -  builds the pooled per-kb marker-accessibility CACHE read by Fig 5
#   panels J (fig5_J_typeaccess.R) and K (fig5_K_examples.R), and the gene-level pooled
#   accessibility heatmap of the consensus meta-cells (Fig5_P3D_marker_access, not a manuscript panel).
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step3_compare/<prefix>.plate.perkb.genes.sparse.rds
#         .../SM2v2_plate/step5_metacell/consensus/<prefix>.consensus.cells.tsv
#         .../SM2v2_plate/step5_metacell/rZ_annotation/consensus/<prefix>.permetacell_{type,marker}_rZ.geom.tsv
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P3D_marker_access.{pdf,png}, _pooled.tsv (THE CACHE), _genes.tsv, _diag.tsv
# Run     Rscript analysis/fig5_biological_impact/fig5_J_access_cache.R   (repo root; ~4 min cold, ~20 s cached)
#         Run BEFORE fig5_J_typeaccess.R and fig5_K_examples.R. ACCESS_RECOMPUTE=1 forces a fresh pool.
# =============================================================================
#
# Part 3 internal lettering kept in titles and cross-references: this is the gene-level "D"
# panel behind the type-level "B | D" pair (manuscript J); the rZ heatmap is "B".
#
# WHAT IS PLOTTED, exactly:
#   1. `step3_compare/<prefix>.plate.perkb.genes.sparse.rds` -- a genes x CELLS dgCMatrix of
#      per-kb gene accessibility (accessibility normalised by gene length). This is the same
#      RAW matrix the upstream 4_3p annotation consumes; it is NOT the 4_3f-smoothed one,
#      because pooling cells into a meta-cell already smooths.
#   2. Cells are POOLED into their consensus meta-cell (sum), using
#      `consensus/<prefix>.consensus.cells.tsv` -- the partition built from the 100 SEACells
#      seeds via the co-assignment matrix F.
#   3. Each meta-cell is then SCALED TO A COMMON TOTAL over ALL genes (not over the marker
#      subset), mirroring 4_3p -- this removes the depth/size difference between meta-cells,
#      which would otherwise dominate everything.
#   4. Only then is the matrix subset to the marker panel and each gene min-max scaled to 0-1.
#
# THE 0-1 SCALING IS PER GENE AND IT MANUFACTURES CONTRAST. Every gene is forced to span the
#   full colour range, including genes with no real variation across meta-cells -- a flat gene
#   and a sharply specific gene look equally colourful. That is the standard marker-heatmap
#   convention and it is what makes the block structure visible at all, but it means COLOUR
#   INTENSITY IS NOT ACCESSIBILITY. The unscaled pooled values are kept in the cache TSV; the
#   per-gene range actually used is exported alongside them.
# The min-max range is taken over the genome's THREE STAGES TOGETHER, not per panel, so the
#   three panels of a row are on one colour scale and comparable. Scaling each panel
#   separately would make every stage look identical.
# Meta-cell partitions differ per stage, so a column position means different cells in each
#   panel. x reuses the panel-B ordering rule (grouped by call, then descending top_rZ) so B
#   and D align, but never read a column across stages as "the same meta-cell".
# Gene sets are intersected across a genome's three stages (At/nd is missing AT2G08610,
#   At:phloem) -- a per-stage gene set on a shared y axis would not be comparable.
# The marker panel is WIDER than the scored type set (maize 221 markers over 39 labels, only
#   20 types scored; At 275 over 20 labels, 18 scored). Genes of unscored types are kept in a
#   trailing block and can never form a diagonal box, because no meta-cell can be called
#   their type.
# COST: the six perkb matrices are 11-69 MB each (~230 MB total). The pooled marker matrix is
#   cached, so re-runs are instant. Force a fresh compute with ACCESS_RECOMPUTE=1.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2); library(patchwork)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source("analysis/_helpers/fig5_part2_helpers.R")   # SOC, PLATE, CFG, STAGES

STEM_D <- "Fig5_P3D_marker_access"
MCD    <- file.path(PLATE, "step5_metacell")
RZDIR  <- file.path(MCD, "rZ_annotation", "consensus")
CONS   <- file.path(MCD, "consensus")
PERKB  <- file.path(PLATE, "step3_compare")
METRIC <- "geom"
RECOMP <- nzchar(Sys.getenv("ACCESS_RECOMPUTE"))

G_LEV <- vapply(CFG, `[[`, character(1), "label")
S_LEV <- names(STAGES)
FIXED <- c("metacell", "seed", "cluster", "n_cells", "top_type", "top_rZ", "margin")

# accessibility ramp: deliberately DIFFERENT from the rZ panels' mako, so the
# two are never confused. White -> blue -> near-black.
AC_PAL   <- colorRampPalette(c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6",
                               "#08519C", "#08214A"))(256)
BOXCOL   <- "#D94801"   # call-vs-own-markers box: orange on a blue ramp
CACHE    <- file.path(OUTDIR, paste0(STEM_D, "_pooled.tsv"))

prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)

# =============================================================================
# 1. the calls + x ordering (same rule as panel B) and the marker panel
# =============================================================================
read_calls <- function(g, s) {
  f <- file.path(RZDIR, sprintf("%s.permetacell_type_rZ.%s.tsv", prefix_for(g, s), METRIC))
  if (!file.exists(f)) stop("missing type rZ table: ", f)
  rz <- fread(f)
  tcol <- setdiff(names(rz), FIXED)
  list(d = rz[, .(metacell, n_cells, top_type, top_rZ)], types = sort(tcol))
}
read_markers <- function(g, s) {
  f <- file.path(RZDIR, sprintf("%s.permetacell_marker_rZ.%s.tsv", prefix_for(g, s), METRIC))
  if (!file.exists(f)) stop("missing marker table: ", f)
  fread(f, select = c("geneID", "name", "type_label"))
}

calls <- list(); types_by_g <- list(); markers <- list()
for (g in names(CFG)) {
  for (s in S_LEV) {
    cc <- read_calls(g, s)
    if (is.null(types_by_g[[g]])) types_by_g[[g]] <- cc$types
    else if (!identical(types_by_g[[g]], cc$types))
      stop(g, ": the scored type set differs between stages")
    tl <- types_by_g[[g]]
    d  <- cc$d
    d[, ypos_t := match(top_type, tl)]
    if (anyNA(d$ypos_t)) stop(prefix_for(g, s), ": a call is not in the scored type set")
    setorder(d, ypos_t, -top_rZ)              # <- panel B's ordering rule
    d[, `:=`(xpos = seq_len(.N), genome = g, stage = s)]
    calls[[prefix_for(g, s)]] <- d
    markers[[prefix_for(g, s)]] <- read_markers(g, s)[, `:=`(genome = g, stage = s)]
  }
}
calls_all <- rbindlist(calls)

# gene set shared by a genome's three stages
gene_keep <- list(); lost_note <- character()
for (g in names(CFG)) {
  sets   <- lapply(S_LEV, function(s) markers[[prefix_for(g, s)]]$geneID)
  shared <- Reduce(intersect, sets)
  lost   <- setdiff(Reduce(union, sets), shared)
  gene_keep[[g]] <- shared
  if (length(lost) > 0) {
    mk  <- rbindlist(markers)[genome == g & geneID %in% lost]
    lab <- unique(mk[, paste0(geneID, " (", type_label, ")")])
    lost_note <- c(lost_note, sprintf("%s is missing %s", G_LEV[[g]],
                                      paste(lab, collapse = ", ")))
    message(sprintf("  [genes] %s: %d absent from >=1 stage, dropped: %s",
                    g, length(lost), paste(lab, collapse = ", ")))
  }
}

# =============================================================================
# 2. POOL the per-kb accessibility into consensus meta-cells   (cached)
# =============================================================================
pool_arm <- function(g, s) {
  pf <- prefix_for(g, s)
  frds <- file.path(PERKB, sprintf("%s.plate.perkb.genes.sparse.rds", pf))
  fcel <- file.path(CONS,  sprintf("%s.consensus.cells.tsv", pf))
  for (f in c(frds, fcel)) if (!file.exists(f)) stop("missing input: ", f)

  cl <- fread(fcel)[, .(cellID, consensus_mc)]
  t0 <- proc.time()[3]
  X  <- readRDS(frds)                     # genes x cells, per-kb accessibility
  if (!inherits(X, "dgCMatrix")) stop(basename(frds), ": not a dgCMatrix")
  if (is.null(rownames(X)) || is.null(colnames(X)))
    stop(basename(frds), ": matrix has no dimnames")

  miss <- setdiff(cl$cellID, colnames(X))
  if (length(miss) > 0)
    stop(pf, ": ", length(miss), " consensus cells absent from the perkb matrix")
  X <- X[, cl$cellID, drop = FALSE]       # order columns to match the map

  # cells -> meta-cells: sum. Group order fixed to cMC-<n> by integer id.
  mc  <- factor(paste0("cMC-", cl$consensus_mc),
                levels = paste0("cMC-", sort(unique(cl$consensus_mc))))
  ind <- sparse.model.matrix(~ 0 + mc)
  colnames(ind) <- levels(mc)
  P <- X %*% ind                          # genes x meta-cells, pooled

  # scale each meta-cell to a COMMON TOTAL over ALL genes (as 4_3p does),
  # BEFORE subsetting to markers -- otherwise the "total" is marker-only and
  # the normalisation silently changes meaning.
  tot <- Matrix::colSums(P)
  if (any(tot <= 0)) stop(pf, ": a meta-cell pooled to zero total accessibility")
  P <- P %*% Diagonal(x = mean(tot) / tot)
  dimnames(P) <- list(rownames(X), levels(mc))

  keep <- intersect(gene_keep[[g]], rownames(P))
  if (length(keep) != length(gene_keep[[g]]))
    stop(pf, ": ", length(gene_keep[[g]]) - length(keep),
         " marker gene(s) absent from the perkb matrix")
  P <- as.matrix(P[keep, , drop = FALSE])
  rm(X); gc(verbose = FALSE)

  message(sprintf("  [pool] %-28s %d genes x %d meta-cells  (%.1f s)",
                  pf, nrow(P), ncol(P), proc.time()[3] - t0))
  data.table(genome = g, stage = s, geneID = rep(rownames(P), ncol(P)),
             metacell = rep(colnames(P), each = nrow(P)), access = as.vector(P))
}

if (!RECOMP && file.exists(CACHE)) {
  pooled <- fread(CACHE)
  message("  [cache] reusing ", CACHE, " (ACCESS_RECOMPUTE=1 to force a recompute)")
  # the cache must match the CURRENT consensus + gene sets, or it is stale
  for (g in names(CFG)) {
    if (!setequal(pooled[genome == g, unique(geneID)], gene_keep[[g]]))
      stop("cached pooled table does not match the current gene set for ", g,
           " -- rerun with ACCESS_RECOMPUTE=1")
    for (s in S_LEV) {
      if (!setequal(pooled[genome == g & stage == s, unique(metacell)],
                    calls[[prefix_for(g, s)]]$metacell))
        stop("cached pooled table does not match the current meta-cells for ",
             prefix_for(g, s), " -- rerun with ACCESS_RECOMPUTE=1")
    }
  }
} else {
  message("  computing pooled accessibility from the perkb matrices (~230 MB to read)")
  pooled <- rbindlist(lapply(names(CFG), function(g)
    rbindlist(lapply(S_LEV, function(s) pool_arm(g, s)))))
  fwrite(pooled, CACHE, sep = "\t")
}

# =============================================================================
# 3. scale each gene to 0-1 over the genome's THREE stages together
# =============================================================================
rngs <- pooled[, .(lo = min(access), hi = max(access)), by = .(genome, geneID)]
rngs[, flat := hi <= lo]
if (rngs[, sum(flat)] > 0)
  message(sprintf("  [scale] %d gene(s) are flat across meta-cells -> set to 0",
                  rngs[, sum(flat)]))
pooled <- merge(pooled, rngs, by = c("genome", "geneID"))
pooled[, scaled := fifelse(flat, 0, (access - lo) / (hi - lo))]

# gene layout: type group (panel B's type order, unscored last), then by the
# gene's overall pooled level so the strongest markers head each group
gi_all <- unique(rbindlist(markers)[, .(genome, geneID, name, type_label)])
gi_all <- gi_all[mapply(function(gg, id) id %in% gene_keep[[gg]], genome, geneID)]
if (anyDuplicated(gi_all[, .(genome, geneID)]) > 0)
  stop("a gene carries more than one type_label within a genome")
str_ <- pooled[, .(strength = mean(access)), by = .(genome, geneID)]
gi_all <- merge(gi_all, str_, by = c("genome", "geneID"))

gord <- list(); ggrp <- list()
for (g in names(CFG)) {
  tl <- types_by_g[[g]]
  gi <- gi_all[genome == g]
  gi[, scored := type_label %in% tl]
  gi[, tord := fifelse(scored, match(type_label, tl), length(tl) + 1L)]
  setorder(gi, tord, type_label, -strength)
  gi[, ypos := seq_len(.N)]
  gord[[g]] <- gi
  ggrp[[g]] <- gi[, .(y0 = min(ypos), y1 = max(ypos), n = .N, scored = scored[1]),
                  by = .(grp = fifelse(scored, type_label,
                                       "(types with too few markers to score)"))]
}

# =============================================================================
# 4. PANELS
# =============================================================================
mk_panel <- function(g, s, ylab_on) {
  pf <- prefix_for(g, s)
  gi <- gord[[g]]; grp <- ggrp[[g]]; d <- calls[[pf]]
  L <- merge(pooled[genome == g & stage == s], gi[, .(geneID, ypos)], by = "geneID")
  L <- merge(L, d[, .(metacell, xpos)], by = "metacell")
  if (nrow(L) != nrow(gi) * nrow(d)) stop(pf, ": pooled matrix is not complete")
  n <- nrow(d); ng <- nrow(gi)

  cg    <- d[, .(x0 = min(xpos), x1 = max(xpos)), by = top_type]
  diagb <- merge(cg, grp[scored == TRUE], by.x = "top_type", by.y = "grp")

  ggplot(L, aes(x = xpos, y = ypos, fill = scaled)) +
    geom_raster() +
    geom_hline(yintercept = grp$y1[-nrow(grp)] + 0.5, colour = "grey75", linewidth = 0.12) +
    geom_vline(xintercept = cg[order(x1), x1][-nrow(cg)] + 0.5,
               colour = "grey75", linewidth = 0.12) +
    geom_rect(data = diagb, inherit.aes = FALSE,
              aes(xmin = x0 - 0.5, xmax = x1 + 0.5, ymin = y0 - 0.5, ymax = y1 + 0.5),
              fill = NA, colour = BOXCOL, linewidth = 0.32) +
    scale_fill_gradientn(colours = AC_PAL, limits = c(0, 1), guide = "none") +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, ng + 0.5), expand = c(0, 0),
                       breaks = grp[, (y0 + y1) / 2],
                       labels = sub("^[A-Za-z]+:", "", grp$grp)) +
    labs(x = sprintf("%s meta-cells, grouped by call (panel B order)", n),
         y = NULL, title = sprintf("%s - %s", G_LEV[[g]], s),
         subtitle = sprintf("%d markers shared by all 3 stages | %d in a scored type",
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

cbar <- ggplot() +
  annotation_raster(as.raster(matrix(rev(AC_PAL), ncol = 1)),
                    xmin = 0, xmax = 1, ymin = 0, ymax = 1, interpolate = TRUE) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = NA,
           colour = "grey25", linewidth = 0.3) +
  annotate("segment", x = 1, xend = 1.35, y = c(0, .25, .5, .75, 1),
           yend = c(0, .25, .5, .75, 1), linewidth = 0.3, colour = "grey25") +
  annotate("text", x = 1.6, y = c(0, .25, .5, .75, 1),
           label = c("0", "0.25", "0.50", "0.75", "1"), hjust = 0,
           size = 2.1, colour = "grey20") +
  coord_cartesian(xlim = c(-0.2, 5.2), ylim = c(-0.03, 1.03), expand = FALSE) +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 5.6, colour = "grey20", hjust = 0),
        plot.margin = margin(30, 2, 34, 2)) +
  labs(title = "pooled per-kb\naccessibility,\nmin-max scaled\nper gene")

P <- list()
for (g in names(CFG)) for (s in S_LEV)
  P[[prefix_for(g, s)]] <- mk_panel(g, s, ylab_on = (s == "PreClean"))

# =============================================================================
# REPORT -- including the on-diagonal contrast, so the picture is not read by eye
# =============================================================================
cat("=== FIG 5 PART 3, gene-level D | pooled gene ACCESSIBILITY per consensus meta-cell ===\n\n")
cat("perkb source : ", PERKB, "\n", sep = "")
cat("partition    : consensus meta-cells (100-seed F consensus)\n\n")

diag_stats <- rbindlist(lapply(names(CFG), function(g) rbindlist(lapply(S_LEV, function(s) {
  gi <- gord[[g]][scored == TRUE, .(geneID, type_label)]
  L  <- merge(pooled[genome == g & stage == s], gi, by = "geneID")
  L  <- merge(L, calls[[prefix_for(g, s)]][, .(metacell, top_type)], by = "metacell")
  L[, on := type_label == top_type]
  r  <- rank(L$scaled); ni <- sum(L$on); no <- sum(!L$on)
  data.table(genome = g, stage = s,
             mean_in  = round(L[on == TRUE,  mean(scaled)], 3),
             mean_out = round(L[on == FALSE, mean(scaled)], 3),
             AUC      = round((sum(r[L$on]) - ni * (ni + 1) / 2) / (ni * no), 3))
}))))
print(diag_stats)
cat("\n  AUC = P(a random on-diagonal gene x meta-cell value exceeds a random\n",
    "  off-diagonal one). 0.5 is no discrimination. Read AUC, not the ratio of\n",
    "  means: the off-diagonal mass sits near zero, which inflates any ratio.\n",
    "  NOTE: COLOUR IS NOT ACCESSIBILITY -- each gene is min-max scaled to 0-1, so a\n",
    "  flat gene looks as colourful as a specific one. Unscaled pooled values are\n",
    "  in ", basename(CACHE), ".\n", sep = "")

comp <- (P[["SM2_B73_B73v5"]] | P[["Clean.SM2v2wd_B73_B73v5"]] | P[["Clean.SM2v2_B73_B73v5"]] | cbar) /
        (P[["SM2_At_TAIR10"]] | P[["Clean.SM2v2wd_At_TAIR10"]] | P[["Clean.SM2v2_At_TAIR10"]] | plot_spacer()) +
  plot_layout(widths = c(1, 1, 1, 0.17)) +
  plot_annotation(
    title = "D. Pooled gene accessibility of the consensus meta-cells (per-kb, min-max scaled 0-1)",
    subtitle = paste0(
      "The ACCESSIBILITY ITSELF, not a z-score. Cells are pooled into their consensus meta-cell (the partition built from the 100 SEACells seeds via the\n",
      "co-assignment matrix F), each meta-cell is scaled to a common total OVER ALL GENES so depth cannot drive the picture, and only then is the matrix\n",
      "subset to the marker panel. y = markers grouped by the type they mark; x = meta-cells grouped by their call (panel B's ordering). ORANGE BOXES mark\n",
      "each call group against its own markers - where a correct call should put its signal.\n",
      "*** COLOUR IS NOT ACCESSIBILITY. Each gene is min-max scaled to 0-1 across the genome's three stages, so every gene spans the full range whether or\n",
      "not it varies - a flat gene looks as colourful as a sharply specific one. That is the standard marker-heatmap convention and it is what makes the block\n",
      "structure visible, but intensity must never be read as an accessibility level. The unscaled pooled values and each gene's range are exported. ***\n",
      "The range is taken over the three stages TOGETHER so the panels of a row are comparable; scaling per panel would make every stage look alike.\n",
      "The marker panel is wider than the scored type set (maize 221 markers / 39 labels / 20 scored; At 275 / 20 / 18) - unscored genes sit in the trailing\n",
      "block and can never form a diagonal box. ",
      if (length(lost_note)) paste0("Marker sets differ across stages (", paste(lost_note, collapse = "; "),
                                    "), so each row is intersected to its shared genes. ")
      else "",
      "Partitions differ per stage: a column is not the same meta-cell across panels."),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 6.3, colour = "grey35",
                                               lineheight = 1.28)))

# =============================================================================
# EXPORT + VERIFY
# =============================================================================
written <- CACHE
for (ext in c("pdf", "png")) {
  f <- file.path(OUTDIR, sprintf("%s.%s", STEM_D, ext))
  # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
  if (ext == "pdf") ggsave(f, comp, width = 13.2, height = 11.0, bg = "white")
  else              ggsave(f, comp, width = 13.2, height = 11.0, dpi = 300, bg = "white")
  written <- c(written, f)
}
f <- file.path(OUTDIR, paste0(STEM_D, "_genes.tsv"))
fwrite(merge(rbindlist(gord), rngs, by = c("genome", "geneID"))[
  order(genome, ypos), .(genome, geneID, name, type_label, scored, ypos,
                         mean_access = round(strength, 5),
                         scale_lo = round(lo, 5), scale_hi = round(hi, 5), flat)],
  f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, paste0(STEM_D, "_diag.tsv")); fwrite(diag_stats, f, sep = "\t")
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
cat("\n[done] ", STEM_D, ".{pdf,png} + cache + 2 TSVs -> ", OUTDIR, "/\n", sep = "")
