#!/usr/bin/env Rscript
# =============================================================================
# fig5_I_consensus.R  -  Fig 5 panel I: the 100-seed SEACells co-assignment matrix F for all six
#   arms (2 genomes x PreClean / WD / ND) with the consensus meta-cell blocks outlined on the
#   diagonal; At shows the complete matrix, maize a fixed 20-block window at native resolution.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/consensus/
#           <prefix>.consensus.{F.rds,cells.tsv,summary.tsv,ksweep.tsv,metacells.tsv}
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P3A_consensus_Fmatrix.{pdf,png} + Fig5_P3A_consensus_Fmatrix_blocks.tsv
# Run     Rscript analysis/fig5_biological_impact/fig5_I_consensus.R    (from the repo root, ~2 min)
# Part 3 internal lettering kept in titles and cross-references: A = this panel = manuscript I.
# =============================================================================
#
# F[i,j] = fraction of the 100 SEACells runs in which cells i and j share a meta-cell
# (upstream: sparse tcrossprod of the seed indicator / 100). The consensus = average-linkage
# hclust on (1 - F), k PINNED to the SEACells target (At 14, maize 286).
#
# DISPLAY: At = the COMPLETE matrix; maize = a ~20-block WINDOW at native resolution (the full
# k = 286 matrix is too big to fit the final plot: a block is ~0.8 pt at print size and the
# panel is ~98% empty). The window is a FIXED RULE, not a curated sample: contiguous blocks
# centred on the midpoint of the block ordering. Upstream previews sample ~30 meta-cells; this
# is deterministic.
#
# ORDERING: the PUBLISHED consensus blocks drive the layout, dendrograms are used only to
# ARRANGE them. A naive full-tree recompute (average linkage on 1-F) is tie-sensitive: measured,
# it reproduces the upstream blocks exactly on 4/6 arms but drifts by 1 cell on maize wd and nd
# (F is quantised to 1/100, merge heights tie, and hclust's tie-break depends on input order).
# The published labels are the figure's ground truth, so instead:
#   between blocks: average linkage on (1 - mean F between block members)
#   within a block: average linkage on the block's own (1 - F) submatrix
# Blocks are contiguous BY CONSTRUCTION (asserted anyway); a tie can only nudge aesthetics,
# never the delimitation.
#
# PDF SIZE: each matrix is embedded as a COMPRESSED RASTER inside an otherwise vector PDF
# (annotation_raster, the same idea as ComplexHeatmap's use_raster); text and boxes stay
# vector. Vector rectangles per cell pair would be ~208M objects for one full maize matrix.
# With At complete (~1k cells) and maize windowed (~0.7-1.5k cells) every panel is native
# resolution, one pixel per cell pair, no binning anywhere, and the whole figure stays a few MB.
#
# TRAPS:
#   The maize panels show a WINDOW, not the whole matrix. The stamps (gate, split-half, n, k)
#     always refer to the FULL consensus, and the complete block layout is exported to the
#     blocks TSV (`shown` marks the window).
#   AGREEMENT / CONSISTENCY, never "accuracy" -- no ground truth anywhere.
#   IN-SAMPLE: F is built from the same 100 seeds that define the consensus; the non-circular
#     check is the split-half stability, stamped per panel.
#   The stored diagonal of F is 0 (self-pairs excluded upstream): it reads as a hairline
#     through dark blocks. Faithful to the input, not a bug.
#   Larger blocks are mechanically less coherent (~ -0.14 per doubling, stage-independent):
#     block AREA shows size, so the confound is visible.
#   wd and nd are INDEPENDENT treatments of the same raw input, NOT a chain.
#   Meta-cell IDs are NOT comparable across stages; boxes are per panel.
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

STEM <- "Fig5_P3A_consensus_Fmatrix"
CONS <- file.path(PLATE, "step5_metacell", "consensus")
N_SEED <- 100L
WIN_K  <- 20L          # arms with more blocks than this show a WIN_K-block window

G_LEV <- vapply(CFG, `[[`, character(1), "label")
S_LEV <- names(STAGES)
STAGE_UP <- c(PreClean = "Pre", wd = "wd", nd = "nd")

# white -> orange -> dark maroon, matching the upstream previews
PAL    <- colorRampPalette(c("#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33",
                             "#7F0000", "#250000"))(256)
BOXCOL <- "#2166AC"    # consensus-block outline: blue, visible on white + maroon

# F values -> raster; raster row 1 = top, so flip rows: position 1 = bottom-left
to_ras <- function(m) {
  idx <- matrix(1L + as.integer(round(pmin(pmax(m, 0), 1) * 255)),
                nrow(m), ncol(m))
  as.raster(matrix(PAL[idx], nrow(m), ncol(m))[nrow(m):1, , drop = FALSE])
}

prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)
cons_file  <- function(g, s, what)
  file.path(CONS, sprintf("%s.consensus.%s.tsv", prefix_for(g, s), what))

# --- stamps: gate + split-half at pinned k ---------------------------------------
read_stamp <- function(g, s) {
  summ <- fread(cons_file(g, s, "summary"))
  ksw  <- fread(cons_file(g, s, "ksweep"))
  if (nrow(summ) != 1L || summ$stage != STAGE_UP[[s]])
    stop(prefix_for(g, s), ": summary stage mismatch")
  if (!isTRUE(summ$gate_pass))
    stop(prefix_for(g, s), ": upstream acceptance gate FAILED -- do not plot this arm")
  if (!summ$k_pinned %in% ksw$k)
    stop(prefix_for(g, s), ": ksweep has no row at the pinned k")
  cbind(summ, splithalf_mean = ksw[k == summ$k_pinned, mean],
              splithalf_sd   = ksw[k == summ$k_pinned, sd])
}

# =============================================================================
# per arm: load F, recompute + validate the consensus tree, build the raster
# =============================================================================
panels <- list(); blocks_out <- list(); report <- list()

for (g in names(CFG)) for (s in S_LEV) {
  pf <- prefix_for(g, s)
  frds <- file.path(CONS, paste0(pf, ".consensus.F.rds"))
  if (!file.exists(frds)) stop("missing upstream consensus F: ", frds)
  x <- readRDS(frds)
  n <- length(x$cells); k <- x$k
  stamp <- read_stamp(g, s)

  # ---- upstream data contract ----------------------------------------------
  if (!identical(dim(x$F), c(n, n)) || length(x$consensus) != n)
    stop(pf, ": F / cells / consensus dimensions disagree")
  if (n != stamp$n_cells || k != stamp$k_pinned)
    stop(pf, ": F.rds n/k disagree with summary.tsv")
  if (uniqueN(x$consensus) != k)
    stop(pf, ": consensus label count != pinned k")
  ctsv <- fread(cons_file(g, s, "cells"))
  m <- merge(data.table(cellID = x$cells, mc_rds = x$consensus),
             ctsv[, .(cellID, mc_tsv = consensus_mc)], by = "cellID")
  if (nrow(m) != n || m[mc_rds != mc_tsv, .N] > 0)
    stop(pf, ": consensus labels in F.rds disagree with cells.tsv")

  # ---- ordering: published blocks first, dendrograms only for arrangement --
  t0 <- proc.time()[3]
  lab <- x$consensus
  idx_by_block <- split(seq_len(n), lab)
  bsz     <- lengths(idx_by_block)
  blk_ids <- as.integer(names(idx_by_block))
  tr0 <- as.data.table(summary(x$F))[x > 0]      # stored triangle (i <= j)

  # between-block mean F from the sparse triplets (off-diagonal mirrored)
  bb <- copy(tr0)[, `:=`(bi = lab[i], bj = lab[j])]
  bb <- rbindlist(list(bb[, .(bi, bj, x)],
                       bb[bi != bj | i != j][, .(bi = bj, bj = bi, x)]))[
    , .(sx = sum(x)), by = .(bi, bj)]
  B <- matrix(0, k, k)
  pi_ <- match(bb$bi, blk_ids); pj_ <- match(bb$bj, blk_ids)
  B[cbind(pi_, pj_)] <- bb$sx / (bsz[pi_] * bsz[pj_])
  border <- if (k > 2L) blk_ids[hclust(as.dist(1 - B), "average")$order] else blk_ids

  # within-block leaf order on the block's dense submatrix (blocks <= ~370)
  ord <- unlist(lapply(border, function(b) {
    ix <- idx_by_block[[as.character(b)]]
    if (length(ix) <= 2L) return(ix)
    ix[hclust(as.dist(1 - as.matrix(x$F[ix, ix])), "average")$order]
  }), use.names = FALSE)
  if (length(ord) != n || anyDuplicated(ord) > 0L)
    stop(pf, ": ordering lost or duplicated cells")
  r <- rle(lab[ord])
  if (length(r$lengths) != k)
    stop(pf, ": consensus blocks not contiguous after ordering")
  t_tree <- proc.time()[3] - t0

  # block boundaries in display positions (1-based, inclusive)
  ends   <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  mcs    <- fread(cons_file(g, s, "metacells"))
  blk <- data.table(genome = g, stage = s, consensus_mc = r$values,
                    ord_start = starts, ord_end = ends, n_cells = r$lengths)
  blk <- merge(blk, mcs[, .(consensus_mc, conf_meanF)],
               by = "consensus_mc", sort = FALSE)
  if (blk[, sum(n_cells)] != n) stop(pf, ": block sizes do not sum to n")
  blocks_out[[pf]] <- blk

  # ---- display: complete matrix (At) or a WIN_K-block native window (maize) -
  # The window is DETERMINISTIC, not curated: contiguous blocks centred on the
  # block holding the midpoint of the display ordering. The full ordering of
  # every block stays in the blocks TSV (`shown` marks the plotted window),
  # and the stamps (gate, split-half, n, k) always refer to the FULL consensus.
  windowed <- k > WIN_K
  if (!windowed) {
    lo <- 1L; hi <- k
  } else {
    mid <- which(blk$ord_start <= ceiling(n / 2) & blk$ord_end >= ceiling(n / 2))[1]
    lo <- hi <- mid
    while (hi - lo + 1L < WIN_K && !(lo == 1L && hi == k)) {
      if (lo > 1L) lo <- lo - 1L
      if (hi - lo + 1L >= WIN_K) break
      if (hi < k) hi <- hi + 1L
    }
  }
  s_win <- blk$ord_start[lo]; e_win <- blk$ord_end[hi]
  nshow <- e_win - s_win + 1L
  blk[, shown := seq_len(.N) %in% lo:hi]
  ras <- to_ras(as.matrix(x$F[ord[s_win:e_win], ord[s_win:e_win]]))
  blk_show <- blk[lo:hi][, .(bx0 = ord_start - s_win,   # 0-based window coords
                             bx1 = ord_end - s_win + 1L)]

  report[[pf]] <- data.table(
    genome = g, stage = s, cells = n, k = k, singletons = stamp$n_singletons,
    order_s = round(t_tree, 1),
    display = if (windowed) sprintf("window %d/%d blocks (%s cells), native",
                                    hi - lo + 1L, k,
                                    trimws(format(nshow, big.mark = ",")))
              else "native (full)",
    gate = sprintf("%.2f > %.2f", stamp$ari_cons_seed, stamp$ari_seed_seed),
    splithalf = sprintf("%.2f (sd %.2f)", stamp$splithalf_mean, stamp$splithalf_sd))

  # ---- panel ---------------------------------------------------------------
  ttl <- sprintf("%s - %s", G_LEV[[g]], s)
  sub <- paste0(
    sprintf("full consensus: %s cells, %d meta-cells%s",
            trimws(format(n, big.mark = ",")), k,
            fifelse(stamp$n_singletons > 0,
                    sprintf(" (%d singletons)", stamp$n_singletons), "")),
    "\n",
    sprintf("gate ARI %.2f > %.2f | split-half %.2f (sd %.2f)",
            stamp$ari_cons_seed, stamp$ari_seed_seed,
            stamp$splithalf_mean, stamp$splithalf_sd),
    "\n",
    if (windowed) sprintf("shown: %d meta-cells (%s cells) around the ordering midpoint, native resolution",
                          hi - lo + 1L, trimws(format(nshow, big.mark = ",")))
    else "shown: the complete matrix, native resolution (one pixel per cell pair)")
  panels[[pf]] <- ggplot() +
    annotation_raster(ras, xmin = 0, xmax = nshow, ymin = 0, ymax = nshow,
                      interpolate = FALSE) +
    geom_rect(data = blk_show, aes(xmin = bx0, xmax = bx1, ymin = bx0, ymax = bx1),
              fill = NA, colour = BOXCOL, linewidth = 0.35) +
    annotate("rect", xmin = 0, xmax = nshow, ymin = 0, ymax = nshow,
             fill = NA, colour = "grey25", linewidth = 0.3) +
    coord_fixed(xlim = c(0, nshow), ylim = c(0, nshow), expand = FALSE) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 8.6, hjust = 0.5),
          plot.subtitle = element_text(size = 5.6, colour = "grey35",
                                       hjust = 0.5, lineheight = 1.25),
          plot.margin = margin(2, 6, 4, 6)) +
    labs(title = ttl, subtitle = sub)
  message(sprintf("  [%s] n=%d k=%d order %.1fs %s",
                  pf, n, k, t_tree, report[[pf]]$display))
}

# =============================================================================
# REPORT
# =============================================================================
cat("=== FIG 5 PANEL I (Part 3 A) | consensus co-assignment matrix, ", N_SEED,
    " seeds ===\n\n", sep = "")
cat("consensus source: ", CONS, "\n\n", sep = "")
rp <- rbindlist(report)[order(genome, factor(stage, levels = S_LEV))]
print(rp)
cat("\n  ordering: PUBLISHED consensus blocks, arranged by average linkage on the\n",
    "  between-block mean F; cells within a block by its own sub-dendrogram. Blocks\n",
    "  contiguous by construction (asserted). A naive full-tree recompute reproduces\n",
    "  the published blocks exactly on 4/6 arms but drifts by 1 cell on maize wd/nd\n",
    "  (F is quantised to 1/100 -> hclust tie-breaks depend on input order), which is\n",
    "  why the published labels drive the layout rather than a recomputed tree.\n",
    "  The stored diagonal of F is 0 (self-pairs excluded upstream) -> hairline diagonal.\n",
    "  display: arms with more than ", WIN_K, " blocks show a ", WIN_K, "-block window around\n",
    "  the ordering midpoint (the full k = 286 matrix cannot fit a printed panel);\n",
    "  stamps refer to the FULL consensus; blocks TSV `shown` marks the window.\n",
    sep = "")

# =============================================================================
# ASSEMBLY: maize row | At row, shared colourbar column
# =============================================================================
cbar <- ggplot() +
  annotation_raster(as.raster(matrix(rev(PAL), ncol = 1)),
                    xmin = 0, xmax = 1, ymin = 0, ymax = 1, interpolate = TRUE) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = NA, colour = "grey25", linewidth = 0.3) +
  annotate("segment", x = 1, xend = 1.35, y = c(0, .25, .5, .75, 1),
           yend = c(0, .25, .5, .75, 1), linewidth = 0.3, colour = "grey25") +
  annotate("text", x = 1.65, y = c(0, .25, .5, .75, 1),
           label = c("0", "0.25", "0.50", "0.75", "1"),
           hjust = 0, size = 2.3, colour = "grey20") +
  coord_cartesian(xlim = c(-0.2, 4.2), ylim = c(-0.03, 1.03), expand = FALSE) +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 6.4, colour = "grey20", hjust = 0),
        plot.margin = margin(28, 2, 30, 2)) +
  labs(title = "co-assignment\nF (fraction of\n100 seeds)")

ordkey <- function(g, s) prefix_for(g, s)
comp <- (panels[[ordkey("B73v5", "PreClean")]] | panels[[ordkey("B73v5", "wd")]] |
         panels[[ordkey("B73v5", "nd")]] | cbar) /
        (panels[[ordkey("TAIR10", "PreClean")]] | panels[[ordkey("TAIR10", "wd")]] |
         panels[[ordkey("TAIR10", "nd")]] | plot_spacer()) +
  plot_layout(widths = c(1, 1, 1, 0.14)) +
  plot_annotation(
    title = sprintf("A. Consensus meta-cells from %d SEACells seeds - the co-assignment matrix", N_SEED),
    subtitle = paste0(
      "F[i,j] = fraction of the ", N_SEED, " independent SEACells runs in which cells i and j land in the same meta-cell. ",
      "Cells are ordered so the PUBLISHED consensus blocks are contiguous:\n",
      "blocks arranged by average linkage on the between-block mean F, cells within a block by its own sub-dendrogram. ",
      "Blue boxes = the consensus meta-cells. At panels show the\nCOMPLETE matrix at native resolution. ",
      "Maize panels show a WINDOW of 20 of the 286 meta-cells, native resolution - the full 14.4k-cell matrix cannot fit a printed panel\n",
      "(a block is ~0.8 pt wide). The window is a FIXED RULE, not curated: contiguous blocks around the midpoint of the block ordering. ",
      "Panel stamps (gate, split-half, n, k) always refer\nto the FULL consensus, and the complete block layout is in the blocks TSV (column `shown` marks the plotted window). ",
      "Dark diagonal blocks = cells that co-assign in most runs;\nthe hairline on the diagonal is the stored zero self-pair. ",
      "Block area shows meta-cell size; larger blocks are mechanically less coherent (about -0.14 per doubling, stage-independent) - ",
      "do not read block darkness alone as\nstage quality. ",
      "F is measured against the same seeds that built the consensus (in-sample); the split-half number in each stamp is the non-circular stability check. ",
      "AGREEMENT, not\naccuracy - nothing is compared to a ground truth. ",
      "wd and nd are independent treatments of the same raw input, not a chain. Meta-cell IDs and orderings are per panel - never\ncompare across panels."),
    theme = theme(plot.title = element_text(face = "bold", size = 12.5),
                  plot.subtitle = element_text(size = 7.2, colour = "grey35",
                                               lineheight = 1.3)))

# =============================================================================
# EXPORT + VERIFY
# =============================================================================
written <- character()
for (ext in c("pdf", "png")) {
  f <- file.path(OUTDIR, sprintf("%s.%s", STEM, ext))
  # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
  if (ext == "pdf") ggsave(f, comp, width = 12.8, height = 9.6, bg = "white")
  else              ggsave(f, comp, width = 12.8, height = 9.6, dpi = 300, bg = "white")
  written <- c(written, f)
}
f <- file.path(OUTDIR, paste0(STEM, "_blocks.tsv"))
fwrite(rbindlist(blocks_out)[, .(genome, stage, consensus_mc, ord_start, ord_end,
                                 n_cells, conf_meanF, shown)], f, sep = "\t")
written <- c(written, f)

cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz <- if (file.exists(f)) file.size(f) else NA_integer_
  good <- !is.na(sz) && sz > (if (grepl("\\.tsv$", f)) 50 else 1000)
  ok <- ok && good
  cat(sprintf("  %-4s %-44s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
cat("\n[done] ", STEM, ".{pdf,png} + blocks TSV -> ", OUTDIR, "/\n", sep = "")
