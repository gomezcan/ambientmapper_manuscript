#!/usr/bin/env Rscript
## =================================================================================================
## 4_3z_consensus_rZ_join.R -- one tidy heatmap-input table per object: marker-type rZ scored on
## the CONSENSUS meta-cells (4_3p output) joined with the consensus confidence attributes (4_3u).
##
## WHY. The Fig-5 heatmap script is visualization-only (figure scripts load precomputed
## results), and a two-file join left to the figure side is a join-key bug waiting to
## happen. This emits ONE long-format table per object with everything the heatmap needs:
##   metacell x type rZ (geom PRIMARY, euclid secondary) + top_type/top_rZ/margin
##   + n_cells, conf_meanF, stab, loyalty_sd, loyalty_min, frac_F_ge0.9, dominant_cluster, purity.
##
## READING RULES (carry into any figure):
##   * DESCRIPTIVE until the meta-cell gene null runs -- top_type is a provisional label, not
##     annotation. Neither existing null randomizes genes.
##   * geom is the PRIMARY metric; at k=14 the Zi cap is (N-1)/sqrt(N)=3.47 (At), so euclid is
##     cap-limited -- 4_3p prints the cap, read it before trusting euclid.
##   * WEIGHT by confidence (conf_meanF / loyalty), never filter; if a figure draws a line it must
##     report the retained fraction per stage (min.c=50 failure mode).
##   * n_cells is stage-entangled for At (nd lost 34% of cells) -- condition on it in any
##     cross-stage per-meta-cell comparison.
##
## RECONCILIATION GATE (hard stop): per-meta-cell n_cells recomputed by 4_3p from the adapter
## membership must EQUAL n_cells in consensus.metacells.tsv, and the meta-cell ID sets must match.
##
## Usage:
##   Rscript 4_3z_consensus_rZ_join.R <rz_dir> <consensus_dir> <prefix> <stage> <genome>
## Reads : <rz_dir>/<prefix>.permetacell_type_rZ.{geom,euclid}.tsv   (4_3p on consensus membership)
##         <consensus_dir>/<prefix>.consensus.metacells.tsv          (4_3u)
## Writes: <rz_dir>/<prefix>.consensus_rZ.heatmap_input.tsv
## =================================================================================================

args <- commandArgs(TRUE)
if (length(args) != 5)
  stop("Usage: Rscript 4_3z_consensus_rZ_join.R <rz_dir> <consensus_dir> <prefix> <stage> <genome>")
RZD <- args[1]; CONS <- args[2]; PREFIX <- args[3]; STAGE <- args[4]; GENOME <- args[5]

g <- read.delim(file.path(RZD, paste0(PREFIX, ".permetacell_type_rZ.geom.tsv")),   check.names = FALSE)
e <- read.delim(file.path(RZD, paste0(PREFIX, ".permetacell_type_rZ.euclid.tsv")), check.names = FALSE)
m <- read.delim(file.path(CONS, paste0(PREFIX, ".consensus.metacells.tsv")))

meta_cols <- c("metacell", "seed", "cluster", "n_cells", "top_type", "top_rZ", "margin")
stopifnot(all(meta_cols %in% names(g)), identical(names(g), names(e)),
          identical(g$metacell, e$metacell))
types <- setdiff(names(g), meta_cols)

gm <- g[, meta_cols]
gm$metacell <- as.character(gm$metacell)
m$metacell  <- paste0("cMC-", m$consensus_mc)
stopifnot(setequal(gm$metacell, m$metacell))
names(m)[names(m) == "n_cells"] <- "n_cells_consensus"
j <- merge(gm, m[, c("metacell", "n_cells_consensus", "conf_meanF", "stab", "loyalty_sd",
                     "loyalty_min", "frac_F_ge0.9", "dominant_cluster", "purity")],
           by = "metacell")
if (!all(j$n_cells == j$n_cells_consensus))
  stop(PREFIX, ": RECONCILIATION FAIL -- n_cells from 4_3p vs consensus.metacells.tsv differ")
j$n_cells_consensus <- NULL

long <- do.call(rbind, lapply(types, function(tp)
  data.frame(metacell = as.character(g$metacell), type = tp,
             rZ_geom = g[[tp]], rZ_euclid = e[[tp]])))
out <- merge(long, j, by = "metacell")
out$stage <- STAGE; out$genome <- GENOME; out$prefix <- PREFIX
out$type <- factor(out$type, levels = types)                 # preserve panel order for plotting
out <- out[order(out$cluster, as.numeric(sub("^cMC-", "", out$metacell)), out$type), ]

f <- file.path(RZD, paste0(PREFIX, ".consensus_rZ.heatmap_input.tsv"))
write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[4_3z] %s: %d meta-cells x %d types -> %s", PREFIX,
                length(unique(out$metacell)), length(types), f))
message(sprintf("[4_3z] top_type spread: %d distinct types across %d meta-cells; median margin %.3f",
                length(unique(j$top_type)), nrow(j), median(j$margin)))
