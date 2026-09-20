#!/usr/bin/env Rscript
###############################################################################
## 4_3p_metacell_rZ.R
##
## Reciprocal ("bidirectional") marker z-score at the META-CELL level -- the 4_3i idea
## applied to SEACells meta-cells instead of single cells.
##
##   Zi[g,mc] = z of gene g ACROSS META-CELLS (row z)
##   Zj[g,mc] = z of meta-cell ACROSS GENES   (col z)   <- requires a perkb matrix
##   rZ_euclid = sqrt( max0(Zi)^2 + max0(Zj)^2 )     (OR-ish, formula as written)
##   rZ_geom   = sqrt( max0(Zi) *  max0(Zj) )        (AND -- the project's PRIMARY metric)
##
## WHY META-CELLS: Test A (4_3n) aggregates over a WHOLE cluster, but maize clusters are
## 0.33-0.46 dominant-label mixtures, so a 30% cell type is diluted ~3x. Meta-cells are
## sub-cluster units (measured purity 0.96-1.00 here), so identity that is real but
## subset-carried should survive aggregation at this granularity when it does not at cluster
## granularity.
##
## NON-CIRCULARITY: meta-cell MEMBERSHIP comes from ACR space (Step 5.1/5.3); the markers
## enter only here, in gene space. Markers never touch construction.
##
## ------------------------------------------------------------------------------------
## THREE DESIGN POINTS, ALL OF WHICH BIT SOMETHING EARLIER IN THIS PROJECT
##
## (1) AGGREGATE THE RAW perkb MATRIX, NEVER THE 4_3f-SMOOTHED ONE. Pooling cells into a
##     meta-cell IS smoothing; scoring the kNN-diffused matrix on top of it is double
##     smoothing. (Also why this runs locally: sparse perkb is 10-70 MB, smoothed maize 3.7 GB.)
##
## (2) Zi CAPS AT (N-1)/sqrt(N), where N = number of meta-cells. This is the failure mode
##     that broke the per-CLUSTER rZ form (N=5 -> cap 1.79 -> Zj dominates; see 4_3i header).
##     Meta-cells sit between that and per-cell (N~1e3 -> 32.6):
##         At 1 seed  N=14  -> 3.47      At 5 seeds N=70 -> 8.25
##         maize      N=286 -> 16.85
##     The cap is PRINTED and written to the results md -- read it before trusting euclid.
##     geom is the more robust combiner here (a product is limited by a capped Zi rather
##     than swamped by an uncapped Zj).
##
## (3) DEPTH. Step 5.3's RAREFY_TO equalised the ACR matrix. The gene matrix is a SEPARATE
##     count matrix built from the same fragments, so it carries the full unrarefied stage
##     difference (At nd was +55% deeper). perkb values are FRACTIONAL (counts / kb), so a
##     hypergeometric rarefaction cannot be applied here and the integer gene matrix was
##     never built for the plate arm. What this script does instead is scale every meta-cell
##     to a COMMON TOTAL, which removes the SCALE difference exactly. It does NOT remove the
##     PRECISION difference (a deeper meta-cell has less shot noise, so its profile is
##     genuinely peakier, which inflates Zi for that stage). REPORT THAT RESIDUAL; removing
##     it needs the gene matrices rebuilt in integer counts.
##     Because every meta-cell is scaled to the same total, the constant is arbitrary --
##     z-scores are invariant to a global rescale -- so no cross-stage coordination is
##     needed here (unlike RAREFY_TO, which had to be one shared value).
##
## Usage:
##   Rscript 4_3p_metacell_rZ.R <out_dir> <prefix> <perkb_sparse_rds> <c2s_csv_list> \
##           <markers_bed> <stage> [topN=6]
##     <c2s_csv_list> = comma-separated cell_to_seacell.tsv paths. More than one => the
##     seeds are POOLED (meta-cell IDs prefixed s<k>:) which raises the Zi cap.
##     Pooled seeds are PSEUDO-REPLICATES: they re-partition the SAME cells, so N is not
##     an independent sample size. Descriptive annotation only -- never feed it to a null.
##
## Outputs (under <out_dir>):
##   <prefix>.metacell_rZ.cluster_mean.tsv     marker x cluster mean rZ, both metrics
##   <prefix>.best_markers_by_cluster.tsv      topN specific markers per cluster, each metric
##   <prefix>.permetacell_type_rZ.<metric>.tsv META-CELL x TYPE mean rZ  <- annotation basis
##   <prefix>.metacell_table.tsv               per meta-cell: cluster, purity, n_cells, depth
##   <prefix>.4_3p.results.md                  manifest + Zi cap + headline calls
###############################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("Usage: Rscript 4_3p_metacell_rZ.R <out_dir> <prefix> <perkb_sparse_rds> <c2s_csv_list> <markers_bed> <stage> [topN]")
out_dir  <- args[1]; prefix <- args[2]; perkb_rds <- args[3]
c2s_list <- strsplit(args[4], ",")[[1]]
markers_bed <- args[5]; stage <- args[6]
topN <- if (length(args) >= 7) as.integer(args[7]) else 6L

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
op <- function(s) file.path(out_dir, paste0(prefix, s))
message(" - 4_3p meta-cell reciprocal z | prefix=", prefix, " stage=", stage,
        " | ", length(c2s_list), " assignment file(s)")

## ------------------------------------------------------------------ (1) markers ----------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")
markers$type_label <- paste0(markers$species, ":", markers$type)
message(" - markers (deduped): ", nrow(markers))

## ------------------------------------------------- (2) meta-cell membership --------------
## Pooling seeds: meta-cell IDs are namespaced s<k>: so identically-named SEACell-0 from two
## seeds do not silently merge into one column.
mem <- do.call(rbind, lapply(seq_along(c2s_list), function(k) {
  f <- c2s_list[k]
  if (!file.exists(f)) stop("missing assignment file: ", f)
  d <- read.table(f, header = TRUE, sep = "\t", quote = "", comment.char = "")
  colnames(d)[1] <- "cellID"
  if (!"SEACell" %in% colnames(d)) stop("no SEACell column in ", f)
  seed_tag <- sub(".*seed([0-9]+).*", "\\1", f)
  if (identical(seed_tag, f)) seed_tag <- as.character(k)
  data.frame(cellID = d$cellID,
             mc = if (length(c2s_list) > 1) paste0("s", seed_tag, ":", d$SEACell) else d$SEACell,
             cluster = as.character(d$LouvainClusters), seed = seed_tag)
}))
message(" - assignments: ", nrow(mem), " cell-rows | ", length(unique(mem$mc)), " meta-cells")

## ------------------------------------------------- (3) perkb matrix, aggregate -----------
sm <- readRDS(perkb_rds)                                  # genes x cells (sparse, perkb)
keep <- intersect(colnames(sm), mem$cellID)
if (!length(keep)) stop("no overlap between matrix colnames and assignment cellIDs")
message(" - matrix ", nrow(sm), " genes x ", ncol(sm), " barcodes | ",
        length(keep), " assigned cells matched")
mem <- mem[mem$cellID %in% keep, ]
sm  <- sm[, unique(mem$cellID), drop = FALSE]

## indicator (cells x meta-cells); agg = sm %*% Ind sums member cells per meta-cell.
mc_lv <- unique(mem$mc)
Ind <- sparseMatrix(i = match(mem$cellID, colnames(sm)), j = match(mem$mc, mc_lv),
                    x = 1, dims = c(ncol(sm), length(mc_lv)),
                    dimnames = list(colnames(sm), mc_lv))
agg <- as.matrix(sm %*% Ind)                              # genes x meta-cells (small, dense)
n_genes <- nrow(agg); n_mc <- ncol(agg)

## per-meta-cell bookkeeping BEFORE normalisation (raw pooled perkb mass = the depth proxy)
mc_tot <- colSums(agg)
mc_n   <- as.integer(table(mem$mc)[mc_lv])
mc_cl  <- tapply(mem$cluster, mem$mc, function(v) names(sort(table(v), decreasing = TRUE))[1])[mc_lv]
mc_pur <- tapply(mem$cluster, mem$mc, function(v) max(table(v)) / length(v))[mc_lv]
mc_sd  <- tapply(mem$seed,    mem$mc, function(v) v[1])[mc_lv]

## DEPTH MATCH: common total per meta-cell. Fixes scale exactly; leaves the precision
## residual (see header point 3). The constant cancels in both z-scores.
agg <- sweep(agg, 2, pmax(mc_tot, .Machine$double.eps), "/") * 1e4
message(" - aggregated -> ", n_genes, " genes x ", n_mc, " meta-cells",
        " | pooled perkb mass median ", round(median(mc_tot)),
        " (range ", round(min(mc_tot)), "-", round(max(mc_tot)), ")")

## ------------------------------------------------- (4) the rZ math (as 4_3i) -------------
zi_cap <- (n_mc - 1) / sqrt(n_mc)
message(" - Zi cap at N=", n_mc, " meta-cells: ", round(zi_cap, 2),
        if (zi_cap < 5) "  *** LOW -- prefer geom; euclid is Zj-dominated here ***" else "")

cmean <- colMeans(agg)
csd   <- apply(agg, 2, sd); csd[csd == 0 | is.na(csd)] <- 1

mk_genes <- intersect(rownames(agg), markers$geneID)
message(" - markers present in matrix: ", length(mk_genes), " / ", nrow(markers))
if (length(mk_genes) < 3) stop("fewer than 3 markers present -- check the panel/genome pairing")
smk <- agg[mk_genes, , drop = FALSE]

rmean <- rowMeans(smk)
rsd   <- apply(smk, 1, sd); rsd[rsd == 0 | is.na(rsd)] <- 1

Zi <- sweep(sweep(smk, 1, rmean, "-"), 1, rsd, "/")       # gene z ACROSS meta-cells
Zj <- sweep(sweep(smk, 2, cmean, "-"), 2, csd,  "/")      # meta-cell z ACROSS genes
Zi[Zi < 0] <- 0; Zj[Zj < 0] <- 0
rZe <- sqrt(Zi^2 + Zj^2)
rZg <- sqrt(Zi *  Zj)
rm(Zi, Zj); gc()

## ------------------------------------------------- (5) cluster-level summaries -----------
cls <- unique(mc_cl)
cls <- if (suppressWarnings(all(!is.na(as.numeric(cls))))) cls[order(as.numeric(cls))] else sort(cls)
mc_by_cluster <- split(mc_lv, mc_cl)

summarize <- function(rZ) {
  mrz <- sapply(cls, function(cl) rowMeans(rZ[, mc_by_cluster[[cl]], drop = FALSE]))
  if (is.null(dim(mrz))) mrz <- matrix(mrz, nrow = length(mk_genes), dimnames = list(mk_genes, cls))
  colnames(mrz) <- cls; rownames(mrz) <- mk_genes
  pk  <- apply(mrz, 1, which.max)
  pv  <- mrz[cbind(seq_len(nrow(mrz)), pk)]
  sec <- apply(mrz, 1, function(v) if (length(v) >= 2) sort(v, decreasing = TRUE)[2] else 0)
  list(mean = mrz, peak = colnames(mrz)[pk], peak_val = pv, spec = pv - sec)
}
Se <- summarize(rZe); Sg <- summarize(rZg)

cm <- data.frame(stage = stage, geneID = mk_genes,
                 name = markers[mk_genes, "name"], type_label = markers[mk_genes, "type_label"],
                 peak_cl_euclid = Se$peak, peak_rZ_euclid = round(Se$peak_val, 4), spec_euclid = round(Se$spec, 4),
                 peak_cl_geom   = Sg$peak, peak_rZ_geom   = round(Sg$peak_val, 4), spec_geom   = round(Sg$spec, 4))
write.table(cm[order(cm$peak_cl_geom, -cm$spec_geom), ], op(".metacell_rZ.cluster_mean.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pick_best <- function(S, tag) {
  do.call(rbind, lapply(cls, function(cl) {
    idx <- which(S$peak == cl); if (!length(idx)) return(NULL)
    d <- data.frame(geneID = mk_genes[idx], name = markers[mk_genes[idx], "name"],
                    type_label = markers[mk_genes[idx], "type_label"], peak_cluster = cl,
                    peak_rZ = round(S$peak_val[idx], 4), specificity = round(S$spec[idx], 4), metric = tag)
    head(d[order(-d$specificity), ], topN)
  }))
}
best_e <- pick_best(Se, "euclid"); best_g <- pick_best(Sg, "geom")
write.table(rbind(best_e, best_g), op(".best_markers_by_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------- (6) META-CELL x TYPE -- the annotation basis -------------
## Same rule as 4_3i: cell-cycle 'dividing' excluded (cross-tissue confound), >=3 markers/type.
tl_tab <- table(markers$type_label[markers$geneID %in% mk_genes])
tl_tab <- tl_tab[!grepl("dividing", names(tl_tab), ignore.case = TRUE)]
type_levels <- names(tl_tab)[tl_tab >= 3]
if (!length(type_levels)) stop("no marker type has >=3 genes present -- panel/genome mismatch?")

agg_types <- function(rZmat) {
  m <- sapply(type_levels, function(tl) {
    g <- intersect(markers$geneID[markers$type_label == tl], rownames(rZmat))
    if (!length(g)) rep(NA_real_, ncol(rZmat)) else colMeans(rZmat[g, , drop = FALSE])
  })
  rownames(m) <- colnames(rZmat); round(m, 4)
}
calls <- list()
for (mn in c("geom", "euclid")) {
  pct <- agg_types(if (mn == "geom") rZg else rZe)
  top <- colnames(pct)[apply(pct, 1, which.max)]
  mx  <- apply(pct, 1, max)
  sec <- apply(pct, 1, function(v) sort(v, decreasing = TRUE)[2])
  calls[[mn]] <- data.frame(metacell = rownames(pct), top_type = top,
                            top_rZ = round(mx, 4), margin = round(mx - sec, 4))
  write.table(data.frame(metacell = rownames(pct), seed = mc_sd[rownames(pct)],
                         cluster = mc_cl[rownames(pct)], n_cells = mc_n[match(rownames(pct), mc_lv)],
                         top_type = top, top_rZ = round(mx, 4), margin = round(mx - sec, 4),
                         pct, check.names = FALSE),
              op(paste0(".permetacell_type_rZ.", mn, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
message(" - meta-cell x type matrices written (", length(type_levels), " types, dividing excluded)")

write.table(data.frame(metacell = mc_lv, seed = mc_sd, cluster = mc_cl,
                       purity = round(mc_pur, 4), n_cells = mc_n,
                       pooled_perkb_mass = round(mc_tot, 1)),
            op(".metacell_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## MARKER x META-CELL rZ. Needed for anything that compares meta-cell PROFILES rather than
## their winning label -- e.g. within- vs between-cluster rZ correlation (4_3q). Small
## (markers x meta-cells), so both metrics are written.
for (mn in c("geom", "euclid")) {
  M <- if (mn == "geom") rZg else rZe
  write.table(data.frame(geneID = rownames(M), name = markers[rownames(M), "name"],
                         type_label = markers[rownames(M), "type_label"],
                         round(M, 4), check.names = FALSE),
              op(paste0(".permetacell_marker_rZ.", mn, ".tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
}

## ------------------------------------------------- (7) results md ------------------------
tally <- table(calls$geom$top_type)
tally <- sort(tally, decreasing = TRUE)
md <- c(
  paste0("# 4_3p meta-cell reciprocal z-score — ", prefix, " (", stage, ")"), "",
  paste0("- perkb matrix: `", perkb_rds, "`  (RAW, not 4_3f-smoothed — meta-cell pooling is already smoothing)"),
  paste0("- assignments: ", length(c2s_list), " file(s)",
         if (length(c2s_list) > 1) "  NOTE: SEEDS POOLED — pseudo-replicates, descriptive only, never an N for a null" else ""),
  paste0("- markers: `", markers_bed, "`  (", length(mk_genes), " of ", nrow(markers), " present)"),
  paste0("- meta-cells: ", n_mc, " | genes: ", n_genes, " | clusters: ", paste(cls, collapse = ", ")),
  paste0("- cells per meta-cell: median ", round(median(mc_n)), " (", min(mc_n), "-", max(mc_n), ")"),
  paste0("- meta-cell cluster purity: median ", round(median(mc_pur), 3)), "",
  "## NOTE: Zi cap",
  paste0("Zi is a z across the ", n_mc, " meta-cells, so it cannot exceed (N-1)/sqrt(N) = **",
         round(zi_cap, 2), "**."),
  paste0("Zj is a z across ", n_genes, " genes and is effectively uncapped. ",
         if (zi_cap < 5) "**This cap is LOW: `euclid` is Zj-dominated — read `geom`.**"
         else "Both axes are usable; `geom` remains primary."), "",
  "## NOTE: Depth residual",
  paste0("Meta-cells are scaled to a common total, which removes the SCALE difference exactly. ",
         "It does NOT remove the PRECISION difference: Step 5.3's rarefaction equalised the ACR ",
         "matrix, not this gene matrix, so a deeper stage still yields intrinsically peakier ",
         "profiles and a mildly inflated Zi. Pooled perkb mass here: median ", round(median(mc_tot)),
         " (", round(min(mc_tot)), "-", round(max(mc_tot)), ")."), "",
  "## Meta-cell type calls — geom (AND, primary)",
  "| type | n meta-cells |", "|---|---|",
  paste0("| ", names(tally), " | ", as.integer(tally), " |"), "",
  "## Best marker per cluster — geom",
  "| cluster | best | peak_rZ | specificity | type |", "|---|---|---|---|---|",
  apply(best_g[!duplicated(best_g$peak_cluster), ], 1, function(r)
    paste0("| ", r["peak_cluster"], " | ", r["name"], " | ", r["peak_rZ"], " | ",
           r["specificity"], " | ", r["type_label"], " |")), "")
writeLines(md, op(".4_3p.results.md"))
message(" - wrote tables + results md -> ", out_dir)
message(" - DONE")
