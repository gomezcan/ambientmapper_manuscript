#!/usr/bin/env Rscript
###############################################################################
## 4_3i_reciprocal_zscore_percell.R
##
## Per-cell bidirectional ("reciprocal") marker z-score on a SMOOTHED PERKB
## gene-activity matrix (genes x cells; CP10k + kNN-Markov diffused by 4_3f).
##
##   Zi[g,cell] = z of gene g ACROSS CELLS  (row z): (x - mean_over_cells) / sd_over_cells   [per gene]
##   Zj[g,cell] = z of cell   ACROSS GENES  (col z): (x - mean_over_genes) / sd_over_genes   [per cell]
##
## Two combiners, computed side by side for comparison:
##   rZ_euclid = sqrt( max(0,Zi)^2 + max(0,Zj)^2 )     -- the formula as written (OR-ish)
##   rZ_geom   = sqrt( max(0,Zi) *  max(0,Zj)     )     -- AND (a gene must be high on BOTH axes)
##
## Both axes are large-N (cells ~1e4, genes ~3e4) so the Euclidean combine is balanced here (unlike
## the 5-cluster pseudobulk form, where Zi caps at 1.79). Zj compares genes within a cell, so the
## matrix MUST be perkb (gene-length-normalized).
##
## Memory-safe: the full matrix is touched ONLY for the per-cell mean/sd vectors (Zj background),
## computed in gene-blocks to avoid a full sm^2 copy. Zi/Zj/rZ are materialized for MARKERS only.
##
## Usage:
##   Rscript 4_3i_reciprocal_zscore_percell.R <out_dir> <prefix> <smoothed_rds> <meta_txt> \
##           <markers_bed> <stage> [cluster_col=LouvainClusters] [topN=6]
##
## Outputs (all under <out_dir>):
##   <prefix>.percell_rZ.cluster_mean.tsv    marker x cluster mean rZ (peak+specificity, both metrics)
##   <prefix>.best_markers_by_cluster.tsv    topN specific markers per cluster, each metric (violin cols)
##   <prefix>.percell_rZ.best.tsv            long: best-marker x cell, rZ_euclid + rZ_geom (+cluster) -> violin
##   <prefix>.4_3i.results.md                manifest + headline
###############################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("Usage: Rscript 4_3i_reciprocal_zscore_percell.R <out_dir> <prefix> <smoothed_rds> <meta_txt> <markers_bed> <stage> [cluster_col] [topN]")
out_dir     <- args[1]
prefix      <- args[2]
smoothed_rds<- args[3]
meta_txt    <- args[4]
markers_bed <- args[5]
stage       <- args[6]
cluster_col <- if (length(args) >= 7) args[7] else "LouvainClusters"
topN        <- if (length(args) >= 8) as.integer(args[8]) else 6L

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(out_dir, paste0(prefix, suffix))
message(" - 4_3i per-cell reciprocal z | prefix=", prefix, " stage=", stage,
        " | cluster_col=", cluster_col, " topN=", topN)

## ----------------------------------------------------------------- (1) markers ------------------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")
markers$type_label <- paste0(markers$species, ":", markers$type)
message(" - markers (deduped): ", nrow(markers))

## ----------------------------------------------------------------- (2) metadata -----------------
b <- read.table(meta_txt, header = TRUE, sep = "\t", quote = "", comment.char = "")
if ("cellID" %in% colnames(b)) rownames(b) <- b$cellID
if (!cluster_col %in% colnames(b))
  stop("cluster column '", cluster_col, "' not in metadata; have: ", paste(colnames(b), collapse = ", "))
b$.cluster <- as.character(b[[cluster_col]])

## ----------------------------------------------------------------- (3) smoothed matrix ----------
sm <- readRDS(smoothed_rds)                       # genes x cells (dense, CP10k+diffused)
cells <- intersect(colnames(sm), rownames(b))
if (length(cells) == 0) stop("no overlap between matrix colnames and metadata cellIDs")
sm <- sm[, cells, drop = FALSE]
b  <- b[cells, , drop = FALSE]
n_genes <- nrow(sm); n_cells <- ncol(sm)
message(" - matrix genes x cells: ", n_genes, " x ", n_cells)

## --------------------------------- Zj background: per-cell mean/sd OVER ALL GENES ----------------
## chunked sum-of-squares over gene blocks (no full sm^2 temp); sample sd via n/(n-1)
cmean <- Matrix::colMeans(sm)
csumsq <- numeric(n_cells); blk <- 2000L
for (s in seq(1L, n_genes, by = blk)) {
  e <- min(s + blk - 1L, n_genes)
  csumsq <- csumsq + Matrix::colSums(sm[s:e, , drop = FALSE]^2)
}
cvar <- (csumsq - n_genes * cmean^2) / (n_genes - 1); cvar[cvar < 0] <- 0
csd  <- sqrt(cvar); csd[csd == 0 | is.na(csd)] <- 1
message(" - per-cell background done (median cell sd = ", round(median(csd), 4), ")")

## --------------------------------- subset to markers; per-cell Zi / Zj --------------------------
mk_genes <- intersect(rownames(sm), markers$geneID)
message(" - markers present in matrix: ", length(mk_genes), " / ", nrow(markers))
smk <- as.matrix(sm[mk_genes, , drop = FALSE])    # markers x cells (small)

rmean <- rowMeans(smk)
rsd   <- sqrt((rowSums(smk^2) - n_cells * rmean^2) / (n_cells - 1)); rsd[rsd == 0 | is.na(rsd)] <- 1

Zi <- sweep(sweep(smk, 1, rmean, "-"), 1, rsd, "/")      # gene z across CELLS
Zj <- sweep(sweep(smk, 2, cmean, "-"), 2, csd, "/")      # cell z across GENES
Zi[Zi < 0] <- 0; Zj[Zj < 0] <- 0                          # max(0, .)
rZe <- sqrt(Zi^2 + Zj^2)                                  # EUCLID (formula as-is)
rZg <- sqrt(Zi *  Zj)                                     # GEOM   (AND)
rm(Zi, Zj); gc()

## ----------------------------------------------------------------- cluster order ----------------
cls <- unique(b$.cluster)
cls <- if (suppressWarnings(all(!is.na(as.numeric(cls))))) cls[order(as.numeric(cls))] else sort(cls)
cell_by_cluster <- split(rownames(b), b$.cluster)

## helper: cluster-mean (markers x clusters) + peak/spec for one per-cell rZ matrix ---------------
summarize <- function(rZ) {
  mrz <- sapply(cls, function(cl) rowMeans(rZ[, cell_by_cluster[[cl]], drop = FALSE]))
  if (is.null(dim(mrz))) mrz <- matrix(mrz, nrow = length(mk_genes), dimnames = list(mk_genes, cls))
  colnames(mrz) <- cls; rownames(mrz) <- mk_genes
  pk  <- apply(mrz, 1, which.max)
  pv  <- mrz[cbind(seq_len(nrow(mrz)), pk)]
  sec <- apply(mrz, 1, function(v) if (length(v) >= 2) sort(v, decreasing = TRUE)[2] else 0)
  list(mean = mrz, peak = colnames(mrz)[pk], peak_val = pv, spec = pv - sec)
}
Se <- summarize(rZe); Sg <- summarize(rZg)

## (A) cluster-mean rZ per marker, both metrics -------------------------------------------------
cm <- data.frame(stage = stage, geneID = mk_genes,
                 name = markers[mk_genes, "name"], type_label = markers[mk_genes, "type_label"],
                 peak_cl_euclid = Se$peak, peak_rZ_euclid = round(Se$peak_val, 4), spec_euclid = round(Se$spec, 4),
                 peak_cl_geom   = Sg$peak, peak_rZ_geom   = round(Sg$peak_val, 4), spec_geom   = round(Sg$spec, 4))
cm <- cm[order(cm$peak_cl_euclid, -cm$spec_euclid), ]
write.table(cm, op(".percell_rZ.cluster_mean.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (B) best markers per cluster (topN by specificity) for each metric ---------------------------
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
write.table(rbind(best_e, best_g), op(".best_markers_by_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (C) per-cell rZ (both metrics) for the UNION of selected best markers -> violin input ---------
bg <- union(best_e$geneID, best_g$geneID)
long <- data.frame(stage = stage,
                   geneID    = rep(bg, times = n_cells),
                   cell      = rep(colnames(rZe), each = length(bg)),
                   rZ_euclid = round(as.vector(rZe[bg, , drop = FALSE]), 4),
                   rZ_geom   = round(as.vector(rZg[bg, , drop = FALSE]), 4))
long$name       <- markers[long$geneID, "name"]
long$type_label <- markers[long$geneID, "type_label"]
long$cluster    <- b[long$cell, ".cluster"]
write.table(long, op(".percell_rZ.best.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (D) per-cell x TYPE mean rZ -- the CLUSTER-FREE annotation basis (cell-cycle excluded, >=3 markers/type)
tl_tab <- table(markers$type_label[markers$geneID %in% mk_genes])
tl_tab <- tl_tab[!grepl("dividing", names(tl_tab), ignore.case = TRUE)]
type_levels <- names(tl_tab)[tl_tab >= 3]
agg_types <- function(rZmat) {                       # -> cells x types (mean rZ over each type's markers)
  m <- sapply(type_levels, function(tl) {
    g <- intersect(markers$geneID[markers$type_label == tl], rownames(rZmat))
    if (!length(g)) rep(NA_real_, ncol(rZmat)) else colMeans(rZmat[g, , drop = FALSE])
  })
  rownames(m) <- colnames(rZmat); round(m, 4)
}
for (mn in c("geom", "euclid")) {
  pct <- agg_types(if (mn == "geom") rZg else rZe)
  write.table(data.frame(cellID = rownames(pct), pct, check.names = FALSE),
              op(paste0(".percell_type_rZ.", mn, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
message(" - per-cell x type matrix written (", length(type_levels), " types, dividing excluded)")

## ----------------------------------------------------------------- results md -------------------
rk1 <- function(best, tag) {
  b1 <- best[!duplicated(best$peak_cluster), ]
  c(paste0("## Best marker per cluster — ", tag, " (rank 1 by specificity)"),
    "| cluster | best | peak_rZ | specificity | type |", "|---|---|---|---|---|",
    apply(b1, 1, function(r) paste0("| ", r["peak_cluster"], " | ", r["name"], " | ", r["peak_rZ"],
                                    " | ", r["specificity"], " | ", r["type_label"], " |")), "")
}
md <- c(
  paste0("# 4_3i per-cell reciprocal z-score - ", prefix, " (", stage, ")"),
  "", "rZ[g,cell] from Zi (gene z across CELLS) and Zj (cell z across GENES).",
  paste0("- Smoothed matrix (perkb): `", smoothed_rds, "`"),
  paste0("- Markers: `", markers_bed, "`  (", length(mk_genes), " in matrix of ", nrow(markers), ")"),
  paste0("- Cells: ", n_cells, " | genes: ", n_genes, " | clusters: ", paste(cls, collapse = ", "),
         " | topN/cluster: ", topN), "",
  rk1(best_e, "EUCLID (formula as-is)"), rk1(best_g, "GEOM (AND)")
)
writeLines(md, op(".4_3i.results.md"))
message(" - wrote tables + results md -> ", out_dir)
message(" - DONE")
