#!/usr/bin/env Rscript
###############################################################################
## 4_3h_reciprocal_zscore.R
##
## Reciprocal (bidirectional) marker z-score per gene x cluster, on a PERKB
## (gene-length-normalized) pseudobulk matrix.
##
##   Zi[g,c] = t(scale(t(CPM)))  -- per-gene z ACROSS CLUSTERS  ("which cluster does g mark")
##                                  == the z that 4_3 already computes; length-IMMUNE (cancels)
##   Zj[g,c] = scale(CPM)        -- per-cluster z ACROSS GENES   ("is g high WITHIN c vs other genes")
##                                  length-SENSITIVE -> this is why the matrix MUST be perkb
##   rZ[g,c] = sqrt( max(0,Zi)^2 + max(0,Zj)^2 )   -- high only if g marks c AND stands out within c
##
## Pseudobulk + CPM are built identically to 4_3_marker_accessibility.R so Zi here == 4_3's z.
## Computed over ALL genes genome-wide, then subset to markers (Zj needs the full gene background).
##
## Usage:
##   Rscript 4_3h_reciprocal_zscore.R <out_dir> <prefix> <perkb_matrix_rds> <meta_txt> \
##           <markers_bed> <stage> [cluster_col=LouvainClusters] [topN=5] [cluster_annotation.tsv]
##
## Outputs (all under <out_dir>):
##   <prefix>.reciprocal_zscore.tsv     long: marker x cluster with Zi, Zj, rZ
##   <prefix>.marker_rZ_summary.tsv     per marker: peak cluster, max rZ/Zi/Zj (+concordance if ann given)
##   <prefix>.best_markers_by_type.tsv  topN markers per cell type by max rZ  (violin input)
##   <prefix>.cluster_top_markers.tsv   per cluster: its topN markers by rZ (what each cluster "is")
##   <prefix>.4_3h.results.md           manifest + headline
###############################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("Usage: Rscript 4_3h_reciprocal_zscore.R <out_dir> <prefix> <perkb_matrix_rds> <meta_txt> <markers_bed> <stage> [cluster_col] [topN] [cluster_annotation.tsv]")
out_dir     <- args[1]
prefix      <- args[2]
matrix_rds  <- args[3]
meta_txt    <- args[4]
markers_bed <- args[5]
stage       <- args[6]
cluster_col <- if (length(args) >= 7) args[7] else "LouvainClusters"
topN        <- if (length(args) >= 8) as.integer(args[8]) else 5L
annp        <- if (length(args) >= 9) args[9] else NA_character_

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(out_dir, paste0(prefix, suffix))
message(" - stage=", stage, " | prefix=", prefix, " | cluster_col=", cluster_col, " | topN=", topN)

## ----------------------------------------------------------------- (1) markers ------------------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")
markers$type_label <- paste0(markers$species, ":", markers$type)
message(" - markers (deduped): ", nrow(markers),
        "  (Zm=", sum(markers$species == "Zm"), " At=", sum(markers$species == "At"), ")")

## ----------------------------------------------------------------- (2) metadata -----------------
b <- read.table(meta_txt, header = TRUE, sep = "\t", quote = "", comment.char = "")
if ("cellID" %in% colnames(b)) rownames(b) <- b$cellID
if (!cluster_col %in% colnames(b))
  stop("cluster column '", cluster_col, "' not in metadata; have: ", paste(colnames(b), collapse = ", "))
b$.cluster <- as.character(b[[cluster_col]])
message(" - metadata cells: ", nrow(b))

## ----------------------------------------------------------------- (3) matrix + align -----------
m <- readRDS(matrix_rds)
cells <- intersect(colnames(m), rownames(b))
if (length(cells) == 0) stop("no overlap between matrix colnames and metadata cellIDs")
m <- m[, cells, drop = FALSE]
b <- b[cells, , drop = FALSE]
message(" - aligned cells (matrix ∩ meta): ", length(cells))

## ----------------------------------------------------------------- cluster order ----------------
cls <- unique(b$.cluster)
cls <- if (suppressWarnings(all(!is.na(as.numeric(cls))))) cls[order(as.numeric(cls))] else sort(cls)
cell_by_cluster <- split(rownames(b), b$.cluster)
message(" - clusters (", length(cls), "): ", paste(cls, collapse = ", "))

## ------------------------------------------- pseudobulk (all genes) -> CPM  (matches 4_3) --------
pb <- sapply(cls, function(cl) Matrix::rowSums(m[, cell_by_cluster[[cl]], drop = FALSE]))
colnames(pb) <- cls
libsize <- colSums(pb); libsize[libsize == 0] <- 1
cpm <- sweep(pb, 2, libsize, "/") * 1e6
cpm <- cpm[rowSums(cpm) > 0, , drop = FALSE]
message(" - genes in pseudobulk (nonzero): ", nrow(cpm))

## ------------------------------------------- Zi, Zj, rZ (all genes) ------------------------------
Zi <- t(scale(t(cpm)))   # per-gene z across clusters  (row z)   -- length-immune
Zj <- scale(cpm)         # per-cluster z across genes  (col z)   -- length-sensitive (needs perkb)
Zi[is.na(Zi)] <- 0
Zj[is.na(Zj)] <- 0
Zi <- as.matrix(Zi); Zj <- as.matrix(Zj)
## clamp negatives by subassignment (pmax(0, matrix) would drop the dim attribute)
Zi_pos <- Zi; Zi_pos[Zi_pos < 0] <- 0
Zj_pos <- Zj; Zj_pos[Zj_pos < 0] <- 0
rZ <- sqrt(Zi_pos^2 + Zj_pos^2)   # inherits dimnames from the matrices

## ------------------------------------------- subset to markers ----------------------------------
mk_genes <- intersect(rownames(rZ), markers$geneID)
message(" - markers present in matrix: ", length(mk_genes), " / ", nrow(markers))
Zi.m <- Zi[mk_genes, , drop = FALSE]; Zj.m <- Zj[mk_genes, , drop = FALSE]; rZ.m <- rZ[mk_genes, , drop = FALSE]

## (A) long table: marker x cluster with all three scores --------------------------------------
long <- data.frame(
  stage      = stage,
  geneID     = rep(mk_genes, times = ncol(rZ.m)),
  cluster    = rep(colnames(rZ.m), each = length(mk_genes)),
  Zi         = round(as.vector(Zi.m), 4),
  Zj         = round(as.vector(Zj.m), 4),
  rZ         = round(as.vector(rZ.m), 4))
long$name       <- markers[long$geneID, "name"]
long$species    <- markers[long$geneID, "species"]
long$type       <- markers[long$geneID, "type"]
long$type_label <- markers[long$geneID, "type_label"]
long <- long[, c("stage","geneID","name","species","type","type_label","cluster","Zi","Zj","rZ")]
write.table(long, op(".reciprocal_zscore.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (B) per-marker summary: peak cluster (argmax rZ), max scores ---------------------------------
peak_idx <- apply(rZ.m, 1, which.max)
summ <- data.frame(
  stage       = stage,
  geneID      = mk_genes,
  name        = markers[mk_genes, "name"],
  species     = markers[mk_genes, "species"],
  type_label  = markers[mk_genes, "type_label"],
  peak_cluster= colnames(rZ.m)[peak_idx],
  max_rZ      = round(rZ.m[cbind(seq_along(mk_genes), peak_idx)], 4),
  Zi_at_peak  = round(Zi.m[cbind(seq_along(mk_genes), peak_idx)], 4),
  Zj_at_peak  = round(Zj.m[cbind(seq_along(mk_genes), peak_idx)], 4))
## optional concordance vs a cluster annotation (peak cluster annotated as the marker's own type?)
if (!is.na(annp) && file.exists(annp)) {
  ann <- read.table(annp, header = TRUE, sep = "\t", quote = "", comment.char = "")
  top_by_cl <- setNames(ann$top_type, as.character(ann$cluster))
  summ$peak_cluster_top_type <- top_by_cl[as.character(summ$peak_cluster)]
  summ$peak_concordant <- !is.na(summ$peak_cluster_top_type) & summ$peak_cluster_top_type == summ$type_label
}
summ <- summ[order(summ$type_label, -summ$max_rZ), ]
write.table(summ, op(".marker_rZ_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (C) best markers per cell type (topN by max_rZ) -- violin input ------------------------------
best <- do.call(rbind, lapply(split(summ, summ$type_label), function(d) head(d[order(-d$max_rZ), ], topN)))
best <- best[order(best$type_label, -best$max_rZ), ]
write.table(best, op(".best_markers_by_type.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (D) per cluster: its topN markers by rZ (what does each cluster "look like"?) -----------------
cltop <- do.call(rbind, lapply(cls, function(cl) {
  v <- sort(rZ.m[, cl], decreasing = TRUE)
  g <- names(v)[seq_len(min(topN, length(v)))]
  data.frame(stage = stage, cluster = cl, rank = seq_along(g), geneID = g,
             name = markers[g, "name"], type_label = markers[g, "type_label"],
             rZ = round(v[g], 4), Zi = round(Zi.m[g, cl], 4), Zj = round(Zj.m[g, cl], 4))
}))
write.table(cltop, op(".cluster_top_markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ----------------------------------------------------------------- results md -------------------
conc <- if ("peak_concordant" %in% colnames(summ))
  sprintf("%d/%d (%.0f%%)", sum(summ$peak_concordant, na.rm = TRUE), nrow(summ),
          100 * mean(summ$peak_concordant, na.rm = TRUE)) else "NA (no annotation given)"
md <- c(
  paste0("# 4_3h reciprocal z-score - ", prefix, " (", stage, ")"),
  "",
  "rZ[g,c] = sqrt( max(0,Zi)^2 + max(0,Zj)^2 );  Zi = per-gene z across clusters, Zj = per-cluster z across genes.",
  paste0("- Matrix (perkb): `", matrix_rds, "`"),
  paste0("- Metadata: `", meta_txt, "`  (cluster `", cluster_col, "`)"),
  paste0("- Markers: `", markers_bed, "`  (", length(mk_genes), " in matrix of ", nrow(markers), ")"),
  paste0("- Cells: ", length(cells), " | clusters: ", paste(cls, collapse = ", ")),
  paste0("- Peak-cluster concordance (marker peaks in a cluster annotated as its own type): ", conc),
  "",
  "## Best marker per cell type (rank 1 by max rZ)",
  "| type | best | max_rZ | Zi@peak | Zj@peak | peak_cl |",
  "|---|---|---|---|---|---|",
  {
    b1 <- best[!duplicated(best$type_label), ]
    apply(b1, 1, function(r) paste0("| ", r["type_label"], " | ", r["name"], " | ", r["max_rZ"],
                                    " | ", r["Zi_at_peak"], " | ", r["Zj_at_peak"], " | ", r["peak_cluster"], " |"))
  }
)
writeLines(md, op(".4_3h.results.md"))
message(" - wrote tables + results md -> ", out_dir)
message(" - DONE")
