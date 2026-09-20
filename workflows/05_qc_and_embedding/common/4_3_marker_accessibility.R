###################################################################################################
## 4_3_marker_accessibility.R
##
## Pseudobulk marker-accessibility annotation for the SM2v2 minc_50 clusters (run once per stage,
## PRE and POST). This script's job is to emit the *tables* that feed the publication figure (built
## later in analysis/), plus quick EXPLORATORY plots for our own eyes here.
##
## Adapted -- dependency-light -- from the maize_282 reference
## (reference annotation code: PlotClusterZscore.R, plot_marker_accessibility.R):
##   - pseudobulk per cluster -> CPM -> per-gene z-score across clusters     (PlotClusterZscore)
##   - per-cluster mean activity + proportion of cells accessible            (plot.new.markers)
## We deliberately DROP the per-cluster 2-Gaussian Mclust thresholding + Markov-affinity smoothing
## of the reference: those are fragile on small clusters (SM2v2's At island is ~500-900 cells) and
## drag in heavy/non-CRAN deps (varistran, mclust, RANN, ...). Pseudobulk z-scores are deterministic
## and robust on SM2v2's ~12 clusters. The smoothed per-marker UMAP panels are deferred to 4_3b.
##
## Two species share one matrix + one panel: genes split by ID prefix (Zm00001eb = maize B73v5,
## AT.G.... = Arabidopsis TAIR); cell types are species-tagged ("Zm:epidermis" vs "At:epidermis")
## because the two atlases reuse type names. CellIDs join matrix<->metadata directly (no remap).
##
## Usage:
##   Rscript 4_3_marker_accessibility.R <out_dir> <prefix> <matrix_rds> <meta_txt> <markers_bed> <stage> [cluster_col] [min_markers]
##     min_markers : a cell type needs >= this many markers-in-matrix to be an eligible top call
##                   (default 3; guards against 1-2-marker types spiking the mean-z call on an
##                   imbalanced/uncapped panel -- see the informative-pruned panel from 4_3e/4_3d)
##
## Outputs (all under <out_dir>; tidy + stage-tagged so the figure scripts can join PRE+POST):
##   <prefix>.marker_overlap_qc.tsv    panel vs matrix overlap per species             (table E / QC)
##   <prefix>.marker_zscore.tsv        per-marker pseudobulk z-score x cluster          (table C)
##   <prefix>.marker_dotplot.tsv       per-marker mean activity + % cells x cluster     (table A)
##   <prefix>.celltype_score.tsv       cluster x cell-type mean z-score (heatmap input) (table B)
##   <prefix>.cluster_annotation.tsv   per-cluster top cell-type call + confidence      (table D)
##   <prefix>.4_3.results.md           manifest + headline numbers (read by the figure scripts)
##   plots/<prefix>.celltype_heatmap.pdf, plots/<prefix>.marker_dotplot.pdf   (EXPLORATORY only)
###################################################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript 4_3_marker_accessibility.R <out_dir> <prefix> <matrix_rds> <meta_txt> <markers_bed> <stage> [cluster_col]")
}
out_dir     <- args[1]
prefix      <- args[2]
matrix_rds  <- args[3]
meta_txt    <- args[4]
markers_bed <- args[5]
stage       <- args[6]
cluster_col <- if (length(args) >= 7) args[7] else "LouvainClusters"
min_markers <- if (length(args) >= 8) as.integer(args[8]) else 3L   # min markers-in-matrix for a callable top type

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(out_dir, paste0(prefix, suffix))   # output-path helper

message(" - stage=", stage, " | prefix=", prefix, " | cluster_col=", cluster_col)

## ----------------------------------------------------------------- (1) markers ------------------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")   # maize v5 vs At TAIR
markers$type_label <- paste0(markers$species, ":", markers$type)         # avoid cross-species name clash
message(" - markers (deduped by geneID): ", nrow(markers),
        "  (Zm=", sum(markers$species == "Zm"), " At=", sum(markers$species == "At"), ")")

## ----------------------------------------------------------------- (2) metadata -----------------
b <- read.table(meta_txt, header = TRUE, sep = "\t", quote = "", comment.char = "")
if ("cellID" %in% colnames(b)) rownames(b) <- b$cellID
if (!cluster_col %in% colnames(b)) {
  stop("cluster column '", cluster_col, "' not in metadata; have: ", paste(colnames(b), collapse = ", "))
}
b$.cluster <- as.character(b[[cluster_col]])
if (!"species" %in% colnames(b)) b$species <- ifelse(grepl("_At_", rownames(b)), "At", "B73")
message(" - metadata cells: ", nrow(b))

## ----------------------------------------------------------------- (3) matrix + align -----------
m <- readRDS(matrix_rds)
cells <- intersect(colnames(m), rownames(b))
if (length(cells) == 0) stop("no overlap between matrix colnames and metadata cellIDs (check naming)")
m <- m[, cells, drop = FALSE]
b <- b[cells, , drop = FALSE]
message(" - aligned cells (matrix ∩ meta): ", length(cells), " / meta ", nrow(b))

## ----------------------------------------------------------------- (E) overlap QC ---------------
genes_mat <- rownames(m)
qc <- do.call(rbind, lapply(c("Zm", "At"), function(sp) {
  pin   <- markers$geneID[markers$species == sp]
  inmat <- intersect(pin, genes_mat)
  data.frame(stage = stage, species = sp, n_panel = length(pin), n_in_matrix = length(inmat),
             frac = if (length(pin) > 0) round(length(inmat) / length(pin), 4) else NA)
}))
write.table(qc, op(".marker_overlap_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - marker overlap: Zm ", qc$n_in_matrix[qc$species == "Zm"], "/", qc$n_panel[qc$species == "Zm"],
        " | At ", qc$n_in_matrix[qc$species == "At"], "/", qc$n_panel[qc$species == "At"])

mk       <- markers[markers$geneID %in% genes_mat, ]   # markers present in the matrix
mk_genes <- mk$geneID

## ----------------------------------------------------------------- cluster order ----------------
cls <- unique(b$.cluster)
cls <- if (suppressWarnings(all(!is.na(as.numeric(cls))))) cls[order(as.numeric(cls))] else sort(cls)
cell_by_cluster <- split(rownames(b), b$.cluster)
message(" - clusters (", length(cls), "): ", paste(cls, collapse = ", "))

## ------------------------------------------- pseudobulk (all genes) -> CPM -> per-gene z ---------
pb <- sapply(cls, function(cl) Matrix::rowSums(m[, cell_by_cluster[[cl]], drop = FALSE]))
colnames(pb) <- cls                                    # genes x clusters (summed counts)
libsize <- colSums(pb); libsize[libsize == 0] <- 1     # guard empty clusters (avoid Inf in scale())
cpm <- sweep(pb, 2, libsize, "/") * 1e6
cpm <- cpm[rowSums(cpm) > 0, , drop = FALSE]
zall <- t(scale(t(cpm)))                               # per-gene z across clusters
zall[is.na(zall)] <- 0

## ----------------------------------------------------------------- (C) per-marker z -------------
zmk <- zall[rownames(zall) %in% mk_genes, , drop = FALSE]
zmk_df <- data.frame(stage   = stage,
                     geneID  = rep(rownames(zmk), times = ncol(zmk)),
                     cluster = rep(colnames(zmk), each = nrow(zmk)),
                     zscore  = as.vector(zmk))          # as.vector is column-major -> matches reps
zmk_df$name       <- mk[zmk_df$geneID, "name"]
zmk_df$species    <- mk[zmk_df$geneID, "species"]
zmk_df$type       <- mk[zmk_df$geneID, "type"]
zmk_df$type_label <- mk[zmk_df$geneID, "type_label"]
zmk_df <- zmk_df[, c("stage", "geneID", "name", "species", "type", "type_label", "cluster", "zscore")]
write.table(zmk_df, op(".marker_zscore.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ----------------------------------------------------------------- (A) dotplot data -------------
libcell <- Matrix::colSums(m); libcell[libcell == 0] <- 1
mk_mat  <- m[rownames(m) %in% mk_genes, , drop = FALSE]
norm_cell <- mk_mat %*% Matrix::Diagonal(x = 1 / libcell) * 1e4   # CP10k per cell
colnames(norm_cell) <- colnames(mk_mat)
dot <- do.call(rbind, lapply(cls, function(cl) {
  ci <- cell_by_cluster[[cl]]
  data.frame(stage = stage, cluster = cl, geneID = rownames(mk_mat),
             mean_activity   = as.numeric(Matrix::rowMeans(norm_cell[, ci, drop = FALSE])),
             frac_accessible = as.numeric(Matrix::rowMeans(mk_mat[, ci, drop = FALSE] > 0)),
             n_cells = length(ci))
}))
dot$name       <- mk[dot$geneID, "name"]
dot$species    <- mk[dot$geneID, "species"]
dot$type       <- mk[dot$geneID, "type"]
dot$type_label <- mk[dot$geneID, "type_label"]
dot <- dot[, c("stage", "cluster", "geneID", "name", "species", "type", "type_label",
               "mean_activity", "frac_accessible", "n_cells")]
write.table(dot, op(".marker_dotplot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ----------------------------------------------- (B) cluster x cell-type mean z-score -----------
mk$.in_z    <- mk$geneID %in% rownames(zmk)
type_levels <- sort(unique(mk$type_label[mk$.in_z]))
ct <- sapply(cls, function(cl) {
  sapply(type_levels, function(tl) {
    g <- mk$geneID[mk$type_label == tl & mk$.in_z]
    if (length(g) == 0) NA_real_ else mean(zmk[g, cl], na.rm = TRUE)
  })
})                                                     # type_label x cluster
if (is.null(dim(ct))) ct <- matrix(ct, nrow = length(type_levels),
                                   dimnames = list(type_levels, cls))
nmk <- sapply(type_levels, function(tl) sum(mk$type_label == tl & mk$.in_z))
ct_df <- data.frame(stage = stage,
                    cluster    = rep(colnames(ct), each = nrow(ct)),
                    type_label = rep(rownames(ct), times = ncol(ct)),
                    mean_zscore = as.vector(ct))
ct_df$species        <- sub(":.*", "", ct_df$type_label)
ct_df$n_markers_used <- nmk[ct_df$type_label]
ct_df <- ct_df[, c("stage", "cluster", "type_label", "species", "mean_zscore", "n_markers_used")]
write.table(ct_df, op(".celltype_score.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ----------------------------------------------------------------- (D) annotation calls ---------
ann <- do.call(rbind, lapply(cls, function(cl) {
  v <- ct[, cl]; v <- v[!is.na(v)]
  v <- v[nmk[names(v)] >= min_markers]           # exclude under-powered types (thin-type mean-z spikes)
  if (length(v) == 0) return(NULL)
  o <- order(v, decreasing = TRUE)
  top_s <- v[o[1]]
  sec   <- if (length(v) >= 2) names(v)[o[2]] else NA
  sec_s <- if (length(v) >= 2) v[o[2]] else NA_real_
  p <- pmax(v, 0)
  ent <- if (sum(p) > 0) { pn <- p[p > 0] / sum(p); -sum(pn * log2(pn)) } else NA_real_
  ci <- cell_by_cluster[[cl]]
  sp_tab <- sort(table(b[ci, "species"]), decreasing = TRUE)
  data.frame(stage = stage, cluster = cl,
             top_type = names(v)[o[1]], top_score = round(top_s, 3),
             second_type = sec, second_score = round(as.numeric(sec_s), 3),
             margin = round(top_s - as.numeric(sec_s), 3),
             shannon_entropy = round(ent, 3),
             n_cells = length(ci), dominant_species = names(sp_tab)[1])
}))
write.table(ann, op(".cluster_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - annotation calls:")
print(ann[, c("cluster", "top_type", "top_score", "margin", "shannon_entropy", "n_cells", "dominant_species")])

## ----------------------------------------------------- exploratory plots (never block tables) ---
tryCatch({
  suppressMessages(library(pheatmap))
  hm <- ct; hm[is.na(hm)] <- 0
  pheatmap(hm, cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 5, fontsize_col = 8,
           color = colorRampPalette(c("dodgerblue4", "white", "firebrick4"))(100),
           main = paste0(prefix, " (", stage, ") cell-type z-score"),
           filename = file.path(out_dir, "plots", paste0(prefix, ".celltype_heatmap.pdf")),
           width = 8, height = 12)
}, error = function(e) message(" ! exploratory heatmap skipped: ", conditionMessage(e)))

tryCatch({
  suppressMessages(library(ggplot2))
  topg <- unique(unlist(lapply(cls, function(cl) {
    s <- sort(zmk[, cl], decreasing = TRUE); names(s)[seq_len(min(3, length(s)))]
  })))
  gg <- dot[dot$geneID %in% topg, ]
  gg$label <- paste0(gg$name, " (", gg$type_label, ")")
  p <- ggplot(gg, aes(x = factor(cluster, levels = cls), y = label)) +
    geom_point(aes(size = frac_accessible, color = mean_activity)) +
    scale_color_viridis_c() + theme_bw() +
    labs(x = "cluster", y = NULL, title = paste0(prefix, " (", stage, ") top-marker dotplot")) +
    theme(axis.text.y = element_text(size = 5))
  ggsave(file.path(out_dir, "plots", paste0(prefix, ".marker_dotplot.pdf")), p, width = 9, height = 12)
}, error = function(e) message(" ! exploratory dotplot skipped: ", conditionMessage(e)))

## ----------------------------------------------------------------- results md -------------------
md <- c(
  paste0("# 4_3 marker-accessibility results - ", prefix, " (", stage, ")"),
  "",
  "Generated by `0_scripts/common/4_3_marker_accessibility.R` (pseudobulk track).",
  paste0("- Matrix: `", matrix_rds, "`"),
  paste0("- Metadata: `", meta_txt, "`  (cluster column `", cluster_col, "`)"),
  paste0("- Markers: `", markers_bed, "`"),
  paste0("- Cells aligned (matrix ∩ meta): ", length(cells)),
  paste0("- Clusters (", length(cls), "): ", paste(cls, collapse = ", ")),
  paste0("- Marker overlap - Zm: ", qc$n_in_matrix[qc$species == "Zm"], "/", qc$n_panel[qc$species == "Zm"],
         ", At: ", qc$n_in_matrix[qc$species == "At"], "/", qc$n_panel[qc$species == "At"]),
  "",
  "## Tables (tidy, stage-tagged)",
  "| file | contents |",
  "|------|----------|",
  paste0("| `", prefix, ".marker_overlap_qc.tsv` | panel ∩ matrix per species |"),
  paste0("| `", prefix, ".marker_zscore.tsv` | per-marker pseudobulk z x cluster |"),
  paste0("| `", prefix, ".marker_dotplot.tsv` | mean activity + %cells x cluster |"),
  paste0("| `", prefix, ".celltype_score.tsv` | cluster x cell-type mean z (heatmap input) |"),
  paste0("| `", prefix, ".cluster_annotation.tsv` | per-cluster top cell-type call |"),
  "",
  "## Cluster annotation (top call)",
  "| cluster | top_type | top_score | margin | entropy | n_cells | dom_species |",
  "|---|---|---|---|---|---|---|",
  apply(ann, 1, function(r) paste0("| ", r["cluster"], " | ", r["top_type"], " | ", r["top_score"],
                                   " | ", r["margin"], " | ", r["shannon_entropy"], " | ",
                                   r["n_cells"], " | ", r["dominant_species"], " |")),
  "",
  "_Plots in `plots/` are EXPLORATORY, not publication figures - the final figure is built in `analysis/` from these tables._"
)
writeLines(md, op(".4_3.results.md"))
message(" - wrote ", length(list.files(out_dir, pattern = "\\.tsv$")), " tables + results md -> ", out_dir)
message(" - DONE")
