###################################################################################################
## 4_3b_marker_smoothing.R
##
## Markov-affinity graph smoothing (imputation) of MARKER gene activity for the marker-on-UMAP
## panels of the figure. scATAC gene activity is ~0/1 per cell, so a per-cell UMAP painted with raw
## marker activity is pure noise; diffusing each marker over each cell's neighbourhood in PC space
## reveals the spatial (cell-type) pattern. This is the deferred companion to 4_3 (which does the
## cluster-level tables and needs NO smoothing -- pseudobulk aggregation is its own denoiser).
##
## Adapted from the maize_282 reference annotation code (functions.plot_marker_accessibility.R::smooth.data),
## but lighter and memory-safe:
##   * smooth ONLY the marker genes (the reference already uses smooth.markers=T),
##   * use the EXISTING reduced dimensions (PC_3-20) as the manifold -- no LSI recompute,
##   * skip the fragile per-cluster Mclust normalization (simple per-cell CP10k instead),
##   * DIFFUSE THE DATA ITERATIVELY (A %*% X, `step` times) instead of forming A^(2^step):
##     the cell graph stays sparse and we never build a dense cells x cells power matrix.
##
## Usage:
##   Rscript 4_3b_marker_smoothing.R <out_dir> <prefix> <matrix_rds> <meta_txt> <rd_txt> <markers_bed> <stage> [k] [step] [max_panels]
##
## Outputs (under <out_dir>; consumed by the figure scripts for the marker-on-UMAP panels):
##   <prefix>.marker_smoothed.rds          smoothed marker x cell matrix (genes x cells)
##   <prefix>.marker_smoothed.manifest.tsv  geneID/name/species/type/type_label for the RDS rows
##   <prefix>.cells_umap.tsv                cellID, umap1, umap2, cluster, species, stage (join key)
##   plots/<prefix>.marker_umaps.png        EXPLORATORY grid (bounded) -- not the final figure
###################################################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript 4_3b_marker_smoothing.R <out_dir> <prefix> <matrix_rds> <meta_txt> <rd_txt> <markers_bed> <stage> [k] [step] [max_panels]")
}
out_dir     <- args[1]
prefix      <- args[2]
matrix_rds  <- args[3]
meta_txt    <- args[4]
rd_txt      <- args[5]
markers_bed <- args[6]
stage       <- args[7]
k           <- if (length(args) >= 8)  as.integer(args[8])  else 25L
step        <- if (length(args) >= 9)  as.integer(args[9])  else 3L
max_panels  <- if (length(args) >= 10) as.integer(args[10]) else 48L
cluster_col <- "LouvainClusters"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(out_dir, paste0(prefix, suffix))
message(" - stage=", stage, " | prefix=", prefix, " | k=", k, " step=", step)

## KNN helper: returns a k-column index matrix INCLUDING self as column 1 (RANN convention).
## Uses RANN/FNN (kd-tree, fast) if available; otherwise a chunked base-R exact KNN so the script
## never hard-fails on a missing package (slower, but memory-bounded for ~20-30k cells).
get_knn_idx <- function(X, k) {
  if (requireNamespace("RANN", quietly = TRUE)) return(RANN::nn2(X, k = k)$nn.idx)
  if (requireNamespace("FNN",  quietly = TRUE)) return(cbind(seq_len(nrow(X)), FNN::get.knn(X, k = k - 1)$nn.index))
  message("   * RANN/FNN absent -> base-R exact KNN fallback (slower)")
  n <- nrow(X); sq <- rowSums(X^2); idx <- matrix(0L, n, k); chunk <- 1024L
  for (s in seq(1L, n, by = chunk)) {
    e  <- min(s + chunk - 1L, n)
    d2 <- outer(sq[s:e], sq, "+") - 2 * (X[s:e, , drop = FALSE] %*% t(X))   # (chunk x n) squared dist
    for (r in seq_len(nrow(d2))) idx[s + r - 1L, ] <- order(d2[r, ])[seq_len(k)]
  }
  idx
}

## ----------------------------------------------------------------- inputs -----------------------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")
markers$type_label <- paste0(markers$species, ":", markers$type)

b <- read.table(meta_txt, header = TRUE, sep = "\t", quote = "", comment.char = "")
if ("cellID" %in% colnames(b)) rownames(b) <- b$cellID
if (!"species" %in% colnames(b)) b$species <- ifelse(grepl("_At_", rownames(b)), "At", "B73")
if (!cluster_col %in% colnames(b)) b[[cluster_col]] <- NA

pcs <- read.table(rd_txt, header = TRUE, sep = "\t", quote = "", comment.char = "", row.names = 1)
m   <- readRDS(matrix_rds)

## ----------------------------------------------------------------- align cells/genes ------------
cells <- Reduce(intersect, list(colnames(m), rownames(b), rownames(pcs)))
if (length(cells) == 0) stop("no common cells across matrix / metadata / reduced_dims")
b   <- b[cells, , drop = FALSE]
pcs <- as.matrix(pcs[cells, , drop = FALSE])
mk_genes <- intersect(markers$geneID, rownames(m))
mk_mat   <- m[mk_genes, cells, drop = FALSE]
message(" - cells (matrix ∩ meta ∩ rd): ", length(cells), " | marker genes: ", length(mk_genes),
        " | PC dims: ", ncol(pcs))

## ----------------------------------------------------------------- per-cell CP10k norm ----------
libcell <- Matrix::colSums(m); libcell[libcell == 0] <- 1
mk_norm <- mk_mat %*% Matrix::Diagonal(x = 1 / libcell[cells]) * 1e4   # genes x cells, lib-normalized
mk_norm <- as.matrix(mk_norm)
rownames(mk_norm) <- mk_genes; colnames(mk_norm) <- cells

## ----------------------------------------------------------------- build Markov graph -----------
message(" - building KNN graph (k=", k, ") ...")
nn <- get_knn_idx(pcs, k)
n  <- nrow(pcs)
A  <- sparseMatrix(i = rep(seq_len(n), each = ncol(nn)), j = as.vector(t(nn)), x = 1, dims = c(n, n))
A  <- A + Matrix::t(A)                                  # symmetrize
A  <- Matrix::Diagonal(x = 1 / Matrix::rowSums(A)) %*% A # row-normalize -> Markov transition

## ----------------------------------------------------------------- iterative diffusion ----------
message(" - diffusing ", length(mk_genes), " markers x ", n, " cells, step=", step, " ...")
X <- t(mk_norm)                                         # cells x genes (dense, small)
for (s in seq_len(step)) X <- as.matrix(A %*% X)        # one hop per iteration; never form A^step
smoothed <- t(X)                                        # genes x cells
rownames(smoothed) <- mk_genes; colnames(smoothed) <- cells

## ----------------------------------------------------------------- save -------------------------
saveRDS(smoothed, op(".marker_smoothed.rds"))
man <- markers[mk_genes, c("geneID", "name", "species", "type", "type_label")]
man$stage <- stage
write.table(man, op(".marker_smoothed.manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

coords <- data.frame(cellID = cells,
                     umap1 = if ("umap1" %in% colnames(b)) b$umap1 else NA,
                     umap2 = if ("umap2" %in% colnames(b)) b$umap2 else NA,
                     cluster = b[[cluster_col]], species = b$species, stage = stage)
write.table(coords, op(".cells_umap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote smoothed RDS (", nrow(smoothed), " x ", ncol(smoothed), "), manifest, coords")

## ----------------------------------------------- exploratory marker-on-UMAP grid (bounded) ------
tryCatch({
  if (!all(c("umap1", "umap2") %in% colnames(b))) stop("no umap1/umap2 in metadata")
  ## pick a bounded, type-spread set: the highest-mean-smoothed marker per type_label, capped
  ord_types <- unique(man$type_label)
  pick <- unlist(lapply(ord_types, function(tl) {
    g <- man$geneID[man$type_label == tl]
    g[which.max(Matrix::rowMeans(smoothed[g, , drop = FALSE]))]
  }))
  if (length(pick) > max_panels) {
    message(" - capping exploratory panels at ", max_panels, " of ", length(pick),
            " type-representative markers (RDS still holds all ", length(mk_genes), ")")
    pick <- pick[seq_len(max_panels)]
  }
  ncol_g <- 6; nrow_g <- ceiling(length(pick) / ncol_g)
  png(file.path(out_dir, "plots", paste0(prefix, ".marker_umaps.png")),
      width = ncol_g * 2, height = nrow_g * 2, units = "in", res = 200, type = "cairo")
  par(mfrow = c(nrow_g, ncol_g), mar = c(1, 1, 2, 1))
  pal <- colorRampPalette(c("grey85", "grey80", "goldenrod2", "firebrick3"))(100)
  uo <- order(b$umap1)  # static
  for (g in pick) {
    v <- smoothed[g, ]; v[is.na(v)] <- 0
    hi <- as.numeric(quantile(v, 0.99)); if (hi <= 0) hi <- max(v) + 1e-9
    v[v > hi] <- hi
    o <- order(v)                              # plot high-activity cells on top
    cidx <- if (hi > 0) pal[cut(v[o], breaks = seq(0, hi, length.out = 101), include.lowest = TRUE)] else pal[1]
    plot(b$umap1[o], b$umap2[o], col = cidx, pch = 16, cex = 0.15,
         main = paste0(man[g, "name"], "\n", man[g, "type_label"]),
         xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", cex.main = 0.6)
  }
  dev.off()
  message(" - exploratory grid: ", length(pick), " panels -> plots/", prefix, ".marker_umaps.png")
}, error = function(e) message(" ! exploratory marker UMAP grid skipped: ", conditionMessage(e)))

message(" - DONE")
