###################################################################################################
## 4_3g_cell_annotation.R   --   STEP 3: per-cell cell-type annotation on SMOOTHED gene activity
##
## Consumes the 4_3f smoothed gene-activity matrix and implements the Marand et al. 2021
## (Cell 184:3041) enrichment classifier (approaches 1-2; NOT the glmnet logistic of approach 3):
##   per cell, per cell type:
##     z = ( mean smoothed activity of the type's markers  -  background mean )
##         / ( background sd / sqrt(n_markers) )
##   background = a random sample of non-marker genes drawn from the smoothed matrix. This z is the
##   closed form of the paper's "1,000 random gene sets of the same size" permutation (the mean of n
##   random genes has sd = pop_sd / sqrt(n)). Then SCALE 0-1 per cell by the max across types (paper).
##   CALL per cell: top type if z_top >= z_thresh AND z_top >= ratio * z_second, else "unknown".
##   Per cluster: majority type if one type > 50%; >=2 majorities -> "mixed"; else "unknown".
##
## Re-runnable at will to tune z_thresh / ratio WITHOUT re-smoothing (that's the point of the split).
##
## Usage:
##   Rscript 4_3g_cell_annotation.R <out_dir> <prefix> <smoothed_rds> <meta_txt> <markers_bed> <stage> \
##           [z_thresh=2] [ratio=1.5] [min_markers=3] [n_bg=2000] [seed=1] [cluster_col=LouvainClusters]
## Outputs (prefix-tagged): .cell_annotation.tsv .cluster_majority.tsv .celltype_meanZ_by_cluster.tsv
##   plots/<prefix>.annotation_umap.png
###################################################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("Usage: Rscript 4_3g_cell_annotation.R <out_dir> <prefix> <smoothed_rds> <meta_txt> <markers_bed> <stage> [z_thresh ratio min_markers n_bg seed cluster_col]")
out_dir <- args[1]; prefix <- args[2]; smoothed_rds <- args[3]; meta_txt <- args[4]
markers_bed <- args[5]; stage <- args[6]
z_thresh    <- if (length(args) >= 7)  as.numeric(args[7])  else 2.0
ratio       <- if (length(args) >= 8)  as.numeric(args[8])  else 1.5
min_markers <- if (length(args) >= 9)  as.integer(args[9])  else 3L
n_bg        <- if (length(args) >= 10) as.integer(args[10]) else 2000L
seed        <- if (length(args) >= 11) as.integer(args[11]) else 1L
cluster_col <- if (length(args) >= 12) args[12] else "LouvainClusters"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
op <- function(s) file.path(out_dir, paste0(prefix, s))
set.seed(seed)
message(" - 4_3g annotate | prefix=", prefix, " stage=", stage,
        " | z_thresh=", z_thresh, " ratio=", ratio, " min_markers=", min_markers, " n_bg=", n_bg)

## ---- inputs -------------------------------------------------------------------------------------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]; rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")
markers$type_label <- paste0(markers$species, ":", markers$type)

b <- read.table(meta_txt, header = TRUE, sep = "\t", quote = "", comment.char = "")
if ("cellID" %in% colnames(b)) rownames(b) <- b$cellID
if (!cluster_col %in% colnames(b)) stop("cluster column '", cluster_col, "' not in metadata")
b$.cluster <- as.character(b[[cluster_col]])

sm <- readRDS(smoothed_rds)                          # genes x cells (dense, smoothed CP10k)
cells <- intersect(colnames(sm), rownames(b))
if (length(cells) == 0) stop("no common cells between smoothed matrix and metadata")
sm <- sm[, cells, drop = FALSE]; b <- b[cells, , drop = FALSE]
message(" - cells (smoothed ∩ meta): ", length(cells), " | genes in smoothed: ", nrow(sm))

## ---- marker + background gene sets --------------------------------------------------------------
mk_in <- intersect(markers$geneID, rownames(sm))
if (length(mk_in) == 0) stop("no marker genes in the smoothed matrix (check gene IDs)")
bg_pool  <- setdiff(rownames(sm), mk_in)
bg_genes <- if (length(bg_pool) > n_bg) sample(bg_pool, n_bg) else bg_pool
mk <- markers[mk_in, ]
types <- sort(unique(mk$type_label)); n_t <- sapply(types, function(tl) sum(mk$type_label == tl))
types <- types[n_t >= min_markers]; n_t <- n_t[types]
if (length(types) == 0) stop("no cell type has >= ", min_markers, " markers present")
message(" - markers in matrix: ", length(mk_in), " | background genes: ", length(bg_genes),
        " | callable types (>=", min_markers, " markers): ", length(types))

## ---- per-cell background + per-type enrichment z -----------------------------------------------
bg <- sm[bg_genes, , drop = FALSE]
bg_mean <- colMeans(bg); bg_sd <- apply(bg, 2, sd); bg_sd[bg_sd == 0 | is.na(bg_sd)] <- .Machine$double.eps
tmean <- t(sapply(types, function(tl) colMeans(sm[mk$geneID[mk$type_label == tl], , drop = FALSE])))  # types x cells
Z <- (tmean - matrix(bg_mean, length(types), length(cells), byrow = TRUE)) /
     (matrix(bg_sd,   length(types), length(cells), byrow = TRUE) / sqrt(n_t))
rownames(Z) <- types; colnames(Z) <- cells

## ---- per-cell call ------------------------------------------------------------------------------
ord      <- apply(Z, 2, order, decreasing = TRUE)
top_ty   <- rownames(Z)[ord[1, ]]; sec_ty <- rownames(Z)[ord[2, ]]
top_z    <- Z[cbind(ord[1, ], seq_len(ncol(Z)))]; sec_z <- Z[cbind(ord[2, ], seq_len(ncol(Z)))]
zmax     <- pmax(apply(Z, 2, max), .Machine$double.eps)
top_scaled <- top_z / zmax
call <- ifelse(top_z >= z_thresh & (sec_z <= 0 | top_z >= ratio * sec_z), top_ty, "unknown")

cell_ann <- data.frame(cellID = cells, cluster = b$.cluster, stage = stage,
                       top_type = top_ty, top_z = round(top_z, 3), top_scaled = round(top_scaled, 3),
                       second_type = sec_ty, second_z = round(sec_z, 3), call = call)
zt <- as.data.frame(t(round(Z, 3))); colnames(zt) <- paste0("z.", gsub("[^A-Za-z0-9]+", "_", rownames(Z)))
cell_ann <- cbind(cell_ann, zt)
write.table(cell_ann, op(".cell_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- per-cluster majority -----------------------------------------------------------------------
cl_levels <- unique(b$.cluster)
cl_levels <- if (suppressWarnings(all(!is.na(as.numeric(cl_levels))))) cl_levels[order(as.numeric(cl_levels))] else sort(cl_levels)
maj <- do.call(rbind, lapply(cl_levels, function(cl) {
  cc <- call[b$.cluster == cl]; nn <- length(cc)
  tab <- sort(table(cc[cc != "unknown"]), decreasing = TRUE); frac <- if (length(tab)) tab / nn else numeric(0)
  major <- names(frac)[frac > 0.5]
  label <- if (length(major) == 1) major else if (length(major) >= 2) "mixed" else "unknown"
  data.frame(stage = stage, cluster = cl, n_cells = nn, majority = label,
             top_call = if (length(tab)) names(tab)[1] else "unknown",
             top_frac = if (length(tab)) round(as.numeric(frac[1]), 3) else 0,
             second_call = if (length(tab) >= 2) names(tab)[2] else NA,
             second_frac = if (length(tab) >= 2) round(as.numeric(frac[2]), 3) else 0,
             frac_unknown = round(mean(cc == "unknown"), 3))
}))
write.table(maj, op(".cluster_majority.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - per-cluster majority:"); print(maj[, c("cluster","n_cells","majority","top_call","top_frac","frac_unknown")], row.names = FALSE)

## ---- cluster x type mean per-cell z (smoothed analog of 4_3 table B) ----------------------------
meanZ <- sapply(cl_levels, function(cl) rowMeans(Z[, b$.cluster == cl, drop = FALSE]))
if (is.null(dim(meanZ))) meanZ <- matrix(meanZ, nrow = length(types), dimnames = list(types, cl_levels))
mz_df <- data.frame(stage = stage, cluster = rep(colnames(meanZ), each = nrow(meanZ)),
                    type_label = rep(rownames(meanZ), times = ncol(meanZ)), mean_cell_z = round(as.vector(meanZ), 3))
write.table(mz_df, op(".celltype_meanZ_by_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- exploratory UMAP: per-cell call + top-type z ----------------------------------------------
tryCatch({
  if (!all(c("umap1", "umap2") %in% colnames(b))) stop("no umap1/umap2 in metadata")
  png(file.path(out_dir, "plots", paste0(prefix, ".annotation_umap.png")), width = 12, height = 5.5, units = "in", res = 200, type = "cairo")
  par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
  labs <- factor(call); pal <- grDevices::rainbow(nlevels(labs))
  plot(b$umap1, b$umap2, col = pal[as.integer(labs)], pch = 16, cex = 0.3,
       main = paste0(prefix, " (", stage, ") per-cell call"), xlab = "umap1", ylab = "umap2")
  legend("topright", legend = levels(labs), col = pal, pch = 16, cex = 0.5, ncol = 2, bty = "n")
  v <- top_z; hi <- as.numeric(quantile(v, 0.99)); v[v > hi] <- hi
  ramp <- colorRampPalette(c("grey85", "goldenrod2", "firebrick3"))(100); o <- order(v)
  plot(b$umap1[o], b$umap2[o], col = ramp[cut(v, 100, labels = FALSE)][o], pch = 16, cex = 0.3,
       main = paste0(prefix, " top-type enrichment z"), xlab = "umap1", ylab = "umap2")
  dev.off(); message(" - wrote annotation UMAP png")
}, error = function(e) message(" ! annotation UMAP skipped: ", conditionMessage(e)))

message(" - DONE ", prefix)
