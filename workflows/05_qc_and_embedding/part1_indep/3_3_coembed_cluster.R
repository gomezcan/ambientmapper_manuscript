#!/usr/bin/env Rscript
###############################################################################
## 3_3_coembed_cluster.R  --  STAGE 3.3: final clustering of the CO-EMBED -> Fig 1F/G/H
##
## Clusters the joint (individual-mapping) co-embedding built by 3_1 at the config chosen from
## the 3_2 grid, and emits <prefix>.updated_metadata_v7.<tag>.txt carrying umap1/umap2 +
## LouvainClusters + Genome -- the exact columns the Fig 1 script (analysis/fig1_ambient_contamination)
## needs for F/G/H.
##
## Copied from 3_0_0b_cluster_pergenome.R (== 2_3_cluster.R). TWO co-embed deviations vs that
## copy, both required and both local to this file (2_3 is untouched):
##   (a) the constant-Genome stamp is now GUARDED -- the co-embed's 2-level plate-of-origin
##       label is preserved (it is the Fig-1F colour + the genome_mixing label);
##   (b) outputs land directly in <outdir> (the driver's <STAGE>_cluster/), not <outdir>/step2_cluster.
## min.t stays 0.001, matching 3_2_coembed_gridscan.R:213, so the grid-chosen config transfers.
##
## Per-genome fork of 3_0_0_Normalization_clustering.R for the split-genome objects.
## Two reasons 3_0_0 cannot be reused directly here:
##   1. res is HARDCODED to 0.5 in 3_0_0's callClusters. Step 2 must apply the resolution
##      FROZEN by 1b (3_0_1), so `resolution` is a REQUIRED arg here.
##   2. 3_0_0's species-normalization + At/B73 cross-tab + genome_mixing assume a `species`
##      column and TWO genomes; the per-genome v6 metadata has NEITHER, so that block
##      crashes (character(0) -> N-row assign). Single genome -> those diagnostics are
##      degenerate and dropped.
##
## Core pipeline is IDENTICAL to 3_0_0: align to v6 -> cleanData(min.c) -> tfidf ->
## reduceDims(SVD, n.pcs) -> projectUMAP(k,min_dist) -> detectDoublets(tryCatch) ->
## Socrates::callClusters(res, cl.method=4, m.clst=50) -> save v7 rds + metadata + reduced dims.
## Apply the SAME frozen (pcs,k_near,min_dist,min_c,res) to BOTH Pre and Post-wd of a genome
## (plan decision B) so any Pre-vs-Post difference is cleaning, not tuning.
##
## Usage:
##   Rscript 3_3_coembed_cluster.R <soc_rds> <meta_v6_tsv> <outdir> <prefix> \
##           <pcs> <k_near> <min_dist> <min_c> <resolution> [seed=1] [genome_label] [m_clst=50]
##     min_c: numeric feature/cell floor (from the frozen config), or "NA" -> data-driven 250 floor.
##     m_clst: min cluster size for callClusters (from the frozen config). Default 50 (unchanged
##             for all existing callers). Non-default values are appended to the output tag so
##             they never collide with a m.clst=50 object of the same pcs/k/res.
###############################################################################

suppressMessages(library(Socrates))
suppressMessages(library(igraph))
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(FNN))
suppressMessages(library(DelayedArray))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript 3_3_coembed_cluster.R <soc_rds> <meta_v6_tsv> <outdir> <prefix> <pcs> <k_near> <min_dist> <min_c> <resolution> [seed=1] [genome_label] [m_clst=50]")
}
obj_path   <- args[1]
meta_path  <- args[2]
outdir     <- args[3]
out_prefix <- args[4]
pcs        <- as.integer(args[5])
k_near     <- as.integer(args[6])
min_dis    <- as.numeric(args[7])
min_c_arg  <- as.numeric(args[8])                                   # "NA" -> data-driven floor
resolution <- as.numeric(args[9])
seed       <- if (length(args) >= 10) as.integer(args[10]) else 1L
genome_lab <- if (length(args) >= 11) args[11] else sub("^(SM2|Clean\\.SM2v2wd|Clean\\.SM2v2)_", "", out_prefix)
m_clst     <- if (length(args) >= 12) as.integer(args[12]) else 50L   # min cluster size; 50 = pipeline default

# ---- TFIDF (identical to 3_0_0) ---------------------------------------------
tfidf <- function(obj, frequencies = TRUE, log_scale_tf = TRUE, scale_factor = 10000,
                  doL2 = FALSE, slotName = "residuals") {
  bmat <- obj$counts
  .safe_tfidf <- function(tf, idf, block_size = 2000e6) {
    tryCatch({ tf * idf }, error = function(e) {
      options(DelayedArray.block.size = block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      tf <- DelayedArray(tf); idf <- as.matrix(idf); tf * idf
    })
  }
  if (frequencies) tf <- t(t(bmat) / Matrix::colSums(bmat)) else tf <- bmat
  if (log_scale_tf) tf@x <- log1p(tf@x * (if (frequencies) scale_factor else 1))
  idf <- log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  tf_idf_counts <- .safe_tfidf(tf, idf)
  if (doL2) {
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x = 1/colNorm)
  }
  rownames(tf_idf_counts) <- rownames(bmat); colnames(tf_idf_counts) <- colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse = TRUE); obj$norm_method <- "tfidf"
  obj
}

# ---- plotUMAP2 (verbatim from 3_0_0) ----------------------------------------
plotUMAP2 <- function(obj, column = "LouvainClusters", cluster_slotName = "Clusters",
                      filter_column = NULL, filter_value = NULL, main = "", xlab = "umap1",
                      ylab = "umap2", colors = NULL, point_size = 0.3, alpha = 1,
                      show_legend = TRUE, device = c("screen", "png", "pdf"), file = NULL,
                      width = 7, height = 7, dpi = 300, rasterize_points = FALSE) {
  device <- match.arg(device)
  if (is.null(obj[[cluster_slotName]])) cluster_slotName <- "meta"
  b <- obj[[cluster_slotName]]
  if (!is.null(filter_column) && !is.null(filter_value)) {
    if (!filter_column %in% colnames(b)) stop(" - ERROR: filter_column '", filter_column, "' missing")
    b <- b[b[[filter_column]] == filter_value, ]
  }
  b <- b[complete.cases(b$umap1) & complete.cases(b$umap2), ]
  if (!column %in% colnames(b)) stop(" - ERROR: column '", column, "' missing from ", cluster_slotName)
  col_data <- b[[column]]
  is_disc  <- is.factor(col_data) || is.character(col_data)
  if (is_disc && is.character(col_data)) b[[column]] <- factor(col_data)
  library(ggplot2)
  aes_map <- aes(x = umap1, y = umap2, color = .data[[column]])
  use_raster <- FALSE
  if (rasterize_points && requireNamespace("ggrastr", quietly = TRUE)) {
    use_raster <- TRUE; geom_fun <- ggrastr::geom_point_rast
  } else geom_fun <- geom_point
  p <- ggplot(b, aes_map) + geom_fun(size = point_size, alpha = alpha) +
    labs(title = main, x = xlab, y = ylab, color = column) + theme_bw()
  if (is_disc) {
    if (is.null(colors)) {
      levs <- levels(b[[column]]); if (is.null(levs)) levs <- sort(unique(b[[column]]))
      pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levs))
    } else pal <- colors
    p <- p + scale_color_manual(values = pal)
  } else if (is.numeric(col_data)) {
    if (is.null(colors)) p <- p + viridis::scale_color_viridis(option = "D")
    else p <- p + scale_color_gradientn(colours = colors)
  } else stop(" - ERROR: column '", column, "' must be factor/character or numeric.")
  if (!show_legend) p <- p + theme(legend.position = "none")
  if (device == "screen") print(p) else {
    if (is.null(file) || !nzchar(file)) {
      base_name <- if (nzchar(main)) main else column
      base_name <- gsub("[^A-Za-z0-9_]+", "_", base_name)
      file <- paste0(base_name, "_UMAP.", if (device == "png") "png" else "pdf")
    }
    if (device == "png") ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi, device = "png")
    else ggplot2::ggsave(file, p, width = width, height = height, device = cairo_pdf)
  }
  invisible(p)
}

# ---- (1) load + align to v6 -------------------------------------------------
obj <- readRDS(obj_path)
meta.data <- read.table(meta_path, header = TRUE, sep = "\t")
meta.data <- meta.data[, !colnames(meta.data) == "...1"]
rownames(meta.data) <- meta.data$cellID

tag <- paste0(".pcs_", pcs, ".k_near_", k_near, ".min_dis_", min_dis,
              ".minc_", ifelse(is.na(min_c_arg), "auto", as.integer(min_c_arg)),
              ".res_", resolution,
              if (m_clst != 50L) paste0(".mclst_", m_clst) else "")   # tag only when non-default -> no collision, indep tags unchanged
out_dir_plots <- file.path(outdir, "plots")
dir.create(out_dir_plots, showWarnings = FALSE, recursive = TRUE)
out_dir_obj <- outdir                    # co-embed: v7 outputs land directly in <STAGE>_cluster/
dir.create(out_dir_obj, showWarnings = FALSE, recursive = TRUE)

obj$counts <- obj$counts[, colnames(obj$counts) %in% rownames(meta.data)]
obj$meta   <- meta.data[colnames(obj$counts), ]

## Genome label -- DO NOT blindly stamp the constant here.
## Per-genome v6 metadata carries no Genome column, so the constant label is correct there
## (original 3_0_0b behaviour, kept for that path). The CO-EMBED metadata built by 3_1 carries a
## real 2-level plate-of-origin label, and that label is exactly what Fig 1F colours by and what
## genome_mixing is computed on -- overwriting it makes the panel unmakeable and mixing degenerate.
## Mirrors the guard in 3_2_coembed_gridscan.R:194-197 so clustering and the grid agree.
if (!"Genome" %in% colnames(obj$meta) || length(unique(na.omit(obj$meta$Genome))) < 2) {
  obj$meta$Genome <- genome_lab
  message(" - Genome: stamped constant '", genome_lab, "'")
} else {
  gt <- table(obj$meta$Genome, useNA = "ifany")
  message(" - Genome: PRESERVED from metadata (", paste(names(gt), gt, sep = "=", collapse = "  "), ")")
}
message(" - ", out_prefix, " (", genome_lab, "): ", ncol(obj$counts), " cells after v6 intersect")

# ---- (2) min.c threshold (data-driven floor 250, or fixed from frozen config) ----
cell.counts   <- log10(Matrix::colSums(obj$counts))
cell.counts.z <- as.numeric(scale(cell.counts))
mask <- cell.counts.z[(cell.counts.z) < -0.5]
cell.counts.threshold <- if (is.na(min_c_arg)) max(c((10^mask), 250)) else min_c_arg
message(" - cleanData min.c (features/cell) = ", signif(cell.counts.threshold, 6))

set.seed(seed)

# ---- (3) clean -> (4) tfidf -> (5) SVD -> (6) UMAP --------------------------
soc.obj <- cleanData(obj, min.c = cell.counts.threshold, min.t = 0.001, max.t = 0, verbose = TRUE)
soc.obj <- tfidf(soc.obj, doL2 = TRUE)
number.sites <- ceiling(nrow(soc.obj$counts) * 0.5)
soc.obj <- reduceDims(soc.obj, method = "SVD", n.pcs = pcs, cor.max = 0.6, num.var = number.sites,
                      verbose = TRUE, scaleVar = TRUE, doSTD = FALSE, doL1 = FALSE, doL2 = TRUE,
                      refit_residuals = FALSE)
soc.obj <- projectUMAP(soc.obj, verbose = TRUE, k.near = k_near, m.dist = min_dis)

plotUMAP2(soc.obj, cluster_slotName = "meta", column = "log10nSites", point_size = 0.2, alpha = 0.5,
          width = 6, height = 5, device = "png",
          file = paste0(out_dir_plots, "/", out_prefix, "_scidiATAC_umap", tag, "_log10nSites.png"))

# ---- (7a) doublets (best-effort; QC-only, never blocks clustering) ----------
ok_dbl <- tryCatch({ soc.obj <- detectDoublets(soc.obj, threads = 6); TRUE },
                   error = function(e) { message(" - detectDoublets SKIPPED: ", conditionMessage(e)); FALSE })
if (ok_dbl) {
  plotUMAP2(soc.obj, cluster_slotName = "meta", column = "doubletscore", point_size = 0.3, alpha = 0.5,
            width = 6, height = 5, device = "png",
            file = paste0(out_dir_plots, "/", out_prefix, "_scidiATAC_umap", tag, "_doubletscore.png"))
}

# ---- (7b) clusters at the FROZEN resolution (else identical to 3_0_0) -------
soc.obj <- Socrates::callClusters(soc.obj, res = resolution, k.near = k_near, verbose = TRUE,
                                  cleanCluster = FALSE, cl.method = 4, e.thresh = 3,
                                  threshold = 3, m.clst = m_clst)

plotUMAP2(soc.obj, cluster_slotName = "Clusters", column = "LouvainClusters", point_size = 0.2, alpha = 0.5,
          width = 6, height = 5, device = "png",
          file = paste0(out_dir_plots, "/", out_prefix, "_scidiATAC_umap", tag, "_LouvainClusters.png"))

# ---- (8) save v7 object + metadata + reduced dims ---------------------------
saveRDS(soc.obj, file = file.path(out_dir_obj, paste0(out_prefix, ".full.SocObj_v7", tag, ".rds")))
write.table(soc.obj$Clusters, file = file.path(out_dir_obj, paste0(out_prefix, ".updated_metadata_v7", tag, ".txt")),
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(soc.obj$PCA[rownames(soc.obj$Clusters), ],
            file = file.path(out_dir_obj, paste0(out_prefix, ".reduced_dimensions_v7", tag, ".txt")),
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

n_clust <- length(unique(soc.obj$Clusters$LouvainClusters))
message(" - ", out_prefix, ": ", n_clust, " clusters at res=", resolution, " m.clst=", m_clst,
        " | ", nrow(soc.obj$Clusters), " clustered cells")
message(" - DONE. v7 object + metadata + reduced dims written to ", out_dir_obj)
