#!/usr/bin/env Rscript
###############################################################################
## 3_0_1b_umap_panels.R  --  UMAP visualization companion to the 3_0_1 res scan
##                           (Part 2, Step 1b; MARKER-FREE)
##
## WHY THIS EXISTS
##   3_0_1_resolution_scan.R measures cluster stability/purity/count in PC space and
##   plots only those metric CURVES -- no UMAP. The by-eye call the plan requires for
##   At ("continuum vs real islands; pick bins by eye") cannot be made from curves, and
##   a single pcs is not a decision -- pcs must be swept (like the 2_0_2b / 3_0_4 grids).
##
##   This renders the SAME embedding Step 2 (3_0_0) would produce -- verbatim tfidf,
##   reduceDims + cor.max, projectUMAP -- across a pcs x resolution GRID, drawn:
##     (a) colored by Leiden cluster at each (pcs,res)  -> granularity vs pcs & res
##     (b) colored by QC (log10 depth, log10 nSites, FRiP, pTSS, pOrg, dif), faceted by pcs
##   Panel (b) is the discriminator: if At's split just tracks a depth/sites gradient,
##   that is the CONTINUUM signature (root gradient), not two biological islands.
##
##   reduceDims is recomputed PER pcs (truncated SVD at exactly that pcs) -- identical to
##   the 3_0_1 scan and to what 3_0_0 would freeze -- so these UMAPs correspond exactly to
##   the stability numbers in <prefix>.resolution_scan.tsv. Each pcs embedding is rescaled
##   (center + shared-scalar divide, aspect preserved) to a common [-1,1] frame so panels
##   are directly comparable in facet_grid.
##
##   NA handling: cells callClusters leaves unlabeled are NOT dropped -- they are counted
##   and reported per (pcs,res) in <prefix>.umap_grid_summary.tsv and drawn grey ("NA").
##
## Usage:
##   Rscript 3_0_1b_umap_panels.R <soc_rds> <meta_tsv> <outdir> <prefix> \
##           <pcs_list> <k_near> <min_dist> [min_c=NA] [res_grid] [seed=1]
##     pcs_list: comma-separated, e.g. "5,8,10,15,20,30"
##     res_grid: comma-separated, e.g. "0.1,0.2,0.3,0.5,0.8,1.0"
##     min_c: numeric cleanData floor, or "NA" -> data-driven 250 floor (matches scan).
###############################################################################

suppressMessages({
  library(Socrates)
  library(Matrix)
  library(FNN)
  library(data.table)
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript 3_0_1b_umap_panels.R <soc_rds> <meta_tsv> <outdir> <prefix> <pcs_list> <k_near> <min_dist> [min_c=NA] [res_grid] [seed=1]")
}
soc_rds   <- args[1]
meta_tsv  <- args[2]
outdir    <- args[3]
prefix    <- args[4]
pcs_vec   <- as.integer(strsplit(args[5], ",")[[1]])
k_near    <- as.integer(args[6])
min_dist  <- as.numeric(args[7])
min_c_arg <- if (length(args) >= 8) as.numeric(args[8]) else NA_real_
RES_GRID  <- if (length(args) >= 9 && nzchar(args[9])) as.numeric(strsplit(args[9], ",")[[1]]) else c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0)
seed      <- if (length(args) >= 10) as.integer(args[10]) else 1L

dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(outdir, paste0(prefix, suffix))

message(" - 3_0_1b UMAP grid | prefix=", prefix, " | pcs={", paste(pcs_vec, collapse = ","), "}",
        " k_near=", k_near, " min_dist=", min_dist,
        " | min_c=", ifelse(is.na(min_c_arg), "NA(data-driven)", min_c_arg),
        " | res={", paste(RES_GRID, collapse = ","), "}")

# -------------------------------------------------------------------------
# helpers -- embedding path VERBATIM from 3_0_1_resolution_scan.R (keep in sync)
# -------------------------------------------------------------------------
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
  rownames(tf_idf_counts) <- rownames(bmat)
  colnames(tf_idf_counts) <- colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse = TRUE)
  obj$norm_method <- "tfidf"
  obj
}

# cluster with the EXACT 3_0_0 call, return a named label vector over cellIDs
cluster_labels <- function(socobj, res, k_near) {
  cc <- Socrates::callClusters(socobj, res = res, k.near = k_near, verbose = FALSE,
                               cleanCluster = FALSE, cl.method = 4,
                               e.thresh = 3, threshold = 3, m.clst = 50)
  lab <- as.character(cc$Clusters$LouvainClusters)
  names(lab) <- rownames(cc$Clusters)
  lab
}

# winsorize for stable continuous color scales
winsorize <- function(x, p = 0.02) {
  qs <- quantile(x, c(p, 1 - p), na.rm = TRUE); pmin(pmax(x, qs[1]), qs[2])
}

# center + shared-scalar rescale to ~[-1,1], aspect ratio preserved (UMAP coords arbitrary)
rescale_xy <- function(u1, u2) {
  u1c <- u1 - mean(u1); u2c <- u2 - mean(u2)
  s <- max(abs(c(u1c, u2c))); if (!is.finite(s) || s == 0) s <- 1
  list(x = u1c / s, y = u2c / s)
}

# -------------------------------------------------------------------------
# Load + align (mirror 3_0_1 / 3_0_0 sec 1)
# -------------------------------------------------------------------------
obj  <- readRDS(soc_rds)
meta <- read.table(meta_tsv, header = TRUE, sep = "\t", check.names = FALSE)
if ("...1" %in% colnames(meta)) meta <- meta[, colnames(meta) != "...1", drop = FALSE]
stopifnot("cellID" %in% colnames(meta))
rownames(meta) <- meta$cellID

shared <- intersect(colnames(obj$counts), rownames(meta))
if (length(shared) == 0) stop("No shared cells between obj$counts and metadata.")
obj$counts <- obj$counts[, shared, drop = FALSE]
obj$meta   <- meta[shared, , drop = FALSE]
message(" - aligned cells (counts n v6 meta): ", length(shared))

# -------------------------------------------------------------------------
# cleanData + tfidf ONCE (shared); reduceDims/UMAP/cluster PER pcs
# -------------------------------------------------------------------------
set.seed(seed)
cell.counts   <- log10(Matrix::colSums(obj$counts))
cell.counts.z <- as.numeric(scale(cell.counts))
mask  <- cell.counts.z[cell.counts.z < -0.5]
min_c <- if (is.na(min_c_arg)) max(c(10^mask, 250), na.rm = TRUE) else min_c_arg
message(" - cleanData min.c = ", signif(min_c, 4))

soc_norm     <- cleanData(obj, min.c = min_c, min.t = 0.001, max.t = 0, verbose = TRUE)
soc_norm     <- tfidf(soc_norm, doL2 = TRUE)
number.sites <- ceiling(nrow(soc_norm$counts) * 0.5)
message(" - normalized: ", ncol(soc_norm$counts), " cells x ", nrow(soc_norm$counts),
        " sites; num.var=", number.sites)

coords_list  <- list()   # (cellID, pcs_req, pcs_act, plab, x, y) rescaled coords per pcs
clus_list    <- list()   # (cellID, pcs_req, res, cluster)
summ_list    <- list()   # (pcs_req, pcs_act, res, n_clusters, n_unlabeled, min_size, median_size)

for (pcs_i in pcs_vec) {
  message(" ===== pcs = ", pcs_i, " =====")
  soc <- reduceDims(soc_norm, method = "SVD", n.pcs = pcs_i, cor.max = 0.6, num.var = number.sites,
                    verbose = TRUE, scaleVar = TRUE, doSTD = FALSE, doL1 = FALSE, doL2 = TRUE,
                    refit_residuals = FALSE)
  soc <- projectUMAP(soc, verbose = FALSE, k.near = k_near, m.dist = min_dist)
  pcs_act <- ncol(soc$PCA)
  if (pcs_act < pcs_i) message("   ! reduceDims returned ", pcs_act, " PCs (< requested ", pcs_i, ") after cor.max drop")
  plab <- sprintf("pcs=%d%s", pcs_i, ifelse(pcs_act < pcs_i, sprintf(" (r%d)", pcs_act), ""))

  um  <- as.data.frame(soc$UMAP)
  cid <- rownames(um)
  rc  <- rescale_xy(um[[1]], um[[2]])
  coords_list[[as.character(pcs_i)]] <- data.frame(
    cellID = cid, pcs_req = pcs_i, pcs_act = pcs_act, plab = plab,
    x = rc$x, y = rc$y, stringsAsFactors = FALSE)

  for (res in RES_GRID) {
    lab    <- cluster_labels(soc, res, k_near)
    lab_al <- lab[cid]                                # align to embedding; NA where unlabeled
    n_unlab <- sum(is.na(lab_al))
    sizes  <- as.integer(table(lab_al[!is.na(lab_al)]))
    n_clust <- length(sizes)
    clus_list[[length(clus_list) + 1L]] <- data.frame(
      cellID = cid, pcs_req = pcs_i, res = res,
      cluster = ifelse(is.na(lab_al), "NA", lab_al), stringsAsFactors = FALSE)
    summ_list[[length(summ_list) + 1L]] <- data.frame(
      prefix = prefix, pcs_req = pcs_i, pcs_act = pcs_act, res = res,
      n_clusters = n_clust, n_unlabeled = n_unlab,
      min_size = if (n_clust) min(sizes) else NA_integer_,
      median_size = if (n_clust) as.numeric(median(sizes)) else NA_real_,
      stringsAsFactors = FALSE)
    message(sprintf("   pcs=%d(act %d) res=%.2f -> %d clusters (%d unlabeled)",
                    pcs_i, pcs_act, res, n_clust, n_unlab))
  }
  rm(soc); gc()
}

coords <- do.call(rbind, coords_list)
clus   <- do.call(rbind, clus_list)
summ   <- do.call(rbind, summ_list)

# ordered facet labels
plab_levels <- sprintf("pcs=%d%s", pcs_vec,
                       vapply(pcs_vec, function(p) {
                         pa <- coords$pcs_act[match(p, coords$pcs_req)]
                         ifelse(pa < p, sprintf(" (r%d)", pa), "") }, character(1)))
coords$plab <- factor(coords$plab, levels = plab_levels)

# -------------------------------------------------------------------------
# Write decision tables: long summary + wide n_clusters matrix (+ n_unlabeled)
# -------------------------------------------------------------------------
summ_path <- op(".umap_grid_summary.tsv")
write.table(summ, summ_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote grid summary (n_clusters + n_unlabeled per pcs x res): ", summ_path)

wide_nclust <- reshape(summ[, c("pcs_req", "res", "n_clusters")],
                       idvar = "pcs_req", timevar = "res", direction = "wide")
colnames(wide_nclust) <- sub("^n_clusters\\.", "res_", colnames(wide_nclust))
message(" - n_clusters grid (rows=pcs, cols=res):")
print(wide_nclust, row.names = FALSE)
if (any(summ$n_unlabeled > 0)) {
  message(" - NOTE: unlabeled (callClusters-dropped) cells present at:")
  print(summ[summ$n_unlabeled > 0, c("pcs_req", "res", "n_clusters", "n_unlabeled")], row.names = FALSE)
} else {
  message(" - NA check: callClusters labeled ALL cells at every (pcs,res) -- zero unlabeled.")
}

# -------------------------------------------------------------------------
# Plotting frames
# -------------------------------------------------------------------------
clus_plot <- merge(clus, coords[, c("cellID", "pcs_req", "plab", "x", "y")],
                   by = c("cellID", "pcs_req"), sort = FALSE)
clus_plot$plab <- factor(clus_plot$plab, levels = plab_levels)
clus_plot$res_label <- factor(sprintf("res=%.2f", clus_plot$res),
                              levels = sprintf("res=%.2f", sort(unique(clus_plot$res))))
# n_clusters annotation per facet cell (top-left)
ann <- summ; ann$plab <- factor(sprintf("pcs=%d%s", ann$pcs_req, ifelse(ann$pcs_act < ann$pcs_req, sprintf(" (r%d)", ann$pcs_act), "")),
                                levels = plab_levels)
ann$res_label <- factor(sprintf("res=%.2f", ann$res), levels = levels(clus_plot$res_label))

# QC overlays (per cell; present columns only)
qc_defs <- list(
  log10_depth  = function(m) log10(m$total + 1),
  log10_nSites = function(m) if ("log10nSites" %in% names(m)) m$log10nSites else log10(m$nSites + 1),
  FRiP = function(m) m$FRiP, pTSS = function(m) m$pTSS,
  pOrg = function(m) m$pOrg, dif  = function(m) m$dif
)
mm <- meta[unique(coords$cellID), , drop = FALSE]
qc_have <- names(qc_defs)[vapply(qc_defs, function(f) !is.null(tryCatch(f(mm), error = function(e) NULL)), logical(1))]
qc_tab <- data.frame(cellID = rownames(mm), stringsAsFactors = FALSE)
for (nm in qc_have) qc_tab[[nm]] <- as.numeric(qc_defs[[nm]](mm))
message(" - QC overlays available: ", paste(qc_have, collapse = ", "))

pt_size  <- if (nrow(mm) > 6000) 0.30 else 0.6
pt_alpha <- if (nrow(mm) > 6000) 0.5 else 0.8
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.key.size = unit(0.35, "cm"),
        axis.text = element_blank(), axis.ticks = element_blank())

npcs <- length(pcs_vec); nres <- length(RES_GRID)
pdf_path <- file.path(outdir, "plots", paste0(prefix, ".umap_panels.pdf"))
pdf(pdf_path, width = max(8, 1.9 * nres + 1.5), height = max(6, 1.9 * npcs + 1))

## Page 1: THE GRID -- pcs (rows) x resolution (cols), colored by cluster
print(
  ggplot(clus_plot, aes(x, y, color = cluster)) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    geom_text(data = ann, aes(x = -0.95, y = 0.95, label = n_clusters),
              inherit.aes = FALSE, size = 3, hjust = 0, vjust = 1, fontface = "bold") +
    facet_grid(plab ~ res_label) +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    scale_color_discrete(guide = "none") +
    base_theme +
    labs(title = paste0(prefix, " -- UMAP cluster granularity across pcs x resolution"),
         subtitle = paste0("k_near=", k_near, " min_dist=", min_dist, " min.c=", signif(min_c, 4),
                           "  |  number top-left = n_clusters; colors per-facet (not comparable across facets)"),
         x = NULL, y = NULL)
)

## QC pages: one per metric, faceted by pcs (does the split track a technical gradient?)
qc_merge_base <- coords[, c("cellID", "pcs_req", "plab", "x", "y")]
for (nm in qc_have) {
  qv <- data.frame(cellID = qc_tab$cellID, value = winsorize(qc_tab[[nm]]), stringsAsFactors = FALSE)
  dd <- merge(qc_merge_base, qv, by = "cellID", sort = FALSE)
  dd$plab <- factor(dd$plab, levels = plab_levels)
  print(
    ggplot(dd, aes(x, y, color = value)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      facet_wrap(~plab, nrow = 1) +
      coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
      scale_color_viridis_c(name = nm, option = "C") +
      base_theme +
      labs(title = paste0(prefix, " -- UMAP colored by ", nm, " (across pcs)"),
           subtitle = "winsorized 2-98%; a cluster split that tracks this = technical continuum, not biology",
           x = NULL, y = NULL)
  )
}

dev.off()
message(" - wrote PDF: ", pdf_path)
message(" - DONE. Page 1 = pcs x res cluster grid; QC pages = continuum check. Decision table: ", summ_path)
