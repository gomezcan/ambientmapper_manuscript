#!/usr/bin/env Rscript
###############################################################################
## 3_0_1c_umap_kgrid.R  --  UMAP + cluster-count grid over pcs x k_near x res
##                          at a chosen m.clst  (Part 2, Step 1b; MARKER-FREE)
##
## WHY THIS EXISTS (extends 3_0_1b)
##   3_0_1b swept pcs x res with k_near, min_dist, min_c, m.clst FIXED. For a SMALL cell
##   set (plate At ~1090), k_near (graph neighbours) and m.clst (min cluster size in
##   callClusters) are large relative to n and can THEMSELVES force the merge to 2 clusters
##   -- so "At = 2" could be a k/m.clst artifact rather than biology. This grid sweeps
##   k_near and lowers m.clst to test that, per pcs, at every resolution.
##
##   k_near enters BOTH projectUMAP (layout) and callClusters (graph), so each (pcs,k_near)
##   has its OWN embedding -- faithful to how Step 2 (3_0_0) would freeze it. Each embedding
##   is rescaled (center + shared-scalar divide, aspect preserved) to a common [-1,1] frame.
##   reduceDims is k-independent -> computed once per pcs and reused across k_near.
##
##   NA handling: cells callClusters leaves unlabeled are counted (n_unlabeled), reported in
##   the summary + an n_unlabeled heatmap, and drawn grey. Lowering m.clst should reduce them.
##
## Usage:
##   Rscript 3_0_1c_umap_kgrid.R <soc_rds> <meta_tsv> <outdir> <prefix> \
##           <pcs_list> <k_list> <min_dist> <min_c> <res_grid> <m_clst> [seed=1]
##     pcs_list "5,8,10"   k_list "10,15,20,30"   res_grid "0.1,0.2,0.3,0.5,0.8,1.0"
##     min_c: numeric cleanData floor, or "NA" -> data-driven 250 floor.   m_clst e.g. 40
## Output: <outdir>/plots/<prefix>.umap_panels.kgrid.pdf  (+ heatmaps)
##         <outdir>/<prefix>.umap_grid_summary.kgrid.tsv  (n_clusters + n_unlabeled per pcs x k x res)
###############################################################################

suppressMessages({
  library(Socrates); library(Matrix); library(FNN)
  library(data.table); library(ggplot2); library(Seurat); library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: Rscript 3_0_1c_umap_kgrid.R <soc_rds> <meta_tsv> <outdir> <prefix> <pcs_list> <k_list> <min_dist> <min_c> <res_grid> <m_clst> [seed=1]")
}
soc_rds  <- args[1]; meta_tsv <- args[2]; outdir <- args[3]; prefix <- args[4]
pcs_vec  <- as.integer(strsplit(args[5], ",")[[1]])
k_vec    <- as.integer(strsplit(args[6], ",")[[1]])
min_dist <- as.numeric(args[7])
min_c_arg<- if (identical(args[8], "NA")) NA_real_ else as.numeric(args[8])
res_grid <- as.numeric(strsplit(args[9], ",")[[1]])
m_clst   <- as.integer(args[10])
seed     <- if (length(args) >= 11) as.integer(args[11]) else 1L

dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(outdir, paste0(prefix, suffix))
message(" - 3_0_1c kgrid | ", prefix, " | pcs={", paste(pcs_vec, collapse=","), "}",
        " k_near={", paste(k_vec, collapse=","), "} res={", paste(res_grid, collapse=","),
        "} m.clst=", m_clst, " min_dist=", min_dist,
        " min_c=", ifelse(is.na(min_c_arg), "NA(data-driven)", min_c_arg))

# ---- helpers (embedding path VERBATIM from 3_0_1 scan; cluster_labels takes m.clst) ----
tfidf <- function(obj, frequencies = TRUE, log_scale_tf = TRUE, scale_factor = 10000,
                  doL2 = FALSE, slotName = "residuals") {
  bmat <- obj$counts
  .safe_tfidf <- function(tf, idf, block_size = 2000e6) {
    tryCatch({ tf * idf }, error = function(e) {
      options(DelayedArray.block.size = block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      tf <- DelayedArray(tf); idf <- as.matrix(idf); tf * idf })
  }
  if (frequencies) tf <- t(t(bmat) / Matrix::colSums(bmat)) else tf <- bmat
  if (log_scale_tf) tf@x <- log1p(tf@x * (if (frequencies) scale_factor else 1))
  idf <- log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  tf_idf_counts <- .safe_tfidf(tf, idf)
  if (doL2) { colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2)); tf_idf_counts <- tf_idf_counts %*% Diagonal(x = 1/colNorm) }
  rownames(tf_idf_counts) <- rownames(bmat); colnames(tf_idf_counts) <- colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse = TRUE); obj$norm_method <- "tfidf"; obj
}
cluster_labels <- function(socobj, res, k_near, m_clst) {
  cc <- Socrates::callClusters(socobj, res = res, k.near = k_near, verbose = FALSE,
                               cleanCluster = FALSE, cl.method = 4,
                               e.thresh = 3, threshold = 3, m.clst = m_clst)
  lab <- as.character(cc$Clusters$LouvainClusters); names(lab) <- rownames(cc$Clusters); lab
}
winsorize  <- function(x, p = 0.02) { qs <- quantile(x, c(p, 1 - p), na.rm = TRUE); pmin(pmax(x, qs[1]), qs[2]) }
rescale_xy <- function(u1, u2) { u1c <- u1 - mean(u1); u2c <- u2 - mean(u2)
  s <- max(abs(c(u1c, u2c))); if (!is.finite(s) || s == 0) s <- 1; list(x = u1c/s, y = u2c/s) }

# ---- load + align ----
obj  <- readRDS(soc_rds)
meta <- read.table(meta_tsv, header = TRUE, sep = "\t", check.names = FALSE)
if ("...1" %in% colnames(meta)) meta <- meta[, colnames(meta) != "...1", drop = FALSE]
stopifnot("cellID" %in% colnames(meta)); rownames(meta) <- meta$cellID
shared <- intersect(colnames(obj$counts), rownames(meta))
if (length(shared) == 0) stop("No shared cells between obj$counts and metadata.")
obj$counts <- obj$counts[, shared, drop = FALSE]; obj$meta <- meta[shared, , drop = FALSE]
message(" - aligned cells: ", length(shared))

# ---- cleanData + tfidf ONCE ----
set.seed(seed)
cell.counts <- log10(Matrix::colSums(obj$counts)); cell.counts.z <- as.numeric(scale(cell.counts))
mask <- cell.counts.z[cell.counts.z < -0.5]
min_c <- if (is.na(min_c_arg)) max(c(10^mask, 250), na.rm = TRUE) else min_c_arg
message(" - cleanData min.c = ", signif(min_c, 4))
soc_norm <- cleanData(obj, min.c = min_c, min.t = 0.001, max.t = 0, verbose = TRUE)
soc_norm <- tfidf(soc_norm, doL2 = TRUE)
number.sites <- ceiling(nrow(soc_norm$counts) * 0.5)
message(" - normalized: ", ncol(soc_norm$counts), " cells x ", nrow(soc_norm$counts), " sites")

# ---- sweep pcs (reduceDims) x k_near (projectUMAP) x res (callClusters) ----
coords_list <- list(); clus_list <- list(); summ_list <- list()
for (pcs_i in pcs_vec) {
  message(" ===== pcs = ", pcs_i, " =====")
  soc_pca <- reduceDims(soc_norm, method = "SVD", n.pcs = pcs_i, cor.max = 0.6, num.var = number.sites,
                        verbose = TRUE, scaleVar = TRUE, doSTD = FALSE, doL1 = FALSE, doL2 = TRUE,
                        refit_residuals = FALSE)
  pcs_act <- ncol(soc_pca$PCA)
  if (pcs_act < pcs_i) message("   ! reduceDims returned ", pcs_act, " PCs (< ", pcs_i, ")")
  for (k in k_vec) {
    set.seed(seed)
    soc <- projectUMAP(soc_pca, verbose = FALSE, k.near = k, m.dist = min_dist)
    um  <- as.data.frame(soc$UMAP); cid <- rownames(um); rc <- rescale_xy(um[[1]], um[[2]])
    coords_list[[paste0(pcs_i, "_", k)]] <- data.frame(
      cellID = cid, pcs_req = pcs_i, pcs_act = pcs_act, k_near = k, x = rc$x, y = rc$y, stringsAsFactors = FALSE)
    for (res in res_grid) {
      lab_al  <- cluster_labels(soc, res, k, m_clst)[cid]
      n_unlab <- sum(is.na(lab_al)); sizes <- as.integer(table(lab_al[!is.na(lab_al)])); n_clust <- length(sizes)
      clus_list[[length(clus_list)+1L]] <- data.frame(
        cellID = cid, pcs_req = pcs_i, k_near = k, res = res,
        cluster = ifelse(is.na(lab_al), "NA", lab_al), stringsAsFactors = FALSE)
      summ_list[[length(summ_list)+1L]] <- data.frame(
        prefix = prefix, pcs_req = pcs_i, pcs_act = pcs_act, k_near = k, res = res, m_clst = m_clst,
        n_clusters = n_clust, n_unlabeled = n_unlab,
        min_size = if (n_clust) min(sizes) else NA_integer_,
        median_size = if (n_clust) as.numeric(median(sizes)) else NA_real_, stringsAsFactors = FALSE)
      message(sprintf("   pcs=%d k=%d res=%.2f -> %d clusters (%d unlabeled)", pcs_i, k, res, n_clust, n_unlab))
    }
    rm(soc); gc()
  }
  rm(soc_pca); gc()
}
coords <- do.call(rbind, coords_list); clus <- do.call(rbind, clus_list); summ <- do.call(rbind, summ_list)

# ---- decision table ----
summ_path <- op(".umap_grid_summary.kgrid.tsv")
write.table(summ, summ_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote grid summary: ", summ_path)
if (any(summ$n_unlabeled > 0)) {
  message(" - unlabeled cells present (max ", max(summ$n_unlabeled), " at some pcs/k/res); see n_unlabeled heatmap.")
} else {
  message(" - NA check: every (pcs,k,res) labeled ALL cells -- zero unlabeled at m.clst=", m_clst)
}

# ---- plotting frames ----
pcs_labs <- paste0("pcs=", pcs_vec)
summ$pcs_lab <- factor(paste0("pcs=", summ$pcs_req), levels = pcs_labs)
summ$k_f     <- factor(summ$k_near, levels = sort(unique(summ$k_near)))
summ$res_f   <- factor(sprintf("res=%.2f", summ$res), levels = sprintf("res=%.2f", sort(unique(summ$res))))
# contrast-aware text color for heatmaps
txtcol <- function(v) { fr <- (v - min(v)) / (max(v) - min(v) + 1e-9); ifelse(fr > 0.55, "black", "white") }

clus_plot <- merge(clus, coords[, c("cellID","pcs_req","k_near","x","y")], by = c("cellID","pcs_req","k_near"), sort = FALSE)
clus_plot$k_f   <- factor(clus_plot$k_near, levels = sort(unique(clus_plot$k_near)))
clus_plot$res_f <- factor(sprintf("res=%.2f", clus_plot$res), levels = levels(summ$res_f))
ann <- summ  # carries n_clusters per (pcs,k,res)

# QC (depth) overlay values
qcv <- data.frame(cellID = rownames(meta), log10_depth = log10(meta$total + 1), stringsAsFactors = FALSE)

pt_size <- 0.55; pt_alpha <- 0.8
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.key.size = unit(0.35,"cm"),
        axis.text = element_blank(), axis.ticks = element_blank())

nres <- length(res_grid); nk <- length(k_vec)
pdf_path <- file.path(outdir, "plots", paste0(prefix, ".umap_panels.kgrid.pdf"))
pdf(pdf_path, width = max(9, 1.7*nres + 2), height = max(7, 1.7*nk + 1.5))

## Page 1: n_clusters heatmap (res x k_near, faceted by pcs) -- the money view
summ$tc_nc <- txtcol(summ$n_clusters)
print(
  ggplot(summ, aes(res_f, k_f, fill = n_clusters)) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = n_clusters, color = tc_nc), size = 3.4, fontface = "bold") +
    facet_wrap(~pcs_lab, nrow = 1) +
    scale_fill_viridis_c(name = "n_clusters", option = "D") +
    scale_color_identity() +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(prefix, " -- n_clusters over k_near x resolution (m.clst=", m_clst, ")"),
         subtitle = "does lowering k_near / m.clst surface >2 clusters, or does At stay at 2?",
         x = "resolution", y = "k_near")
)

## Page 2: n_unlabeled heatmap
summ$tc_nu <- txtcol(summ$n_unlabeled)
print(
  ggplot(summ, aes(res_f, k_f, fill = n_unlabeled)) +
    geom_tile(color = "grey85") +
    geom_text(aes(label = n_unlabeled, color = tc_nu), size = 3.0) +
    facet_wrap(~pcs_lab, nrow = 1) +
    scale_fill_viridis_c(name = "n_unlabeled", option = "A") +
    scale_color_identity() +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(prefix, " -- n_unlabeled (callClusters-dropped) over k_near x resolution (m.clst=", m_clst, ")"),
         subtitle = "high = algorithm abandoning cells = over-resolving a continuum", x = "resolution", y = "k_near")
)

## Pages 3..: per pcs, UMAP grid facet_grid(k_near ~ res), colored by cluster
for (pcs_i in pcs_vec) {
  d  <- clus_plot[clus_plot$pcs_req == pcs_i, , drop = FALSE]
  aa <- ann[ann$pcs_req == pcs_i, , drop = FALSE]
  print(
    ggplot(d, aes(x, y, color = cluster)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      geom_text(data = aa, aes(x = -0.95, y = 0.95, label = n_clusters), inherit.aes = FALSE,
                size = 3, hjust = 0, vjust = 1, fontface = "bold") +
      facet_grid(k_f ~ res_f) +
      coord_fixed(xlim = c(-1,1), ylim = c(-1,1)) +
      scale_color_discrete(guide = "none") + base_theme +
      labs(title = paste0(prefix, " -- pcs=", pcs_i, " : UMAP by cluster, k_near (rows) x resolution (cols), m.clst=", m_clst),
           subtitle = "number top-left = n_clusters; colors per-facet; each k_near row is its own embedding", x = NULL, y = NULL)
  )
}

## Pages: per pcs, UMAP colored by log10_depth, faceted by k_near (continuum check under each k)
for (pcs_i in pcs_vec) {
  cc <- coords[coords$pcs_req == pcs_i, c("cellID","k_near","x","y")]
  dd <- merge(cc, qcv, by = "cellID", sort = FALSE); dd$value <- winsorize(dd$log10_depth)
  dd$k_f <- factor(dd$k_near, levels = sort(unique(dd$k_near)))
  print(
    ggplot(dd, aes(x, y, color = value)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      facet_wrap(~k_f, nrow = 1, labeller = labeller(k_f = function(z) paste0("k=", z))) +
      coord_fixed(xlim = c(-1,1), ylim = c(-1,1)) +
      scale_color_viridis_c(name = "log10_depth", option = "C") + base_theme +
      labs(title = paste0(prefix, " -- pcs=", pcs_i, " : UMAP colored by log10_depth across k_near"),
           subtitle = "if the split tracks this gradient = technical continuum, not biology", x = NULL, y = NULL)
  )
}

dev.off()
message(" - wrote PDF: ", pdf_path)
message(" - DONE. Page 1 = n_clusters heatmap (k x res x pcs); page 2 = n_unlabeled; then per-pcs UMAP grids + depth overlays.")
