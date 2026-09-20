#!/usr/bin/env Rscript
###############################################################################
## 2_2_resolution_scan.R  --  per-genome pcs x resolution stability scan
##                              (Part 2, Step 1b; MARKER-FREE)
##
## Step 1a (2_0_2b) could not principledly pick pcs (knn_preservation is mechanically
## monotone in pcs for B73v5 and flat for TAIR10), so 1b ARBITRATES pcs jointly with
## resolution using the non-circular objective: bootstrap CLUSTER STABILITY.
##
## Fix k_near + min_dist (min_dist is UMAP-viz-only -- clustering is on the SVD graph).
## Scan pcs x resolution. For each (pcs, res) score by:
##   - granularity structure : n_clusters, min/median cluster size, small-cluster frac
##   - kNN cluster purity     : mean frac of a cell's PCA-space kNN sharing its label
##   - stability (ARI)        : mean adjusted Rand vs full over B subsample re-clusters
##                              (embedding FIXED per pcs, cells resampled)
## Marker coherence is intentionally OUT (design decision: 1b marker-free).
##
## cleanData + tfidf are computed ONCE (shared); reduceDims + projectUMAP + the whole
## resolution scan are redone PER pcs, so each pcs uses exactly the Step-2 (3_0_0) path.
## Clustering call is IDENTICAL to 3_0_0 step (7):
##   Socrates::callClusters(cl.method=4, m.clst=50, e.thresh=3, threshold=3, cleanCluster=F).
## Pick the (pcs, resolution) at the stability knee with sane granularity, then freeze the
## full (pcs, k_near, min_dist, resolution, min.c) and apply it IDENTICALLY to Pre + Post-wd.
##
## Usage:
##   Rscript 2_2_resolution_scan.R <soc_rds> <meta_tsv> <outdir> <prefix> \
##           <pcs_list> <k_near> <min_dist> [min_c=NA] [seed=1] [n_boot=10] [boot_frac=0.8]
##     pcs_list: comma-separated, e.g. "20,30,50"
##     min_c: numeric feature/cell floor for cleanData, or "NA" -> data-driven 250 floor.
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
  stop("Usage: Rscript 2_2_resolution_scan.R <soc_rds> <meta_tsv> <outdir> <prefix> <pcs_list> <k_near> <min_dist> [min_c=NA] [seed=1] [n_boot=10] [boot_frac=0.8] [res_grid]")
}
soc_rds   <- args[1]
meta_tsv  <- args[2]
outdir    <- args[3]
prefix    <- args[4]
pcs_vec   <- as.integer(strsplit(args[5], ",")[[1]])
k_near    <- as.integer(args[6])
min_dist  <- as.numeric(args[7])
min_c_arg <- if (length(args) >= 8)  as.numeric(args[8])  else NA_real_
seed      <- if (length(args) >= 9)  as.integer(args[9])  else 1L
n_boot    <- if (length(args) >= 10) as.integer(args[10]) else 10L
boot_frac <- if (length(args) >= 11) as.numeric(args[11]) else 0.8

# plan-specified default grid; arg 12 overrides. Small-n arms (e.g. plate At ~800 cells,
# which fragments into 94 clusters at res=1) need to sweep BELOW 0.3.
RES_GRID <- if (length(args) >= 12) as.numeric(strsplit(args[12], ",")[[1]]) else c(0.3, 0.5, 0.8, 1.0, 1.5)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(outdir, paste0(prefix, suffix))

message(" - 3_0_1 pcs x resolution scan | prefix=", prefix,
        " | pcs={", paste(pcs_vec, collapse = ","), "} k_near=", k_near, " min_dist=", min_dist,
        " | min_c=", ifelse(is.na(min_c_arg), "NA(data-driven)", min_c_arg),
        " | n_boot=", n_boot, " boot_frac=", boot_frac)

# -------------------------
# helpers
# -------------------------
tfidf <- function(obj,
                  frequencies = TRUE,
                  log_scale_tf = TRUE,
                  scale_factor = 10000,
                  doL2 = FALSE,
                  slotName = "residuals") {
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

# base-R adjusted Rand index (no mclust/aricode dependency)
adj_rand <- function(a, b) {
  ct <- table(a, b)
  csum <- function(x) sum(choose(x, 2))
  n  <- length(a)
  ii <- csum(as.vector(ct))
  ai <- csum(rowSums(ct))
  bj <- csum(colSums(ct))
  expected <- ai * bj / choose(n, 2)
  maxi     <- (ai + bj) / 2
  denom    <- maxi - expected
  if (denom == 0) return(1)      # both partitions trivial -> perfectly "agree"
  (ii - expected) / denom
}

# subset a Socrates object to a cell set, keeping all slot dims consistent
subset_soc <- function(soc, cells) {
  s <- soc
  s$counts <- soc$counts[, cells, drop = FALSE]
  s$PCA    <- soc$PCA[cells, , drop = FALSE]
  if (!is.null(soc$UMAP))      s$UMAP      <- soc$UMAP[cells, , drop = FALSE]
  if (!is.null(soc$residuals)) s$residuals <- soc$residuals[, cells, drop = FALSE]
  s$meta <- soc$meta[cells, , drop = FALSE]
  s
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

# mean fraction of a cell's PCA-space kNN that share its cluster label
knn_cluster_purity <- function(pca, labels, k) {
  nn  <- FNN::get.knn(pca, k = k)$nn.index
  lab <- labels[rownames(pca)]                 # align to pca row order
  mean(vapply(seq_len(nrow(nn)), function(i) mean(lab[nn[i, ]] == lab[i]), numeric(1)))
}

# ---------------------------------------------------------------------------
# PER-CLUSTER CONTAMINATION PROFILE  (indep arm)
#
# On the indep arm every object maps BOTH plate libraries onto ONE reference, so the
# cell's physical origin is recoverable from the cellID suffix and the contamination
# question can be asked of each cluster WITHOUT any marker or reference to Fig 3:
#     SM2_B73v5  -> B73-plate genuine, At-plate  contaminant
#     SM2_TAIR10 -> At-plate  genuine, B73-plate contaminant
# Reproduced at the frozen config and matches Fig 5 Part 1 panel E exactly
# (B73v5 Pre cl5 = 80.6% contaminant vs a 34.5% object baseline; post-clean max 25.0%).
#
# REPORT ENRICHMENT, NOT THE RAW FRACTION, when comparing stages. The object baseline
#   moves enormously between Pre and Post (34.5% -> 2.9% on B73v5), so a raw per-cluster
#   percentage is not comparable across stages: post-clean's worst cluster (25.0%) is
#   8.6x ITS baseline while Pre's worst (80.6%) is only 2.3x its own. Same composition
#   trap as the 2p(1-p) one on the mixing panel. `contam_baseline` is emitted alongside
#   so every ratio here is auditable.
#
# DIAGNOSTIC, NOT AN OBJECTIVE. Like lib_mixing in 2_0_2b, these columns are reported
#   and deliberately kept OUT of the recommendation gate below. Selecting the config that
#   best isolates contamination would tune the partition to the answer it is meant to
#   measure. Use them to check that the verdict is stable across (pcs, res).
#
# `contam_captured_frac` is the isolation measure that matters: of ALL contaminant cells,
# what fraction sits in clusters that are enriched >= enrich_mult x baseline? High = the
# configuration concentrates contamination into identifiable clusters; low = it is smeared
# uniformly and no partition will separate it.
# ---------------------------------------------------------------------------
plate_library <- function(cell_ids) {
  lib <- rep(NA_character_, length(cell_ids))
  lib[grepl("-SM2_At_",  cell_ids, fixed = TRUE)] <- "At"
  lib[grepl("-SM2_B73_", cell_ids, fixed = TRUE)] <- "B73"
  lib
}

cluster_contamination <- function(labels, cell_ids, genuine, min_size = 50, enrich_mult = 2) {
  na_out <- data.frame(contam_baseline = NA_real_, contam_cluster_max = NA_real_,
                       contam_enrich_max = NA_real_, n_clusters_contam_enriched = NA_integer_,
                       contam_captured_frac = NA_real_)
  if (is.na(genuine)) return(na_out)
  lib <- plate_library(cell_ids)
  ok  <- !is.na(lib)
  # single-library object (the plate arm) -> nothing to measure, emit NA and move on
  if (sum(ok) < 2L || length(unique(lib[ok])) < 2L) return(na_out)

  is_con   <- (lib != genuine)[ok]
  lab      <- labels[ok]
  baseline <- mean(is_con)
  n_by <- tapply(is_con, lab, length)
  c_by <- tapply(is_con, lab, sum)
  frac <- c_by / n_by
  big  <- n_by >= min_size
  if (!any(big) || !is.finite(baseline) || baseline <= 0) {
    na_out$contam_baseline <- round(baseline, 4)
    return(na_out)
  }
  cl_max <- max(frac[big])

  # The enrichment THRESHOLD saturates at high baseline: a fraction cannot exceed 1,
  #   so `enrich_mult * baseline >= 1` makes the threshold unreachable by construction.
  #   Measured case: TAIR10 baseline is 80.0% (Pre) / 68.7% (wd), so at enrich_mult=2 the
  #   bar is 160% / 137% and the count is vacuously 0 -- which reads as "nothing is
  #   enriched" when the truth is "the question is undefined here". Emit NA instead; the
  #   RATIO (contam_enrich_max) stays valid and carries the signal on its own
  #   (TAIR10 measured 1.11x/1.14x = essentially uniform contamination, no separable cluster).
  reachable <- (enrich_mult * baseline) < 1
  if (reachable) {
    enr        <- big & (frac >= enrich_mult * baseline)
    n_enriched <- as.integer(sum(enr))
    captured   <- if (sum(c_by) > 0) round(sum(c_by[enr]) / sum(c_by), 4) else NA_real_
  } else {
    n_enriched <- NA_integer_
    captured   <- NA_real_
  }

  data.frame(
    contam_baseline            = round(baseline, 4),
    contam_cluster_max         = round(cl_max, 4),
    contam_enrich_max          = round(cl_max / baseline, 3),
    n_clusters_contam_enriched = n_enriched,
    contam_captured_frac       = captured
  )
}

# -------------------------
# Load + align (mirror 3_0_0 sec 1)
# -------------------------
obj <- readRDS(soc_rds)
meta <- read.table(meta_tsv, header = TRUE, sep = "\t", check.names = FALSE)
if ("...1" %in% colnames(meta)) meta <- meta[, colnames(meta) != "...1", drop = FALSE]
stopifnot("cellID" %in% colnames(meta))
rownames(meta) <- meta$cellID

shared <- intersect(colnames(obj$counts), rownames(meta))
if (length(shared) == 0) stop("No shared cells between obj$counts and metadata.")
obj$counts <- obj$counts[, shared, drop = FALSE]
obj$meta   <- meta[shared, , drop = FALSE]
message(" - aligned cells (counts n v6 meta): ", length(shared))

# Which plate library is GENUINE in this object, from the reference in the prefix.
# NA -> the per-cluster contamination columns are emitted as NA (plate arm / unknown ref),
# which is a silent no-op rather than an error.
genuine_lib <- if (grepl("TAIR10", prefix)) "At" else if (grepl("B73v5", prefix)) "B73" else NA_character_
if (is.na(genuine_lib)) {
  message(" - genuine plate library: UNKNOWN for prefix '", prefix, "' -> contamination columns will be NA")
} else {
  message(" - genuine plate library: SM2_", genuine_lib,
          " (everything else in this object is contaminant)")
}

# -------------------------
# cleanData + tfidf ONCE (shared across pcs); reduceDims/UMAP/scan PER pcs
# -------------------------
set.seed(seed)
cell.counts   <- log10(Matrix::colSums(obj$counts))
cell.counts.z <- as.numeric(scale(cell.counts))
mask <- cell.counts.z[cell.counts.z < -0.5]
min_c <- if (is.na(min_c_arg)) max(c(10^mask, 250), na.rm = TRUE) else min_c_arg
message(" - cleanData min.c = ", signif(min_c, 4))

soc_norm <- cleanData(obj, min.c = min_c, min.t = 0.001, max.t = 0, verbose = TRUE)
soc_norm <- tfidf(soc_norm, doL2 = TRUE)
number.sites <- ceiling(nrow(soc_norm$counts) * 0.5)
message(" - normalized once: ", ncol(soc_norm$counts), " cells x ", nrow(soc_norm$counts), " sites; num.var=", number.sites)

# -------------------------
# Scan pcs x resolution
# -------------------------
rows <- list()
for (pcs_i in pcs_vec) {
  message(" ===== pcs = ", pcs_i, " =====")
  soc <- reduceDims(soc_norm, method = "SVD", n.pcs = pcs_i, cor.max = 0.6, num.var = number.sites,
                    verbose = TRUE, scaleVar = TRUE, doSTD = FALSE, doL1 = FALSE, doL2 = TRUE,
                    refit_residuals = FALSE)
  soc <- projectUMAP(soc, verbose = FALSE, k.near = k_near, m.dist = min_dist)
  pca     <- soc$PCA[rownames(soc$meta), , drop = FALSE]
  pcs_act <- ncol(pca)                              # cor.max may drop some (this is the Step-2 behaviour too)
  cellIDs <- rownames(pca)
  N       <- length(cellIDs)
  if (pcs_act < pcs_i) message("   ! reduceDims returned ", pcs_act, " PCs (< requested ", pcs_i, ") after cor.max drop")

  for (j in seq_along(RES_GRID)) {
    res <- RES_GRID[j]
    lab_full <- cluster_labels(soc, res, k_near)
    lab_full <- lab_full[cellIDs]                   # align to embedding order
    keep     <- cellIDs[!is.na(lab_full)]           # guard: any cell callClusters left unlabeled
    if (length(keep) < length(cellIDs))
      message("   ! ", length(cellIDs) - length(keep), " cells unlabeled at pcs=", pcs_i, " res=", res)
    lab_k    <- lab_full[keep]
    sizes    <- as.integer(table(lab_k))
    n_clust  <- length(sizes)

    purity <- knn_cluster_purity(pca[keep, , drop = FALSE], lab_k, k_near)

    # stability: B subsample re-clusters on the FIXED (this-pcs) embedding, ARI vs full
    aris <- rep(NA_real_, n_boot)
    n_take <- max(50L, floor(boot_frac * length(keep)))
    for (b in seq_len(n_boot)) {
      set.seed(seed * 100000L + pcs_i * 1000L + j * 100L + b)   # deterministic, distinct
      cells_b <- sample(keep, n_take)
      lab_b <- tryCatch(cluster_labels(subset_soc(soc, cells_b), res, k_near),
                        error = function(e) { message("     ! boot ", b, " skipped: ", conditionMessage(e)); NULL })
      if (is.null(lab_b)) next
      common <- intersect(names(lab_b), keep)
      if (length(common) < 50) next
      aris[b] <- adj_rand(lab_k[common], lab_b[common])
    }
    n_ok <- sum(is.finite(aris))

    # per-cluster contamination profile (diagnostic; NOT part of the recommendation gate)
    contam <- cluster_contamination(lab_k, keep, genuine_lib)

    rows[[length(rows) + 1L]] <- cbind(data.frame(
      prefix = prefix, pcs = pcs_act, pcs_requested = pcs_i, resolution = res,
      n_clusters = n_clust, min_size = min(sizes), median_size = as.numeric(median(sizes)),
      n_singletons = sum(sizes == 1), n_small_clusters = sum(sizes < 50),
      frac_cells_small = round(sum(sizes[sizes < 50]) / sum(sizes), 4),
      knn_cluster_purity = round(purity, 4),
      stability_ARI_mean = round(mean(aris, na.rm = TRUE), 4),
      stability_ARI_sd   = round(sd(aris,   na.rm = TRUE), 4),
      n_boot_ok = n_ok, n_cells = N, stringsAsFactors = FALSE
    ), contam)

    message(sprintf("   pcs=%d res=%.2f | clusters=%d min=%d medsize=%.0f | purity=%.3f | ARI=%.3f+/-%.3f (n=%d)%s",
                    pcs_act, res, n_clust, min(sizes), median(sizes), purity,
                    mean(aris, na.rm = TRUE), sd(aris, na.rm = TRUE), n_ok,
                    if (is.finite(contam$contam_enrich_max))
                      sprintf(" | contam base=%.1f%% max=%.1f%% (%.1fx) captured=%.0f%% in %d cl",
                              100 * contam$contam_baseline, 100 * contam$contam_cluster_max,
                              contam$contam_enrich_max, 100 * contam$contam_captured_frac,
                              contam$n_clusters_contam_enriched)
                    else ""))
  }
  rm(soc); gc()
}

scan <- do.call(rbind, rows)
scan_tsv <- op(".resolution_scan.tsv")
write.table(scan, scan_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote scan table: ", scan_tsv)

# -------------------------
# Recommend (pcs, resolution)  [DRAFT -- confirm by eye before freezing]
# Sane gate: >=3 clusters, min_size>=20, <10% of cells in <50-cell clusters.
# Among sane rows within 0.02 ARI of the peak, take the FINEST (most clusters) --
# i.e. the most granular partition that is still near-peak stable across pcs & res;
# tie-break toward FEWER pcs (parsimony).
# -------------------------
sane <- scan[scan$n_clusters >= 3 & scan$min_size >= 20 & scan$frac_cells_small < 0.10, , drop = FALSE]
if (nrow(sane) == 0) { warning("No config passed the sanity gate; falling back to all rows."); sane <- scan }
best_ari <- max(sane$stability_ARI_mean, na.rm = TRUE)
near     <- sane[is.finite(sane$stability_ARI_mean) & sane$stability_ARI_mean >= best_ari - 0.02, , drop = FALSE]
near     <- near[order(-near$n_clusters, near$pcs), , drop = FALSE]
rec      <- near[1, , drop = FALSE]

# also: the single most-stable config per pcs (for the user to eyeball the pcs axis)
per_pcs_best <- do.call(rbind, lapply(split(sane, sane$pcs), function(d) d[which.max(d$stability_ARI_mean), , drop = FALSE]))

message(" - RECOMMENDED (draft): pcs=", rec$pcs, " res=", rec$resolution,
        " (", rec$n_clusters, " clusters, ARI=", rec$stability_ARI_mean, ", purity=", rec$knn_cluster_purity, ")")
message(" - most-stable per pcs:")
print(per_pcs_best[, c("pcs", "resolution", "n_clusters", "stability_ARI_mean", "knn_cluster_purity")])

rec_txt <- op(".resolution_scan.recommendation.txt")
writeLines(c(
  paste0("# 3_0_1 pcs x resolution recommendation (DRAFT) -- ", prefix),
  paste0("# fixed: k_near=", k_near, " min_dist=", min_dist, " min_c=", signif(min_c, 4)),
  "# heuristic: finest (most clusters) config within 0.02 ARI of peak stability, sane granularity, fewer pcs on ties.",
  "# CONFIRM by eye (scan table + PDF) before freezing into <genome>.cluster_config.tsv.",
  paste0("recommended_pcs\t",        rec$pcs),
  paste0("recommended_resolution\t", rec$resolution),
  paste0("n_clusters\t",             rec$n_clusters),
  paste0("stability_ARI_mean\t",     rec$stability_ARI_mean),
  paste0("knn_cluster_purity\t",     rec$knn_cluster_purity),
  "#",
  "# most-stable resolution per pcs:",
  paste0("# pcs=", per_pcs_best$pcs, " res=", per_pcs_best$resolution,
         " nclust=", per_pcs_best$n_clusters, " ARI=", per_pcs_best$stability_ARI_mean)
), rec_txt)

genome <- sub("^(SM2|Clean\\.SM2v2wd|Clean\\.SM2v2)_", "", prefix)   # SM2_B73v5 -> B73v5
cfg <- data.frame(
  genome = genome, prefix = prefix,
  pcs = rec$pcs, k_near = k_near, min_dist = min_dist, min_c = signif(min_c, 6),
  resolution = rec$resolution, n_clusters = rec$n_clusters,
  stability_ARI_mean = rec$stability_ARI_mean, knn_cluster_purity = rec$knn_cluster_purity,
  note = "DRAFT from 3_0_1 (pcs arbitrated); confirm before Step 2", stringsAsFactors = FALSE
)
cfg_tsv <- op(".cluster_config.draft.tsv")
write.table(cfg, cfg_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote draft config: ", cfg_tsv)

# -------------------------
# Diagnostic PDF: stability + granularity vs resolution, one line per pcs
# -------------------------
scan$pcs_f <- factor(scan$pcs)
rec$pcs_f  <- factor(rec$pcs, levels = levels(scan$pcs_f))   # rec was extracted before pcs_f existed; without this the geom_point(data=rec) highlight layers inherit color=pcs_f and error -> empty PDF
pdf_path <- file.path(outdir, "plots", paste0(prefix, ".resolution_scan.pdf"))
pdf(pdf_path, width = 6.5, height = 4.5)
print(
  ggplot(scan, aes(resolution, stability_ARI_mean, color = pcs_f, group = pcs_f)) +
    geom_line() +
    geom_errorbar(aes(ymin = stability_ARI_mean - stability_ARI_sd,
                      ymax = stability_ARI_mean + stability_ARI_sd), width = 0.03, alpha = 0.5) +
    geom_point(size = 2.5) +
    geom_point(data = rec, shape = 21, size = 5, stroke = 1.3, fill = NA, color = "red") +
    scale_color_viridis_d(name = "pcs", end = 0.85) +
    theme_bw(base_size = 11) +
    labs(title = paste0(prefix, " -- cluster stability vs (pcs, resolution)"),
         subtitle = "red ring = recommended (draft); error bars = +/-1 sd bootstrap ARI",
         x = "Leiden resolution", y = "stability (mean bootstrap ARI)")
)
print(
  ggplot(scan, aes(resolution, n_clusters, color = pcs_f, group = pcs_f)) +
    geom_line() + geom_point(size = 2.5) +
    geom_point(data = rec, shape = 21, size = 5, stroke = 1.3, fill = NA, color = "red") +
    scale_color_viridis_d(name = "pcs", end = 0.85) +
    theme_bw(base_size = 11) +
    labs(title = paste0(prefix, " -- cluster count vs (pcs, resolution)"),
         x = "resolution", y = "n_clusters")
)
print(
  ggplot(scan, aes(resolution, knn_cluster_purity, color = pcs_f, group = pcs_f)) +
    geom_line() + geom_point(size = 2.5) +
    scale_color_viridis_d(name = "pcs", end = 0.85) +
    theme_bw(base_size = 11) +
    labs(title = paste0(prefix, " -- kNN cluster purity vs (pcs, resolution)"),
         x = "resolution", y = "kNN cluster purity")
)
dev.off()
message(" - wrote PDF: ", pdf_path)
message(" - DONE")
