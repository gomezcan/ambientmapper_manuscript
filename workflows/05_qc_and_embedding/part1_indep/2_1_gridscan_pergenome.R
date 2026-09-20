#!/usr/bin/env Rscript
###############################################################################
## 2_1_gridscan_pergenome.R  --  per-genome UMAP parameter scan (Part 2, Step 1a)
##
## Fork of 2_0_2_umap_grid_scan.R for the SPLIT-genome (independent-mapping)
## objects, where each Socrates object holds a SINGLE REFERENCE (B73v5 OR TAIR10).
## Two reasons the parent script cannot be used as-is here:
##   1. The cross-species "genome mixing" objective is DEGENERATE with one genome.
##   2. The per-genome minDepth200 v6 metadata has NO Genome/species column, so the
##      parent would hard-stop() at its `Genome` requirement.
## This variant therefore:
##   - drops the Genome requirement + all genome_mixing / genome_mixing_centered terms
##   - scores combos by:   score = knn_preservation - 0.5 * qc_abs_cor_max
##   - colours the diagnostic UMAP by log10nSites (a QC covariate) instead of Genome
##
## ---------------------------------------------------------------------------
## MIXING RESTORED AS A DIAGNOSTIC (indep arm)
##
## Reason (1) above is TRUE for the PLATE arm and FALSE for the INDEP arm, and this
## script is used by both. An indep object maps EVERY barcode to one reference, so it
## holds BOTH plate libraries at once -- mixing is not degenerate here, it is the
## genuine-vs-contaminant question:
##     SM2_B73v5  : B73-plate genuine (16,768) vs At-plate contaminant (8,775)
##     SM2_TAIR10 : At-plate  genuine   (597)  vs B73-plate contaminant (2,340)
## The label comes from the cellID suffix (`-SM2_At_*` / `-SM2_B73_*`), verified
## to be present on 100% of barcodes in both v6 files. When no cellID
## carries a suffix (the plate arm, which is single-library by construction), the
## columns are emitted as NA and nothing else changes -- so this is safe for both arms.
##
## REPORTED, NEVER SCORED -- and that is deliberate. `score` is left exactly as it
##   was: knn_preservation - 0.5 * qc_abs_cor_max. Library mixing is the RESULT this
##   arm reports (Fig 5 Part 1 panel D), so folding it into the objective would tune
##   the embedding to the answer -- the same failure mode as letting treatments shape
##   the space they are evaluated in (cf. the Pre-primary consensus peak set, and the
##   `nla` trap). Its purpose here is ROBUSTNESS: if the mixing verdict swings across
##   the grid, the verdict is a parameter artifact and must be reported as such.
##   NOTE the parent 2_0_2 DOES score on mixing (`mix_penalty`); do not copy that here.
##
## obs/exp, never raw. The random-intermingling baseline is 2p(1-p) at minority
##   fraction p, and cleaning changes p enormously (it removes the minority by
##   construction), so a raw drop is mostly the baseline collapsing mechanically.
##   Measured on the combined arm: an 8.7x raw "improvement" was 4.4x composition x
##   2.0x structure. obs/exp = raw / 2p(1-p); 1.0 = indistinguishable from random.
##   Both are emitted (`lib_mixing_raw`, `lib_mixing_obsexp`) plus `p_minor_library`
##   so the correction is auditable rather than implicit.
##   Formula verified against the Fig-5 Part-1 numbers: B73v5 Pre
##   raw 0.391, p=0.3449 -> obs/exp 0.865, matching the reported 0.865 exactly.
##
## SPACE: mixing is computed in the UMAP embedding (like the parent 2_0_2), because
##   what this grid scans IS the UMAP. Fig 5 Part 1 reports mixing in the SVD space
##   (`reduced_dimensions_v7`). Same metric, different space -- the values are NOT
##   interchangeable and should not be quoted against each other. Use these for
##   across-grid stability; quote Fig 5's for the paper.
## Everything else is kept byte-identical to 2_0_2 (fixed TF-IDF + SVD once, the
## pcs x k_near x min_dist grid, the min_c arg + PCA-cap logic) so the chosen
## embedding transfers directly into 3_0_1 / 3_0_0.
##
## Optimize on the PRE (superset) object; the winning (pcs, k_near, min_dist) is
## frozen and applied identically to Pre AND Post-wd downstream (plan decision B).
##
## Usage:
##   Rscript 2_1_gridscan_pergenome.R <soc_rds> <meta_tsv> <outdir> [seed] [min_c]
##     min_c: numeric feature/cell floor for cleanData, or "NA" -> data-driven 250 floor.
###############################################################################

suppressPackageStartupMessages({
  library(Socrates)
  library(Matrix)
  library(FNN)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 5) {
  stop("Usage: Rscript 2_1_gridscan_pergenome.R <soc_rds> <meta_tsv> <outdir> [seed] [min_c]")
}

soc_rds  <- args[1]
meta_tsv <- args[2]
outdir   <- args[3]
seed      <- if (length(args) >= 4) as.integer(args[4]) else 1L
min_c_arg <- if (length(args) >= 5) as.numeric(args[5]) else NA_real_  # NA -> data-driven 250 floor (mirror 3_0_0)
mc_tag    <- if (is.na(min_c_arg)) "" else paste0(".minc_", as.integer(min_c_arg))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

# -------------------------
# TFIDF (identical to 2_0_2 / 3_0_0)
# -------------------------
tfidf <- function(obj,
                  frequencies = TRUE,
                  log_scale_tf = TRUE,
                  scale_factor = 10000,
                  doL2 = FALSE,
                  slotName = "residuals") {
  bmat <- obj$counts
  .safe_tfidf <- function(tf, idf, block_size = 2000e6) {
    tryCatch({
      tf * idf
    }, error = function(e) {
      options(DelayedArray.block.size = block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      tf <- DelayedArray(tf)
      idf <- as.matrix(idf)
      tf * idf
    })
  }
  if (frequencies) {
    tf <- t(t(bmat) / Matrix::colSums(bmat))
  } else {
    tf <- bmat
  }
  if (log_scale_tf) {
    tf@x <- log1p(tf@x * (if (frequencies) scale_factor else 1))
  }
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
  return(obj)
}

# -------------------------
# Diagnostics
# -------------------------

# Plate library from the cellID suffix. This is the cell's PHYSICAL ORIGIN, and on the
# indep arm it is NOT the reference: an At-plate barcode mapped to B73v5 is a
# contaminant. Returns NA for any cellID without a recognisable suffix, which is how
# the plate arm (single-library objects) silently opts out.
plate_library <- function(cell_ids) {
  lib <- rep(NA_character_, length(cell_ids))
  lib[grepl("-SM2_At_",  cell_ids, fixed = TRUE)] <- "At"
  lib[grepl("-SM2_B73_", cell_ids, fixed = TRUE)] <- "B73"
  lib
}

# kNN library mixing, raw and composition-corrected.
#   raw    = mean fraction of a cell's k nearest neighbours carrying the OTHER library
#   obsexp = raw / 2p(1-p), p = minority library fraction; 1.0 = random intermingling
# Returns all-NA when the object is single-library or too small to support k neighbours.
lib_mixing <- function(emb_umap, lib, k = 30) {
  ok <- !is.na(lib) & is.finite(emb_umap[, 1]) & is.finite(emb_umap[, 2])
  if (sum(ok) < k + 1L || length(unique(lib[ok])) < 2L) {
    return(list(raw = NA_real_, obsexp = NA_real_, p_minor = NA_real_, n_labelled = sum(ok)))
  }
  emb <- emb_umap[ok, , drop = FALSE]
  l   <- lib[ok]
  nn  <- FNN::get.knn(emb, k = k)$nn.index
  raw <- mean(vapply(seq_len(nrow(nn)), function(i) mean(l[nn[i, ]] != l[i]), numeric(1)))
  p   <- min(table(l)) / length(l)
  exp_rand <- 2 * p * (1 - p)
  list(raw = raw,
       obsexp = if (is.finite(exp_rand) && exp_rand > 0) raw / exp_rand else NA_real_,
       p_minor = as.numeric(p),
       n_labelled = sum(ok))
}

knn_preservation <- function(emb_pca, emb_umap, k = 30) {
  nn_pca  <- FNN::get.knn(emb_pca,  k = k)$nn.index
  nn_umap <- FNN::get.knn(emb_umap, k = k)$nn.index
  mean(vapply(seq_len(nrow(nn_pca)), function(i) {
    length(intersect(nn_pca[i, ], nn_umap[i, ])) / k
  }, numeric(1)))
}

qc_cor_summary <- function(emb_umap, meta, qc_cols = c("log10nSites", "pOrg")) {
  out <- list()
  for (qc in qc_cols) {
    nm1 <- paste0("rho_umap1_", qc)
    nm2 <- paste0("rho_umap2_", qc)
    if (!qc %in% colnames(meta)) {
      out[[nm1]] <- NA_real_
      out[[nm2]] <- NA_real_
      next
    }
    x <- meta[[qc]]
    ok <- is.finite(x) & is.finite(emb_umap[, 1]) & is.finite(emb_umap[, 2])
    if (sum(ok) < 50) {
      out[[nm1]] <- NA_real_
      out[[nm2]] <- NA_real_
    } else {
      out[[nm1]] <- suppressWarnings(cor(emb_umap[ok, 1], x[ok], method = "spearman"))
      out[[nm2]] <- suppressWarnings(cor(emb_umap[ok, 2], x[ok], method = "spearman"))
    }
  }
  as.data.frame(out)
}

# QC-coloured UMAP panel (replaces the Genome-coloured panel of 2_0_2)
plot_umap_qc <- function(umap, meta, title, subtitle,
                         qc_col = "log10nSites", point_size = 0.12, alpha = 0.6) {
  df <- data.frame(
    umap1 = umap[, 1],
    umap2 = umap[, 2],
    qc    = if (qc_col %in% colnames(meta)) meta[[qc_col]] else NA_real_
  )
  ggplot(df, aes(umap1, umap2, color = qc)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_viridis_c(name = qc_col, na.value = "grey70") +
    theme_bw(base_size = 11) +
    labs(title = title, subtitle = subtitle, x = "umap1", y = "umap2") +
    theme(legend.position = "right")
}

# -------------------------
# Load + align
# -------------------------
message("Reading Socrates object: ", soc_rds)
obj <- readRDS(soc_rds)

message("Reading metadata: ", meta_tsv)
meta <- read.table(meta_tsv, sep = "\t", header = TRUE, check.names = FALSE)
if ("...1" %in% colnames(meta)) meta <- meta[, colnames(meta) != "...1", drop = FALSE]
stopifnot("cellID" %in% colnames(meta))
rownames(meta) <- meta$cellID

shared_cells <- intersect(colnames(obj$counts), rownames(meta))
if (length(shared_cells) == 0) stop("No shared cells between obj$counts and metadata rownames.")
obj$counts <- obj$counts[, shared_cells, drop = FALSE]
obj$meta   <- meta[shared_cells, , drop = FALSE]
message("Aligned cells (counts n meta): ", length(shared_cells))

if (!"log10nSites" %in% colnames(obj$meta) && "nSites" %in% colnames(obj$meta)) {
  obj$meta$log10nSites <- log10(obj$meta$nSites + 1)
}

# -------------------------
# Fix TF-IDF + SVD once (identical to 2_0_2)
# -------------------------
set.seed(seed)

cell.counts <- log10(Matrix::colSums(obj$counts))
cell.counts.z <- as.numeric(scale(cell.counts))
mask <- cell.counts.z[cell.counts.z < -0.5]
cell.counts.threshold <- if (is.na(min_c_arg)) max(c(10^mask, 250), na.rm = TRUE) else min_c_arg

min_t <- 0.001
message("cleanData once: min.c=", signif(cell.counts.threshold, 4), " min.t=", min_t)
soc <- cleanData(obj, min.c = cell.counts.threshold, min.t = min_t, max.t = 0, verbose = TRUE)

message("TFIDF once")
soc <- tfidf(soc, doL2 = TRUE)

pcs_grid <- c(20, 25, 30, 40, 50)
max_pcs  <- max(pcs_grid) + 1

number.sites <- ceiling(nrow(soc$counts) * 0.50)
message("reduceDims once (SVD): n.pcs=", max_pcs, " num.var=", number.sites)

soc <- reduceDims(soc,
                  method = "SVD",
                  n.pcs  = max_pcs,
                  cor.max = 0.6,
                  num.var = number.sites,
                  verbose = TRUE,
                  scaleVar = TRUE,
                  doSTD = FALSE,
                  doL1  = FALSE,
                  doL2  = TRUE,
                  refit_residuals = FALSE)

# IMPORTANT: make sure PCA rows align to soc$meta rows (cells)
pca_full <- soc$PCA[rownames(soc$meta), , drop = FALSE]

# -------------------------
# Grid to scan (small) -- cap pcs to what reduceDims actually returned
# -------------------------
avail_pcs <- ncol(pca_full)
pcs_grid  <- pcs_grid[pcs_grid <= avail_pcs]
if (length(pcs_grid) == 0) pcs_grid <- avail_pcs
message("Available PCs after reduceDims: ", avail_pcs, " -> scanning pcs = ", paste(pcs_grid, collapse = ", "))

k_grid        <- c(15, 20, 30)
min_dist_grid <- c(0.05, 0.15, 0.30)

grid <- expand.grid(
  pcs      = pcs_grid,
  k_near   = k_grid,
  min_dist = min_dist_grid,
  stringsAsFactors = FALSE
)
message("Grid size: ", nrow(grid), " combinations")

# Diagnostics k must be <= smallest k_near
diag_k <- min(k_grid)

# -------------------------
# Run scan
# -------------------------
metrics_out <- vector("list", nrow(grid))

plots_pdf <- file.path(outdir, paste0("plots/UMAP_grid_scan", mc_tag, ".QC.pdf"))
pdf(plots_pdf, width = 6, height = 5)

for (i in seq_len(nrow(grid))) {
  pcs_i      <- grid$pcs[i]
  k_near_i   <- grid$k_near[i]
  min_dist_i <- grid$min_dist[i]

  tag <- paste0("pcs=", pcs_i, " | k_near=", k_near_i, " | min_dist=", min_dist_i)

  soc_i <- soc
  soc_i$PCA <- pca_full[, seq_len(pcs_i), drop = FALSE]

  soc_i <- projectUMAP(
    soc_i,
    m.dist = min_dist_i,
    k.near = k_near_i,
    metric = "cosine",
    svd_slotName  = "PCA",
    umap_slotName = "UMAP",
    verbose = FALSE,
    seed = seed
  )

  umap <- as.matrix(soc_i$UMAP[rownames(soc_i$meta), c("umap1", "umap2"), drop = FALSE])

  pres   <- knn_preservation(emb_pca = soc_i$PCA, emb_umap = umap, k = diag_k)
  qc_cor <- qc_cor_summary(umap, soc_i$meta, qc_cols = c("log10nSites", "pOrg"))
  cor_abs_max <- max(abs(unlist(qc_cor)), na.rm = TRUE)

  # DIAGNOSTIC ONLY -- deliberately absent from `score` (see header).
  mixd <- lib_mixing(umap, plate_library(rownames(soc_i$meta)), k = diag_k)

  metrics_out[[i]] <- cbind(
    data.frame(
      pcs = pcs_i,
      k_near = k_near_i,
      min_dist = min_dist_i,
      diag_k = diag_k,
      knn_preservation = pres,
      qc_abs_cor_max = cor_abs_max,
      lib_mixing_raw = mixd$raw,
      lib_mixing_obsexp = mixd$obsexp,
      p_minor_library = mixd$p_minor,
      n_labelled_library = mixd$n_labelled,
      n_cells = nrow(umap),
      stringsAsFactors = FALSE
    ),
    qc_cor
  )

  subtitle <- paste0(
    "kNN_pres=", sprintf("%.3f", pres),
    " | max|rho|=", sprintf("%.3f", cor_abs_max),
    if (is.finite(mixd$obsexp)) paste0(" | libMix obs/exp=", sprintf("%.3f", mixd$obsexp),
                                       " (raw ", sprintf("%.3f", mixd$raw),
                                       ", p=", sprintf("%.3f", mixd$p_minor), ")") else "",
    " | diag_k=", diag_k
  )

  print(plot_umap_qc(umap, soc_i$meta, title = tag, subtitle = subtitle))

  message("Done [", i, "/", nrow(grid), "]: ", tag,
          " | pres=", sprintf("%.3f", pres),
          " max|rho|=", sprintf("%.3f", cor_abs_max),
          if (is.finite(mixd$obsexp)) paste0(" libMix_obsexp=", sprintf("%.3f", mixd$obsexp)) else "")

  rm(soc_i)
}

dev.off()

metrics_df <- rbindlist(metrics_out, use.names = TRUE, fill = TRUE)
metrics_tsv <- file.path(outdir, paste0("UMAP_grid_scan", mc_tag, ".metrics.tsv"))
fwrite(metrics_df, metrics_tsv, sep = "\t", quote = FALSE)

# -------------------------
# Rank / choose "best"
# -------------------------
# Single-genome objective: reward neighbourhood preservation, penalise QC leakage.
# `lib_mixing_*` is INTENTIONALLY NOT IN THE SCORE -- it is the reported result on the
# indep arm, so optimising against it would tune the embedding to its own answer.
# Use it to check the verdict's STABILITY across the grid, not to pick the winner.
metrics_df <- as.data.frame(fread(metrics_tsv, sep = "\t")) %>%
  mutate(
    score = knn_preservation - 0.5 * qc_abs_cor_max
  ) %>%
  arrange(desc(score))

best_row <- metrics_df[1, , drop = FALSE]
best_tsv <- file.path(outdir, paste0("UMAP_grid_scan", mc_tag, ".best.tsv"))
fwrite(best_row, best_tsv, sep = "\t", quote = FALSE)

message("Wrote metrics:   ", metrics_tsv)
message("Wrote UMAP PDF:  ", plots_pdf)
message("Best combo TSV:  ", best_tsv)
message("Best parameters:")
print(best_row)

# knn_preservation vs QC-leakage diagnostic scatter (replaces the genome-mixing scatter)
best <- metrics_df[1, ]
bad  <- metrics_df %>% arrange(knn_preservation) %>% slice(1)

scatter_pdf <- file.path(outdir, paste0("plots/knn_preservation_grid_scan", mc_tag, ".QCleak.pdf"))
pdf(scatter_pdf, width = 9, height = 4)
print(
  ggplot(metrics_df, aes(x = knn_preservation, y = qc_abs_cor_max)) +
    geom_point(aes(color = score), size = 3) +
    scale_color_viridis_c(name = "score") +
    geom_point(data = best, shape = 21, size = 5, stroke = 1.2, fill = "white", color = "black") +
    geom_point(data = bad,  shape = 21, size = 5, stroke = 1.2, fill = "white", color = "red") +
    geom_text(data = best,
              aes(label = paste0("Selected\npcs=", pcs, ", k=", k_near, ", min_dist=", min_dist)),
              hjust = -0.1, vjust = -0.5, size = 3) +
    theme_bw() +
    labs(x = "kNN preservation (PCA -> UMAP)  [higher better]",
         y = "max |QC corr| (UMAP vs log10nSites, pOrg)  [lower better]",
         title = "Per-genome UMAP robustness across parameter grid",
         subtitle = "black = selected (max score); red = worst kNN preservation")
)
dev.off()
message("Wrote scatter:   ", scatter_pdf)

# ---------------------------------------------------------------------------
# Library-mixing companion scatter -- same layout as the combined arm's Supp Fig 3
# (mixing on y, kNN preservation on x, coloured by QC leakage). Emitted only when the
# object actually carries two plate libraries, i.e. the indep arm.
#
# y is obs/exp, NOT raw. 1.0 = indistinguishable from random intermingling. The raw
#   value is composition-driven and a Pre-vs-Post drop in it is mostly the 2p(1-p)
#   baseline collapsing as the contaminant becomes rare -- see the header.
# Mixing is NOT in `score`; the selected point is still chosen on
#   knn_preservation - 0.5*qc_leak. This panel shows whether the mixing verdict is
#   STABLE across the grid. A tight horizontal band = the verdict is a property of the
#   data; a wide vertical spread = it is a property of the parameters, and must be
#   reported that way.
# ---------------------------------------------------------------------------
if (any(is.finite(metrics_df$lib_mixing_obsexp))) {
  mix_pdf <- file.path(outdir, paste0("knn_preservation_grid_scan", mc_tag, ".LibMixing.pdf"))
  mix_pdf <- file.path(outdir, "plots", basename(mix_pdf))
  rng <- range(metrics_df$lib_mixing_obsexp, na.rm = TRUE)
  message(sprintf("Library mixing obs/exp across grid: min=%.3f max=%.3f spread=%.3f",
                  rng[1], rng[2], diff(rng)))
  pdf(mix_pdf, width = 9, height = 4)
  print(
    ggplot(metrics_df, aes(x = knn_preservation, y = lib_mixing_obsexp)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
      annotate("text", x = min(metrics_df$knn_preservation, na.rm = TRUE), y = 1,
               label = "random intermingling", hjust = 0, vjust = -0.6, size = 3, colour = "grey30") +
      geom_point(aes(color = qc_abs_cor_max), size = 3) +
      scale_color_viridis_c(name = "max |QC corr|") +
      geom_point(data = best, shape = 21, size = 5, stroke = 1.2, fill = "white", color = "black") +
      geom_text(data = best,
                aes(label = paste0("Selected\npcs=", pcs, ", k=", k_near, ", min_dist=", min_dist)),
                hjust = -0.1, vjust = -0.5, size = 3) +
      theme_bw() +
      labs(x = "kNN preservation (PCA -> UMAP)  [higher better]",
           y = "Plate-library mixing, observed / expected",
           title = "Library mixing across parameter grid (indep arm)",
           subtitle = paste0("obs/exp; 1.0 = random intermingling. Mixing is NOT in the selection score. ",
                             "Grid spread = ", sprintf("%.3f", diff(rng))))
  )
  dev.off()
  message("Wrote mixing scatter: ", mix_pdf)
} else {
  message("Library mixing: not emitted (single-library object -- expected on the plate arm).")
}
