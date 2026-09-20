#!/usr/bin/env Rscript
# Fig S3 (panels A to C): UMAP parameter robustness of the At/B73 co-projection (concatenated reference).
# TF-IDF and SVD are fixed once; a grid over pcs x k_near (UMAP n_neighbors) x min_dist (45 combinations)
# is scanned and each combination scored for kNN preservation (SVD vs UMAP neighbourhoods, same k),
# genome mixing (raw and composition-centred) and QC leakage (max |Spearman rho| of UMAP1/2 vs log10nSites, pOrg).
# B, C: genome-coloured UMAPs, one page per combination (plots/UMAP_grid_scan[.minc_N].Genome.pdf).
# A: drawn as shipped by figS3A_replot.R from the metrics table this script writes (UMAP_grid_scan[.minc_N].metrics.tsv).
# Inputs: data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2/step1_integrate/SM2.full.SocObj.rds and
#   .../SM2/step2_metaqc/SM2.FRiP0.2.FULL.minDepth200.stagePreClean.updated_metadata_v4.txt (written by figS2.R)
# Usage:  Rscript analysis/supplementary/figS3.R [soc_rds] [meta_tsv] [outdir] [seed] [min_c]   (HPC, about 30 GB; see figS3.sh)
# Output: <outdir>/UMAP_grid_scan[.minc_N].{metrics,best}.tsv and <outdir>/plots/*.pdf (default outdir figures/supplementary/figS3)

suppressPackageStartupMessages({
  library(Socrates)
  library(Matrix)
  library(FNN)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

# ---- CONFIG ----------------------------------------------------------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis/socrates"
OUTDIR <- "figures/supplementary/figS3"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 5) {
  stop("Usage: Rscript figS3.R [soc_rds] [meta_tsv] [outdir] [seed] [min_c]")
}
soc_rds   <- if (length(args) >= 1) args[1] else file.path(DATA, "SM2/step1_integrate/SM2.full.SocObj.rds")
meta_tsv  <- if (length(args) >= 2) args[2] else file.path(DATA, "SM2/step2_metaqc/SM2.FRiP0.2.FULL.minDepth200.stagePreClean.updated_metadata_v4.txt")
outdir    <- if (length(args) >= 3) args[3] else OUTDIR
seed      <- if (length(args) >= 4) as.integer(args[4]) else 1L
# min_c: minimum reads per cell for cleanData; NA = data-driven floor with a minimum of 250 (as in the clustering pipeline)
min_c_arg <- if (length(args) >= 5) as.numeric(args[5]) else NA_real_
mc_tag    <- if (is.na(min_c_arg)) "" else paste0(".minc_", as.integer(min_c_arg))

dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

# -------------------------
# TFIDF (identical to the clustering pipeline)
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


cleanObj <- function(raw, conf) {
  # Identify shared cells between raw meta and conf object clusters
  shared <- intersect(rownames(raw$meta), rownames(conf$Clusters))
  raw$counts <- raw$counts[, !colnames(raw$counts) %in% shared]
  shared.sites <- intersect(rownames(raw$counts), rownames(conf$residuals))
  raw$counts <- raw$counts[shared.sites, ]
  raw$counts <- raw$counts[, Matrix::colSums(raw$counts) > 0]
  raw$counts <- raw$counts[Matrix::rowSums(raw$counts) > 0, ]
  raw$meta <- raw$meta[colnames(raw$counts), ]
  return(raw)
}

# -------------------------
# Diagnostics
# -------------------------
# kNN preservation: mean fraction of a cell's k nearest neighbours in SVD space that are also
# among its k nearest neighbours in UMAP space
knn_preservation <- function(emb_pca, emb_umap, k = 30) {
  nn_pca  <- FNN::get.knn(emb_pca,  k = k)$nn.index
  nn_umap <- FNN::get.knn(emb_umap, k = k)$nn.index
  mean(vapply(seq_len(nrow(nn_pca)), function(i) {
    length(intersect(nn_pca[i, ], nn_umap[i, ])) / k
  }, numeric(1)))
}

# genome mixing: mean fraction of a cell's k nearest UMAP neighbours that carry the other genome label
genome_mixing <- function(emb_umap, genome, k = 30) {
  nn <- FNN::get.knn(emb_umap, k = k)$nn.index
  genome <- as.character(genome)
  mean(vapply(seq_len(nrow(nn)), function(i) {
    mean(genome[nn[i, ]] != genome[i])
  }, numeric(1)))
}

# centred genome mixing: observed minus the fraction expected from the global genome composition
genome_mixing_centered <- function(emb_umap, genome, k = 30) {
  nn <- FNN::get.knn(emb_umap, k = k)$nn.index
  genome <- as.character(genome)
  p_global <- table(genome) / length(genome)

  mix <- vapply(seq_len(nrow(nn)), function(i) {
    g <- genome[i]
    expected_other <- 1 - as.numeric(p_global[g])
    observed_other <- mean(genome[nn[i, ]] != g)
    observed_other - expected_other
  }, numeric(1))

  mean(mix)  # 0 means as mixed as the global expectation
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

plot_umap_genome <- function(umap, meta, title, subtitle,
                             point_size = 0.12, alpha = 0.6) {
  df <- data.frame(
    umap1  = umap[, 1],
    umap2  = umap[, 2],
    Genome = meta$Genome
  )
  ggplot(df, aes(umap1, umap2, color = Genome)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = c("At" = "#AB82FF", "B73" = "#FFB90F"), na.value = "grey70") +
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

if (!"Genome" %in% colnames(obj$meta)) {
  if ("species" %in% colnames(obj$meta)) obj$meta$Genome <- obj$meta$species
}
if (!"Genome" %in% colnames(obj$meta)) stop("Metadata must contain 'species' or 'Genome' column.")

if (!"log10nSites" %in% colnames(obj$meta) && "nSites" %in% colnames(obj$meta)) {
  obj$meta$log10nSites <- log10(obj$meta$nSites + 1)
}

# -------------------------
# Fix TF-IDF + SVD once
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
max_pcs  <- max(pcs_grid) +1

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
# Grid to scan (small)
# -------------------------
# reduceDims(cor.max) can return fewer PCs than requested (on low-coverage sets more depth-correlated
# components are dropped), so the pcs grid is capped to what is actually available.
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

# Diagnostics k must be <= smallest k_near (otherwise UMAP kNN not meaningful)
diag_k <- min(k_grid)

# -------------------------
# Run scan (panels B, C: one genome-coloured UMAP page per combination)
# -------------------------
metrics_out <- vector("list", nrow(grid))

plots_pdf <- file.path(outdir, paste0("plots/UMAP_grid_scan", mc_tag, ".Genome.pdf"))
pdf(plots_pdf, width = 6, height = 5)

for (i in seq_len(nrow(grid))) {
  pcs_i      <- grid$pcs[i]
  k_near_i   <- grid$k_near[i]
  min_dist_i <- grid$min_dist[i]

  tag <- paste0("pcs=", pcs_i, " | k_near=", k_near_i, " | min_dist=", min_dist_i)

  # projectUMAP reads from soc_i$PCA and writes soc_i$UMAP + meta umap1/umap2
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

  pres <- knn_preservation(emb_pca = soc_i$PCA, emb_umap = umap, k = diag_k)

  mix_raw <- genome_mixing(emb_umap = umap, genome = soc_i$meta$Genome, k = diag_k)
  mix_ctr <- genome_mixing_centered(emb_umap = umap, genome = soc_i$meta$Genome, k = diag_k)

  qc_cor <- qc_cor_summary(umap, soc_i$meta, qc_cols = c("log10nSites", "pOrg"))
  cor_abs_max <- max(abs(unlist(qc_cor)), na.rm = TRUE)

  metrics_out[[i]] <- cbind(
    data.frame(
      pcs = pcs_i,
      k_near = k_near_i,
      min_dist = min_dist_i,
      diag_k = diag_k,
      knn_preservation = pres,
      genome_mixing = mix_raw,
      genome_mixing_centered = mix_ctr,
      qc_abs_cor_max = cor_abs_max,
      n_cells = nrow(umap),
      stringsAsFactors = FALSE
    ),
    qc_cor
  )

  subtitle <- paste0(
    "kNN_pres=", sprintf("%.3f", pres),
    " | mix=", sprintf("%.3f", mix_raw),
    " | mix_ctr=", sprintf("%.3f", mix_ctr),
    " | max|rho|=", sprintf("%.3f", cor_abs_max),
    " | diag_k=", diag_k
  )

  print(plot_umap_genome(umap, soc_i$meta, title = tag, subtitle = subtitle))

  message("Done [", i, "/", nrow(grid), "]: ", tag,
          " | pres=", sprintf("%.3f", pres),
          " mix=", sprintf("%.3f", mix_raw),
          " mix_ctr=", sprintf("%.3f", mix_ctr),
          " max|rho|=", sprintf("%.3f", cor_abs_max))

  rm(soc_i)
}

dev.off()

metrics_df <- rbindlist(metrics_out, use.names = TRUE, fill = TRUE)
metrics_tsv <- file.path(outdir, paste0("UMAP_grid_scan", mc_tag, ".metrics.tsv"))
fwrite(metrics_df, metrics_tsv, sep = "\t", quote = FALSE)

metrics_df <- read.table(metrics_tsv,  sep = "\t", header = TRUE)

# -------------------------
# Rank / choose the best-scoring combination
# -------------------------
# Use centered mixing as the objective (0 is ideal); raw mixing depends on class proportions.
# score = knn_preservation - 2 |centred mixing| - 0.5 max |QC corr|
metrics_df <- as.data.frame(metrics_df) %>%
  mutate(
    mix_penalty = abs(genome_mixing_centered),   # closer to 0 is better
    score = knn_preservation - 2 * mix_penalty - 0.5 * qc_abs_cor_max
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

# -------------------------
# Panel A draft: kNN preservation vs genome mixing across the grid.
# Panel A as shipped is redrawn from the metrics TSV by figS3A_replot.R.
# -------------------------
best <- metrics_df[1, ]
bad <- metrics_df %>%
  arrange(knn_preservation) %>%
  slice(1)

plots_pdf <- file.path(outdir, paste0("plots/knn_preservation_grid_scan", mc_tag, ".GenomeMixing.pdf"))
pdf(plots_pdf, width = 9, height = 3)

ggplot(metrics_df, aes(x = knn_preservation, y = genome_mixing)) +
  geom_point(aes(color = qc_abs_cor_max), size = 3) +
  scale_color_viridis_c(name = "max |QC corr|") +

  scale_y_continuous() +

  # Highlight the best-scoring combination
  geom_point(
    data = best,
    aes(x = knn_preservation, y = genome_mixing),
    shape = 21, size = 5, stroke = 1.2, fill = "white", color = "black"
  ) +

  geom_point(
    data = bad,
    aes(x = knn_preservation, y = genome_mixing),
    shape = 21, size = 5, stroke = 1.2, fill = "white", color = "red"
  ) +

  geom_text(
    data = best,
    aes(label = paste0("Selected\npcs=", pcs,
                       ", k=", k_near,
                       ", min_dist=", min_dist)),
    hjust = -0.1, vjust = -0.5, size = 3
  ) +

  geom_text(
    data = bad,
    aes(label = paste0("Selected\npcs=", pcs,
                       ", k=", k_near,
                       ", min_dist=", min_dist)),
    hjust = -0.1, vjust = -0.5, size = 3
  ) +

  theme_bw() +
  labs(
    x = "kNN preservation (PCA to UMAP)",
    y = "Genome mixing",
    title = "UMAP robustness across parameter grid",
    subtitle = "Cross-species co-projection is stable and not driven by QC covariates"
  )

dev.off()
