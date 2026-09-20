#!/usr/bin/env Rscript
###############################################################################
## 3_2_coembed_gridscan.R -- UMAP parameter scan + quantitative diagnostics for the At/B73 co-projection
##
## Fix TF-IDF + SVD once; then scan a small grid over:
##   pcs, k_near (UMAP n_neighbors), min_dist
##
## For each combo compute:
##   1) kNN preservation (PCA vs UMAP neighborhoods; same k = diag_k)
##   2) genome mixing (optionally baseline-centered)
##   3) QC correlation leakage (max |Spearman rho| of UMAP1/2 vs log10nSites, pOrg)
##
## Writes:
##   - metrics TSV (one row per combo)
##   - PDF of Genome-colored UMAPs labeled with diagnostics
##   - best-row TSV (by a composite score)
##
## Usage:
##   Rscript 3_2_coembed_gridscan.R <soc_rds> <meta_tsv> <outdir> [seed] [min_c]
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
  stop("Usage: Rscript 3_2_coembed_gridscan.R <soc_rds> <meta_tsv> <outdir> [seed] [min_c]")
}

soc_rds  <- args[1]
meta_tsv <- args[2]
outdir   <- args[3]
seed      <- if (length(args) >= 4) as.integer(args[4]) else 1L
min_c_arg <- if (length(args) >= 5) as.numeric(args[5]) else NA_real_  # NA -> data-driven 250 floor (mirror 3_0_0)
mc_tag    <- if (is.na(min_c_arg)) "" else paste0(".minc_", as.integer(min_c_arg))


dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# TFIDF (identical to the combined-arm pipeline)
# -------------------------
# TFIDF normalization
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
knn_preservation <- function(emb_pca, emb_umap, k = 30) {
  nn_pca  <- FNN::get.knn(emb_pca,  k = k)$nn.index
  nn_umap <- FNN::get.knn(emb_umap, k = k)$nn.index
  mean(vapply(seq_len(nrow(nn_pca)), function(i) {
    length(intersect(nn_pca[i, ], nn_umap[i, ])) / k
  }, numeric(1)))
}

genome_mixing <- function(emb_umap, genome, k = 30) {
  nn <- FNN::get.knn(emb_umap, k = k)$nn.index
  genome <- as.character(genome)
  mean(vapply(seq_len(nrow(nn)), function(i) {
    mean(genome[nn[i, ]] != genome[i])
  }, numeric(1)))
}

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
  
  mean(mix)  # 0 means “as mixed as global expectation”
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
pca_full[1:2,]

# -------------------------
# Grid to scan (small)
# -------------------------
# reduceDims(cor.max) can return fewer PCs than requested — on the low-coverage min.c=50
# set more depth-correlated components are dropped — so cap the pcs grid to what is actually
# available, else pca_full[, seq_len(pcs_i)] goes subscript-out-of-bounds (crash at pcs=50).
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
diag_k <- min(k_grid)  # safest default; you can set diag_k=15 explicitly

# -------------------------
# Run scan
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

 metrics_df <- read.table(metrics_tsv,  sep = "\t", header = T)
 
# -------------------------
# Rank / choose “best”
# -------------------------
# Use centered mixing as the objective (0 is ideal); raw mixing depends on class proportions.
metrics_df <- as.data.frame(metrics_df) %>%
  mutate(
    mix_penalty = abs(genome_mixing_centered),   # closer to 0 is better
    score = knn_preservation - 2 * mix_penalty - 0.5 * qc_abs_cor_max
  ) %>%
  arrange(desc(score))

best_row <- metrics_df[1, , drop = FALSE]
best_tsv <- file.path(outdir, paste0("UMAP_grid_scan", mc_tag, ".best.tsv"))
fwrite(best_row, best_tsv, sep = "\t", quote = FALSE)
# (removed stray read.ftable() here — it clobbered the scored metrics_df and aborted the
#  final knn_preservation/GenomeMixing scatter plot below.)

message("Wrote metrics:   ", metrics_tsv)
message("Wrote UMAP PDF:  ", plots_pdf)
message("Best combo TSV:  ", best_tsv)
message("Best parameters:")
print(best_row)

# or manually select the one you used
library(ggplot2)
best <- metrics_df[1, ]
bad <- metrics_df %>%
  arrange(knn_preservation) %>%
  slice(1)

plots_pdf <- file.path(outdir, paste0("plots/knn_preservation_grid_scan", mc_tag, ".GenomeMixing.pdf"))
pdf(plots_pdf, width = 9, height = 3)

ggplot(metrics_df, aes(x = knn_preservation, y = genome_mixing)) +
  geom_point(aes(color = qc_abs_cor_max), size = 3) +
  scale_color_viridis_c(name = "max |QC corr|") +
  
  # Expected mixing
  #geom_hline(yintercept = 0.4, linetype = "dashed", color = "grey40") +
  scale_y_continuous() +   # was limits=c(0.2,0.5); clipped post mixing (~0.015). Auto-scale to the data.
  
  # Highlight chosen UMAP
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
  
  # Interpretive annotation
  annotate(
    "text",
    x = min(metrics_df$knn_preservation),
    y = 0.4,
    hjust = 0,
    label = "Higher kNN preservation\n→ stable genome mixing\n→ low QC dependence",
    size = 3.5
  ) +
  
  theme_bw() +
  labs(
    x = "kNN preservation (PCA → UMAP)",
    y = "Genome mixing",
    title = "UMAP robustness across parameter grid",
    subtitle = "Cross-species co-projection is stable and not driven by QC covariates"
  )

dev.off()
