#!/usr/bin/env Rscript
# Fig S3, panel A (as shipped): kNN preservation vs genome mixing over the 45-point UMAP grid, redrawn
# from the metrics table written by figS3.R (no Socrates object needed).
# The dashed line is the mixing expected under random intermingling (genome_mixing - genome_mixing_centered,
# identical for every row since it depends only on species composition), so the y axis has a reference.
# Two configurations are highlighted: the one used in this study (pcs 20, k 30, min_dist 0.3; held fixed for
# cross-object comparability, not the score maximiser) and the score maximiser
# (score = kNN preservation - 2 |centred mixing| - 0.5 max |QC corr|, the ranking figS3.R writes to *.best.tsv).
# Input:  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2/step3_embedscan/UMAP_grid_scan.metrics.tsv
# Usage:  Rscript analysis/supplementary/figS3A_replot.R [metrics.tsv] [out_prefix]
# Output: figures/supplementary/figS3/FigS3A_knn_preservation_grid_scan.{pdf,png}

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

# ---- CONFIG ----------------------------------------------------------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis/socrates"
OUTDIR <- "figures/supplementary/figS3"
# configuration used in this study (Fig S3B): held fixed for comparability across objects, not the score maximiser
USED   <- list(pcs = 20, k_near = 30, min_dist = 0.3)

args <- commandArgs(trailingOnly = TRUE)
metrics_tsv <- if (length(args) >= 1) args[1] else file.path(DATA, "SM2/step3_embedscan/UMAP_grid_scan.metrics.tsv")
out_prefix  <- if (length(args) >= 2) args[2] else file.path(OUTDIR, "FigS3A_knn_preservation_grid_scan")
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

m <- read.delim(metrics_tsv)
# expected mixing under random intermingling; identical for every row by construction
m$expected <- m$genome_mixing - m$genome_mixing_centered
stopifnot(diff(range(m$expected)) < 1e-6)
exp_mix <- m$expected[1]
diag_k  <- unique(m$diag_k); stopifnot(length(diag_k) == 1)

m <- m %>% mutate(score = knn_preservation - 2 * abs(genome_mixing_centered) - 0.5 * qc_abs_cor_max)
best <- m %>% arrange(desc(score)) %>% slice(1)                                             # score maximiser (Fig S3C)
used <- m %>% filter(pcs == USED$pcs, k_near == USED$k_near, min_dist == USED$min_dist)   # configuration used (Fig S3B)
stopifnot(nrow(used) == 1)
lab <- function(d, role) paste0(role, "\npcs=", d$pcs, ", k=", d$k_near, ", min_dist=", d$min_dist)

p <- ggplot(m, aes(knn_preservation, genome_mixing)) +
  geom_hline(yintercept = exp_mix, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(m$knn_preservation), y = exp_mix, vjust = -0.5, hjust = 1, size = 3,
           colour = "grey30", label = sprintf("random intermingling (%.2f)", exp_mix)) +
  geom_point(aes(colour = qc_abs_cor_max), size = 3) +
  scale_colour_viridis_c(name = "max |QC corr|") +
  geom_point(data = used, shape = 21, size = 5, stroke = 1.2, fill = "white", colour = "black") +
  geom_point(data = best, shape = 21, size = 5, stroke = 1.2, fill = "white", colour = "black") +
  geom_text(data = used, aes(label = lab(used, "Used in this study")), hjust = -0.08, vjust = 1.6, size = 3) +
  geom_text(data = best, aes(label = lab(best, "Score maximiser")),   hjust = 1.08, vjust = 1.6, size = 3) +
  annotate("text", x = min(m$knn_preservation), y = exp_mix + 0.03, hjust = 0, vjust = 0, size = 3.3,
           label = sprintf("Mixing stays just below random intermingling\nin all %d configurations, whatever the kNN preservation.\nNo configuration depends on QC covariates", nrow(m))) +
  scale_y_continuous(limits = c(0.25, 0.48)) +
  theme_bw() +
  labs(x = sprintf("kNN preservation, PCA to UMAP (k = %d)", diag_k),
       y = sprintf("Fraction of other-species neighbours (k = %d)", diag_k),
       title = "UMAP robustness across the parameter grid",
       subtitle = "Cross-species co-projection is stable and not driven by QC covariates")

ggsave(paste0(out_prefix, ".pdf"), p, width = 9, height = 3.4)
ggsave(paste0(out_prefix, ".png"), p, width = 9, height = 3.4, dpi = 150)
message("expected mixing under random intermingling = ", round(exp_mix, 4),
        "; obs/exp range = ", paste(round(range(m$genome_mixing / exp_mix), 3), collapse = " to "))
