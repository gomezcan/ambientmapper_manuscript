#!/usr/bin/env Rscript
# Fig 4, panels H to K (AmbientMapper C0 vs Souporcell v2.1, both on the raw, pre-cleaning BAMs).
#   H: B73/Mo17 (rep1 + rep2 pooled), AmbientMapper call x Souporcell status
#   I: multi-genotype (7 genotypes), AmbientMapper call x Souporcell status
#   J: B73/Mo17 (rep1 + rep2 pooled), AmbientMapper genome_1 x Souporcell genotype
#      (Souporcell singlets only; cluster IDs mapped to genotypes via Genotype_ID_key.v2.txt)
#   K: multi-genotype (7 genotypes), same as J
# Inputs (DATA = data/processed/zhang2024): <s>/genotyping_runs/<AM_RUN_TAG>/C0/<s>_cells_calls.tsv.gz,
#   <s>/souporcell/supervised/<s>.min500/{clusters.tsv, Genotype_ID_key.v2.txt}
# Sources analysis/_helpers/fig4_helpers.R. Run from the repo root: Rscript analysis/fig4_robustness/fig4_part3.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(data.table)
})

source("analysis/_helpers/fig4_helpers.R")

# -------------------------
# 0) CONFIG
# -------------------------
DATA   <- "data/processed/zhang2024"
OUTDIR <- "figures/main/fig4"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

DATASETS <- c("B73Mo17_rep1", "B73Mo17_rep2", "multiGenotypes_rep1")

SAMPLE_GROUP <- c(
  "B73Mo17_rep1"        = "B73/Mo17 (rep1+rep2)",
  "B73Mo17_rep2"        = "B73/Mo17 (rep1+rep2)",
  "multiGenotypes_rep1" = "Multi (7 genotypes)"
)
GROUP_LEVELS <- c("B73/Mo17 (rep1+rep2)", "Multi (7 genotypes)")

# `low_reads` excluded by design: Souporcell .min500 filter means joined
# barcodes always have >=500 reads, so low_reads is structurally 0.
AM_CALL_LEVELS <- c("single_clean", "dirty_singlet", "weak_doublet",
                    "doublet", "ambiguous", "empty")

SOUP_STATUS_LEVELS <- c("singlet", "doublet", "unassigned")

AM_RUN_TAG <- "4cfg_2026-05-01"   # AmbientMapper genotyping run tag (C0 configuration)

AM_FILE <- function(sample) {
  file.path(DATA, sample, "genotyping_runs", AM_RUN_TAG, "C0",
            paste0(sample, "_cells_calls.tsv.gz"))
}
SOUP_DIR <- function(sample) {
  file.path(DATA, sample, "souporcell", "supervised",
            paste0(sample, ".min500"))
}

# Extract 26-char barcode prefix (everything before the first dash)
strip_bc <- function(x) sub("-.*$", "", x)

# -------------------------
# 1) LOAD + JOIN per dataset
# -------------------------
load_pair <- function(sample) {
  cat("Loading", sample, "...\n")

  am <- fread(AM_FILE(sample),
              select = c("barcode", "call", "genome_1"))
  am[, bc := strip_bc(barcode)]
  cat("  AM:    ", format(nrow(am), big.mark = ","), "BCs\n")

  soup <- fread(file.path(SOUP_DIR(sample), "clusters.tsv"),
                select = c("barcode", "status", "assignment"))
  soup[, bc := strip_bc(barcode)]
  cat("  Soup:  ", format(nrow(soup), big.mark = ","), "BCs\n")

  joined <- merge(am[,   .(bc, am_call = call, am_genome_1 = genome_1)],
                  soup[, .(bc, soup_status = status,
                           soup_assignment = assignment)],
                  by = "bc")
  cat("  Joined:", format(nrow(joined), big.mark = ","), "BCs\n\n")

  joined[, sample := sample]
  joined
}

dat <- rbindlist(lapply(DATASETS, load_pair))
dat[, am_call     := factor(am_call,     levels = AM_CALL_LEVELS)]
dat[, soup_status := factor(soup_status, levels = SOUP_STATUS_LEVELS)]
dat[, group       := factor(SAMPLE_GROUP[sample], levels = GROUP_LEVELS)]

cat("Pooled BCs per group:\n")
print(dat[, .N, by = group])

# -------------------------
# 1b) LOAD cluster -> genotype mapping per dataset
# -------------------------
load_cluster_map <- function(sample) {
  m <- fread(file.path(SOUP_DIR(sample), "Genotype_ID_key.v2.txt"))
  setnames(m, c("Genotype_ID", "Cluster_ID"), c("genotype", "cluster"))
  m[, cluster := as.character(cluster)]
  m
}

# Both B73Mo17 reps share the same B73/Mo17 plate design and Souporcell prior,
# so cluster 0 = B73, cluster 1 = Mo17 in both reps. Use rep1's mapping for the
# pooled B73/Mo17 panel.
cmap_b73mo17 <- load_cluster_map("B73Mo17_rep1")
cmap_multi   <- load_cluster_map("multiGenotypes_rep1")

cat("\nCluster -> Genotype maps:\n")
cat("  B73Mo17:\n"); print(cmap_b73mo17)
cat("  Multi:  \n"); print(cmap_multi)

# Sanity check: B73Mo17_rep2 mapping should match rep1
cmap_b73mo17_rep2 <- load_cluster_map("B73Mo17_rep2")
if (!isTRUE(all.equal(
      cmap_b73mo17[order(genotype), .(genotype, cluster)],
      cmap_b73mo17_rep2[order(genotype), .(genotype, cluster)]))) {
  warning("B73Mo17 rep1 vs rep2 cluster->genotype mapping differs! ",
          "Pooling assumes they match.")
}

# =============================================================================
# 2) PANEL BUILDER H/I: AM call x Soup status confusion matrix
# =============================================================================
build_call_confmat <- function(df, title) {
  conf <- df[, .(N = .N), by = .(am_call, soup_status)]

  full_grid <- CJ(am_call     = AM_CALL_LEVELS,
                  soup_status = SOUP_STATUS_LEVELS)
  conf <- merge(full_grid, conf, by = c("am_call", "soup_status"), all.x = TRUE)
  conf[is.na(N), N := 0L]

  conf[, am_call     := factor(am_call,     levels = rev(AM_CALL_LEVELS))]
  conf[, soup_status := factor(soup_status, levels = SOUP_STATUS_LEVELS)]

  conf[, row_total := sum(N), by = am_call]
  conf[, pct       := ifelse(row_total > 0, 100 * N / row_total, NA_real_)]

  conf[, lbl := ifelse(is.na(pct),
                       scales::comma(N),
                       sprintf("%s\n(%.0f%%)", scales::comma(N), pct))]

  ggplot(conf, aes(x = soup_status, y = am_call, fill = N)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = lbl,
                  color = ifelse(N > max(N, na.rm = TRUE) * 0.45, "w", "b")),
              size = 2.7, lineheight = 0.9) +
    scale_fill_gradient(low = "#f0f0f0", high = "#2c3e50",
                        trans = "log1p", guide = "none") +
    scale_color_manual(values = c("b" = "grey15", "w" = "white"),
                       guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = "Souporcell status", y = "AmbientMapper call",
         title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title   = element_text(face = "bold"),
          panel.grid   = element_blank(),
          axis.title.x = element_text(margin = margin(b = 4)),
          axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 10))
}

# =============================================================================
# 3) PANEL BUILDER J/K: AM genome_1 x Soup genotype confusion matrix
#     Filtered to Soup status == "singlet" so the cluster ID maps to a single
#     genotype. AM call type is NOT restricted (AM may call any type;
#     genome_1 is the top genome regardless).
# =============================================================================
build_genotype_confmat <- function(df, cmap, title) {
  s <- df[soup_status == "singlet"]
  cmap_vec <- setNames(cmap$genotype, cmap$cluster)
  s[, soup_genotype := cmap_vec[as.character(soup_assignment)]]

  # Rows: AM genome_1 (each named genome + None, if any)
  geno_order <- c(cmap$genotype, "None")
  row_levels <- intersect(geno_order, unique(s$am_genome_1))
  col_levels <- cmap$genotype

  conf <- s[, .(N = .N), by = .(am_genome_1, soup_genotype)]
  full <- CJ(am_genome_1   = row_levels,
             soup_genotype = col_levels)
  conf <- merge(full, conf, by = c("am_genome_1", "soup_genotype"), all.x = TRUE)
  conf[is.na(N), N := 0L]

  conf[, am_genome_1   := factor(am_genome_1,   levels = rev(row_levels))]
  conf[, soup_genotype := factor(soup_genotype, levels = col_levels)]

  # Column-normalized percentage (within each Soup genotype, what fraction
  # does AM call as each genome). Diagonal = correct calls.
  conf[, col_total := sum(N), by = soup_genotype]
  conf[, pct       := ifelse(col_total > 0, 100 * N / col_total, NA_real_)]

  conf[, lbl := ifelse(is.na(pct),
                       scales::comma(N),
                       sprintf("%s\n(%.0f%%)", scales::comma(N), pct))]

  ggplot(conf, aes(x = soup_genotype, y = am_genome_1, fill = N)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = lbl,
                  color = ifelse(N > max(N, na.rm = TRUE) * 0.45, "w", "b")),
              size = 2.4, lineheight = 0.9) +
    scale_fill_gradient(low = "#f0f0f0", high = "#1f4e79",
                        trans = "log1p", guide = "none") +
    scale_color_manual(values = c("b" = "grey15", "w" = "white"),
                       guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = "Souporcell genotype (singlets)",
         y = "AmbientMapper genome_1",
         title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title   = element_text(face = "bold"),
          panel.grid   = element_blank(),
          axis.title.x = element_text(margin = margin(b = 4)),
          axis.text.x  = element_text(size = 9, angle = 0),
          axis.text.y  = element_text(size = 9))
}

# =============================================================================
# 4) BUILD PANELS
# =============================================================================
pH <- build_call_confmat(dat[group == "B73/Mo17 (rep1+rep2)"],
                         "B73/Mo17 --methods")
pI <- build_call_confmat(dat[group == "Multi (7 genotypes)"],
                         "Multi --methods")
pJ <- build_genotype_confmat(dat[group == "B73/Mo17 (rep1+rep2)"],
                             cmap_b73mo17, "B73/Mo17 --genotype")
pK <- build_genotype_confmat(dat[group == "Multi (7 genotypes)"],
                             cmap_multi,  "Multi --genotype")

# =============================================================================
# 5) ASSEMBLE & EXPORT
# =============================================================================
# Layout (2 rows x 2 cols):
#   H | I
#   J | K
# Multi panels (I, K) get slightly more width because K is a 7-genotype matrix.
fig <- (pH | pI) / (pJ | pK) +
  plot_layout(widths = c(1, 1.3), heights = c(1, 0.95)) +
  plot_annotation(tag_levels = list(c("H", "I", "J", "K"))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(OUTDIR, "Fig4_part3_HtoK.pdf"), fig,
       width = 13, height = 11)
ggsave(file.path(OUTDIR, "Fig4_part3_HtoK.png"), fig,
       width = 13, height = 11, dpi = 300)

# Individual panels
ggsave(file.path(OUTDIR, "Fig4_H_B73Mo17_method.pdf"),   pH, width = 5.5, height = 5.5)
ggsave(file.path(OUTDIR, "Fig4_I_multi_method.pdf"),     pI, width = 6.5, height = 5.5)
ggsave(file.path(OUTDIR, "Fig4_J_B73Mo17_genotype.pdf"), pJ, width = 4.5, height = 4.5)
ggsave(file.path(OUTDIR, "Fig4_K_multi_genotype.pdf"),   pK, width = 7,   height = 6.5)

cat("\nDone. Outputs in:", OUTDIR, "\n")
