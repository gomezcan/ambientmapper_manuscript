#!/usr/bin/env Rscript
# Fig 4, panels L to N (effect of AmbientMapper cleaning on allele purity, WASP-corrected alignments).
#   L: allele purity (WASP-corrected p_top1) before vs after cleaning, per barcode (hexbin),
#      B73/Mo17 rep1 + rep2 pooled and multi-genotype, faceted by AmbientMapper call class
#   M: contaminant reads removed per cell at informative sites, by call class x depth bin
#   N: selectivity = contaminant-removal rate / on-target-removal rate
#   (M and N are built as one patchwork object, pM, and lettered separately in the assembled figure)
# Inputs (DATA = data/processed/zhang2024): <s>/diagnostics/06_48_barcode_purity/<s>_barcode_purity.tsv.gz
#   (barcode x 1-Mb-block purity on the WASP-corrected raw vs cleaned BAMs, workflows/03b_variant_based_comparison/)
# Sources analysis/_helpers/fig4_helpers.R. Run from the repo root: Rscript analysis/fig4_robustness/fig4_part4.R

# Purity is ABSOLUTE (WASP removes the single-reference mapping bias).
# weak_doublet -> singlet everywhere: AmbientMapper's weak_doublet calls are
# genetically clean singlets mislabelled by the k>=2 doublet model (pairfrac ~0.93,
# top1or2 ~1.0, p_top1 ~ singlet; the evidence is Table S4). Presentational only.
# Read and contaminant accounting is at INFORMATIVE marker sites: T1 = ref + alt
# reads at {top1-hom AND informative} panel sites, M1 = reads matching the assigned
# genome, off = T1 - M1. Informative reads are ~42 to 87% of the total.

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
DATA    <- "data/processed/zhang2024"
OUTDIR  <- "figures/main/fig4"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SAMPLES <- c("B73Mo17_rep1", "B73Mo17_rep2", "multiGenotypes_rep1")
# B73/Mo17 rep1+rep2 pooled into one group; multi kept separate.
GROUP_OF <- c(B73Mo17_rep1 = "B73Mo17", B73Mo17_rep2 = "B73Mo17",
              multiGenotypes_rep1 = "multiGenotypes_rep1")
GROUP_LABELS <- c(B73Mo17 = "B73/Mo17 (rep1+rep2)",
                  multiGenotypes_rep1 = "Multi (7 genotypes)")
GROUP_ORDER  <- c("B73Mo17", "multiGenotypes_rep1")
GROUP_LABELS_SHORT <- c(B73Mo17 = "B73/Mo17", multiGenotypes_rep1 = "Multi (7 geno.)")

DEPTH_LEVELS <- c("200-500", "500-1000", "1000-2000",
                  "2000-5000", "5000-10000", ">10000")

CLASS_COLORS <- c(singlet = "#264653", doublet = "#E76F51")  # same singlet/doublet palette as part 3
GROUP_LTY    <- c("B73/Mo17" = "solid", "Multi (7 geno.)" = "22")

# Facets with fewer barcodes than this are dropped: removes the multi-genotype
# singlet facet (n = 960), too few barcodes to interpret.
MIN_BC_FACET <- 1000

genodiag <- function(s, sub, file)
  file.path(DATA, s, "diagnostics", sub, file)

hex_layer <- function(bins = 45) {
  if (requireNamespace("hexbin", quietly = TRUE)) geom_hex(bins = bins)
  else geom_bin2d(bins = bins)
}

# =============================================================================
# 1) LOAD 06_48 per-barcode purity (once; feeds L, M and N)
# =============================================================================
load_0648 <- function(s) {
  f <- genodiag(s, "06_48_barcode_purity", paste0(s, "_barcode_purity.tsv.gz"))
  if (!file.exists(f)) stop("Missing 06_48: ", f)
  dt <- fread(f, sep = "\t",
              select = c("barcode", "mode", "T1", "M1", "p_top1",
                         "class", "top1", "depth_bin"))
  dt[, `:=`(group = GROUP_OF[s], sample = s)]
  dt
}
cat("Loading 06_48 barcode purity for", length(SAMPLES), "samples...\n")
bp <- rbindlist(lapply(SAMPLES, load_0648))
bp[class == "weak_doublet", class := "singlet"]        # weak_doublet folded into singlet
bp <- bp[class %in% c("singlet", "doublet")]
bp[, class := factor(class, levels = names(CLASS_COLORS))]

# =============================================================================
# 2) PANEL L: absolute allele purity, raw vs clean (hexbin, by class)
# =============================================================================
# `sample` in the key: identical barcode sequences recur across rep1/rep2.
wide <- dcast(bp, group + sample + barcode + class + top1 + depth_bin ~ mode,
              value.var = "p_top1")
setnames(wide, c("raw", "clean"), c("p_raw", "p_clean"))
wide <- wide[is.finite(p_raw) & is.finite(p_clean)]

facet_n <- wide[, .N, by = .(group, class)]
wide <- wide[facet_n[N >= MIN_BC_FACET], on = .(group, class)]
wide[, facet := paste0(GROUP_LABELS[group], "\n", as.character(class))]
facet_order <- unlist(lapply(GROUP_ORDER, function(g)
  paste0(GROUP_LABELS[g], "\n",
         intersect(c("singlet", "doublet"),
                   as.character(unique(wide[group == g]$class))))))
wide[, facet := factor(facet, levels = facet_order)]

med <- wide[, .(p_raw_med = median(p_raw), p_clean_med = median(p_clean),
                n = .N), by = .(facet)]
med[, lab := sprintf("median %.3f -> %.3f\n(n=%s)",
                     p_raw_med, p_clean_med, format(n, big.mark = ","))]

pL <- ggplot(wide, aes(p_raw, p_clean)) +
  hex_layer(bins = 45) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  geom_text(data = med, aes(x = 0.02, y = 0.98, label = lab),
            hjust = 0, vjust = 1, size = 2.4, color = "grey15", inherit.aes = FALSE) +
  scale_fill_viridis_c(trans = "log10", name = "barcodes") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  facet_wrap(~ facet, nrow = 1) +
  labs(x = "Allele purity, raw  (WASP-corrected p_top1)",
       y = "Allele purity, cleaned",
       title = "L. Allele purity increases after AM cleaning (mass above diagonal)",
       subtitle = "WASP-corrected absolute purity; per-barcode; AM weak_doublet folded into singlet (06_48 merge)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5),
        legend.position = "right",
        strip.text = element_text(size = 7.5))

# =============================================================================
# 3) PANELS M and N: AM selectively removes contaminant reads
# =============================================================================
w <- dcast(bp, group + sample + barcode + class + depth_bin ~ mode,
           value.var = c("T1", "M1"))
w <- w[!is.na(T1_raw) & !is.na(T1_clean)]
w[, off_raw        := T1_raw   - M1_raw]        # contaminant marker reads, raw
w[, contam_dropped := off_raw - (T1_clean - M1_clean)]
w[, match_dropped  := M1_raw  - M1_clean]

magg <- w[, .(
  n_bc            = .N,
  contam_per_cell = mean(contam_dropped),                 # panel M
  contam_rm_rate  = sum(contam_dropped) / sum(off_raw),   # panel N numerator
  match_rm_rate   = sum(match_dropped)  / sum(M1_raw)     # panel N denominator
), by = .(group, class, depth_bin)]
magg[, selectivity := contam_rm_rate / match_rm_rate]
magg[, depth_bin := factor(depth_bin, levels = DEPTH_LEVELS)]
magg[, gl := factor(GROUP_LABELS_SHORT[group], levels = GROUP_LABELS_SHORT[GROUP_ORDER])]

m_theme <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        legend.position = "right", legend.key.height = unit(0.8, "lines"))

pM1 <- ggplot(magg, aes(depth_bin, contam_per_cell, color = class,
                        group = interaction(gl, class), linetype = gl)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.9) +
  scale_color_manual(values = CLASS_COLORS, name = "AM call class") +
  scale_linetype_manual(values = GROUP_LTY, name = "Dataset") +
  scale_y_log10(labels = label_number(accuracy = 1)) +
  labs(x = "Per-cell depth (total reads)",
       y = "Contaminant reads removed / cell",
       title = "M. AM selectively removes contaminant reads",
       subtitle = "Reads removed per cell (informative sites)") +
  m_theme

pM2 <- ggplot(magg, aes(depth_bin, selectivity, color = class,
                        group = interaction(gl, class), linetype = gl)) +
  geom_hline(yintercept = 1, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.9) +
  scale_color_manual(values = CLASS_COLORS, name = "AM call class") +
  scale_linetype_manual(values = GROUP_LTY, name = "Dataset") +
  scale_y_continuous(labels = label_number(accuracy = 1, suffix = "×"),
                     limits = c(0, NA)) +
  labs(x = "Per-cell depth (total reads)",
       y = "Selectivity (contaminant / on-target removal rate)",
       title = " ",
       subtitle = "AM removes contaminants at 6-42x higher rate than genuine reads") +
  m_theme

pM <- pM1 + pM2 + plot_layout(guides = "collect") & theme(legend.position = "right")

# =============================================================================
# 4) ASSEMBLY & EXPORT
# =============================================================================
cap <- paste(
  "Purity and read counts at informative marker sites (WASP-corrected).",
  "weak_doublet folded into singlet (06_48 merge, Table S4)."
)

fig <- (pL / pM) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(caption = cap) &
  theme(plot.caption = element_text(size = 7.5, hjust = 0))

ggsave(file.path(OUTDIR, "Fig4_part4_LtoN.pdf"), fig,
       width = 12, height = 7.5, limitsize = FALSE)
ggsave(file.path(OUTDIR, "Fig4_part4_LtoN.png"), fig,
       width = 12, height = 7.5, dpi = 200, limitsize = FALSE)

# Individual panels (M and N are the two halves of pM)
ggsave(file.path(OUTDIR, "Fig4_L_purity_abs.pdf"), pL, width = 9.5, height = 3.4)
ggsave(file.path(OUTDIR, "Fig4_MN_reads_removed.pdf"), pM, width = 10, height = 3.6)

cat("\nDone. Outputs in:", OUTDIR, "\n")
cat("Panel L medians:\n"); print(med[order(facet), .(facet, p_raw_med, p_clean_med, n)])
