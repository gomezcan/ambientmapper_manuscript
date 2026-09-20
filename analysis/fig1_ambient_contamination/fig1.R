#!/usr/bin/env Rscript
# fig1.R -- Figure 1, panels B to H (ambient contamination in the interspecies barnyard).
#   B barnyard cloud (pre-clean winner reads), C off-target asymmetry by expected genome, D binned
#   off-target profile with cumulative capture, E the ratio trap (pseudocounted dominance ratio vs
#   depth), F to H the co-embedding UMAP by plate of origin, by Leiden cluster, and per-cluster composition.
# Inputs : data/processed/scifiATAC_B73_Arabidopsis/SM2/decontam_with_design_alpha05_v2/  (B to E)
#          data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_indep/coembed/Pre_cluster/  (F to H)
# Sources: analysis/_helpers/plotting.R
# Run    : Rscript analysis/fig1_ambient_contamination/fig1.R   (from the repository root)
#
# Inputs for B to E are AmbientMapper decontam outputs:
#   <SAMPLE>_barcode_policy.tsv.gz            per-barcode policy: expected_genome (design), allowed_set
#   <SAMPLE>_pre_barcode_genome_counts.tsv.gz  winner-evidence read counts per barcode x genome (PRE-clean)
#   <SAMPLE>_pre_barcode_composition.tsv.gz    per-barcode composition from winner counts (PRE-clean);
#                                              with a design: expected_genome, expected_frac,
#                                              contamination_rate = 1 - expected_frac
#   <SAMPLE>_cells_calls.decontam.tsv.gz       genotyping calls + metrics (ratio_top1_top2, purity, n_reads)
#   <SAMPLE>_barcode_postclean.tsv.gz          post-clean gate (keep_postclean) and post-clean read metrics
#
# Nomenclature:
#   - "Reads" are winner-evidence reads from the assignment files unless stated as n_reads from
#     cells_calls (barcode total).
#   - "Off-target fraction" = 1 - expected_frac (design-aware).
#   - "Dominance ratio" in panel E is the count-based pseudocounted (k1 + 1) / (k2 + 1), not
#     AmbientMapper's ratio_top1_top2 (see the panel E notes).

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(data.table)
  library(RColorBrewer)
})

source("analysis/_helpers/plotting.R")   # cols_expected_fill, cols_expected_line, cols_species

# -------------------------
# 0) CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SAMPLE <- "SM2"
INDIR  <- file.path(DATA, "SM2", "decontam_with_design_alpha05_v2")   # AmbientMapper run read by B to E

# Panels F to H: the independent (multi-reference) co-embedding of the pre-clean SM2v2 objects.
SOC_DIR   <- file.path(DATA, "socrates")
COEMB_CFG <- "pcs_20.k_near_30.min_dis_0.3.minc_50.res_0.5"

OUTDIR <- "figures/main/fig1"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Threshold for the binned quantification panel (D)
MIN_TOTAL_READS_ZOOM <- 10

# Bin edges for off-target fraction (fraction of reads NOT from expected genome)
breaks_custom <- c(-Inf, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, Inf)
labels_custom <- c("0–1%", "1–2%", "2–4%", "4–8%", "8–16%", "16–32%", "32–64%", ">64%")

# -------------------------
# 1) LOAD TABLES
# -------------------------
policy <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_barcode_policy.tsv.gz")),
                   show_col_types = FALSE) %>%
  mutate(
    barcode = as.character(barcode),
    expected_genome = na_if(as.character(expected_genome), "")
  )

pre_counts <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_pre_barcode_genome_counts.tsv.gz")),
                       show_col_types = FALSE) %>%
  mutate(
    barcode = as.character(barcode),
    genome  = as.character(genome),
    n_winner_reads = as.numeric(n_winner_reads)
  )

pre_comp <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_pre_barcode_composition.tsv.gz")),
                     show_col_types = FALSE) %>%
  mutate(barcode = as.character(barcode))

calls <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_cells_calls.decontam.tsv.gz")),
                  show_col_types = FALSE) %>%
  mutate(barcode = as.character(barcode))

post_gate <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_barcode_postclean.tsv.gz")),
                      show_col_types = FALSE) %>%
  mutate(
    barcode = as.character(barcode),
    keep_postclean = as.logical(keep_postclean)
  )

# =============================================================================
# DATA PREP: Socrates co-embedding metadata (panels F to H)
# =============================================================================
# Panels F to H use the INDEPENDENT (multi-reference) co-embedding, so they share the
# multi-reference design of the barnyard panels B to E and of Fig 3. The concatenated-reference
# (ZmATcombined) object is the Fig S1 robustness check only.
#
# The configuration (COEMB_CFG) is the same pcs / k / min_dist the concatenated Fig 1 was frozen
# at, held fixed for cross-object comparability (it is not the score maximiser), so this panel and
# Fig S1 differ only in mapping strategy. `Genome` is the plate-of-origin ground truth (At / B73)
# and matches cols_species. The `[,-c(1)]` drops the R-exported row-name column.
metaSoc <- fread(file.path(SOC_DIR, "SM2v2_indep/coembed/Pre_cluster",
                           paste0("SM2v2_coembed_Pre.updated_metadata_v7.", COEMB_CFG, ".txt")))[,-c(1)]

# `LouvainClusters` is the column name Socrates writes; the clustering is Leiden
# (Socrates::callClusters(cl.method = 4) forwards to Seurat::FindClusters(algorithm = 4)).

# Filter low-count clusters: clusters with < 100 cells are too noisy to plot
cluster_counts <- table(metaSoc$LouvainClusters)
valid_clusters <- names(cluster_counts[cluster_counts >= 100])

metaSocReduce <- metaSoc %>%
  filter(LouvainClusters %in% valid_clusters) %>%
  mutate(LouvainClusters = as.factor(LouvainClusters))

# -------------------------
# 2) STANDARDIZE / MERGE KEYS
# -------------------------
# Ensure expected_genome is available everywhere
pre_comp <- pre_comp %>%
  left_join(policy %>% select(barcode, expected_genome), by = "barcode") %>%
  mutate(expected_genome = coalesce(expected_genome.x, expected_genome.y)) %>%
  select(-expected_genome.x, -expected_genome.y)

calls <- calls %>%
  left_join(policy %>% select(barcode, expected_genome), by = "barcode") %>%
  left_join(post_gate %>% select(barcode, keep_postclean), by = "barcode") %>%
  mutate(
    n_reads = suppressWarnings(as.numeric(n_reads)),
    purity  = suppressWarnings(as.numeric(purity)),
    ratio_top1_top2 = suppressWarnings(as.numeric(ratio_top1_top2))
  )

# For composition, ensure contamination_rate exists (design-aware)
# If missing, derive it from expected_frac.
if (!("contamination_rate" %in% names(pre_comp))) {
  if (!("expected_frac" %in% names(pre_comp))) {
    stop("pre_comp missing contamination_rate and expected_frac; cannot compute off-target fraction.")
  }
  pre_comp <- pre_comp %>%
    mutate(expected_frac = suppressWarnings(as.numeric(expected_frac)),
           contamination_rate = 1 - expected_frac)
} else {
  pre_comp <- pre_comp %>%
    mutate(contamination_rate = suppressWarnings(as.numeric(contamination_rate)))
}

pre_comp <- pre_comp %>%
  mutate(total_reads = suppressWarnings(as.numeric(total_reads))) %>%
  filter(!is.na(total_reads), total_reads >= 0, !is.na(contamination_rate))

# -------------------------
# 3) PICK BARNYARD AXES (top 2 genomes by winner evidence)
# -------------------------
top2 <- pre_counts %>%
  group_by(genome) %>%
  summarise(total = sum(n_winner_reads, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1:2) %>%
  pull(genome)

GENOME_X <- top2[1]
GENOME_Y <- top2[2]
message("[Fig1] Barnyard axes: ", GENOME_X, " vs ", GENOME_Y)

# =============================================================================
# PANEL B: Barnyard Cloud (pre-clean)
# =============================================================================

barnyard <- pre_counts %>%
  filter(genome %in% c(GENOME_X, GENOME_Y)) %>%
  mutate(winner_reads_log10 = log10(n_winner_reads + 1)) %>%
  select(barcode, genome, winner_reads_log10) %>%
  pivot_wider(names_from = genome, values_from = winner_reads_log10, values_fill = 0) %>%
  left_join(policy %>% select(barcode, expected_genome), by = "barcode") %>%
  mutate(expected_genome = replace_na(expected_genome, "unknown"))

# Nomenclature:
# - Axes are log10(winner_reads + 1), not raw reads.
# - Color scale is hex bin count (density), also shown on log scale.
pB <- ggplot(barnyard, aes(x = .data[[GENOME_X]], y = .data[[GENOME_Y]])) +
  geom_hex(bins = 120) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Barcodes\n(log10)") +
  theme_minimal(base_size = 11) +
  labs(
    title = "B. The Barnyard Cloud (pre-clean)",
    x = paste0(GENOME_X, " winner reads (log10+1)"),
    y = paste0(GENOME_Y, " winner reads (log10+1)")
  )

# =============================================================================
# PANEL C: Asymmetric contamination ("Soup")
# =============================================================================
# Nomenclature:
# - This is the core asymmetry claim: off-target fraction differs by expected genome.
# - Use the two dominant expected groups (usually Arabidopsis vs B73/maize).
top_expected <- pre_comp %>%
  filter(!is.na(expected_genome), expected_genome != "") %>%
  count(expected_genome, sort = TRUE) %>%
  slice(1:2) %>%
  pull(expected_genome)

E1 <- top_expected[1]
E2 <- top_expected[2]
message("[Fig1] Top expected groups: ", E1, " and ", E2)

# Standardize labels for cleaner legends
comp_FC <- pre_comp %>%
  mutate(
    Plate_Origin = case_when(
      str_detect(expected_genome, regex("At", ignore_case = TRUE)) ~ "Expected: Arabidopsis",
      str_detect(expected_genome, regex("B73|maize|zea", ignore_case = TRUE)) ~ "Expected: Maize",
      TRUE ~ paste0("Expected: ", expected_genome)
    )
  )

pC <- comp_FC %>%
  filter(expected_genome %in% c(E1, E2),
         !is.na(contamination_rate), total_reads > 0) %>%
  ggplot(aes(x = Plate_Origin, y = contamination_rate, fill = Plate_Origin)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = cols_expected_fill) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "C. Asymmetric contamination (“Soup”)",
    subtitle = "Off-target fraction reflects exogenous\nbackground dominated by the higher-biomass genome",
    x = "Expected genome (design)",
    y = "Off-target fraction (1 − expected fraction)"
  )


# =============================================================================
# PANEL D: Contamination Profile (Zoom) — binned + cumulative capture
# =============================================================================
# This panel quantifies C by:
#   - Binning off-target fraction
#   - Counting barcodes per bin per expected group
#   - Computing cumulative percent of barcodes captured up to each bin
#
# IMPORTANT: for stable bin order, we force contam_bin to a factor with fixed levels.
comp_FC <- comp_FC %>%
  mutate(
    contam_bin = cut(
      contamination_rate,
      breaks = breaks_custom,
      labels = labels_custom,
      include.lowest = TRUE,
      right = FALSE
    ),
    contam_bin = factor(contam_bin, levels = labels_custom)
  )


summary_stats <- comp_FC %>%
  filter(total_reads > MIN_TOTAL_READS_ZOOM, !is.na(contam_bin)) %>%
  count(Plate_Origin, contam_bin, name = "count") %>%
  group_by(Plate_Origin) %>%
  arrange(contam_bin, .by_group = TRUE) %>%
  mutate(
    cum_count   = cumsum(count),
    total_cells = sum(count),
    cum_pct     = cum_count / total_cells
  ) %>%
  ungroup()

pD <- ggplot(summary_stats, aes(x = contam_bin, group = Plate_Origin)) +
  geom_col(
    aes(y = count, fill = Plate_Origin),
    position = position_dodge(width = 0.8),
    width = 0.7, alpha = 0.75, color = "black", linewidth = 0.2
  ) +
  geom_line(
    aes(y = cum_count, color = Plate_Origin),
    position = position_dodge(width = 0.8),
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = cum_count, color = Plate_Origin),
    position = position_dodge(width = 0.8),
    size = 2, shape = 21, fill = "white", stroke = 0.3
  ) +
  geom_text(
    aes(y = cum_count, label = scales::percent(cum_pct, accuracy = 1)),
    position = position_dodge(width = 0.8),
    vjust = -0.6,
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_y_log10(labels = comma,
                name = paste0("Barcodes (log10); total_reads > ", MIN_TOTAL_READS_ZOOM)) +
  scale_fill_manual(values = cols_expected_fill) +
  scale_color_manual(values = cols_expected_line) +
  theme_minimal(base_size = 11) +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "D. Contamination profile (zoom)",
    subtitle = "Bars = counts per bin;\nline/points = cumulative count;\nlabels = cumulative % within group",
    x = "Off-target fraction (binned)",
    y = NULL,
    fill = "Physical origin",
    color = "Physical origin"
  )


# =============================================================================
# PANEL E: The Ratio Trap (pre-clean)
# =============================================================================
# X: dominance ratio (log10) · Y: barcode read depth (log10) · Fill: hex bin count
#
# WHY NOT cells_calls$ratio_top1_top2:
#   That column is p_top1/p_top2, where the p_* are MODEL-ESTIMATED mixing
#   proportions. p_top2 == 0 for 130,673 barcodes (78.0%) -- barcodes with no
#   read on any second genome -- so the ratio is Inf and ggplot silently dropped
#   all of them as non-finite. The dropped set is not noise: it is the SHALLOW
#   barcodes (99.8% of those with <=25 reads), i.e. exactly the evidence this
#   panel exists to show.
#
# THE FIX: rebuild the ratio from the integer winner-read counts that panel B
#   already loads, with a Laplace pseudocount:   ratio_pc = (k1 + 1) / (k2 + 1)
#   - defined for 100% of barcodes (nothing dropped)
#   - Spearman 0.981 vs the model ratio where the model ratio exists, so the
#     ordering -- and therefore the published story -- is unchanged
#   - k2 == 0 and p_top2 == 0 select the same barcodes (21 of 167,633 differ,
#     where the model zeroed a 1-2 read count)
#
# CONSEQUENCE, which is the point of the panel: for k2 == 0 the ratio collapses
#   to k1 + 1, so those barcodes lie on the line ratio = reads + 1. That line is
#   the DETECTION CEILING -- the largest dominance ratio observable at a given
#   depth. It is drawn explicitly below; without the label a reader mistakes it
#   for structure in the data. The diagonal is NOT an artifact of the +1: any
#   zero-safe transform (asinh differences included) must place "no contaminant
#   observed" at a depth-dependent value, because that is the real information
#   content -- zero contaminant in 20 reads is weaker evidence than in 5,000.
#
# Panel E therefore uses this count-based pseudocounted ratio rather than
# AmbientMapper's ratio_top1_top2.

# Top-1 / top-2 winner read counts per barcode (k2 = 0 when only one genome present)
ratio_counts <- pre_counts %>%
  filter(!is.na(n_winner_reads)) %>%
  group_by(barcode) %>%
  arrange(desc(n_winner_reads), .by_group = TRUE) %>%
  summarise(
    k1 = first(n_winner_reads),
    k2 = if (n() >= 2) nth(n_winner_reads, 2) else 0,
    .groups = "drop"
  )

ratio_dat <- calls %>%
  select(barcode, n_reads) %>%
  inner_join(ratio_counts, by = "barcode") %>%
  filter(!is.na(n_reads), n_reads > 0) %>%
  mutate(ratio_pc = (k1 + 1) / (k2 + 1))

message("[Fig1] Panel E: ", nrow(ratio_dat), " barcodes plotted; ",
        sum(ratio_dat$k2 == 0), " (",
        round(100 * mean(ratio_dat$k2 == 0), 1),
        "%) have no second-genome read and sit on the detection ceiling")

# Detection ceiling: ratio = reads + 1 (cannot observe more dominance than reads)
ceiling_df <- tibble(x = 10^seq(0, log10(max(ratio_dat$ratio_pc)), length.out = 300)) %>%
  mutate(y = x - 1) %>%
  filter(y >= 1)

pE <- ggplot(ratio_dat, aes(x = ratio_pc, y = n_reads)) +
  geom_hex(bins = 160) +
  scale_fill_viridis_c(option = "viridis", trans = "log10", name = "Barcodes") +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  geom_line(data = ceiling_df, aes(x = x, y = y), inherit.aes = FALSE,
            colour = "grey20", linetype = "22", linewidth = 0.45) +
  annotate("text",
           x = max(ratio_dat$ratio_pc) * 0.75, y = max(ratio_dat$n_reads) * 0.92,
           label = "detection ceiling\nratio = reads + 1",
           hjust = 1, vjust = 1, size = 2.8, colour = "grey20", lineheight = 0.95) +
  theme_minimal(base_size = 11) +
  labs(
    title = "E. The Ratio Trap",
    subtitle = "Barcodes with no second-genome read fall on the ceiling:\nthey look like perfect singlets only because they are shallow",
    x = "Dominance ratio (Top1+1)/(Top2+1)",
    y = "Reads per barcode"
  )


# =============================================================================
# PANEL F: UMAP by Physical Origin (The "Blur")
# =============================================================================
# Goal: Show that physical At and B73 cells overlap in UMAP space due to ambient noise.
# QC gates for the displayed nuclei: total >= 500, pTSS >= 0.2, FRiP >= 0.2 (same in G and H).

F_data <- metaSocReduce %>%
  filter(total >= 500) %>%
  filter(pTSS >= 0.2) %>%
  filter(FRiP >= 0.2)


pF <- F_data %>%
  ggplot(aes(x = umap1, y = umap2, color = Genome)) +
  # Use smaller points and lower alpha to see density overlap
  geom_point(data = F_data %>% filter(Genome == "B73"),
             aes(alpha = 0.25, color = Genome), size = 0.05) +
  geom_point(data = F_data %>% filter(Genome == "At"),
             aes(alpha = 0.25, color = Genome), size = 0.05) +
  scale_color_manual(values = cols_species) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "F. UMAP: Physical Origin (Pre-Clean)",
    subtitle = "Ambient noise causes species to mix in feature space",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Expected Genome"
  )
# =============================================================================
# PANEL G: UMAP by Cluster (The "Structure")
# =============================================================================
# Goal: visualize the Leiden clusters to map them to Panel H.

# Create a palette large enough for all clusters
n_clusters <- length(unique(metaSocReduce$LouvainClusters))
custom_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_clusters)
names(custom_colors) <- as.character(sort(unique(metaSocReduce$LouvainClusters)))

pG <- metaSocReduce %>%
  filter(total >= 500) %>%
  filter(pTSS >= 0.2) %>%
  filter(FRiP >= 0.2) %>%
  ggplot(aes(x = umap1, y = umap2, color = as.factor(LouvainClusters))) +
  geom_point(alpha = 0.5, size = 0.05) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 11) +
  # Override legend dot size for visibility
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "G. UMAP: Leiden Clusters",
    subtitle = "Standard clustering on contaminated data",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cluster"
  )
# =============================================================================
# PANEL H: Cluster Composition (The "Quantification")
# =============================================================================
# Goal: Show that specific clusters are NOT pure (they contain both species).

# 1. Calculate Composition
dataH <- metaSocReduce %>%
  filter(total >= 500) %>%
  filter(pTSS >= 0.2) %>%
  filter(FRiP >= 0.2) %>%
  count(LouvainClusters, Genome) %>%
  group_by(LouvainClusters) %>%
  mutate(
    Total = sum(n),
    FracComposition = n / Total
  ) %>%
  ungroup()

# 2. Plot Stacked Bars
pH <- dataH %>%
  ggplot(aes(y = fct_rev(as.factor(LouvainClusters)), x = FracComposition, fill = Genome)) +
  # Use geom_col for pre-calculated y-values
  # position = "fill" ensures bars span 0-100%, emphasizing purity vs mix
  geom_col(position = "fill", alpha = 0.85, width = 0.8) +

  scale_fill_manual(values = cols_species) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) +

  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "H. Cluster Purity",
    subtitle = "Many clusters contain a mix of species (Artifacts)",
    x = "Fraction of Cells",
    y = "Leiden Cluster",
    fill = "Physical Origin"
  )

# =============================================================================
# EXPORT (individual + assembled)
# =============================================================================
panels <- list(
  list(p = pB, stem = "Fig1B_barnyard_pre",        w_png = 6,   w_pdf = 6),
  list(p = pC, stem = "Fig1C_soup_asymmetry",      w_png = 4.5, w_pdf = 4.5),
  list(p = pD, stem = "Fig1D_contam_profile_zoom", w_png = 5.5, w_pdf = 5.5),
  list(p = pE, stem = "Fig1E_ratio_trap",          w_png = 6,   w_pdf = 6),
  list(p = pF, stem = "Fig1F_umap_mix",            w_png = 4,   w_pdf = 4),
  list(p = pG, stem = "Fig1G_umap_cluster",        w_png = 5,   w_pdf = 5),
  list(p = pH, stem = "Fig1H_cluster_composition", w_png = 5,   w_pdf = 3)
)
for (pn in panels) {
  ggsave(file.path(OUTDIR, paste0(pn$stem, ".png")), pn$p,
         width = pn$w_png, height = 4, dpi = 300)
  ggsave(file.path(OUTDIR, paste0(pn$stem, ".pdf")), pn$p,
         width = pn$w_pdf, height = 4)
}

# Composite preview of the seven panels (the manuscript figure is assembled by hand)
fig <- (pB | pD) / (pC | pE) / (pF | pG) / pH +
  plot_annotation(tag_levels = 'A')
ggsave(file.path(OUTDIR, "Fig1BtoH_preview.png"), fig,
       width = 10, height = 16, dpi = 300)
