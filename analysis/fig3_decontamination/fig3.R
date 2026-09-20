#!/usr/bin/env Rscript
# fig3.R -- Figure 3, panels A, B, C and E (decontamination against the plate-design ground truth).
#   A barnyard restoration (pre vs post winner reads), B off-target fraction pre vs post by call
#   group (the ECDF is the manuscript panel, the density version is exported for comparison),
#   C reads removed vs original reads coloured by cleaning outcome, E specificity of removal
#   (reads removed per genome by well type). Panel D is built by fig3D_rescue_anatomy.R.
# Inputs : data/processed/scifiATAC_B73_Arabidopsis/SM2v2/decontam_with_design_alpha05_v2/
# Sources: analysis/_helpers/plotting.R
# Run    : Rscript analysis/fig3_decontamination/fig3.R   (from the repository root; fig3.sh wraps it for SLURM)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(viridis)
  library(ggpubr)
  library(scattermore)   # panel C: hundreds of thousands of points
})

source("analysis/_helpers/plotting.R")   # status_cols, stage_cols, cols_species, nm_polish

# -------------------------
# 0) CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SAMPLE <- "SM2v2"
INDIR  <- file.path(DATA, "SM2v2", "decontam_with_design_alpha05_v2")   # WD, design-aware run
OUTDIR <- "figures/main/fig3"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_READS_PLOT        <- 10     # barcodes with more pre-clean reads than this are plotted
MIN_READS_POST_CLEAN  <- 100    # rescue class: minimum post-clean reads
MIN_ALLOWED_FRAC_POST <- 0.90   # rescue class: minimum post-clean expected-genome fraction

# -------------------------
# 1) LOAD DATA
# -------------------------
policy <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_barcode_policy.tsv.gz")), show_col_types = FALSE) %>%
  mutate(barcode = as.character(barcode), expected_genome = na_if(as.character(expected_genome), ""))

pre_counts <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_pre_barcode_genome_counts.tsv.gz")), show_col_types = FALSE)
post_counts <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_post_barcode_genome_counts.tsv.gz")), show_col_types = FALSE)

pre_comp <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_pre_barcode_composition.tsv.gz")), show_col_types = FALSE)
post_comp <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_post_barcode_composition.tsv.gz")), show_col_types = FALSE)

calls <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_cells_calls.decontam.tsv.gz")), show_col_types = FALSE)
post_gate <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_barcode_postclean.tsv.gz")), show_col_types = FALSE)

# -------------------------
# 2) DATA PREP & STATUS DEFINITION
# -------------------------

# A. Prepare Calls with Groups (current call vocabulary first, older vocabulary kept for compatibility)
calls_clean <- calls %>%
  mutate(
    call_group = case_when(
      call %in% c("single_clean", "dirty_singlet") ~ "Single",
      call %in% c("doublet", "weak_doublet")       ~ "Doublet",
      call == "ambiguous"                          ~ "Ambiguous",
      call == "low_reads"                          ~ "Low reads",
      call == "doublet_confident"                  ~ "Doublet",
      call == "indistinguishable"                  ~ "Indistinguishable",
      str_detect(call, "^ambiguous")               ~ "Ambiguous",
      TRUE                                         ~ "Ambiguous"
    )
  ) %>%
  mutate(call_group = factor(call_group, levels = c("Single", "Doublet", "Ambiguous", "Low reads", "Indistinguishable")))

# B. Cleaning outcome per barcode (pre vs post comparison)
#    Expected fraction and totals come from the pre and post composition tables.
#    pre_comp carries expected_genome from the design-aware run; a design-free run leaves it empty.
pre_meta <- pre_comp %>%
  filter(!is.na(expected_genome)) %>%
  select(barcode, expected_genome, pre_total=total_reads, pre_expected_frac=expected_frac, pre_contam=contamination_rate)

post_meta <- post_comp %>%
  select(barcode, post_total=total_reads, post_expected_frac=expected_frac)

rescue_data <- pre_meta %>%
  left_join(post_meta, by="barcode") %>%
  left_join(calls_clean %>% select(barcode, call_group), by="barcode") %>%
  mutate(
    post_total = replace_na(post_total, 0),
    post_expected_frac = replace_na(post_expected_frac, 0),
    reads_removed = pmax(0, pre_total - post_total),
    pct_removed = if_else(pre_total > 0, reads_removed / pre_total, 0),

    # --- STATUS DEFINITIONS ---
    # 1. Rescue:
    #    Initially looked like the wrong genome (< 50% expected-genome reads), but post-cleaning is
    #    high quality (>= 90% expected) and has enough reads (>= 100).
    is_rescue = (pre_expected_frac < 0.5) &
      (post_total >= MIN_READS_POST_CLEAN) &
      (post_expected_frac >= MIN_ALLOWED_FRAC_POST),

    # 2. Heavily Contaminated:
    #    Had high contamination (>50%) AND we removed most of it (>50% reads lost).
    is_heavy = (pre_contam >= 0.5) & (pct_removed >= 0.5),

    # 3. Clean (Preserved):
    #    We removed very little (<10%), implying it was already clean.
    is_clean = (pct_removed < 0.1),

    # Assign Labels (status_cols in analysis/_helpers/plotting.R uses the same strings)
    status = case_when(
      is_rescue ~ "Rescue",
      is_heavy  ~ "Heavily Contaminated",
      is_clean  ~ "Clean (Preserved)",
      TRUE      ~ "Mixed"
    )
  ) %>%
  filter(pre_total > MIN_READS_PLOT)


# -------------------------
# PANEL A: Barnyard (Pre vs Post)
# -------------------------
# (Using simple pre/post facets, no call grouping needed here)
top2 <- pre_counts %>%
  group_by(genome) %>% summarise(tot=sum(n_winner_reads)) %>%
  arrange(desc(tot)) %>% slice(1:2) %>% pull(genome)
GX <- top2[1]; GY <- top2[2]

barn_all <- bind_rows(
  pre_counts %>% mutate(Stage="Pre"),
  post_counts %>% mutate(Stage="Post")
) %>%
  filter(genome %in% c(GX, GY)) %>%
  pivot_wider(names_from=genome, values_from=n_winner_reads, values_fill=0) %>%
  mutate(Stage = factor(Stage, levels=c("Pre", "Post")))

pA <- ggplot(barn_all, aes(x = .data[[GX]]+1, y = .data[[GY]]+1)) +
  geom_hex(bins=100) +
  scale_x_log10(labels = trans_format("log10", label_math())) +
  scale_y_log10(labels = trans_format("log10", label_math())) +
  scale_fill_viridis_c(option="magma", trans="log10", name = "Barcodes") +
  facet_wrap(~Stage) +
  theme_pubclean(base_size=11) + nm_polish +
  labs(title="Barnyard Restoration",
       x=paste0(GX, " winner reads (+1)"),
       y=paste0(GY, " winner reads (+1)"))

# -------------------------
# PANEL B: Contamination Profile (Faceted by Call)
# -------------------------
# Off-target fraction pre vs post cleaning, faceted by call group. Low-read and
# indistinguishable barcodes are excluded here and in panels C and E.
comp_long <- bind_rows(
  pre_comp %>% mutate(Stage="Raw (Pre)"),
  post_comp %>% mutate(Stage="Cleaned (Post)")
) %>%
  left_join(calls_clean %>% select(barcode, call_group), by="barcode") %>%
  filter(!is.na(call_group), !call_group %in% c("Indistinguishable", "Low reads")) %>%
  filter(total_reads > MIN_READS_PLOT)

# Density version (exported for comparison, not the manuscript panel)
pB <- ggplot(comp_long, aes(x = 1 - expected_frac, fill = Stage)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.02, alpha = 0.55, position = "identity") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = stage_cols) +
  facet_wrap(~call_group, scales = "free_y") +
  theme_pubclean(base_size = 11) + nm_polish +
  labs(title = "Contamination Shift",
       x = "Off-target fraction", y = "Density")

# ECDF version: this is manuscript panel B
pB_ecdf <- ggplot(comp_long, aes(x = 1 - expected_frac, color = Stage)) +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = stage_cols) +
  facet_wrap(~call_group) +
  theme_pubclean(base_size = 11) + nm_polish +
  labs(title = "Contamination Shift (ECDF)",
       x = "Off-target fraction", y = "Cumulative fraction of barcodes")

# -------------------------
# PANEL C: Reads removed vs original reads (Faceted + Status Color)
# -------------------------
pC <- rescue_data %>%
  filter(!is.na(call_group), !call_group %in% c("Indistinguishable", "Low reads")) %>%
  ggplot(aes(x = pre_total, y = reads_removed)) +

  # Main Points colored by Status
  geom_scattermore(aes(color = status),
                   pixels = c(250, 250),
                   alpha = 0.5, size = 1) +

  geom_point(data= rescue_data %>% filter(status=="Rescue"),
             aes(color = status), alpha = 0.5, size = 0.2) +
  # Identity Line (100% removal)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +

  # Highlight rescue barcodes with black rings
  geom_point(data = subset(rescue_data, is_rescue),
             shape = 21, color = "black", fill = NA, size = 0.3, stroke = 0.1) +

  scale_x_log10(labels = trans_format("log10", label_math())) +
  scale_y_log10(labels = trans_format("log10", label_math())) +
  scale_color_manual(values = status_cols, name = "Status") +
  facet_wrap(~call_group) +
  theme_bw(base_size = 11) + nm_polish +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(title = "Reads Removed vs Original",
       subtitle = "Status shows impact of cleaning per genotype call",
       x = "Original reads", y = "Removed reads")

# -------------------------
# PANEL E: Specificity (Faceted by Call)
# -------------------------
# Reads removed per (barcode x genome), aggregated by call group, well type (expected genome) and
# removed genome. Note: in the design-aware run every barcode carries a single-genome allowed_set,
# so removal of expected-genome reads is zero by construction; the panel reports how many
# off-target reads were removed per well type.
diff_counts <- pre_counts %>%
  select(barcode, genome, pre=n_winner_reads) %>%
  full_join(post_counts %>% select(barcode, genome, post=n_winner_reads), by=c("barcode","genome")) %>%
  mutate(removed = replace_na(pre,0) - replace_na(post,0)) %>%
  left_join(policy %>% select(barcode, expected_genome), by="barcode") %>%
  left_join(calls_clean %>% select(barcode, call_group), by="barcode") %>%
  filter(!is.na(expected_genome), !is.na(call_group), !call_group %in% c("Indistinguishable", "Low reads"))

# Aggregate per Call Group -> Expected Genome -> Removed Genome
spec_stats <- diff_counts %>%
  group_by(call_group, expected_genome, genome) %>%
  summarise(total_removed = sum(removed), .groups="drop")

pE_spec <- ggplot(spec_stats, aes(x = expected_genome, y = total_removed, fill = genome)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.25, alpha = 0.85) +
  scale_fill_manual(values = cols_species, name = "Removed genome") +
  scale_y_continuous(
    trans  = pseudo_log_trans(base = 10),
    breaks = c(0, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7),
    labels = comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(~call_group) +
  theme_pubclean(base_size = 11) + nm_polish +
  labs(title = "Specificity", x = "Well type (expected genome)",
       y = "Reads removed (log10)")



# -------------------------
# Barcode depth distribution (Pre vs Post) per expected genome
# -------------------------
# Constructed for inspection only: not exported and not a manuscript panel.
# Goal:
#   Quantify how many GENOTYPED barcodes fall into read-depth bins
#   (<10, 10–100, 100–500, >500) before vs after cleaning, stratified by expected genome.
#
# Notes:
#   - Uses total_reads from *_barcode_composition.tsv.gz (pre/post)
#   - Restricts to barcodes present in calls_clean (i.e., genotyped/called)
#   - Uses expected_genome carried by the composition tables (design-aware grouping)

depth_bins <- function(x) {
  case_when(
    x < 10 ~ "<10",
    x >= 10  & x < 100 ~ "10–100",
    x >= 100 & x < 500 ~ "100–500",
    x >= 500 ~ ">500",
    TRUE ~ NA_character_
  )
}

comp_depth_long <- bind_rows(
  pre_comp  %>% mutate(Stage = "Raw (Pre)"),
  post_comp %>% mutate(Stage = "Cleaned (Post)")
) %>%
  # keep only barcodes that were genotyped (i.e., appear in calls)
  inner_join(calls_clean %>% select(barcode, call_group), by = "barcode") %>%
  filter(!is.na(expected_genome)) %>%
  mutate(
    read_bin = depth_bins(total_reads),
    read_bin = factor(read_bin, levels = c("<10", "10–100", "100–500", ">500")),
    Stage    = factor(Stage, levels = c("Raw (Pre)", "Cleaned (Post)"))
  ) %>%
  filter(!is.na(read_bin))

depth_stats <- comp_depth_long %>%
  group_by(expected_genome, Stage, read_bin) %>%
  summarise(n_barcodes = n_distinct(barcode), .groups = "drop")

p_depth_bins <- ggplot(depth_stats, aes(x = read_bin, y = n_barcodes, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "black", alpha = 0.85) +
  scale_y_sqrt(labels = comma) +
  scale_fill_manual(values = stage_cols) +
  facet_wrap(~ expected_genome) +
  theme_pubclean(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "Barcode Read-Depth Distribution (Genotyped Barcodes)",
    subtitle = "Counts of genotyped barcodes by total reads before vs after cleaning",
    x = "Total reads per barcode (binned)",
    y = "Barcodes (sqrt)"
  )

# -------------------------
# ASSEMBLY
# -------------------------
# Composite preview only: the manuscript figure is assembled by hand, with panel D from
# fig3D_rescue_anatomy.R. Panel B in the manuscript is the ECDF.
fig3 <- (pA | pB_ecdf) / (pC | pE_spec) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "E")))

ggsave(file.path(OUTDIR, "Fig3_preview_ABCE.png"), fig3, width = 8, height = 7, dpi = 300)
ggsave(file.path(OUTDIR, "Fig3_preview_ABCE.pdf"), fig3, width = 8, height = 7)

# Individual panels (manuscript letters: A, B = the ECDF, C, E)
ggsave(file.path(OUTDIR, "Fig3a_Banyard_restoration.pdf"), pA, width = 3.5, height = 3)
ggsave(file.path(OUTDIR, "Fig3b_ContaminationShift.pdf"), pB, width = 5.5, height = 3)
ggsave(file.path(OUTDIR, "Fig3b_ContaminationShift_ecdf.pdf"), pB_ecdf, width = 5.5, height = 3)
ggsave(file.path(OUTDIR, "Fig3c_ReadsRemoved.pdf"), pC, width = 4, height = 3)
ggsave(file.path(OUTDIR, "Fig3E_Specificity.pdf"), pE_spec, width = 4, height = 3)
