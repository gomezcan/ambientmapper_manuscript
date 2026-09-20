#!/usr/bin/env Rscript
# fig3D_rescue_anatomy.R -- Figure 3, panel D only: read depth of the rescue barcode class before
#   and after cleaning (Fig3d_RescueReadLoss). The exploratory panels R1 to R4 and the ND read-fate
#   panel are NOT in the manuscript; their code is kept but off by default (FIG3D_EXPLORATORY=1).
# Inputs : data/processed/scifiATAC_B73_Arabidopsis/SM2v2/decontam_with_design_alpha05_v2/  (WD)
#          data/processed/scifiATAC_B73_Arabidopsis/SM2v2/decontam_without_design_alpha05_v2/  (ND,
#          read only when the exploratory panels are enabled)
# Sources: analysis/_helpers/plotting.R
# Run    : Rscript analysis/fig3_decontamination/fig3D_rescue_anatomy.R   (from the repository root;
#          fig3D_rescue_anatomy.sh wraps it for SLURM)
#
# Purpose
#   Non-circular evidence for the rescue-barcode class. The obvious panel (pre vs post
#   expected-genome fraction) is CIRCULAR and is deliberately NOT built: post-clean expected_frac
#   is exactly 1.000 for every barcode that retains a read, by construction of the design-aware
#   allow-list. Every threshold that DEFINES the class is drawn on the panel it affects as a
#   labelled dashed line, so no cutoff is hidden.
#
#   Production panel D (ships): theme_pubclean + nm_polish, short noun-phrase title, no subtitle
#   and no caption (n, the median removal and the meaning of the dashed floor live in the
#   manuscript legend), 4 x 3 inches to match the other Fig 3 panels.
#     Fig3d_RescueReadLoss  total depth before vs after cleaning.
#
#   Exploratory panels (do not ship; FIG3D_EXPLORATORY=1 builds them): theme_bw + nm_polish_expl,
#   interpretive titles and explanatory captions.
#     R1 how buried they were - PRE-clean expected-genome fraction of the class
#     R2 what the model alone would have done - call_group of the class
#     R3 recovered depth - post-clean read depth of the class
#     R4 what is lost - absolute reads removed, class vs the rest
#     ND read fate - well-genome reads retained, WD vs ND (cut from the figure; the numbers answer
#        "why is the loss worth paying")
#
# NOTE: `rescue_data` below is rebuilt EXACTLY as in fig3.R (same thresholds, same joins, same
#       status cascade) so the class is identical to the one in panel C of the main figure.
# NOTE: panel objects carry a `_rescue` suffix (pD_rescue, pND_rescue) so they cannot collide with
#       the panel objects of fig3.R if both scripts are sourced into one session.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggpubr)     # theme_pubclean(), for the production panel
})

source("analysis/_helpers/plotting.R")   # status_cols, status_levels, stage_cols, nm_polish

# -------------------------
# 0) CONFIG
# -------------------------
DATA     <- "data/processed/scifiATAC_B73_Arabidopsis"
SAMPLE   <- "SM2v2"
INDIR    <- file.path(DATA, "SM2v2", "decontam_with_design_alpha05_v2")     # WD, design-aware
INDIR_ND <- file.path(DATA, "SM2v2", "decontam_without_design_alpha05_v2")  # ND, design-free (exploratory only)
OUTDIR   <- "figures/main/fig3"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Exploratory panels R1 to R4 and the ND read-fate panel: off by default, they are not in the manuscript.
BUILD_EXPLORATORY <- Sys.getenv("FIG3D_EXPLORATORY", "0") %in% c("1", "TRUE", "true")

# thresholds -- identical to fig3.R
MIN_READS_PLOT        <- 10
MIN_READS_POST_CLEAN  <- 100
MIN_ALLOWED_FRAC_POST <- 0.90
PRE_EXPECTED_FRAC_MAX <- 0.50   # the "looked like the wrong genome" boundary

# a conventional scATAC QC depth floor, shown for reference only (panel R3)
QC_DEPTH_REFERENCE    <- 500

RESCUE_LABEL <- "Rescue"

call_group_cols <- c("Single" = "#4daf4a", "Doublet" = "#e41a1c", "Ambiguous" = "#7fb3d5")

# Extension used ONLY by the exploratory panels R1 to R4, which carry explanatory
# captions. The production panel must not use this: interpretation belongs in
# the manuscript legend, not on a main-figure panel.
nm_polish_expl <- nm_polish + theme(
  plot.subtitle = element_text(size = 8, color = "grey30"),
  plot.caption  = element_text(size = 7, color = "grey35", hjust = 0)
)

# -------------------------
# 1) LOAD DATA (WD)
# -------------------------
pre_comp  <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_pre_barcode_composition.tsv.gz")),
                      show_col_types = FALSE)
post_comp <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_post_barcode_composition.tsv.gz")),
                      show_col_types = FALSE)
calls     <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_cells_calls.decontam.tsv.gz")),
                      show_col_types = FALSE)

# -------------------------
# 2) REBUILD rescue_data EXACTLY AS fig3.R DOES
# -------------------------
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
  mutate(call_group = factor(call_group,
                             levels = c("Single", "Doublet", "Ambiguous",
                                        "Low reads", "Indistinguishable")))

pre_meta <- pre_comp %>%
  filter(!is.na(expected_genome)) %>%
  select(barcode, expected_genome, pre_total = total_reads,
         pre_expected_frac = expected_frac, pre_contam = contamination_rate)

post_meta <- post_comp %>%
  select(barcode, post_total = total_reads, post_expected_frac = expected_frac)

rescue_data <- pre_meta %>%
  left_join(post_meta, by = "barcode") %>%
  left_join(calls_clean %>% select(barcode, call_group), by = "barcode") %>%
  mutate(
    post_total         = replace_na(post_total, 0),
    post_expected_frac = replace_na(post_expected_frac, 0),
    reads_removed      = pmax(0, pre_total - post_total),
    pct_removed        = if_else(pre_total > 0, reads_removed / pre_total, 0),

    is_rescue = (pre_expected_frac < PRE_EXPECTED_FRAC_MAX) &
                (post_total >= MIN_READS_POST_CLEAN) &
                (post_expected_frac >= MIN_ALLOWED_FRAC_POST),
    is_heavy  = (pre_contam >= 0.5) & (pct_removed >= 0.5),
    is_clean  = (pct_removed < 0.1),

    status = case_when(
      is_rescue ~ RESCUE_LABEL,
      is_heavy  ~ "Heavily Contaminated",
      is_clean  ~ "Clean (Preserved)",
      TRUE      ~ "Mixed"
    )
  ) %>%
  filter(pre_total > MIN_READS_PLOT) %>%
  mutate(status = factor(status, levels = status_levels))

rescue <- rescue_data %>% filter(is_rescue)
N_RESCUE <- nrow(rescue)

# hard stop if the class drifts away from the 400 barcodes quoted in the manuscript legend
stopifnot(N_RESCUE == 400)
message("Rescue class: ", N_RESCUE, " barcodes (of ", nrow(rescue_data), " with pre_total > ",
        MIN_READS_PLOT, ")")

# comparator for R1: barcodes that clear the SAME post-clean gate but are not rescue
gate_pass_other <- rescue_data %>%
  filter(!is_rescue,
         post_total >= MIN_READS_POST_CLEAN,
         post_expected_frac >= MIN_ALLOWED_FRAC_POST)
N_OTHER <- nrow(gate_pass_other)
message("Comparator (non-rescue, post-clean gate passed): ", N_OTHER, " barcodes")

# #############################################################################
# EXPLORATORY PANELS R1 to R4 -- not in the manuscript, built only when FIG3D_EXPLORATORY=1
# #############################################################################
if (BUILD_EXPLORATORY) {

  # ===========================================================================
  # PANEL R1 -- how buried they were
  # ===========================================================================
  r1_dat <- bind_rows(
    rescue          %>% mutate(grp = sprintf("%s (n = %s)", RESCUE_LABEL, comma(N_RESCUE))),
    gate_pass_other %>% mutate(grp = sprintf("Other gate-passing (n = %s)", comma(N_OTHER)))
  ) %>%
    mutate(grp = factor(grp, levels = c(sprintf("%s (n = %s)", RESCUE_LABEL, comma(N_RESCUE)),
                                        sprintf("Other gate-passing (n = %s)", comma(N_OTHER)))))

  r1_cols <- setNames(c("#2c7bb6", "grey55"), levels(r1_dat$grp))

  med_pre_rescue <- median(rescue$pre_expected_frac)
  q_pre_rescue   <- quantile(rescue$pre_expected_frac, c(0.25, 0.75))

  pR1 <- ggplot(r1_dat, aes(x = pre_expected_frac, fill = grp)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.02,
                   position = "identity", alpha = 0.6, colour = NA) +
    geom_vline(xintercept = PRE_EXPECTED_FRAC_MAX, linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    geom_vline(xintercept = med_pre_rescue, linetype = "dotted",
               colour = "#08306b", linewidth = 0.4) +
    annotate("text", x = PRE_EXPECTED_FRAC_MAX - 0.015, y = Inf,
             label = "class boundary\npre expected fraction < 0.50",
             hjust = 1, vjust = 1.15, size = 2.5, colour = "black", lineheight = 0.95) +
    annotate("text", x = med_pre_rescue + 0.015, y = Inf,
             label = sprintf("rescue median %.3f\n(IQR %.3f to %.3f)",
                             med_pre_rescue, q_pre_rescue[1], q_pre_rescue[2]),
             hjust = 0, vjust = 3.6, size = 2.5, colour = "#08306b", lineheight = 0.95) +
    # coord_cartesian, not scale limits: scale limits would silently DROP the edge
    # bins instead of clipping them, and this panel must not lose barcodes.
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = percent_format(accuracy = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_fill_manual(values = r1_cols, name = NULL) +
    theme_bw(base_size = 11) + nm_polish_expl +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7.5)) +
    labs(title = "Rescue barcodes were deeply buried, not borderline",
         subtitle = "PRE-clean expected-genome fraction (design well genome)",
         x = "Pre-clean expected-genome fraction",
         y = "Density",
         caption = paste0("The two groups are disjoint at 0.50 by definition of the class. ",
                          "The informative feature is WHERE\ninside (0, 0.50) the rescue ",
                          "barcodes sit. They pile up near 0.1, not just under the cutoff."))

  # ===========================================================================
  # PANEL R2 -- what the model alone would have done
  # ===========================================================================
  r2_dat <- rescue %>%
    filter(!is.na(call_group)) %>%
    count(call_group, .drop = FALSE) %>%
    filter(n > 0 | call_group %in% c("Single", "Doublet", "Ambiguous")) %>%
    filter(call_group %in% c("Single", "Doublet", "Ambiguous")) %>%
    mutate(pct = 100 * n / N_RESCUE)

  n_not_single <- sum(r2_dat$n[r2_dat$call_group != "Single"])
  pct_not_single <- 100 * n_not_single / N_RESCUE

  pR2 <- ggplot(r2_dat, aes(x = call_group, y = n, fill = call_group)) +
    geom_col(width = 0.68, colour = "black", linewidth = 0.25, alpha = 0.9) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
              vjust = -0.25, size = 3, lineheight = 0.95) +
    annotate("segment", x = 1.6, xend = 3.4, y = 232, yend = 232,
             linewidth = 0.4, colour = "grey20") +
    annotate("text", x = 2.5, y = 240,
             label = sprintf("discarded or mislabelled without the design:\n%d of %d barcodes (%.1f%%)",
                             n_not_single, N_RESCUE, pct_not_single),
             vjust = 0, size = 2.8, fontface = "bold", colour = "grey15", lineheight = 0.95) +
    scale_y_continuous(limits = c(0, 300), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = call_group_cols, guide = "none") +
    theme_bw(base_size = 11) + nm_polish_expl +
    labs(title = "Half the class is not a clean singlet to the model",
         subtitle = sprintf("AmbientMapper call for the %d rescue barcodes, design ignored",
                            N_RESCUE),
         x = "Genotype call (model only)", y = "Barcodes",
         caption = paste0("A design-blind pipeline keeps only the Single calls. ",
                          "Doublet and Ambiguous barcodes are\ndropped or given the wrong ",
                          "genotype even though the well genome is known."))

  # ===========================================================================
  # PANEL R3 -- recovered depth
  # ===========================================================================
  med_post   <- median(rescue$post_total)
  n_above500 <- sum(rescue$post_total >= QC_DEPTH_REFERENCE)
  pct_above500 <- 100 * n_above500 / N_RESCUE

  pR3 <- ggplot(rescue, aes(x = post_total)) +
    geom_histogram(bins = 34, fill = "#2c7bb6", colour = "white",
                   linewidth = 0.15, alpha = 0.9) +
    geom_vline(xintercept = MIN_READS_POST_CLEAN, linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    geom_vline(xintercept = QC_DEPTH_REFERENCE, linetype = "dashed",
               colour = "#b2182b", linewidth = 0.4) +
    geom_vline(xintercept = med_post, linetype = "dotted",
               colour = "#08306b", linewidth = 0.4) +
    # all three annotations anchored hjust = 0 to the RIGHT of their line, so none
    # can be clipped by the left panel edge
    annotate("text", x = MIN_READS_POST_CLEAN * 1.04, y = Inf,
             label = sprintf("selection floor\nmin_reads_post_clean = %d", MIN_READS_POST_CLEAN),
             hjust = 0, vjust = 1.15, size = 2.4, lineheight = 0.95) +
    annotate("text", x = QC_DEPTH_REFERENCE * 1.06, y = Inf,
             label = sprintf("conventional QC floor = %d reads\n%d of %d (%.1f%%) clear it",
                             QC_DEPTH_REFERENCE, n_above500, N_RESCUE, pct_above500),
             hjust = 0, vjust = 1.15, size = 2.4, colour = "#b2182b", lineheight = 0.95) +
    # vjust 5.2 clears BOTH lines of the two-line selection-floor label above it
    annotate("text", x = med_post * 1.06, y = Inf,
             label = sprintf("median %.1f reads", med_post),
             hjust = 0, vjust = 5.2, size = 2.4, colour = "#08306b") +
    scale_x_log10(breaks = c(100, 200, 500, 1000, 2000), labels = comma) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    # clip, do not drop: scale limits would remove barcodes from the histogram
    coord_cartesian(xlim = c(80, 4000)) +
    theme_bw(base_size = 11) + nm_polish_expl +
    labs(title = "Recovered depth is real but modest",
         subtitle = sprintf("Post-clean read depth of the %d rescue barcodes", N_RESCUE),
         x = "Post-clean reads per barcode (log10)", y = "Barcodes",
         caption = paste0("Most rescue barcodes sit between the 100-read selection floor and ",
                          "500 reads. They are\nrecovered nuclei of modest depth, not deep ones."))

  # ===========================================================================
  # PANEL R4 -- what is lost
  # ===========================================================================
  r4_stats <- rescue_data %>%
    group_by(status) %>%
    summarise(n = n(),
              med_removed = median(reads_removed),
              med_pct     = median(pct_removed),
              .groups = "drop")

  # n and the medians live in the AXIS labels, not as in-panel text: the "Clean
  # (Preserved)" box sits at y = 1 (median 0 reads removed) and any in-panel label
  # at that height collides with it.
  r4_axis_labels <- setNames(
    sprintf("%s\nn = %s\nmedian %s reads\n(%.0f%% of pre)",
            str_replace(as.character(r4_stats$status), " \\(", "\n("),
            comma(r4_stats$n), comma(r4_stats$med_removed), 100 * r4_stats$med_pct),
    as.character(r4_stats$status))

  pR4 <- ggplot(rescue_data, aes(x = status, y = reads_removed + 1, fill = status)) +
    geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.55) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3,
                 fill = "white", colour = "grey20") +
    scale_y_log10(labels = comma, breaks = c(1, 10, 100, 1e3, 1e4, 1e5),
                  expand = expansion(mult = c(0.05, 0.06))) +
    scale_fill_manual(values = status_cols, guide = "none") +
    scale_x_discrete(labels = r4_axis_labels) +
    theme_bw(base_size = 11) + nm_polish_expl +
    theme(axis.text.x = element_text(size = 7, lineheight = 0.95)) +
    labs(title = "Rescue barcodes lose the most reads in absolute terms",
         subtitle = "Reads removed per barcode, by cleaning outcome",
         x = "Cleaning outcome", y = "Reads removed (+1, log10)",
         caption = paste0("Heavily contaminated barcodes lose a larger FRACTION of their reads ",
                          "but far fewer reads,\nbecause they are shallow. Rescue barcodes are ",
                          "deep and mostly ambient."))

  # exploratory composite + individual panels, at the per-panel size of the 10 x 8 composite
  fig_expl <- (pR1 | pR2) / (pR3 | pR4) +
    plot_annotation(tag_levels = 'A',
                    title = sprintf("Anatomy of the rescue barcode class (%s, design-aware run)", SAMPLE))

  ggsave(file.path(OUTDIR, "FigS_rescue_R1_buried.pdf"),       pR1, width = 5.0, height = 4.0)
  ggsave(file.path(OUTDIR, "FigS_rescue_R2_calls.pdf"),        pR2, width = 5.0, height = 4.0)
  ggsave(file.path(OUTDIR, "FigS_rescue_R3_depth.pdf"),        pR3, width = 5.0, height = 4.0)
  ggsave(file.path(OUTDIR, "FigS_rescue_R4_readsremoved.pdf"), pR4, width = 5.0, height = 4.0)

  ggsave(file.path(OUTDIR, "FigS_rescue_anatomy.pdf"), fig_expl, width = 10, height = 8)
  ggsave(file.path(OUTDIR, "FigS_rescue_anatomy.png"), fig_expl, width = 10, height = 8, dpi = 300)
}

# #############################################################################
# PRODUCTION PANEL D -- Fig3d_RescueReadLoss is the only panel of this script in the manuscript.
#
# Styling deliberately differs from R1 to R4 above: theme_pubclean + nm_polish, short
# noun-phrase title, no subtitle and no explanatory caption. Interpretation goes in the
# manuscript legend.
# #############################################################################

# =============================================================================
# PANEL D -- Fig3d_RescueReadLoss : total depth before and after cleaning
# =============================================================================
d_dat <- rescue %>%
  select(barcode, `Raw (Pre)` = pre_total, `Cleaned (Post)` = post_total) %>%
  pivot_longer(-barcode, names_to = "Stage", values_to = "reads") %>%
  mutate(Stage = factor(Stage, levels = c("Raw (Pre)", "Cleaned (Post)")))

d_med <- d_dat %>% group_by(Stage) %>% summarise(med = median(reads), .groups = "drop")
med_pre_total  <- d_med$med[d_med$Stage == "Raw (Pre)"]
med_post_total <- d_med$med[d_med$Stage == "Cleaned (Post)"]
med_pct_removed <- median(rescue$pct_removed)

# guards: the panel can never silently drift off the numbers quoted in the manuscript legend
stopifnot(med_pre_total == 1607.5,
          med_post_total == 153.5,
          round(100 * med_pct_removed, 1) == 88.1)

pD_rescue <- ggplot(d_dat, aes(x = Stage, y = reads, fill = Stage)) +
  geom_violin(scale = "width", width = 0.8, colour = NA, alpha = 0.75) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3,
               fill = "white", colour = "grey20") +
  geom_hline(yintercept = MIN_READS_POST_CLEAN, linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  # Short on-panel label only. The parameter name it stands for
  # (min_reads_post_clean) is spelled out in the manuscript legend, not here.
  annotate("text", x = 0.45, y = MIN_READS_POST_CLEAN * 0.82, hjust = 0, vjust = 1,
           label = sprintf("min. %d reads", MIN_READS_POST_CLEAN), size = 2.4) +
  # accuracy = 0.1: the medians are 1,607.5 and 153.5, do not round away the half read.
  # x nudged +0.11 to clear the 0.16-wide inner boxplot rather than sit on its edge.
  geom_text(data = d_med,
            aes(x = as.numeric(Stage) + 0.11, y = med, label = comma(med, accuracy = 0.1)),
            inherit.aes = FALSE, hjust = 0, vjust = -0.4, size = 2.6, fontface = "bold") +
  scale_y_log10(labels = comma, breaks = c(100, 1000, 10000),
                expand = expansion(mult = c(0.10, 0.08))) +
  scale_fill_manual(values = stage_cols, guide = "none") +
  theme_pubclean(base_size = 11) + nm_polish +
  # Short labels only. n, the median removal and the meaning of the dashed
  # floor all live in the manuscript legend, so the panel carries no subtitle.
  labs(title = "Rescue Barcode Read Loss", x = NULL, y = "Reads per barcode")

# 4 x 3 inches, matching the other Fig 3 panel exports. No composite is built here:
# the main figure is assembled by hand.
ggsave(file.path(OUTDIR, "Fig3d_RescueReadLoss.pdf"), pD_rescue, width = 4, height = 3)
ggsave(file.path(OUTDIR, "Fig3d_RescueReadLoss.png"), pD_rescue, width = 4, height = 3, dpi = 300)

# #############################################################################
# EXPLORATORY: ND read fate -- well-genome reads retained, WD against ND. Cut from the figure;
# built only when FIG3D_EXPLORATORY=1 (the only consumer of the ND input).
# #############################################################################
if (BUILD_EXPLORATORY) {

  # ND (design-free) post-clean per-genome winner counts. ND leaves `expected_genome` empty in
  # its composition tables, so the well genome comes from the WD side and is joined in.
  nd_post_counts <- read_tsv(file.path(INDIR_ND, paste0(SAMPLE, "_post_barcode_genome_counts.tsv.gz")),
                             show_col_types = FALSE)

  # Mode palette. Deliberately a neutral slate pair so it collides with NEITHER the stage
  # vocabulary (pink/green) NOR the species vocabulary (At blue / B73 red) NOR the status
  # palette. The reference column (WD) is the lighter value so it recedes, and the measured
  # column (ND) is dominant.
  mode_cols <- c("WD" = "#B0BEC5", "ND" = "#37474F")

  # ND retains only reads from its own called genome, so a rescue barcode whose ND
  # call is the other genome keeps ZERO well-genome reads. Absent (barcode, genome)
  # rows in the ND post counts are true zeros, hence replace_na(0).
  nd_expected_reads <- nd_post_counts %>%
    inner_join(rescue %>% select(barcode, expected_genome), by = "barcode") %>%
    filter(genome == expected_genome) %>%
    select(barcode, nd_exp_reads = n_winner_reads)

  nd_dat <- rescue %>%
    select(barcode, WD = post_total) %>%
    left_join(nd_expected_reads, by = "barcode") %>%
    mutate(ND = replace_na(nd_exp_reads, 0)) %>%
    select(barcode, WD, ND) %>%
    pivot_longer(-barcode, names_to = "Mode", values_to = "exp_reads") %>%
    mutate(Mode = factor(Mode, levels = c("WD", "ND")))

  wd_exp_total <- sum(nd_dat$exp_reads[nd_dat$Mode == "WD"])
  nd_exp_total <- sum(nd_dat$exp_reads[nd_dat$Mode == "ND"])
  n_nd_zero    <- sum(nd_dat$exp_reads[nd_dat$Mode == "ND"] == 0)
  pct_nd_lost  <- 100 * (1 - nd_exp_total / wd_exp_total)

  # guards on every number annotated on the panel
  stopifnot(wd_exp_total == 81754, nd_exp_total == 6631,
            n_nd_zero == 369, round(pct_nd_lost, 1) == 91.9,
            median(nd_dat$exp_reads[nd_dat$Mode == "WD"]) == 153.5,
            median(nd_dat$exp_reads[nd_dat$Mode == "ND"]) == 0)

  pND_rescue <- ggplot(nd_dat, aes(x = Mode, y = exp_reads + 1, fill = Mode)) +
    geom_violin(scale = "width", width = 0.8, colour = NA, alpha = 0.85) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3,
                 fill = "white", colour = "grey20") +
    # the measured result, placed above the ND violin
    annotate("text", x = 2, y = 1900, hjust = 0.5, vjust = 0.5, size = 2.5,
             colour = "grey15", lineheight = 0.95,
             label = sprintf("zero well-genome reads\nfor %d of %d barcodes",
                             n_nd_zero, N_RESCUE)) +
    # the aggregate, in the empty lower-left region under the WD violin
    annotate("text", x = 0.45, y = 1.55, hjust = 0, vjust = 0.5, size = 2.5,
             colour = "grey15", lineheight = 0.95,
             label = sprintf("ND discards %.1f%% of well-genome reads\n(%s of %s)",
                             pct_nd_lost, comma(wd_exp_total - nd_exp_total),
                             comma(wd_exp_total))) +
    scale_y_log10(labels = comma, breaks = c(1, 10, 100, 1000),
                  expand = expansion(mult = c(0.10, 0.16))) +
    scale_fill_manual(values = mode_cols, guide = "none") +
    theme_pubclean(base_size = 11) + nm_polish +
    labs(title = "Nuclear Read Fate by Cleaning Mode",
         # wrapped to two lines: as a single line this overruns the 4-inch panel and
         # the clipped remainder is exactly the honesty caveat, which must stay visible
         subtitle = paste0("WD retains all well-genome reads by construction.\n",
                           "The result is the ND column."),
         x = "Cleaning mode", y = "Well-genome reads (+1)")

  # exported under the exploratory prefix, not as a Fig3* production part
  ggsave(file.path(OUTDIR, "FigS_rescue_ND_readfate.pdf"), pND_rescue, width = 4, height = 3)
  ggsave(file.path(OUTDIR, "FigS_rescue_ND_readfate.png"), pND_rescue, width = 4, height = 3, dpi = 300)
}

# -------------------------
# CONSOLE SUMMARY (numbers quoted on the panels)
# -------------------------
cat("\n================ RESCUE CLASS SUMMARY ================\n")
cat(sprintf("N rescue                       : %d\n", N_RESCUE))
cat(sprintf("N comparator (gate, non-rescue) : %d\n", N_OTHER))
cat("\nExpected genome of the rescue class:\n"); print(table(rescue$expected_genome))

cat("\n--------- SHIPPING PANEL: Fig3d_RescueReadLoss ---------\n")
cat(sprintf("D median depth  : Raw (Pre) %s, Cleaned (Post) %s reads\n",
            comma(med_pre_total, accuracy = 0.1), comma(med_post_total, accuracy = 0.1)))
cat(sprintf("D median removed: %.1f%%\n", 100 * med_pct_removed))

if (BUILD_EXPLORATORY) {
  cat("\n--------- EXPLORATORY PANELS (not in the manuscript) ---------\n")
  cat("\nR1 pre-clean expected fraction (rescue):\n"); print(summary(rescue$pre_expected_frac))
  cat("\nR2 call_group:\n"); print(r2_dat)
  cat(sprintf("   not Single: %d of %d (%.1f%%)\n", n_not_single, N_RESCUE, pct_not_single))
  cat("\nR3 post-clean depth (rescue):\n"); print(summary(rescue$post_total))
  cat(sprintf("   >= %d reads: %d (%.1f%%)\n", QC_DEPTH_REFERENCE, n_above500, pct_above500))
  cat("\nR4 reads removed by status:\n"); print(r4_stats)
  cat(sprintf("\nND well-genome reads retained: WD %s, ND %s\n",
              comma(wd_exp_total), comma(nd_exp_total)))
  cat(sprintf("ND zero well-genome reads : %d of %d barcodes\n", n_nd_zero, N_RESCUE))
  cat(sprintf("ND discards               : %.1f%% (%s of %s)\n",
              pct_nd_lost, comma(wd_exp_total - nd_exp_total), comma(wd_exp_total)))
}
cat("======================================================\n")
