#!/usr/bin/env Rscript
# Supplementary Fig S6: genotyping parameter validation on the maize root library (Phase 3,
# Root1 sub1k_B subsample, 26-genome reference, B73-only library).
#   A: configuration parameter table (Phase 3 configs)
#   B: B73 dominance per AmbientMapper call class (19 configs x 3 classes; color = g1 == B73 / n_total)
#   C: call-class share per config (same layout as B; color = class / n_total)
#   D: depth-stratified F1        E: failure mode breakdown
#   F: singlet genome assignment (26-genome breakdown, 4 configs)
# Inputs (DATA = data/processed/marand2021_B73_root): Root1_rep1/eval/phase3_2026-04-09/phase3_summary_metrics.tsv,
#   Root1_rep1/sub1k_B/barcodes_with_bin.tsv, Root1_rep1/sub1k_B/genotyping_runs/<tag>/<cfg>/sub1k_B_cells_calls.tsv.gz
# Sources analysis/_helpers/fig4_helpers.R. Run from the repo root: Rscript analysis/supplementary/figS6.R (sbatch analysis/supplementary/figS6.sh)

# Denominator = total barcodes with total_nuclear_reads >= 500. sub1k_B is shown alone
# because it faithfully reproduces the full Root1 result; sub1k_A was used during
# development but is omitted from the final figure. Companion figure S5 covers the
# synthetic benchmark (Phase 2).

suppressPackageStartupMessages({
  if (requireNamespace("tidyverse", quietly = TRUE)) {
    library(tidyverse)
  } else {
    library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  }
  library(patchwork)
  library(scales)
  library(data.table)
})

source("analysis/_helpers/fig4_helpers.R")

# -------------------------
# 0) CONFIG
# -------------------------
DATA         <- "data/processed/marand2021_B73_root"
SAMPLE       <- "Root1_rep1"
PHASE3_FILE  <- file.path(DATA, SAMPLE, "eval",
                          "phase3_2026-04-09", "phase3_summary_metrics.tsv")
PHASE3_CALLS_BASE <- file.path(DATA, SAMPLE)
FACTORIAL_TAG_P3 <- "factorial_phase3_2026-04-09"

OUTDIR <- "figures/supplementary/figS6"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_READS <- 500
TRACK     <- "sub1k_B"  # faithfully reproduces full Root1; sub1k_A dropped.
GENOME_PANEL_CONFIGS <- c("C0", "C1g_naked", "S08_full_xa0", "S10_full_bic3")
DEPTH_FG_CONFIGS <- c("C0", "C1g_naked", "C3b_mq50", "S08_full_xa0",
                       "S10_full_bic3", "S13_full_eta2")

stopifnot("Phase 3 summary not found" = file.exists(PHASE3_FILE))

P3_CONFIGS <- CONFIG_ORDER[grepl("3", CONFIG_PARAM_TABLE$phases)]

# =============================================================================
# 1) LOAD DATA
# =============================================================================
cat("Loading Phase 3 summary (", TRACK, " only)...\n", sep = "")
p3_raw <- fread(PHASE3_FILE) %>% as_tibble() %>%
  filter(track == TRACK)
stopifnot("Phase 3 summary has no rows for TRACK" = nrow(p3_raw) > 0)

# Pool across depth bins (sum numerators / sum denominators)
p3_pooled <- p3_raw %>%
  group_by(track, method) %>%
  summarise(
    n_total          = sum(n_total),
    n_correct_relaxed = sum(n_correct_relaxed),
    n_wrong_genome   = sum(n_wrong_genome),
    n_called_doublet = sum(n_called_doublet),
    n_called_ambig   = sum(n_called_ambig),
    sens_relaxed     = n_correct_relaxed / n_total,
    prec_relaxed     = if (n_correct_relaxed + n_wrong_genome > 0)
      n_correct_relaxed / (n_correct_relaxed + n_wrong_genome) else NA_real_,
    wrong_genome_rate = n_wrong_genome / n_total,
    doublet_rate     = n_called_doublet / n_total,
    ambig_rate       = n_called_ambig / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    f1_relaxed = if_else(
      sens_relaxed + prec_relaxed > 0,
      2 * sens_relaxed * prec_relaxed / (sens_relaxed + prec_relaxed),
      0
    )
  )

p3_avg <- p3_pooled %>%
  group_by(method) %>%
  summarise(
    across(c(sens_relaxed, prec_relaxed, f1_relaxed,
             wrong_genome_rate, doublet_rate, ambig_rate,
             n_total, n_correct_relaxed, n_called_doublet), mean),
    .groups = "drop"
  )

p3_depth <- p3_raw %>%
  group_by(method, depth_bin) %>%
  summarise(
    across(c(sens_relaxed, prec_relaxed, f1_relaxed), mean),
    .groups = "drop"
  )

# Raw cells_calls for the genome breakdown (Panel F) and the B/C heatmaps (all 19 configs).
# Filter via total_nuclear_reads >= MIN_READS to match eval_phase_factorial.R.
cat("Loading raw cells_calls + barcodes_with_bin for all", length(P3_CONFIGS),
    "configs on", TRACK, "...\n")
calls_list <- list()
for (tk in TRACK) {
  bc_file <- file.path(PHASE3_CALLS_BASE, tk, "barcodes_with_bin.tsv")
  bc <- fread(bc_file, select = c("barcode", "total_nuclear_reads", "depth_bin"))
  for (cfg in P3_CONFIGS) {
    f <- file.path(PHASE3_CALLS_BASE, tk, "genotyping_runs",
                   FACTORIAL_TAG_P3, cfg,
                   paste0(tk, "_cells_calls.tsv.gz"))
    if (!file.exists(f)) { warning("Missing: ", f); next }
    dt <- fread(f, select = c("barcode", "call", "genome_1", "n_reads"))
    # cells_calls barcodes have a "-sub1k_{A,B}" suffix; bc table is bare
    # 16-mer. Strip the suffix before joining.
    dt[, barcode := sub("-.*$", "", barcode)]
    dt <- merge(dt, bc, by = "barcode", all.x = TRUE)
    dt <- dt[!is.na(total_nuclear_reads) & total_nuclear_reads >= MIN_READS]
    dt[, `:=`(track = tk, method = cfg)]
    calls_list[[length(calls_list) + 1]] <- dt
  }
}
calls_raw <- rbindlist(calls_list)
cat("  Loaded", format(nrow(calls_raw), big.mark = ","),
    "post-filter barcodes across", uniqueN(calls_raw$method), "configs\n")

# Subset used by Panel F (4-config genome breakdown, singlets only)
singlets_raw <- calls_raw[call %in% c("single_clean", "dirty_singlet", "weak_doublet") &
                           method %in% GENOME_PANEL_CONFIGS]
cat("  Subset of", nrow(singlets_raw), "singlet rows for Panel F\n")

# B/C joint metric table -- same logic as S5 (denominator = all barcodes
# passing MIN_READS filter; numerator B = in call_class AND g1==B73;
# numerator C = in call_class regardless of g1).
calls_raw[, call_class := fcase(
  call %in% c("single_clean", "dirty_singlet"),  "singlet",
  call %in% c("doublet", "weak_doublet"),         "doublet",
  call == "ambiguous",                             "ambiguous",
  default =                                        "empty"
)]
calls_raw[, is_b73_g1 := !is.na(genome_1) & genome_1 == "B73"]

p3_total <- calls_raw[, .(n_total = .N), by = .(track, method)]
CALL_CLASSES <- c("singlet", "doublet", "ambiguous")
p3_class <- calls_raw[call_class %in% CALL_CLASSES,
                       .(n_in_class     = .N,
                         n_in_class_b73 = sum(is_b73_g1)),
                       by = .(track, method, call_class)]
p3_grid <- CJ(track      = unique(calls_raw$track),
              method     = P3_CONFIGS,
              call_class = CALL_CLASSES, unique = TRUE)
p3_metrics <- merge(p3_grid, p3_class,
                    by = c("track", "method", "call_class"), all.x = TRUE)
p3_metrics[is.na(n_in_class),     n_in_class     := 0L]
p3_metrics[is.na(n_in_class_b73), n_in_class_b73 := 0L]
p3_metrics <- merge(p3_metrics, p3_total,
                    by = c("track", "method"), all.x = TRUE)
p3_metrics[, frac_class     := n_in_class     / n_total]
p3_metrics[, frac_class_b73 := n_in_class_b73 / n_total]

cat("  B/C metrics computed; n_total range across configs =",
    min(p3_total$n_total), "-", max(p3_total$n_total), "\n")

# =============================================================================
# PANEL A: Configuration Parameter Table (Phase 3 configs only)
# =============================================================================
cat("Building Panel A...\n")

param_cols <- c("mapq", "xa", "topk", "wdisc", "friend", "xmap", "eta", "w_amb", "bic")
param_labels <- c("MAPQ", "XA", "Top-K\nreclass", "Winner\ndiscount",
                   "Friend\nrescue", "Cross-map\n(L2)", "Eta\niters",
                   "w_amb", "BIC\nmargin")

p3_param_table <- CONFIG_PARAM_TABLE %>%
  filter(method %in% P3_CONFIGS)

ptable_long <- p3_param_table %>%
  mutate(method = factor(method, levels = rev(P3_CONFIGS)),
         across(all_of(param_cols), as.character)) %>%
  pivot_longer(cols = all_of(param_cols), names_to = "param", values_to = "value") %>%
  mutate(param = factor(param, levels = param_cols, labels = param_labels))

c0_long <- p3_param_table %>%
  filter(method == "C0") %>%
  mutate(across(all_of(param_cols), as.character)) %>%
  pivot_longer(cols = all_of(param_cols), names_to = "param", values_to = "c0_value") %>%
  mutate(param = factor(param, levels = param_cols, labels = param_labels))

ptable_long <- ptable_long %>%
  left_join(c0_long %>% select(param, c0_value), by = "param") %>%
  mutate(is_deviation = as.character(value) != as.character(c0_value))

pA <- ggplot(ptable_long, aes(x = param, y = method)) +
  geom_tile(aes(fill = is_deviation), color = "white", linewidth = 0.4) +
  geom_text(aes(label = value), size = 2.4) +
  scale_fill_manual(values = c("TRUE" = "#fff2cc", "FALSE" = "#f5f5f5"),
                    guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = c(0.5, 0.5))) +
  labs(x = NULL, y = NULL,
       title = "Configuration parameter table (Phase 3 configs)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7, hjust = 1,
                                    color = CONFIG_COLORS[rev(P3_CONFIGS)]),
        axis.text.x = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# =============================================================================
# Config ordering: configs are ordered by pooled relaxed F1 (S10 at top,
# C1g_naked at bottom), matching the manuscript narrative.
# =============================================================================
f1_order <- p3_avg %>%
  filter(method %in% P3_CONFIGS) %>%
  arrange(desc(f1_relaxed)) %>%
  pull(method)

# Plot-ready data: factorize method (top of plot = highest F1) and call_class.
p3m <- p3_metrics %>%
  as_tibble() %>%
  mutate(
    method     = factor(method, levels = rev(f1_order)),
    method_lab = factor(CONFIG_LABELS[as.character(method)],
                         levels = CONFIG_LABELS[rev(f1_order)]),
    call_class = factor(call_class, levels = CALL_CLASSES,
                         labels = c("Singlet",   # single_clean + dirty_singlet
                                    "Doublet",   # doublet + weak_doublet
                                    "Ambiguous"))
  )

# build_heatmap(): rows = config, cols = call class.
build_heatmap <- function(data, fill_col, label_col, title, subtitle,
                          hide_y = FALSE) {
  g <- ggplot(data, aes(x = call_class, y = method_lab,
                         fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = .data[[label_col]]),
              size = 2.2, color = "grey15") +
    HEAT_PALETTE() +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    heat_theme
  if (hide_y) g <- g + theme(axis.text.y = element_blank())
  g
}

# =============================================================================
# PANEL B: B73-dominance per AM call class
#   numerator = barcodes in call_class AND genome_1 == "B73"
#   denominator = barcodes with total_nuclear_reads >= MIN_READS per config
# =============================================================================
cat("Building Panel B (B73 dominance per call class)...\n")
pB <- build_heatmap(p3m, "frac_class_b73", "n_in_class_b73",
                    "B73 dominance per AM call class",
                    "g1 == B73, by call class (singlet | doublet | ambiguous)")

# =============================================================================
# PANEL C: AM call-class share (any g1)
#   numerator = barcodes in call_class (any genome_1)
#   denominator = same as B
# =============================================================================
cat("Building Panel C (call-class share)...\n")
pC <- build_heatmap(p3m, "frac_class", "n_in_class",
                    "Call-class share per config",
                    "fraction routed to each AM call class (any g1)",
                    hide_y = TRUE)

# =============================================================================
# PANEL D: Depth-Stratified F1
# =============================================================================
cat("Building Panel D...\n")

depth_labels <- c(b2 = "250-500", b3 = "500-1K", b4 = "1K-2K",
                  b5 = "2K-5K", b6 = "5K+")
p3_depth_p3 <- p3_depth %>%
  filter(method %in% P3_CONFIGS) %>%
  mutate(depth_label = factor(depth_labels[depth_bin],
                               levels = depth_labels))

bg_data <- p3_depth_p3 %>% filter(!method %in% DEPTH_FG_CONFIGS)
fg_data <- p3_depth_p3 %>% filter(method %in% DEPTH_FG_CONFIGS)

pD <- ggplot() +
  geom_line(data = bg_data,
            aes(x = depth_label, y = f1_relaxed, group = method),
            color = "grey85", linewidth = 0.3, alpha = 0.7) +
  geom_line(data = fg_data,
            aes(x = depth_label, y = f1_relaxed, group = method,
                color = method),
            linewidth = 0.8) +
  geom_point(data = fg_data,
             aes(x = depth_label, y = f1_relaxed, color = method,
                 shape = method),
             size = 2) +
  scale_color_manual(values = CONFIG_COLORS,
                     labels = CONFIG_LABELS[DEPTH_FG_CONFIGS],
                     name = NULL) +
  scale_shape_manual(values = CONFIG_SHAPES,
                     labels = CONFIG_LABELS[DEPTH_FG_CONFIGS],
                     name = NULL) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(x = "Read depth bin", y = "F1 (relaxed)",
       title = "Depth-stratified F1") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6)) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

# =============================================================================
# PANEL E: Failure Modes
# =============================================================================
cat("Building Panel E...\n")

p3_failures <- p3_avg %>%
  filter(method %in% P3_CONFIGS) %>%
  select(method, wrong_genome_rate, doublet_rate, ambig_rate) %>%
  pivot_longer(cols = c(wrong_genome_rate, doublet_rate, ambig_rate),
               names_to = "failure", values_to = "rate") %>%
  mutate(
    method = factor(method, levels = rev(f1_order)),
    failure = factor(failure,
                     levels = c("wrong_genome_rate", "doublet_rate", "ambig_rate"),
                     labels = c("Wrong genome", "Called doublet", "Ambiguous"))
  )

pE <- ggplot(p3_failures, aes(x = method, y = rate, fill = failure)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Wrong genome"  = "#d62728",
                                "Called doublet" = "#ff9896",
                                "Ambiguous"      = "#aec7e8"),
                    name = NULL) +
  scale_y_continuous(labels = percent_format()) +
  coord_flip() +
  labs(x = NULL, y = "Fraction of barcodes",
       title = "Failure mode breakdown") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "bottom")

# =============================================================================
# PANEL F: Singlet Genome Assignment, 26-genome breakdown
# =============================================================================
cat("Building Panel F...\n")

genome_counts <- singlets_raw[, .(N = .N),
                               by = .(method, genome_1)] %>%
  as_tibble() %>%
  group_by(method, genome_1) %>%
  summarise(N = sum(N), .groups = "drop") %>%
  mutate(is_B73 = genome_1 == TRUE_GENOME)

genome_order <- genome_counts %>%
  group_by(genome_1) %>%
  summarise(total = sum(N), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genome_1)

genome_counts <- genome_counts %>%
  mutate(
    genome_1 = factor(genome_1, levels = genome_order),
    method   = factor(method, levels = GENOME_PANEL_CONFIGS),
    method_label = CONFIG_LABELS[as.character(method)]
  )

genome_prec <- genome_counts %>%
  group_by(method, method_label) %>%
  summarise(
    n_B73   = sum(N[is_B73]),
    n_total = sum(N),
    prec    = n_B73 / n_total,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("Precision: %.1f%%  (%s / %s)",
                         100 * prec,
                         format(n_B73, big.mark = ","),
                         format(n_total, big.mark = ",")))

pF <- ggplot(genome_counts, aes(x = genome_1, y = N, fill = is_B73)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = genome_highlight_colors, guide = "none") +
  scale_y_continuous(labels = comma_format()) +
  facet_wrap(~ method_label, ncol = 2, scales = "free_y") +
  geom_text(data = genome_prec,
            aes(x = length(genome_order) / 2, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1.3,
            size = 2.8, color = "grey30") +
  labs(x = NULL, y = "Singlet count",
       title = "Genome assignment for singlets (B73-only library, 26 NAM references)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 8))

# =============================================================================
# ASSEMBLY & EXPORT
# =============================================================================
cat("Assembling figure...\n")

# Layout:
#   A: param table (top, full width)
#   B + C: heatmaps side-by-side, full width (B = B73 dominance; C = call share)
#   D + E: depth-F1 + failure modes (lower-middle row, side-by-side)
#   F: genome breakdown (bottom, full width)
design <- "
AAAAAAAAAA
BBBBBCCCCC
DDDDDEEEEE
FFFFFFFFFF
"

fig <- pA + pB + pC + pD + pE + pF +
  plot_layout(design = design, heights = c(2.5, 6, 4, 5)) +
  plot_annotation(
    tag_levels = "A",
    title = "Supplementary Figure S6: Genotyping parameter validation, Root1 sub1k_B",
    subtitle = sprintf("Phase 3: 26-genome reference (B73-only library), 19 configs on sub1k_B (faithfully reproduces full Root1). Denominator = barcodes with total_nuclear_reads ≥ %d.", MIN_READS),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10))
  )

ggsave(file.path(OUTDIR, "Fig_S6_root1.pdf"), fig,
       width = 16, height = 22, limitsize = FALSE)
ggsave(file.path(OUTDIR, "Fig_S6_root1.png"), fig,
       width = 16, height = 22, dpi = 300, limitsize = FALSE)

ggsave(file.path(OUTDIR, "S6A_param_table.pdf"),       pA, width = 10, height = 5)
ggsave(file.path(OUTDIR, "S6B_b73_dominance.pdf"),     pB, width =  5, height = 7)
ggsave(file.path(OUTDIR, "S6C_call_class.pdf"),        pC, width =  4, height = 7)
ggsave(file.path(OUTDIR, "S6D_depth_f1.pdf"),          pD, width =  7, height = 5)
ggsave(file.path(OUTDIR, "S6E_failure_modes.pdf"),     pE, width =  6, height = 5)
ggsave(file.path(OUTDIR, "S6F_genome_breakdown.pdf"),  pF, width = 14, height = 6)

cat("Done. Outputs in:", OUTDIR, "\n")
