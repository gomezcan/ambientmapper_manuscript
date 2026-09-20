#!/usr/bin/env Rscript
# Supplementary Fig S5: genotyping parameter validation on the synthetic benchmark (Phase 2;
# Track B = all ortholog peaks, Track B-disc = discriminative peaks; 3-genome reference).
#   A: configuration parameter table (Phase 2 configs)
#   B: B73-dominance heatmaps per call class (singlet / doublet / ambiguous; genome_1 == B73)
#   C: call-class share heatmaps (any genome_1); C >= B cell-wise, C - B = wrong-genome-call rate
#   D: read-level diagnostics on true B73 singlets (read fraction per genome, NM, AS, MAPQ)
# Inputs (DATA = data/processed/synthetic): synthetic/eval/phase2_2026-04-09/phase2_summary_metrics.tsv,
#   <track>/barcoded/<ds>/truth_table.tsv, <track>/<ds>/genotyping_runs/<tag>/<cfg>/<ds>_cells_calls.tsv.gz,
#   <track>/<ds>/cell_map_ref_chunks/*filtered.tsv.gz. Sources analysis/_helpers/fig4_helpers.R.
# Run from the repo root: Rscript analysis/supplementary/figS5.R (about 24 GB RAM; sbatch analysis/supplementary/figS5.sh)

# All panels share denominator = true B73 singlets with n_reads >= 500.
# Companion figure S6 covers the Root1 sub1k panels (Phase 3).

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
DATA         <- "data/processed/synthetic"
PHASE2_FILE  <- file.path(DATA, "synthetic", "eval",
                          "phase2_2026-04-09", "phase2_summary_metrics.tsv")
FACTORIAL_TAG_P2 <- "factorial_phase2_2026-04-09"

# Track B and Track B-disc base dirs
P2_TRACK_DIRS <- c("Track B"      = file.path(DATA, "synthetic"),
                   "Track B-disc" = file.path(DATA, "synthetic_disc"))
OUTDIR <- "figures/supplementary/figS5"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

stopifnot("Phase 2 summary not found" = file.exists(PHASE2_FILE))

# Configs that appear in Phase 2
P2_CONFIGS <- CONFIG_ORDER[grepl("2", CONFIG_PARAM_TABLE$phases)]

# =============================================================================
# 1) LOAD DATA
# =============================================================================
cat("Loading Phase 2 summary...\n")
p2_raw <- fread(PHASE2_FILE) %>% as_tibble()

cat("Loading Phase 2 raw cells_calls + truth tables...\n")
P2_DATASETS <- c(
  "alpha_000",
  paste0("alpha_", sprintf("%03d", c(2, 5, 10, 20, 30, 40, 50)), "_Il14H"),
  paste0("alpha_", sprintf("%03d", c(2, 5, 10, 20, 30, 40, 50)), "_Ki11")
)

p2_calls_list <- list()
for (track_name in names(P2_TRACK_DIRS)) {
  base_dir <- P2_TRACK_DIRS[[track_name]]
  for (ds in P2_DATASETS) {
    truth_f <- file.path(base_dir, "barcoded", ds, "truth_table.tsv")
    if (!file.exists(truth_f)) { warning("Missing truth: ", truth_f); next }
    truth <- fread(truth_f, select = c("barcode", "type"))

    for (cfg in P2_CONFIGS) {
      calls_f <- file.path(base_dir, ds, "genotyping_runs",
                           FACTORIAL_TAG_P2, cfg,
                           paste0(ds, "_cells_calls.tsv.gz"))
      if (!file.exists(calls_f)) next
      dt <- fread(calls_f, select = c("barcode", "call", "genome_1", "n_reads"))
      dt[, `:=`(track = track_name, method = cfg, dataset = ds)]
      dt <- merge(dt, truth, by = "barcode", all.x = TRUE)
      p2_calls_list[[length(p2_calls_list) + 1]] <- dt
    }
  }
}
p2_calls_raw <- rbindlist(p2_calls_list, fill = TRUE)
cat("  Loaded", format(nrow(p2_calls_raw), big.mark = ","), "Phase 2 rows\n")

# =============================================================================
# B/C joint metric table.
# Denominator (everywhere): true B73 singlets passing n_reads >= MIN_READS.
# For each (track, config, dataset, call_class in {singlet, doublet, ambiguous}):
#   - n_total       : true B73 singlets passing filter
#   - n_in_class    : same set, classified into call_class (numerator for C)
#   - n_in_class_b73: same as above AND genome_1 == "B73"  (numerator for B)
#   - frac_class      = n_in_class      / n_total   ("call-class share" -> pC)
#   - frac_class_b73  = n_in_class_b73  / n_total   ("B73-dominant share" -> pB)
# call_class mapping:
#   single_clean | dirty_singlet      -> singlet
#   doublet     | weak_doublet         -> doublet
#   ambiguous                          -> ambiguous
#   (empty is residual; not shown)
# Matches eval_phase_factorial.R: MIN_READS = 500.
MIN_READS <- 500

p2_singlets <- p2_calls_raw[type == "singlet" & !is.na(call) & n_reads >= MIN_READS]
p2_singlets[, call_class := fcase(
  call %in% c("single_clean", "dirty_singlet"),  "singlet",
  call %in% c("doublet", "weak_doublet"),         "doublet",
  call == "ambiguous",                             "ambiguous",
  default =                                        "empty"
)]
p2_singlets[, is_b73_g1 := !is.na(genome_1) & genome_1 == "B73"]

# Per-cell denominator (across all classes)
p2_total <- p2_singlets[, .(n_total = .N), by = .(track, method, dataset)]

CALL_CLASSES <- c("singlet", "doublet", "ambiguous")
p2_class <- p2_singlets[call_class %in% CALL_CLASSES,
                         .(n_in_class     = .N,
                           n_in_class_b73 = sum(is_b73_g1)),
                         by = .(track, method, dataset, call_class)]
# Fill in zero-cells (combinations that have no rows)
p2_grid <- CJ(track       = unique(p2_singlets$track),
              method      = unique(p2_singlets$method),
              dataset     = unique(p2_singlets$dataset),
              call_class  = CALL_CLASSES, unique = TRUE)
p2_metrics <- merge(p2_grid, p2_class,
                    by = c("track", "method", "dataset", "call_class"),
                    all.x = TRUE)
p2_metrics[is.na(n_in_class),     n_in_class     := 0L]
p2_metrics[is.na(n_in_class_b73), n_in_class_b73 := 0L]
p2_metrics <- merge(p2_metrics, p2_total,
                    by = c("track", "method", "dataset"), all.x = TRUE)
p2_metrics[, frac_class     := n_in_class     / n_total]
p2_metrics[, frac_class_b73 := n_in_class_b73 / n_total]

p2_metrics[, contam_series := fcase(
  dataset == "alpha_000", "none",
  grepl("Il14H", dataset), "Il14H",
  grepl("Ki11", dataset),  "Ki11"
)]
p2_metrics[, true_alpha := as.numeric(
  sub("_.*", "", sub("alpha_", "", dataset))
) / 100]
p2_metrics[dataset == "alpha_000", true_alpha := 0]

cat("  Phase 2 metrics computed for", uniqueN(p2_metrics$method),
    "configs,", MIN_READS, "min reads, denom range",
    range(p2_total$n_total)[1], "-", range(p2_total$n_total)[2],
    "true singlets per dataset\n")

# =============================================================================
# PANEL A: Configuration Parameter Table (Phase 2 configs only)
# =============================================================================
cat("Building Panel A...\n")

param_cols <- c("mapq", "xa", "topk", "wdisc", "friend", "xmap", "eta", "w_amb", "bic")
param_labels <- c("MAPQ", "XA", "Top-K\nreclass", "Winner\ndiscount",
                   "Friend\nrescue", "Cross-map\n(L2)", "Eta\niters",
                   "w_amb", "BIC\nmargin")

p2_param_table <- CONFIG_PARAM_TABLE %>%
  filter(method %in% P2_CONFIGS)

ptable_long <- p2_param_table %>%
  mutate(method = factor(method, levels = rev(P2_CONFIGS)),
         across(all_of(param_cols), as.character)) %>%
  pivot_longer(cols = all_of(param_cols), names_to = "param", values_to = "value") %>%
  mutate(param = factor(param, levels = param_cols, labels = param_labels))

c0_long <- p2_param_table %>%
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
       title = "Configuration parameter table (Phase 2 configs)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7, hjust = 1,
                                    color = CONFIG_COLORS[rev(P2_CONFIGS)]),
        axis.text.x = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# =============================================================================
# Duplicate alpha_000 across both contaminant facets so the 0% column shows
# under each contaminant heading (same data, plotted in both Il14H + Ki11).
# =============================================================================
p2m_a0   <- p2_metrics[contam_series == "none"] %>% as_tibble() %>%
            select(-contam_series) %>%
            crossing(contam_series = c("Il14H", "Ki11"))
p2m_rest <- p2_metrics[contam_series != "none"] %>% as_tibble()
p2m <- bind_rows(p2m_a0, p2m_rest) %>%
  filter(method %in% P2_CONFIGS) %>%
  mutate(method      = factor(method, levels = P2_CONFIGS),
         method_lab  = factor(CONFIG_LABELS[as.character(method)],
                               levels = CONFIG_LABELS[P2_CONFIGS]),
         call_class  = factor(call_class, levels = CALL_CLASSES),
         alpha_label = factor(sprintf("%g%%", true_alpha * 100),
                               levels = sprintf("%g%%",
                                                 sort(unique(true_alpha)) * 100)))

# -----------------------------------------------------------------------------
# HEAT_PALETTE() and heat_theme live in analysis/_helpers/fig4_helpers.R (shared with S6).
# -----------------------------------------------------------------------------

# Shared geometry for B & C sub-panels: rows = config (method_lab),
# cols = alpha, facets = track x contaminant.
#   fill_col  = fraction (0-1) used for the color
#   label_col = integer count to display in the cell (barcodes in that class)
build_heatmap <- function(data, fill_col, label_col, title, subtitle,
                          hide_y = FALSE) {
  g <- ggplot(data, aes(x = alpha_label, y = method_lab,
                         fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = .data[[label_col]]),
              size = 2, color = "grey15") +
    facet_grid(track ~ contam_series) +
    HEAT_PALETTE() +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    heat_theme
  if (hide_y) g <- g + theme(axis.text.y = element_blank())
  g
}

class_subsets <- list(
  singlet   = p2m %>% filter(call_class == "singlet"),
  doublet   = p2m %>% filter(call_class == "doublet"),
  ambiguous = p2m %>% filter(call_class == "ambiguous")
)

# =============================================================================
# PANEL B: B73-dominance per AM call class
#   numerator = barcodes in call_class AND genome_1 == "B73"
#   denominator = true B73 singlets with n_reads >= MIN_READS
# =============================================================================
cat("Building Panel B (B73-dominance per class)...\n")
pB1 <- build_heatmap(class_subsets$singlet,   "frac_class_b73", "n_in_class_b73",
                     "B73 dominant + singlet call",
                     "single_clean / dirty_singlet AND g1 == B73")
pB2 <- build_heatmap(class_subsets$doublet,   "frac_class_b73", "n_in_class_b73",
                     "B73 dominant + doublet call",
                     "doublet / weak_doublet AND g1 == B73", hide_y = TRUE)
pB3 <- build_heatmap(class_subsets$ambiguous, "frac_class_b73", "n_in_class_b73",
                     "B73 dominant + ambiguous call",
                     "ambiguous AND g1 == B73", hide_y = TRUE)

pB <- pB1 + pB2 + pB3 +
  plot_layout(widths = c(1.3, 1, 1), guides = "collect") &
  theme(legend.position = "bottom")

# =============================================================================
# PANEL C: Call-class share (regardless of dominant genome)
#   numerator = barcodes in call_class
#   denominator = true B73 singlets with n_reads >= MIN_READS
# =============================================================================
cat("Building Panel C (call-class share)...\n")
pC1 <- build_heatmap(class_subsets$singlet,   "frac_class", "n_in_class",
                     "Called singlet (any g1)",
                     "single_clean / dirty_singlet")
pC2 <- build_heatmap(class_subsets$doublet,   "frac_class", "n_in_class",
                     "Called doublet (any g1)",
                     "doublet / weak_doublet", hide_y = TRUE)
pC3 <- build_heatmap(class_subsets$ambiguous, "frac_class", "n_in_class",
                     "Called ambiguous",
                     "ambiguous", hide_y = TRUE)

pC <- pC1 + pC2 + pC3 +
  plot_layout(widths = c(1.3, 1, 1), guides = "collect") &
  theme(legend.position = "bottom")

# =============================================================================
# PANEL D: Read-level diagnostics (cross-mapping & quality across alphas).
# =============================================================================
cat("Building Panel D (read-level diagnostics)...\n")
cat("  Loading chunk files for cross-mapping & quality analysis...\n")
genome_order_phi <- c("B73", "Il14H", "Ki11")

CHUNK_ALPHAS <- c("alpha_000",
  "alpha_002_Il14H", "alpha_005_Il14H", "alpha_010_Il14H",
  "alpha_020_Il14H", "alpha_050_Il14H",
  "alpha_002_Ki11", "alpha_005_Ki11", "alpha_010_Ki11",
  "alpha_020_Ki11", "alpha_050_Ki11")

chunk_list <- list()
for (track_name in names(P2_TRACK_DIRS)) {
  base_dir <- P2_TRACK_DIRS[[track_name]]
  for (ds in CHUNK_ALPHAS) {
    truth_f <- file.path(base_dir, "barcoded", ds, "truth_table.tsv")
    if (!file.exists(truth_f)) next
    singlet_bcs <- fread(truth_f, select = c("barcode", "type"))[
      type == "singlet", barcode]

    chunks_dir <- file.path(base_dir, ds, "cell_map_ref_chunks")
    chunk_files <- list.files(chunks_dir, pattern = "filtered[.]tsv[.]gz$",
                               full.names = TRUE)
    if (length(chunk_files) == 0) next
    dt <- rbindlist(lapply(chunk_files, function(f)
      fread(f, select = c("Read", "BC", "Genome", "AS", "MAPQ", "NM", "XAcount"))))
    dt <- dt[BC %in% singlet_bcs]
    dt[, `:=`(track = track_name, dataset = ds)]
    chunk_list[[length(chunk_list) + 1]] <- dt
  }
}
chunk_reads <- rbindlist(chunk_list)

chunk_reads[, contam_series := fcase(
  dataset == "alpha_000", "none",
  grepl("Il14H", dataset), "Il14H",
  grepl("Ki11", dataset),  "Ki11"
)]
chunk_reads[, true_alpha := as.numeric(
  sub("_.*", "", sub("alpha_", "", dataset))) / 100]
chunk_reads[dataset == "alpha_000", true_alpha := 0]
chunk_reads[, alpha_label := sprintf("%g%%", true_alpha * 100)]

cat("  Loaded", format(nrow(chunk_reads), big.mark = ","),
    "reads from true B73 singlets across", length(CHUNK_ALPHAS), "datasets\n")

chunk_a0 <- chunk_reads[contam_series == "none"]
chunk_a0_il <- copy(chunk_a0)[, contam_series := "Il14H"]
chunk_a0_ki <- copy(chunk_a0)[, contam_series := "Ki11"]
chunk_rest <- chunk_reads[contam_series != "none"]
chunk_all <- rbindlist(list(chunk_a0_il, chunk_a0_ki, chunk_rest))

alpha_order <- sprintf("%g%%", c(0, 2, 5, 10, 20, 50))
chunk_all[, alpha_label := factor(alpha_label, levels = alpha_order)]
chunk_all[, Genome := factor(Genome, levels = genome_order_phi)]

phi_alpha <- chunk_all[, .N, by = .(track, contam_series, alpha_label, Genome)]
phi_alpha[, frac := N / sum(N), by = .(track, contam_series, alpha_label)]

# Shared theme for the D row (slightly different from B/C: y-axis is genome,
# not config).
d_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid       = element_blank(),
        strip.text       = element_text(face = "bold", size = 7),
        plot.title       = element_text(size = 9, face = "bold"),
        legend.position  = "bottom",
        legend.key.width = unit(0.9, "cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 6))

# pD1: read fraction per genome -- same Blues palette as B/C (single sequential
# 0-100%, no midpoint trick).
pD1 <- ggplot(phi_alpha, aes(x = alpha_label, y = Genome, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f", frac * 100)), size = 2.3) +
  facet_grid(track ~ contam_series) +
  scale_fill_distiller(palette = "Blues", direction = 1,
                       limits = c(0, 1), labels = percent_format(),
                       name = "Read fraction", oob = scales::squish) +
  scale_y_discrete(limits = rev(genome_order_phi)) +
  labs(x = "Contamination level", y = NULL,
       title = "Read fraction per genome (true B73 singlets)") +
  d_theme

# pD2-pD4: quality metrics. Unified viridis palette across all three so the
# diagnostic strip reads as a coherent block; each panel keeps its own range
# (intrinsic to that metric) shown in its legend.
nm_alpha <- chunk_all[, .(mean_NM = mean(NM), pct_NM0 = mean(NM == 0) * 100),
                       by = .(track, contam_series, alpha_label, Genome)]

pD2 <- ggplot(nm_alpha, aes(x = alpha_label, y = Genome, fill = mean_NM)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", mean_NM)), size = 2.3, color = "white") +
  facet_grid(track ~ contam_series) +
  scale_fill_viridis_c(option = "C", limits = c(0, 2), oob = scales::squish,
                        name = "Mean NM") +
  scale_y_discrete(limits = rev(genome_order_phi)) +
  labs(x = "Contamination level", y = NULL,
       title = "Mismatches (NM) per read") +
  d_theme

as_alpha <- chunk_all[, .(mean_AS = mean(AS)),
                       by = .(track, contam_series, alpha_label, Genome)]

pD3 <- ggplot(as_alpha, aes(x = alpha_label, y = Genome, fill = mean_AS)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", mean_AS)), size = 2.3, color = "white") +
  facet_grid(track ~ contam_series) +
  scale_fill_viridis_c(option = "C", limits = c(65, 76), oob = scales::squish,
                        name = "Mean AS") +
  scale_y_discrete(limits = rev(genome_order_phi)) +
  labs(x = "Contamination level", y = NULL,
       title = "Alignment score (AS)") +
  d_theme

mapq_alpha <- chunk_all[, .(mean_MAPQ = mean(MAPQ), pct_XA0 = mean(XAcount == 0) * 100),
                          by = .(track, contam_series, alpha_label, Genome)]

pD4 <- ggplot(mapq_alpha, aes(x = alpha_label, y = Genome, fill = mean_MAPQ)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", mean_MAPQ)), size = 2.3, color = "white") +
  facet_grid(track ~ contam_series) +
  scale_fill_viridis_c(option = "C", limits = c(57, 60), oob = scales::squish,
                        name = "Mean MAPQ") +
  scale_y_discrete(limits = rev(genome_order_phi)) +
  labs(x = "Contamination level", y = NULL,
       title = "Mapping quality (MAPQ)") +
  d_theme

pD <- pD1 + pD2 + pD3 + pD4 +
  plot_layout(nrow = 1, widths = c(1, 1, 1, 1)) +
  plot_annotation(
    title = "Read-level diagnostics on true B73 singlets",
    theme = theme(plot.title = element_text(size = 10, face = "bold"))
  )

# =============================================================================
# ASSEMBLY & EXPORT
# =============================================================================
cat("Assembling figure...\n")

# A: param table (top)
# B: 3-panel B73-dominance heatmap (singlet, doublet, ambiguous)
# C: 3-panel call-class share heatmap (singlet, doublet, ambiguous)
# D: 4-panel read-level diagnostics (phi, NM, AS, MAPQ)
fig <- pA / pB / pC / pD +
  plot_layout(heights = c(2.5, 5, 5, 4)) +
  plot_annotation(
    tag_levels = "A",
    title = "Supplementary Figure S5: Genotyping parameter validation, synthetic benchmark",
    subtitle = sprintf("Phase 2: 3-genome synthetic, Track B (all peaks) + Track B-disc (discriminative peaks). Denominator = true B73 singlets with n_reads ≥ %d.",
                        MIN_READS),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10))
  )

ggsave(file.path(OUTDIR, "Fig_S5_synthetic.pdf"), fig,
       width = 16, height = 22, limitsize = FALSE)
ggsave(file.path(OUTDIR, "Fig_S5_synthetic.png"), fig,
       width = 16, height = 22, dpi = 300, limitsize = FALSE)

ggsave(file.path(OUTDIR, "S5A_param_table.pdf"), pA,         width = 10, height = 5)
ggsave(file.path(OUTDIR, "S5B_b73_dominance.pdf"), pB,        width = 16, height = 6)
ggsave(file.path(OUTDIR, "S5C_call_class.pdf"), pC,           width = 16, height = 6)
ggsave(file.path(OUTDIR, "S5D_read_diagnostics.pdf"), pD,     width = 16, height = 5)

cat("Done. Outputs in:", OUTDIR, "\n")
