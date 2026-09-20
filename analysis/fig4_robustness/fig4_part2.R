#!/usr/bin/env Rscript
# Fig 4, panels D to G (maize root library Root1, one genotype mapped against 26 NAM genomes).
#   D: cross-mapping profile phi for the B73 anchor (reported descriptively, see section 2)
#   E: call distribution, C0 baseline vs S10_full_bic3
#   F: dominant genome (genome_1) by call type, C0 vs S10_full_bic3
#   G: wrong-genome rate across the five Phase 4 configurations
# Inputs (DATA = data/processed/marand2021_B73_root): Root1_rep1/genotyping_runs/_archive/eval_xmap_v2_phi_B73.tsv,
#   Root1_rep1/genotyping_runs/<PHASE4_TAG>/{C0,S10_full_bic3}/Root1_rep1_cells_calls.tsv.gz,
#   Root1_rep1/eval/phase4_2026-04-28/phase4_summary_metrics.tsv (from workflows/03_genotyping/eval_phase_factorial.R phase4)
# Sources analysis/_helpers/fig4_helpers.R. Run from the repo root: Rscript analysis/fig4_robustness/fig4_part2.R

# Pure B73 scATAC library mapped against 26 NAM genomes (~95% sequence
# identity). Phase 4 of the genotyping factorial: C0 baseline vs Phase-3
# winner S10_full_bic3 (mq50 + friend OFF + xmap + w_amb=0.5 + bic_margin=3).
# phi is reference-panel-derived and not config-dependent.

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
DATA       <- "data/processed/marand2021_B73_root"
SAMPLE     <- "Root1_rep1"
PHASE4_TAG <- "factorial_phase4_2026-04-28"
GRID_DIR   <- file.path(DATA, SAMPLE, "genotyping_runs")
PHASE4_DIR <- file.path(GRID_DIR, PHASE4_TAG)
EVAL_DIR   <- file.path(GRID_DIR, "_archive")   # phi summary table lives here
OUTDIR     <- "figures/main/fig4"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_READS  <- 250

# Two-config subset for the main figure: baseline vs Phase 4 winner
RUNS         <- c("C0", "S10_full_bic3")
RUN_LABELS   <- c(C0            = "C0 baseline\n(mq20, default)",
                  S10_full_bic3 = "S10 winner\n(mq50 + xmap + w_amb=0.5)")

# =============================================================================
# 1) LOAD DATA
# =============================================================================

# ---- Panel D: phi profile (reference-panel-only, not config-dependent) ----
phi_B73 <- fread(file.path(EVAL_DIR, "eval_xmap_v2_phi_B73.tsv")) %>%
  as_tibble() %>%
  mutate(
    is_true       = target_genome == TRUE_GENOME,
    target_genome = fct_reorder(target_genome, phi_B73, .desc = TRUE)
  )

# ---- Panels E/F/G: cells_calls per config ----
load_phase4_calls <- function(run_name) {
  f <- file.path(PHASE4_DIR, run_name, paste0(SAMPLE, "_cells_calls.tsv.gz"))
  if (!file.exists(f)) stop("Missing: ", f)
  dt <- fread(f, sep = "\t",
              select = c("barcode", "call", "genome_1", "n_reads"))
  dt[, run := run_name]
  dt
}

cat("Loading Root1 Phase 4 cells_calls for", length(RUNS), "configs...\n")
calls_both <- rbindlist(lapply(RUNS, load_phase4_calls))
cat("  Total BCs:", format(nrow(calls_both), big.mark = ","), "\n")

# Drop low_reads and require >= MIN_READS (empty calls are kept regardless of depth)
calls_filt <- calls_both[call != "low_reads" &
                         (call == "empty" | n_reads >= MIN_READS)]
cat("  Filtered (>=", MIN_READS, "reads):",
    format(nrow(calls_filt), big.mark = ","), "BCs\n")

# =============================================================================
# 2) PANEL D: phi profile
# phi is reported descriptively, as a measure of reference-panel redundancy
# (how much B73-anchored read mass maps to each other genome). It is a property
# of the reference panel, not of a configuration: the C0 baseline shown in E
# and F applies no cross-mapping correction, S10_full_bic3 does (xmap on).
# =============================================================================
b73_phi_val <- phi_B73$phi_B73[phi_B73$target_genome == "B73"]

pD <- ggplot(phi_B73, aes(x = target_genome, y = phi_B73, fill = is_true)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = genome_highlight_colors, guide = "none") +
  annotate("text",
           x = which(levels(phi_B73$target_genome) == "B73"),
           y = b73_phi_val + 0.008,
           label = sprintf("B73: %.1f%%", 100 * b73_phi_val),
           hjust = 0.5, size = 3.0, color = "#2ca02c", fontface = "bold") +
  labs(x = NULL,
       y = expression(phi["g|B73"]),
       title = expression("Cross-mapping profile " * phi *
                          " for B73 anchor (26 NAM genomes)")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))

# =============================================================================
# 3) PANEL E: call distribution (C0 vs winner)
# =============================================================================
# Use the n_reads-filtered set for the call distribution. Without this, the
# "low_reads" category (>50% of barcodes for both configs) dominates the panel
# and the proportions of called classes are invisible.
calls_dist <- calls_filt[, .N, by = .(run, call)] %>%
  as_tibble() %>%
  complete(run, call = call_levels, fill = list(N = 0)) %>%
  filter(call %in% call_levels) %>%
  mutate(
    display = factor(RUN_LABELS[run], levels = RUN_LABELS[RUNS]),
    call    = factor(call, levels = call_levels)
  )

totals <- calls_dist %>%
  group_by(display) %>%
  summarise(total = sum(N), .groups = "drop")

pE <- ggplot(calls_dist, aes(x = display, y = N, fill = call)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(data = totals,
            aes(x = display, y = 1.03,
                label = paste0("n=", format(total, big.mark = ","))),
            inherit.aes = FALSE, size = 2.6, fontface = "italic") +
  scale_fill_manual(values = call_colors, name = "Call type") +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Fraction of barcodes",
       title = "Call distribution") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold")) +
  guides(fill = guide_legend(nrow = 1))

# =============================================================================
# 4) PANEL F: dominant genome by call type
# =============================================================================
genome1_by_call <- calls_filt[call != "empty" & genome_1 != "None",
                              .(N = .N), by = .(run, call, genome_1)]
genome1_binary <- genome1_by_call[, .(N = sum(N)),
                                  by = .(run, call,
                                         is_B73 = genome_1 == TRUE_GENOME)] %>%
  as_tibble() %>%
  mutate(
    display = factor(RUN_LABELS[run], levels = RUN_LABELS[RUNS]),
    call    = factor(call, levels = call_levels)
  )

genome1_fracs <- genome1_binary %>%
  group_by(run, display, call) %>%
  summarise(
    n_total  = sum(N),
    n_B73    = sum(N[is_B73]),
    frac_B73 = n_B73 / n_total,
    .groups  = "drop"
  )

pF <- ggplot(genome1_binary, aes(x = call, y = N, fill = is_B73)) +
  geom_col(position = "fill", width = 0.75) +
  geom_text(data = genome1_fracs,
            aes(x = call, y = frac_B73 / 2,
                label = sprintf("%.1f%%", 100 * frac_B73)),
            inherit.aes = FALSE, size = 2.3, color = "white",
            fontface = "bold") +
  geom_text(data = genome1_fracs,
            aes(x = call, y = 1.04,
                label = paste0("n=", format(n_total, big.mark = ","))),
            inherit.aes = FALSE, size = 2.0, fontface = "italic") +
  scale_fill_manual(values = genome_highlight_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.08))) +
  facet_wrap(~ display, ncol = 2) +
  labs(x = NULL, y = "Fraction of barcodes",
       title = expression("Dominant genome (genome"[1]*") by call type")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# =============================================================================
# 5) PANEL G: wrong-genome rate across the five Phase 4 configurations
# (baseline, two single-knob changes, two stacked configurations)
# =============================================================================
phase4_metrics <- fread(file.path(DATA, SAMPLE, "eval", "phase4_2026-04-28",
                                  "phase4_summary_metrics.tsv"))

# Order configs by ascending wrong_genome_rate
metric_order <- phase4_metrics[order(wrong_genome_rate), method]
cfg_short <- c(
  C0                       = "C0",
  C3b_mq50                 = "C3b\nmq50",
  C4a_bic3                 = "C4a\nbic3",
  S04_nofr_mq50_xmap_xaUL  = "S04\nstacked",
  S10_full_bic3            = "S10\nstacked\n(winner)"
)
cfg_group <- c(
  C0                       = "baseline",
  C3b_mq50                 = "single-knob",
  C4a_bic3                 = "single-knob",
  S04_nofr_mq50_xmap_xaUL  = "stacked",
  S10_full_bic3            = "stacked"
)
group_colors <- c(baseline = "#4d4d4d",
                  "single-knob" = "#3182bd",
                  stacked = "#238b45")

phase4_metrics[, label := factor(cfg_short[method],
                                 levels = cfg_short[metric_order])]
phase4_metrics[, group := cfg_group[method]]
phase4_metrics[, fp_pct := 100 * wrong_genome_rate]

pG <- ggplot(phase4_metrics, aes(x = label, y = fp_pct, fill = group)) +
  geom_col(width = 0.7, color = "grey20") +
  geom_text(aes(label = sprintf("%.2f%%", fp_pct)),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%", x),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Wrong-genome rate",
       title = "Wrong-genome rate (Phase 4)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

# =============================================================================
# 6) ASSEMBLY & EXPORT
# =============================================================================
fig <- (pD | pE | pF | pG) +
  plot_layout(widths = c(1.2, 0.55, 0.85, 0.9)) +
  plot_annotation(tag_levels = list(c("D", "E", "F", "G"))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(OUTDIR, "Fig4_part2_DtoG.pdf"), fig,
       width = 16, height = 5.5)
ggsave(file.path(OUTDIR, "Fig4_part2_DtoG.png"), fig,
       width = 16, height = 5.5, dpi = 300)

# Individual panels
ggsave(file.path(OUTDIR, "Fig4_D_phi.pdf"),             pD, width = 7,    height = 3.5)
ggsave(file.path(OUTDIR, "Fig4_E_call_dist.pdf"),       pE, width = 4.5,  height = 4.5)
ggsave(file.path(OUTDIR, "Fig4_F_genome1_by_call.pdf"), pF, width = 5.5,  height = 4)
ggsave(file.path(OUTDIR, "Fig4_G_wrong_genome.pdf"),    pG, width = 6,    height = 4)

cat("\nDone. Outputs in:", OUTDIR, "\n")
