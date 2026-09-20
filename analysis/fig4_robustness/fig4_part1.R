#!/usr/bin/env Rscript
# Fig 4, panels A to C (synthetic benchmark, Track B = all ortholog peaks, Track B-disc = discriminative peaks).
#   A: per-read edit distance (NM) at alpha = 0, Track B vs Track B-disc
#   B: singlet and doublet call counts vs true alpha (pre-decontamination, C0 configuration)
#   C: barnyard of B73 vs non-B73 winner reads per barcode, pre vs post decontamination
# Inputs (DATA = data/processed/synthetic): <track>/alpha_000/cell_map_ref_chunks/*_filtered.tsv.gz,
#   <track>/<ds>/genotyping_runs/<SYN_C0_TAG>/C0/<ds>_cells_calls.tsv.gz,
#   <track>/<ds>/<SYN_DECON_DIR>/<ds>_{pre,post}_barcode_genome_counts.tsv.gz (15 datasets per track)
# Sources analysis/_helpers/fig4_helpers.R. Run from the repo root:
#   Rscript analysis/fig4_robustness/fig4_part1.R   (or sbatch analysis/fig4_robustness/fig4_part1.sh)

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
DATA           <- "data/processed/synthetic"
OUTDIR         <- "figures/main/fig4"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SYN_C0_TAG     <- "factorial_phase2_2026-04-09"          # genotyping run tag (C0 baseline)
SYN_DECON_DIR  <- "decontam_without_design_alpha05_C0"   # design-free decontamination run

# alpha label helper: 0 -> "0", 0.02 -> "0.02", 0.1 -> "0.1"
alpha_lab <- function(x) {
  ifelse(x == 0, "0",
         ifelse(x < 0.1, sprintf("%.2f", x), sprintf("%.1f", x)))
}
ALPHA_LEVELS_CHR <- alpha_lab(SYN_ALPHA_LEVELS)

# =============================================================================
# 1) LOAD SYNTHETIC DATA (Panels A, B, C)
# =============================================================================

# ---- Panel A: per-read NM from alpha_000 chunks ----
load_alpha0_nm <- function(track_dir) {
  chunk_dir <- file.path(DATA, track_dir, "alpha_000", "cell_map_ref_chunks")
  fs <- list.files(chunk_dir, pattern = "_filtered\\.tsv\\.gz$", full.names = TRUE)
  if (length(fs) == 0) {
    warning("No filtered chunks for ", track_dir, "/alpha_000")
    return(data.table())
  }
  rbindlist(lapply(fs, function(f) fread(f, select = "NM")))
}

cat("Loading per-read NM for Panel A...\n")
nm_b    <- load_alpha0_nm("synthetic")[, track := SYN_TRACKS["synthetic"]]
nm_disc <- load_alpha0_nm("synthetic_disc")[, track := SYN_TRACKS["synthetic_disc"]]
nm_all  <- rbind(nm_b, nm_disc)
nm_all[, track := factor(track, levels = unname(SYN_TRACKS))]
cat("  ", format(nrow(nm_all), big.mark = ","), "reads loaded\n")
cat("  median NM by track:\n")
print(nm_all[, .(median_NM = median(NM), mean_NM = round(mean(NM), 3), n = .N), by = track])

# ---- Panel C: pre/post barcode_genome_counts (60 files) ----
load_genome_counts <- function(track_dir, ds, stage) {
  f <- file.path(DATA, track_dir, ds, SYN_DECON_DIR,
                 paste0(ds, "_", stage, "_barcode_genome_counts.tsv.gz"))
  if (!file.exists(f) || file.info(f)$size == 0) {
    warning("Missing genome counts: ", f)
    return(NULL)
  }
  dt <- fread(f)
  agg <- dt[, .(B73     = sum(n_winner_reads[genome == "B73"]),
                non_B73 = sum(n_winner_reads[genome != "B73"])),
            by = barcode]
  agg[, `:=`(track_dir = track_dir, dataset = ds, stage = stage)]
  agg
}

cat("\nLoading pre/post genome counts for Panel C (60 files)...\n")
barn_rows <- list()
for (track_dir in names(SYN_TRACKS)) {
  for (ds in SYN_DATASETS) {
    for (stage in c("pre", "post")) {
      barn_rows[[length(barn_rows) + 1L]] <- load_genome_counts(track_dir, ds, stage)
    }
  }
}
barn <- rbindlist(barn_rows[!sapply(barn_rows, is.null)], fill = TRUE)

# Annotate with parsed alpha and contaminant
parsed <- lapply(barn$dataset, parse_synthetic_dataset)
barn[, true_alpha  := vapply(parsed, `[[`, numeric(1),   "true_alpha")]
barn[, contaminant := vapply(parsed, `[[`, character(1), "contaminant")]
barn[, track := factor(SYN_TRACKS[track_dir], levels = unname(SYN_TRACKS))]
barn[, stage := factor(stage, levels = c("pre", "post"),
                       labels = c("Pre-decontam", "Post-decontam"))]
cat("  ", format(nrow(barn), big.mark = ","), "barcode-rows loaded\n")

# ---- Panel B: singlet/doublet counts from pre-decontam C0 cells_calls ----
load_cells_calls_c0 <- function(track_dir, ds) {
  f <- file.path(DATA, track_dir, ds, "genotyping_runs", SYN_C0_TAG, "C0",
                 paste0(ds, "_cells_calls.tsv.gz"))
  if (!file.exists(f) || file.info(f)$size == 0) {
    warning("Missing cells_calls: ", f)
    return(NULL)
  }
  dt <- fread(f, select = c("barcode", "call"))
  dt[, `:=`(track_dir = track_dir, dataset = ds)]
}

cat("\nLoading C0 cells_calls for Panel B (30 files)...\n")
calls_rows <- list()
for (track_dir in names(SYN_TRACKS)) {
  for (ds in SYN_DATASETS) {
    calls_rows[[length(calls_rows) + 1L]] <- load_cells_calls_c0(track_dir, ds)
  }
}
calls_syn <- rbindlist(calls_rows[!sapply(calls_rows, is.null)], fill = TRUE)
parsed2 <- lapply(calls_syn$dataset, parse_synthetic_dataset)
calls_syn[, true_alpha  := vapply(parsed2, `[[`, numeric(1),   "true_alpha")]
calls_syn[, contaminant := vapply(parsed2, `[[`, character(1), "contaminant")]
calls_syn[, track := factor(SYN_TRACKS[track_dir], levels = unname(SYN_TRACKS))]

# Coarse class: singlet (single_clean + dirty_singlet), doublet (doublet + weak_doublet)
calls_syn[, call_class := fcase(
  call %in% c("single_clean", "dirty_singlet"), "Singlets",
  call %in% c("doublet", "weak_doublet"),       "Doublets",
  default = NA_character_
)]
counts_long <- calls_syn[!is.na(call_class),
                         .(N = .N), by = .(track, dataset, true_alpha,
                                           contaminant, call_class)]

# alpha_000 has no contaminant: duplicate the row under both Il14H and Ki11
# so each contaminant series anchors at the same point at alpha = 0.
alpha0_rows <- counts_long[contaminant == "none"]
counts_anchored <- rbind(
  counts_long[contaminant != "none"],
  alpha0_rows[, .(track, dataset, true_alpha, contaminant = "Il14H", call_class, N)],
  alpha0_rows[, .(track, dataset, true_alpha, contaminant = "Ki11",  call_class, N)]
)
counts_anchored[, contaminant := factor(contaminant, levels = c("Il14H", "Ki11"))]
counts_anchored[, call_class  := factor(call_class,  levels = c("Singlets", "Doublets"))]

# Fill missing (track, contaminant, call_class, true_alpha) combos with N = 0
# (high-alpha datasets often have 0 singlet calls under C0, which leaves gaps in
# the line plot otherwise).
counts_anchored <- as_tibble(counts_anchored) %>%
  complete(track, contaminant, call_class, true_alpha = SYN_ALPHA_LEVELS,
           fill = list(N = 0))

# Discrete-factor alpha axis: gives evenly spaced ticks (0, 0.02, 0.05, 0.1, ...)
# instead of crowding all sub-0.1 values together on a continuous axis.
counts_anchored <- counts_anchored %>%
  mutate(alpha_f = factor(alpha_lab(true_alpha), levels = ALPHA_LEVELS_CHR))

cat("  ", nrow(counts_anchored), "rows for plotting (after complete)\n")
cat("  N range by (track, call_class):\n")
print(as.data.table(counts_anchored)[
  , .(min_N = min(N), max_N = max(N), mean_N = round(mean(N), 1)),
  by = .(track, call_class)])

# =============================================================================
# 2) PANELS A, B, C
# =============================================================================

# ---- Panel A: NM density (full data, y-axis clipped to see both spike and tail) ----
track_colors <- c("Track B" = "#a6cee3", "Track B-disc" = "#1f78b4")

nm_stats <- nm_all[, .(
  mean_NM   = mean(NM),
  median_NM = median(NM),
  frac_zero = mean(NM == 0)
), by = track]

xq995 <- quantile(nm_all$NM, 0.995, na.rm = TRUE)

# Build legend label that includes the mean
nm_all[, track_lab := factor(
  sprintf("%s  (mean=%.2f, %.0f%% NM=0)",
          track, nm_stats$mean_NM[match(track, nm_stats$track)],
          100 * nm_stats$frac_zero[match(track, nm_stats$track)]),
  levels = sprintf("%s  (mean=%.2f, %.0f%% NM=0)",
                   nm_stats$track, nm_stats$mean_NM,
                   100 * nm_stats$frac_zero)
)]

lab_colors <- setNames(track_colors[as.character(nm_stats$track)],
                       levels(nm_all$track_lab))

pA <- ggplot(nm_all, aes(x = NM, fill = track_lab, color = track_lab)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, boundary = -0.5,
                 position = "identity", alpha = 0.5, linewidth = 0.3) +
  scale_color_manual(values = lab_colors, name = NULL) +
  scale_fill_manual(values = lab_colors, name = NULL) +
  scale_x_continuous(limits = c(-0.5, xq995 + 0.5),
                     breaks = scales::breaks_pretty(n = 5)) +
  coord_cartesian(ylim = c(0, 0.45)) +   # clip the NM=0 spike
  labs(x = "Edit distance per read (NM)",
       y = "Density",
       title = expression("Per-read SNP content (" * alpha * "=0)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.97, 0.97),
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# ---- Panel C: Barnyard ----
barn[, alpha_f := factor(alpha_lab(true_alpha), levels = ALPHA_LEVELS_CHR)]

pBarn <- ggplot(barn, aes(x = log10(B73 + 1), y = log10(non_B73 + 1),
                       color = alpha_f)) +
  geom_point(size = 0.6, alpha = 0.65) +
  scale_color_manual(values = SYN_ALPHA_PALETTE,
                     name = expression("True " * alpha),
                     drop = FALSE,
                     na.value = "grey70") +
  facet_grid(stage ~ track) +
  coord_fixed() +
  labs(x = expression(log[10] * "(B73 reads + 1)"),
       y = expression(log[10] * "(Il14H + Ki11 reads + 1)"),
       title = "Read distribution per barcode") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# ---- Panel B: call counts vs alpha (pre-decontam C0) ----
pCounts <- ggplot(counts_anchored,
             aes(x = alpha_f, y = N, color = contaminant,
                 shape = call_class, linetype = call_class,
                 group = interaction(contaminant, call_class))) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  scale_color_manual(values = SYN_CONTAM_COLORS[c("Il14H", "Ki11")],
                     name = "Contaminant") +
  scale_shape_manual(values = c(Singlets = 16L, Doublets = 17L),
                     name = "Call class") +
  scale_linetype_manual(values = c(Singlets = "solid", Doublets = "dashed"),
                        name = "Call class") +
  facet_wrap(~ track, ncol = 2) +
  labs(x = expression("True " * alpha),
       y = "Barcode count",
       title = "Pre-decontam call counts vs contamination") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.box = "horizontal")

# =============================================================================
# 3) ASSEMBLY & EXPORT (synthetic only)
# =============================================================================

# Layout: A density | B counts | C barnyard
fig <- (pA | pCounts | pBarn) +
  plot_layout(widths = c(1.0, 1.3, 1.7)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(OUTDIR, "Fig4_part1_AtoC.pdf"), fig,
       width = 17, height = 5.5)
ggsave(file.path(OUTDIR, "Fig4_part1_AtoC.png"), fig,
       width = 17, height = 5.5, dpi = 300)

# Individual panels
ggsave(file.path(OUTDIR, "Fig4A_synth_density.pdf"),  pA,      width = 4.5, height = 4.5)
ggsave(file.path(OUTDIR, "Fig4B_synth_counts.pdf"),   pCounts, width = 7.0, height = 4.5)
ggsave(file.path(OUTDIR, "Fig4C_synth_barnyard.pdf"), pBarn,   width = 7.5, height = 6)

cat("\nDone. Outputs in:", OUTDIR, "\n")
