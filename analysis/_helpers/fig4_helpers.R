# fig4_helpers.R
# Shared constants and helpers for Fig 4 (parts 1 to 4) and Supplementary Figs S5 and S6.
# Source from the repo root: source("analysis/_helpers/fig4_helpers.R")
# HEAT_PALETTE() and heat_theme need ggplot2, scales and grid installed (called with :: prefixes).

# =============================================================================
# Constants
# =============================================================================
TRUE_GENOME <- "B73"   # ground-truth genotype of the maize root library (Root1)

call_levels <- c("single_clean", "dirty_singlet", "doublet", "weak_doublet",
                 "ambiguous", "empty")

call_colors <- c(
  single_clean   = "#2ca02c",
  dirty_singlet  = "#98df8a",
  doublet        = "#d62728",
  weak_doublet   = "#ff9896",
  ambiguous      = "#aec7e8",
  empty          = "#c7c7c7"
)

genome_highlight_colors <- c("TRUE" = "#2ca02c", "FALSE" = "#7f7f7f")

# =============================================================================
# Unified config encoding for the genotyping validation figures (S5, S6)
# =============================================================================

# Canonical order: baseline, knockouts, filter/threshold, phase-2-only, stacked
CONFIG_ORDER <- c(
  # Baseline
  "C0",
  # Knockouts
  "C1c_nofriend", "C1g_naked",
  # Filter/threshold (shared across phases)
  "C3b_mq50", "C4a_bic3", "C4d_wamb005_wdon",
  # Phase-2-only
  "C2a_xmap", "C2c_xmap_eta0", "C4c_eta0", "C4c_eta5", "C4d_wamb05_wdon",
  # Stacked progressive
  "S01_nofr_mq50", "S02_nofr_mq50_xaUL", "S03_nofr_mq50_xmap",
  "S04_nofr_mq50_xmap_xaUL", "S05_nofr_mq50_wamb",
  "S06_nofr_mq50_wamb_xaUL", "S07_nofr_mq50_xmap_wamb",
  # Stacked full
  "S08_full_xa0", "S09_full_xaUL", "S10_full_bic3",
  # Stacked sanity
  "S11_full_friendon", "S12_full_wamb01", "S13_full_eta2"
)

# Human-readable labels
CONFIG_LABELS <- c(
  C0                       = "C0: Baseline",
  C1c_nofriend             = "C1c: No friend rescue",
  C1g_naked                = "C1g: Minimal (no rescue)",
  C3b_mq50                 = "C3b: MAPQ 50",
  C4a_bic3                 = "C4a: BIC margin 3",
  C4d_wamb005_wdon         = "C4d: w_amb=0.05",
  C2a_xmap                 = "C2a: + xmap",
  C2c_xmap_eta0            = "C2c: + xmap, eta=0",
  C4c_eta0                 = "C4c: No eta learning",
  C4c_eta5                 = "C4c: Boosted eta (5)",
  C4d_wamb05_wdon          = "C4d: w_amb=0.5",
  S01_nofr_mq50            = "S01: MQ50",
  S02_nofr_mq50_xaUL       = "S02: MQ50 + XA unlim",
  S03_nofr_mq50_xmap       = "S03: MQ50 + xmap",
  S04_nofr_mq50_xmap_xaUL  = "S04: MQ50 + xmap + XA unlim",
  S05_nofr_mq50_wamb       = "S05: MQ50 + w_amb=0.5",
  S06_nofr_mq50_wamb_xaUL  = "S06: MQ50 + w_amb=0.5 + XA unlim",
  S07_nofr_mq50_xmap_wamb  = "S07: MQ50 + xmap + w_amb=0.5",
  S08_full_xa0             = "S08: Full stack (XA=0)",
  S09_full_xaUL            = "S09: Full stack (XA unlim)",
  S10_full_bic3            = "S10: Full + BIC 3 *",
  S11_full_friendon        = "S11: Full + friend ON",
  S12_full_wamb01          = "S12: Full, w_amb=0.1",
  S13_full_eta2            = "S13: Full, eta=2"
)

# Group membership
CONFIG_GROUPS <- c(
  C0                       = "baseline",
  C1c_nofriend             = "knockout",
  C1g_naked                = "knockout",
  C3b_mq50                 = "filter",
  C4a_bic3                 = "filter",
  C4d_wamb005_wdon         = "filter",
  C2a_xmap                 = "phase2_only",
  C2c_xmap_eta0            = "phase2_only",
  C4c_eta0                 = "phase2_only",
  C4c_eta5                 = "phase2_only",
  C4d_wamb05_wdon          = "phase2_only",
  S01_nofr_mq50            = "stacked",
  S02_nofr_mq50_xaUL       = "stacked",
  S03_nofr_mq50_xmap       = "stacked",
  S04_nofr_mq50_xmap_xaUL  = "stacked",
  S05_nofr_mq50_wamb       = "stacked",
  S06_nofr_mq50_wamb_xaUL  = "stacked",
  S07_nofr_mq50_xmap_wamb  = "stacked",
  S08_full_xa0             = "stacked_full",
  S09_full_xaUL            = "stacked_full",
  S10_full_bic3            = "stacked_full",
  S11_full_friendon        = "stacked_sanity",
  S12_full_wamb01          = "stacked_sanity",
  S13_full_eta2            = "stacked_sanity"
)

# Colors: one per config, grouped by family
CONFIG_COLORS <- c(
  # Baseline: dark gray
  C0                       = "#4d4d4d",
  # Knockouts: orange
  C1c_nofriend             = "#e6550d",
  C1g_naked                = "#fd8d3c",
  # Filter/threshold: blue
  C3b_mq50                 = "#3182bd",
  C4a_bic3                 = "#6baed6",
  C4d_wamb005_wdon         = "#9ecae1",
  # Phase-2-only: purple
  C2a_xmap                 = "#6a51a3",
  C2c_xmap_eta0            = "#807dba",
  C4c_eta0                 = "#9e9ac8",
  C4c_eta5                 = "#bcbddc",
  C4d_wamb05_wdon          = "#dadaeb",
  # Stacked progressive: green ramp (light to dark)
  S01_nofr_mq50            = "#c7e9c0",
  S02_nofr_mq50_xaUL       = "#a1d99b",
  S03_nofr_mq50_xmap       = "#74c476",
  S04_nofr_mq50_xmap_xaUL  = "#41ab5d",
  S05_nofr_mq50_wamb       = "#238b45",
  S06_nofr_mq50_wamb_xaUL  = "#006d2c",
  S07_nofr_mq50_xmap_wamb  = "#00441b",
  # Stacked full: teal/dark green
  S08_full_xa0             = "#005a32",
  S09_full_xaUL            = "#238443",
  S10_full_bic3            = "#d62728",
  # Stacked sanity: muted green
  S11_full_friendon        = "#78c679",
  S12_full_wamb01          = "#addd8e",
  S13_full_eta2            = "#d9f0a3"
)

# Shapes: unique per group, winner gets star
CONFIG_SHAPES <- c(
  C0                       = 16L,  # filled circle
  C1c_nofriend             = 17L,  # filled triangle
  C1g_naked                = 15L,  # filled square
  C3b_mq50                 = 16L,
  C4a_bic3                 = 17L,
  C4d_wamb005_wdon         = 15L,
  C2a_xmap                 = 16L,
  C2c_xmap_eta0            = 17L,
  C4c_eta0                 = 15L,
  C4c_eta5                 = 18L,  # diamond
  C4d_wamb05_wdon          = 3L,   # plus
  S01_nofr_mq50            = 16L,
  S02_nofr_mq50_xaUL       = 17L,
  S03_nofr_mq50_xmap       = 15L,
  S04_nofr_mq50_xmap_xaUL  = 18L,
  S05_nofr_mq50_wamb       = 16L,
  S06_nofr_mq50_wamb_xaUL  = 17L,
  S07_nofr_mq50_xmap_wamb  = 15L,
  S08_full_xa0             = 16L,
  S09_full_xaUL            = 17L,
  S10_full_bic3            = 8L,   # star (winner)
  S11_full_friendon        = 16L,
  S12_full_wamb01          = 17L,
  S13_full_eta2            = 15L
)

# Parameter table: one row per config, columns = the 9 tunable parameters.
# Values represent DEVIATIONS from the C0 baseline (the C0 row shows the defaults).
CONFIG_PARAM_TABLE <- tibble::tibble(
  method = CONFIG_ORDER,
  label  = CONFIG_LABELS[CONFIG_ORDER],
  group  = CONFIG_GROUPS[CONFIG_ORDER],
  mapq   = c(20, 20, 20, 50, 20, 20, 20, 20, 20, 20, 20,
             50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50),
  xa     = c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0",
             "0", "unlim", "0", "unlim", "0", "unlim", "0",
             "0", "unlim", "0", "0", "0", "0"),
  topk   = c("on", "on", "off", "on", "on", "on", "on", "on", "on", "on", "on",
             "on", "on", "on", "on", "on", "on", "on",
             "on", "on", "on", "on", "on", "on"),
  wdisc  = c("on", "on", "off", "on", "on", "on", "on", "on", "on", "on", "on",
             "on", "on", "on", "on", "on", "on", "on",
             "on", "on", "on", "on", "on", "on"),
  friend = c("on", "off", "off", "on", "on", "on", "on", "on", "on", "on", "on",
             "off", "off", "off", "off", "off", "off", "off",
             "off", "off", "off", "on", "off", "off"),
  xmap   = c("off", "off", "off", "off", "off", "off", "on", "on", "off", "off", "off",
             "off", "off", "on", "on", "off", "off", "on",
             "on", "on", "on", "on", "on", "on"),
  eta    = c(2, 2, 2, 2, 2, 2, 2, 0, 0, 5, 2,
             2, 2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 2),
  w_amb  = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1, 0.1, 0.5,
             0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.5),
  bic    = c(6, 6, 6, 6, 3, 6, 6, 6, 6, 6, 6,
             6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 6, 6, 6),
  phases = c("2,3", "2,3", "3", "2,3", "2,3", "2,3",
             "2", "2", "2", "2", "2",
             "3", "3", "3", "3", "3", "3", "3",
             "3", "3", "3", "3", "3", "3")
)

# =============================================================================
# Synthetic benchmark constants (Fig 4 panels A to C, Fig S5)
# =============================================================================

# Track directory -> display label
SYN_TRACKS <- c("synthetic" = "Track B", "synthetic_disc" = "Track B-disc")

# 15 datasets per track (mirrors DATASETS_SYN in workflows/03_genotyping/eval_phase_factorial.R)
SYN_DATASETS <- c(
  "alpha_000",
  paste0("alpha_", sprintf("%03d", c(2, 5, 10, 20, 30, 40, 50)), "_Il14H"),
  paste0("alpha_", sprintf("%03d", c(2, 5, 10, 20, 30, 40, 50)), "_Ki11")
)

# Numeric alpha levels (8 distinct values, including 0)
SYN_ALPHA_LEVELS <- c(0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

# Contaminant identity colors
SYN_CONTAM_COLORS <- c(Il14H = "#1f77b4", Ki11 = "#d62728", none = "#7f7f7f")

# Alpha-level palette (light -> dark for low -> high contamination)
SYN_ALPHA_PALETTE <- c(
  "0"    = "#fde725",
  "0.02" = "#a0da39",
  "0.05" = "#4ac16d",
  "0.1"  = "#1fa187",
  "0.2"  = "#277f8e",
  "0.3"  = "#365c8d",
  "0.4"  = "#46327e",
  "0.5"  = "#440154"
)

# parse_synthetic_dataset(): given a dataset name like "alpha_005_Il14H",
# return (true_alpha, contaminant). alpha_000 -> (0, "none").
# Encoding: 3-digit integer / 100 (e.g. alpha_005 = 0.05 = 5%).
parse_synthetic_dataset <- function(ds) {
  m <- regmatches(ds, regexec("^alpha_(\\d{3})(?:_(Il14H|Ki11))?$", ds))[[1]]
  if (length(m) < 2L) return(list(true_alpha = NA_real_, contaminant = NA_character_))
  list(
    true_alpha  = as.numeric(m[2]) / 100,
    contaminant = if (length(m) >= 3 && nzchar(m[3])) m[3] else "none"
  )
}

# =============================================================================
# Shared visual style for the validation heatmaps (S5, S6).
# =============================================================================
HEAT_PALETTE <- function(name = "Fraction")
  ggplot2::scale_fill_distiller(palette = "Blues", direction = 1,
                                 limits = c(0, 1),
                                 labels = scales::percent_format(),
                                 name = name, oob = scales::squish)

heat_theme <- ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid       = ggplot2::element_blank(),
                 axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
                 axis.text.y      = ggplot2::element_text(size = 7),
                 strip.text       = ggplot2::element_text(size = 7, face = "bold"),
                 plot.title       = ggplot2::element_text(size = 9, face = "bold"),
                 plot.subtitle    = ggplot2::element_text(size = 8),
                 legend.position  = "bottom",
                 legend.key.width = grid::unit(1.2, "cm"))
