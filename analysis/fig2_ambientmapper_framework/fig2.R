#!/usr/bin/env Rscript
# fig2.R -- Figure 2, panels B to E (noise anatomy of the genotyping model).
#   B empty gate: delta-BIC(empty) against JSD-to-eta with the noise-regime overlay,
#   C singlet-vs-doublet delta-BIC histogram by regime, D delta-BIC(empty) against read depth,
#   E the pre-clean barnyard cloud faceted by genotyping call.
# Inputs : data/processed/scifiATAC_B73_Arabidopsis/SM2v2/decontam_with_design_alpha05_v2/
#          (<SAMPLE>_cells_calls.decontam.tsv.gz for B to E, <SAMPLE>_pre_barcode_genome_counts.tsv.gz for E)
# Sources: analysis/_helpers/plotting.R
# Run    : Rscript analysis/fig2_ambientmapper_framework/fig2.R   (from the repository root;
#          FIG2_SAMPLE=SM2 rebuilds the panels from the run in SM2/ into figures/main/fig2_SM2)
#
# Noise regimes drawn on the panels:
#   1) Structural noise (distributional): "ambient soup" -> low JSD-to-eta + empty model favored
#   2) Statistical noise (insufficient evidence): weak model separation / low evidence (often low reads)
#   3) Biological mixing: doublet-favored with meaningful minor fraction (structured multi-genome signal)
#
# Panel B: Empty Gate (Structure vs Likelihood) with regime overlays
#   x = delta_empty = bic_non_empty - bic_empty  (positive => empty better)
#   y = jsd_to_eta  (low => soup-like; high => structured away from ambient)
# Panel C: Singlet vs Doublet separability with regime overlays + clear sign annotation
#   x = delta_sd_plot = (bic_single - bic_doublet)  (positive => doublet better)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(viridis)
})

source("analysis/_helpers/plotting.R")   # regime_cols, call_cols

# -------------------------
# 0) CONFIG
# -------------------------
DATA <- "data/processed/scifiATAC_B73_Arabidopsis"

# SM2v2 is the run shown in the manuscript. The run in SM2/ was produced by an older caller whose
# call vocabulary (doublet_confident / ambiguous_dirty_singlet / ambiguous_weak_doublet /
# indistinguishable) no longer exists; the schema shim in section 1b maps both vocabularies, so
# either run can be plotted. The two builds write to separate output directories.
SAMPLE <- Sys.getenv("FIG2_SAMPLE", "SM2v2")
if (!SAMPLE %in% c("SM2", "SM2v2")) stop("FIG2_SAMPLE must be SM2 or SM2v2, got: ", SAMPLE)
INDIR  <- file.path(DATA, SAMPLE, "decontam_with_design_alpha05_v2")
OUTDIR <- if (SAMPLE == "SM2v2") "figures/main/fig2" else "figures/main/fig2_SM2"
message("[Fig2] SAMPLE=", SAMPLE, "  INDIR=", INDIR, "  OUTDIR=", OUTDIR)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# plotting filter (keeps hex/density from being dominated by ultra-low reads)
MIN_READS_PLOT <- 10

# Regime-overlay thresholds. They are named after the AmbientMapper gates they illustrate
# (--genotyping-empty-bic-margin, --genotyping-bic-margin, --genotyping-doublet-minor-min) but are
# the figure's own display thresholds for the dashed gate lines and the regime colouring; they
# are not read from the run configuration.
EMPTY_BIC_MARGIN       <- 10
BIC_MARGIN_SD          <- 6
DOUBLET_MINOR_MIN      <- 0.20

# "low evidence" definition (statistical noise)
LOW_EVIDENCE_READS_MAX <- 10


# -------------------------
# Optional learned ambient + empty-gate parameters (from genotyping output)
# -------------------------
# The genotyping run emits eta as <SAMPLE>_eta.tsv.gz and does not write *_empty_jsd_tau.json, so
# both files are absent for SM2v2. `eta_final` is loaded only when present and is not used
# downstream; when tau is absent it is derived from the data in section 2b (the documented fallback).
eta_path <- file.path(DATA, SAMPLE, "final", paste0(SAMPLE, "_eta_final.json"))
tau_path <- file.path(DATA, SAMPLE, "final", paste0(SAMPLE, "_empty_jsd_tau.json"))

eta_final <- if (file.exists(eta_path)) jsonlite::fromJSON(eta_path) else NULL
tau_obj   <- if (file.exists(tau_path)) jsonlite::fromJSON(tau_path) else NULL
if (is.null(tau_obj)) message("[Fig2] NOTE: ", basename(tau_path),
    " absent -> tau will be derived from the data (see section 2b).")

JSD_TAU_MANUAL <- if (is.null(tau_obj)) NA_real_ else as.numeric(tau_obj$tau)
EMPTY_SEED_Q   <- if (is.null(tau_obj)) NA_real_ else as.numeric(tau_obj$quantile)
EMPTY_SEED_BIC <- if (is.null(tau_obj)) NA_real_ else as.numeric(tau_obj$seed_bic_min)

# -------------------------
# 1) LOAD DATA
# -------------------------
calls <- read_tsv(
  file.path(INDIR, paste0(SAMPLE, "_cells_calls.decontam.tsv.gz")),
  show_col_types = FALSE
)

# -------------------------
# 1b) SCHEMA SHIM
# -------------------------
# The SM2v2 cells_calls schema renamed three columns and dropped three derived ones.
# Every mapping below was verified EXACT against the 167,633-row SM2 file, where both
# the old and new quantities are present (0 mismatches at tol 1e-9):
#     ratio12            == ratio_top1_top2      ( == p_top1 / p_top2 )
#     purity_best        -> purity               (direct rename)
#     doublet_minor_frac -> minor                (direct rename)
#     bic_best_non_empty == pmin(bic_single, bic_doublet)
#     delta_empty        == bic_best_non_empty - bic_empty   (positive => empty better)
#     delta_math         == bic_doublet - bic_single         (negative => doublet better)
# NOTE: `purity`/`minor` are the fitted mixture parameters, NOT p_top1/p_top2 -- do not
# "simplify" them to the posterior masses, they differ on ~22% of barcodes.
.ren <- c(ratio12 = "ratio_top1_top2", purity_best = "purity",
          doublet_minor_frac = "minor")
for (from in names(.ren)) {
  to <- .ren[[from]]
  if (from %in% names(calls) && !to %in% names(calls)) {
    names(calls)[match(from, names(calls))] <- to
    message("[Fig2 shim] renamed ", from, " -> ", to)
  }
}
if (!"bic_best_non_empty" %in% names(calls) &&
    all(c("bic_single","bic_doublet") %in% names(calls))) {
  calls$bic_best_non_empty <- pmin(calls$bic_single, calls$bic_doublet)
  message("[Fig2 shim] derived bic_best_non_empty = pmin(bic_single, bic_doublet)")
}
if (!"delta_empty" %in% names(calls) &&
    all(c("bic_best_non_empty","bic_empty") %in% names(calls))) {
  calls$delta_empty <- calls$bic_best_non_empty - calls$bic_empty
  message("[Fig2 shim] derived delta_empty = bic_best_non_empty - bic_empty")
}
if (!"delta_math" %in% names(calls) &&
    all(c("bic_doublet","bic_single") %in% names(calls))) {
  calls$delta_math <- calls$bic_doublet - calls$bic_single
  message("[Fig2 shim] derived delta_math = bic_doublet - bic_single")
}
stopifnot(all(c("delta_empty","delta_math","jsd_to_eta","n_reads","call") %in% names(calls)))

# -------------------------
# 2) PREP DATA (robust parsing)
# -------------------------
num_or_na <- function(x) suppressWarnings(as.numeric(x))

calls <- calls %>%
  mutate(
    barcode = if ("barcode" %in% names(.)) as.character(barcode) else NA_character_,
    n_reads = num_or_na(n_reads),

    # core diagnostics
    jsd_to_eta   = num_or_na(jsd_to_eta),
    delta_empty  = num_or_na(delta_empty),   # bic_best_non_empty - bic_empty (positive => empty better)
    delta_math   = num_or_na(delta_math),    # python: bic_doublet - bic_single (negative => doublet better)

    # convert to "singlet - doublet" for intuitive plotting
    delta_sd_plot = -1 * delta_math,         # (bic_single - bic_doublet); positive => doublet better

    # optional fields
    minor  = if ("minor"  %in% names(.)) num_or_na(minor)  else NA_real_,
    purity = if ("purity" %in% names(.)) num_or_na(purity) else NA_real_,
    ratio_top1_top2 = if ("ratio_top1_top2" %in% names(.)) num_or_na(ratio_top1_top2) else NA_real_,
    p_top1 = if ("p_top1" %in% names(.)) num_or_na(p_top1) else NA_real_,

    call = as.character(call),

    # coarse call group for some aesthetics
    call_plot = case_when(
      str_detect(call, "empty") ~ "empty",
      str_detect(call, "single") ~ "single",
      str_detect(call, "doublet") ~ "doublet",
      str_detect(call, "ambiguous") ~ "ambiguous",
      str_detect(call, "indist") ~ "indistinguishable",
      TRUE ~ "ambiguous"
    )
  ) %>%
  filter(!is.na(n_reads), !is.na(jsd_to_eta), !is.na(delta_empty), !is.na(delta_sd_plot)) %>%
  filter(n_reads >= MIN_READS_PLOT)

# call table (diagnostic); status_flag is absent in the SM2v2 schema
if ("status_flag" %in% names(calls)) print(table(calls[, c("call","status_flag")])) else
  print(table(calls$call))

# -------------------------
# 2b) Derive JSD tau if not set
# -------------------------
# Idea: tau should represent "ambient-like structure".
# We derive a conservative tau from barcodes where empty is strongly favored
# (delta_empty >= EMPTY_BIC_MARGIN), taking a low quantile of their JSD distribution.
# A learned tau, when present, takes precedence (JSD_TAU_MANUAL above).
if (is.na(JSD_TAU_MANUAL)) {
  tau_pool <- calls %>% filter(delta_empty >= EMPTY_BIC_MARGIN)
  if (nrow(tau_pool) >= 100) {
    JSD_TAU <- unname(quantile(tau_pool$jsd_to_eta, probs = 0.90, na.rm = TRUE))
  } else {
    # fallback if there are few "empty-favored" points
    JSD_TAU <- unname(quantile(calls$jsd_to_eta, probs = 0.10, na.rm = TRUE))
  }
} else {
  JSD_TAU <- JSD_TAU_MANUAL
}

# -------------------------
# 2c) Regime classification (3 noise axes)
# -------------------------
# Structural noise: empty-favored AND ambient-like structure (low JSD).
# Statistical noise: weak model separation OR low-evidence tail (often low reads).
# Biological mixing: doublet-favored with meaningful minor fraction.
# Otherwise: structured cell.

calls <- calls %>%
  mutate(
    # "weak separation" flags (statistical uncertainty)
    weak_empty_sep = abs(delta_empty) < EMPTY_BIC_MARGIN,
    weak_sd_sep    = abs(delta_sd_plot) < BIC_MARGIN_SD,

    low_evidence   = n_reads <= LOW_EVIDENCE_READS_MAX,

    # biological mixing gate (requires minor if available; else uses delta alone)
    mixing_gate = case_when(
      !is.na(minor) ~ (delta_sd_plot >= BIC_MARGIN_SD) & (minor >= DOUBLET_MINOR_MIN),
      TRUE          ~ (delta_sd_plot >= BIC_MARGIN_SD)
    ),

    regime = case_when(
      (delta_empty >= EMPTY_BIC_MARGIN) & (jsd_to_eta <= JSD_TAU) ~ "structural_noise_soup",
      mixing_gate ~ "biological_mixing_doublet",
      (low_evidence & (weak_empty_sep | weak_sd_sep)) | (weak_empty_sep & weak_sd_sep) ~ "statistical_noise_low_evidence",
      TRUE ~ "structured_cell"
    ),

    regime = factor(regime, levels = c(
      "structured_cell",
      "biological_mixing_doublet",
      "statistical_noise_low_evidence",
      "structural_noise_soup"
    ))
  )

# -------------------------
# 3) PANEL B: Empty Gate (Structure vs Likelihood)
# -------------------------
# Point layer coloured by regime, restricted to barcodes with > 50 reads and |delta_empty| < 100.
set.seed(1)
overlay_B <- calls %>%
  filter(n_reads>50) %>%
  filter(delta_empty > -100, delta_empty < 100)

pB <- overlay_B %>%
  ggplot(aes(x = delta_empty, y = jsd_to_eta)) +
  # regime overlay
  geom_point(
    data = overlay_B,
    aes(color = regime),
    size = 0.35, alpha = 0.30, inherit.aes = TRUE
  ) +
  scale_color_manual(values = regime_cols, name = "Noise regime") +

  # gates
  geom_vline(xintercept = EMPTY_BIC_MARGIN, linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = JSD_TAU, linetype = "dashed", linewidth = 0.4) +

  # quadrant annotations (short, figure-ready)
  annotate(
    "label",
    x = 0.95 * 100, y = 0.40 * max(calls$jsd_to_eta, na.rm = TRUE),
    label = "Structured but empty-favored\n(rare; check thresholds)",
    label.size = 0.2, size = 3, alpha = 0.85
  ) +
  annotate(
    "label",
    x = 0.40 * 100, y = 0.15 * max(calls$jsd_to_eta, na.rm = TRUE),
    label = "Empty / soup:\nΔBIC≥margin + low JSD",
    label.size = 0.2, size = 3, alpha = 0.85
  ) +
  annotate(
    "label",
    x = -0.40 * 100, y = 0.40 * max(calls$jsd_to_eta, na.rm = TRUE),
    label = "Cells:\nstructured (high JSD)\nempty not favored",
    label.size = 0.2, size = 3, alpha = 0.85
  ) +
  annotate(
    "label",
    x = -0.40 * 100, y = 0.15 * max(calls$jsd_to_eta, na.rm = TRUE),
    label = "Low-evidence tail:\nweak separation",
    label.size = 0.2, size = 3, alpha = 0.85
  ) +

  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "B. Empty gate separates soup-like barcodes from structured cellular signal",
    subtitle = paste0(
      "ΔBIC_empty = bic_non_empty − bic_empty (positive favors Empty).  ",
      "JSD-to-η quantifies structure vs ambient.\n",
      "Dashed lines: ΔBIC margin = ", EMPTY_BIC_MARGIN,
      ",  JSD τ = ", signif(JSD_TAU, 3),
      if (is.na(JSD_TAU_MANUAL)) " (derived from data)" else " (learned)"
    ),
    x = expression(Delta*"BIC (empty advantage)  =  BIC"[non-empty]*" − BIC"[empty]),
    y = "JSD(P_barcode || η_ambient)"
  )

# -------------------------
# 4) PANEL C: Statistical Separability (Singlet vs Doublet)
# -------------------------
# Show ΔBIC (singlet - doublet) with clear sign annotation + mixing regime emphasis.
# Positive values => doublet has lower BIC (preferred).
calls_C <- calls %>%
  filter(call_plot %in% c("single", "doublet", "ambiguous", "indistinguishable")) %>%
  filter(n_reads>50) %>%
  mutate(delta_capped = pmax(-50, pmin(50, delta_sd_plot)))

pC <- ggplot(calls_C, aes(x = delta_capped)) +
  geom_histogram(aes(fill = regime), binwidth = 1, alpha = 0.6, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = regime_cols, name = "Noise regime") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.02))) +

  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.35) +
  geom_vline(xintercept = c(-BIC_MARGIN_SD, BIC_MARGIN_SD), linetype = "dashed", linewidth = 0.35) +

  annotate(
    "label",
    x = -18, y = Inf,
    label = "Left: singlet preferred\n( BIC_single < BIC_doublet )",
    vjust = 1.15, label.size = 0.2, size = 3, alpha = 0.85
  ) +
  annotate(
    "label",
    x = 24, y = Inf,
    label = "Right: doublet preferred\n( BIC_doublet < BIC_single )",
    vjust = 1.15, label.size = 0.2, size = 3, alpha = 0.85
  ) +
  annotate(
    "text",
    x = BIC_MARGIN_SD, y = Inf,
    label = paste0("±", BIC_MARGIN_SD, " = model margin"),
    vjust = 2.2, hjust = -0.05, size = 3.2, fontface = "italic"
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "C. Singlet vs doublet model selection is statistically separable",
    subtitle = "ΔBIC = (BIC_single − BIC_doublet). Positive values favor Doublet; negative values favor Singlet.",
    x = expression(Delta*"BIC (singlet − doublet)"),
    y = "Barcodes"
  )

# -------------------------
# 5) PANEL D: Evidence vs Empty preference (ΔBIC_empty vs reads)
# -------------------------
pD <- overlay_B %>%
  filter(n_reads>5) %>%
  mutate(
    delta_cap = pmax(-40, pmin(40, delta_empty))
  ) %>%
  ggplot(aes(x = delta_cap, y = n_reads)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Barcodes\n(log10)") +
  geom_vline(xintercept = EMPTY_BIC_MARGIN, linetype = "dashed", linewidth = 0.4) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotate("text", x = EMPTY_BIC_MARGIN + 10, y = 30, label = "Empty preferred",
           hjust = 0, size = 3) +
  annotate("text", x = EMPTY_BIC_MARGIN - 10, y = 30, label = "Cell preferred",
           hjust = 1, size = 3) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  geom_hline(yintercept = 10, linetype="dotted") +
  labs(
    title = "D. Evidence vs Empty preference",
    x = expression(Delta*BIC[empty]~"(bic_non_empty - bic_empty; capped)"),
    y = "Reads (log10)"
  )

# -------------------------
# 6) PANEL E: the pre-clean barnyard cloud, faceted by genotyping call
# -------------------------
pre_counts <- read_tsv(file.path(INDIR, paste0(SAMPLE, "_pre_barcode_genome_counts.tsv.gz")),
                       show_col_types = FALSE) %>%
  mutate(
    barcode = as.character(barcode),
    genome  = as.character(genome),
    n_winner_reads = as.numeric(n_winner_reads)
  )

top2 <- pre_counts %>%
  group_by(genome) %>%
  summarise(total = sum(n_winner_reads, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1:2) %>%
  pull(genome)

GENOME_X <- top2[1]
GENOME_Y <- top2[2]
message("[Fig2] Barnyard axes: ", GENOME_X, " vs ", GENOME_Y)

barnyard_calls <- pre_counts %>%
  filter(genome %in% c(GENOME_X, GENOME_Y)) %>%
  mutate(winner_reads_log10 = log10(n_winner_reads + 1)) %>%
  select(barcode, genome, winner_reads_log10) %>%
  pivot_wider(names_from = genome, values_from = winner_reads_log10, values_fill = 0) %>%
  left_join(calls %>% filter(n_reads >50) %>% select(barcode, call), by = "barcode") %>%
  mutate(call = replace_na(call, "unknown")) %>% filter(call != "unknown")

# Call levels are data-driven: known classes keep a fixed display order, any unseen class is
# appended rather than dropped, and nothing may become NA (a hard-coded level set silently
# collapsed this panel when the call vocabulary changed between runs).
.preferred <- c("single_clean", "dirty_singlet", "doublet", "weak_doublet",
                "ambiguous", "low_reads",                      # current vocabulary
                "doublet_confident", "ambiguous_dirty_singlet", # legacy SM2 vocabulary
                "ambiguous_weak_doublet", "ambiguous_low_depth", "indistinguishable")
.present <- unique(as.character(barnyard_calls$call))
.levels  <- c(intersect(.preferred, .present), setdiff(.present, .preferred))
barnyard_calls$call <- factor(barnyard_calls$call, levels = .levels)
stopifnot(!any(is.na(barnyard_calls$call)))
message("[Fig2 panel E] call levels: ", paste(.levels, collapse = ", "))


pE <- barnyard_calls %>%
  filter(!call %in% c("indistinguishable")) %>%
  ggplot(aes(x = .data[[GENOME_X]], y = .data[[GENOME_Y]])) +
  geom_hex(bins = 120) +
  scale_fill_viridis_c(option = "D", trans = "log10", name = "Barcodes\n(log10)") +
  facet_wrap( .~call, ncol=5) +
  theme_minimal(base_size = 11) +
  labs(
    title = "E. The Barnyard Cloud (post-genotyping)",
    x = paste0(GENOME_X, " winner reads (log10+1)"),
    y = paste0(GENOME_Y, " winner reads (log10+1)")
  )

# -------------------------
# 7) EXPORT
# -------------------------
ggsave(file.path(OUTDIR, "Fig2B_EmptyGate_NoiseRegimes.png"), pB, width = 10, height = 4, dpi = 320)
ggsave(file.path(OUTDIR, "Fig2C_DeltaBIC_NoiseRegimes.png"), pC, width = 10, height = 4, dpi = 320)

ggsave(file.path(OUTDIR, "Fig2B_EmptyGate_NoiseRegimes.pdf"), pB, width = 7, height = 5, device = cairo_pdf, dpi = 200)
ggsave(file.path(OUTDIR, "Fig2C_DeltaBIC_NoiseRegimes.pdf"), pC, width = 7, height = 5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig2D_EmptyGate_vs_reads.pdf"), pD, width = 7, height = 4, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Fig2E_Barnyard_per_call.pdf"), pE, width = 8, height = 2.5, device = cairo_pdf)

message("Saved: Fig2B to Fig2E in ", OUTDIR)
