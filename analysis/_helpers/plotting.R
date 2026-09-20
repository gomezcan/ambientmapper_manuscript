# analysis/_helpers/plotting.R
# Shared colour palettes and theme fragments for the AmbientMapper manuscript figure scripts.
#
# Usage: source("analysis/_helpers/plotting.R") from the repository root, before any palette or
# theme fragment below is used. Sourced by fig1.R, fig2.R, fig3.R, fig3D_rescue_anatomy.R and figS1.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# =============================================================================
# Colour palettes
# =============================================================================

# --- Species / plate of origin (Fig 1 F to H, Fig 3E) ---
# Arabidopsis blue, maize (B73) red, used consistently across Fig 1 to 3.
cols_species <- c(
  "At"  = "#377eb8",
  "B73" = "#e41a1c"
)

# --- Figure 1: expected genome from the plate design (interspecies barnyard) ---
cols_expected_fill <- c(
  "Expected: Arabidopsis" = "#377eb8",
  "Expected: Maize"       = "#e41a1c"
)
cols_expected_line <- c(
  "Expected: Arabidopsis" = "#204d70",
  "Expected: Maize"       = "#8c0e0f"
)

# --- Figure 2: AmbientMapper call types ---
call_cols <- c(
  "empty"              = "#777777",
  "ambiguous"          = "#1874CD",
  "indistinguishable"  = "#CDC8B1",
  "doublet"            = "#CD1076",
  "single"             = "#BF3EFF"
)

# --- Figure 2: noise regime classification ---
regime_cols <- c(
  "structural_noise_soup"            = "#6C757D",   # gray
  "statistical_noise_low_evidence"   = "#F4A261",   # orange
  "biological_mixing_doublet"        = "#9B5DE5",   # purple
  "structured_cell"                  = "#2A9D8F"    # green
)

# --- Figure 3: cleaning outcome per barcode (post-decontamination) ---
# "Rescue" is the class the manuscript calls rescue barcodes: barcodes that looked like the wrong
# genome before cleaning and are high quality after it. The label is assigned by the Fig 3 scripts
# themselves (it is not a value read from any AmbientMapper output).
status_cols <- c(
  "Clean (Preserved)"    = "grey80",
  "Mixed"                = "#fdae61",
  "Heavily Contaminated" = "#6959CD",
  "Rescue"               = "#2c7bb6"
)
status_levels <- c("Clean (Preserved)", "Mixed", "Heavily Contaminated", "Rescue")

# --- Figure 3: cleaning stage ---
stage_cols <- c("Raw (Pre)" = "#FF83FA", "Cleaned (Post)" = "#43CD80")

# =============================================================================
# Theme fragments
# =============================================================================

# Polish layered on top of ggpubr::theme_pubclean() / theme_bw() by the Fig 3 panels.
nm_polish <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.25, color = "grey92"),
  strip.text       = element_text(face = "bold"),
  strip.background = element_rect(fill = "grey95", color = NA),
  axis.ticks       = element_line(color = "black", linewidth = 0.3),
  axis.line        = element_line(color = "black", linewidth = 0.3),
  plot.title       = element_text(face = "bold", size = 11),
  plot.subtitle    = element_text(size = 9, color = "grey30")
)

#' Minimal theme for publication figures
#' @param base_size Base font size (default 11)
theme_ambientmapper <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_blank()
    )
}
