#!/usr/bin/env Rscript
# Fig S1 (panels a, b): cross-species contamination under concatenated-reference mapping of the SM2
# scifi-ATAC library (B73v5 + TAIR10 concatenated reference, MAPQ >= 10, no AmbientMapper).
# a: hex density of per-barcode log10 maize vs Arabidopsis reads. b: contamination-rate profile per
# plate origin (bars = barcodes per bin, labels = cumulative fraction of the library captured).
# Inputs: data/processed/scifiATAC_B73_Arabidopsis/bed/SM2_{At,B73}/*ZmATcombined_scifiATAC.mq10.tn5.bed.gz
#   Source files: 3_Mapping/SM2_{At,B73}/bed/SM2_{At,B73}_ZmATcombined_scifiATAC.mq10.tn5.bed.gz
#   (per-library tn5 BEDs of the concatenated-reference mapping; not part of the public deposit).
# Run:    Rscript analysis/supplementary/figS1.R      (from the repository root)
# Output: figures/supplementary/figS1/FigS1.{pdf,png}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(patchwork)
  library(viridis)
})

# ---- CONFIG ----------------------------------------------------------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
BEDDIR <- file.path(DATA, "bed")
OUTDIR <- "figures/supplementary/figS1"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# species of a read = chromosome prefix in the concatenated reference
pattern_arab  <- "^At"
pattern_maize <- "^Zm"
bed_pattern   <- "ZmATcombined_scifiATAC.mq10.tn5.bed.gz$"

# barcodes with more than this many nuclear reads enter the contamination profile (panel b)
MIN_TOTAL_READS <- 10

# contamination-rate bins (fraction of a barcode's reads that map to the other species)
breaks_custom <- c(-Inf, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, Inf)
labels_custom <- c("0-1%", "1-2%", "2-4%", "4-8%", "8-16%", "16-32%", "32-64%", ">64%")

# one directory per plate half (library id and expected species)
inputs <- list(
  list(path = file.path(BEDDIR, "SM2_At"),  lib = "SM2", origin = "Expected: Arabidopsis"),
  list(path = file.path(BEDDIR, "SM2_B73"), lib = "SM2", origin = "Expected: Maize")
)

# ---- Helpers ---------------------------------------------------------------------------
chop <- function(myStr, mySep, myField) {
  choppedString <- sapply(strsplit(myStr, mySep), "[", myField)
  if (length(myField) > 1) {
    choppedString <- apply(choppedString, 2, function(x) paste0(x[!is.na(x)], collapse = mySep))
  }
  return(choppedString)
}

# Count reads per barcode and species from one tn5 BED (chr = column 1, cellID = column 4)
process_bed <- function(dir_path, library_id, expected_species) {
  bed_file <- list.files(dir_path, pattern = bed_pattern, full.names = TRUE)
  if (length(bed_file) == 0) {
    warning(paste("No BED file found in:", dir_path))
    return(NULL)
  } else if (length(bed_file) > 1) {
    warning(paste("Multiple BED files in:", dir_path, "- Using the first one."))
    bed_file <- bed_file[1]
  }
  message(paste("Processing:", library_id, "-", expected_species))
  dt <- fread(bed_file, select = c(1, 4), col.names = c("chr", "cellID"))

  # Classify reads by species
  dt[, species_read := fcase(
    grepl(pattern_arab, chr), "Arabidopsis",
    grepl(pattern_maize, chr), "Maize",
    default = "Other"
  )]
  # Nuclear reads only
  dt <- dt[species_read %in% c("Maize", "Arabidopsis")]

  counts <- dt[, .N, by = .(cellID, species_read)]
  counts_wide <- dcast(counts, cellID ~ species_read, value.var = "N", fill = 0)
  counts_wide[, Library := library_id]
  counts_wide[, Plate_Origin := expected_species]
  return(counts_wide)
}

# ---- Load ------------------------------------------------------------------------------
all_data_list <- lapply(inputs, function(x) process_bed(x$path, x$lib, x$origin))
barnyard_df <- rbindlist(all_data_list, fill = TRUE)

# barcode without the library suffix
barnyard_df$cellID <- chop(barnyard_df$cellID, "[_]", 1)

# ---- Metrics ---------------------------------------------------------------------------
# a file with reads of only one species would leave the other column missing after dcast
barnyard_df[is.na(Maize), Maize := 0]
barnyard_df[is.na(Arabidopsis), Arabidopsis := 0]

barnyard_df[, total_reads := Maize + Arabidopsis]
barnyard_df[, log_maize := log10(Maize + 1)]
barnyard_df[, log_arab := log10(Arabidopsis + 1)]

# contamination rate = reads of the other species / total, per plate origin
barnyard_df[, contamination_rate := ifelse(
  Plate_Origin == "Expected: Arabidopsis",
  Maize / total_reads,       # noise = maize reads in an Arabidopsis well
  Arabidopsis / total_reads  # noise = Arabidopsis reads in a maize well
)]

# ---- Panel a: barnyard hex density -----------------------------------------------------
p_hex <- ggplot(barnyard_df, aes(x = log_maize, y = log_arab)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Cell Count\n(Log10)") +
  labs(title = "Reads Density",
       subtitle = "High density on axes = Singlets;\nDensity in middle = Doublets/Ambient",
       x = "Log10 Maize Reads",
       y = "Log10 Arabidopsis Reads") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

# ---- Panel b: contamination profile ----------------------------------------------------
barnyard_df[, contam_bin := cut(contamination_rate,
                                breaks = breaks_custom,
                                labels = labels_custom,
                                include.lowest = TRUE)]

# per plate origin: barcodes per bin and cumulative fraction of the library captured
summary_stats <- barnyard_df[total_reads > MIN_TOTAL_READS] %>%
  group_by(Plate_Origin, contam_bin) %>%
  summarise(count = n(), .groups = 'drop_last') %>%
  mutate(
    cum_count = cumsum(count),
    total_cells = sum(count),
    cum_pct = cum_count / total_cells
  ) %>%
  ungroup()

p_single <- ggplot(summary_stats, aes(x = contam_bin, group = Plate_Origin)) +
  # bars: barcodes per bin
  geom_col(aes(y = count, fill = Plate_Origin),
           position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.7, color = "black", size = 0.2) +
  # line and points: cumulative count
  geom_line(aes(y = cum_count, color = Plate_Origin),
            position = position_dodge(width = 0.8),
            size = 1) +
  geom_point(aes(y = cum_count, color = Plate_Origin),
             position = position_dodge(width = 0.8),
             size = 2, shape = 21, fill = "white") +
  # labels: cumulative fraction of the library
  geom_text(aes(y = cum_count, label = scales::percent(cum_pct, accuracy = 1)),
            position = position_dodge(width = 0.8),
            vjust = -0.5,
            size = 3,
            fontface = "bold",
            show.legend = FALSE) +
  # log scale keeps the small bars of the contaminated bins visible
  scale_y_log10(labels = scales::comma, name = "Number of Barcodes (Log10)") +
  scale_fill_manual(values = c("Observed: Arabidopsis" = "#377eb8", "Observed: Maize" = "#e41a1c")) +
  scale_color_manual(values = c("Observed: Arabidopsis" = "#204d70", "Observed: Maize" = "#8c0e0f")) +
  labs(title = "Contamination Profile",
       subtitle = "Bars = Count per bin;\nLabels = Cumulative % of library captured",
       x = "Contamination Rate (Bin)",
       fill = "Mapping preferences",
       color = "Mapping preferences") +
  theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ---- Assembly --------------------------------------------------------------------------
final_fig <- (p_hex | p_single) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Cross-Species Contamination in Standard scifi-ATAC",
                  subtitle = "Maize vs. Arabidopsis (Standard Mapping to concatenate genomes)",
                  tag_levels = 'a') & theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "FigS1.pdf"), final_fig, width = 10, height = 7)
ggsave(file.path(OUTDIR, "FigS1.png"), final_fig, width = 10, height = 7, dpi = 350)
