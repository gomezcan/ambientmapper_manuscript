#!/usr/bin/env Rscript
# Fig S2 (panels A to E): Socrates meta-QC of the concatenated-reference SM2 object, PreClean stage.
# A to C: MetaQC cascades (depth, pTSS, FRiP, pOrg, dif, qc_check) for the Full, Arabidopsis and B73
# subsets of the same merged metadata (<name>.FRiP0.2.FULL / .FRiP0.2.At / .FRiP0.4.B73 *.QC_FIGURES.pdf).
# D, E: per-species violins and ECDFs of the depth-filtered (v1) distributions (<name>.<tags>.D_E.pdf).
# Input:  the merged Socrates metadata, e.g.
#   data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2/step1_integrate/SM2.full.metadata_updated.txt
# Usage:  Rscript analysis/supplementary/figS2.R <meta.tsv> <name> <depth_filter_raw> [stage] [outdir]
#   e.g.  Rscript analysis/supplementary/figS2.R <the file above> SM2 200 PreClean figures/supplementary/figS2
# Also writes <outdir>/<name>.<tags>.updated_metadata_v{1..6}.txt; the v4 table of the Full subset is the input of figS3.R.

suppressPackageStartupMessages({
  library(MASS)
  library(viridis)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggpubr)
  library(patchwork)
})

# ---- CONFIG ----------------------------------------------------------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis/socrates"
OUTDIR <- "figures/supplementary/figS2"
# Statistical-test annotations (Wilcoxon on the violins, KS labels on the ECDFs) are computed but
# not drawn by default. KS annotations not shown in the manuscript.
SHOW_TEST_ANNOTATIONS <- FALSE

# Required columns: cellID, species, stage, total, pTSS, FRiP, pOrg, tss_z, acr_z, log10nSites
# Optional: dif, unique (if missing, unique = total)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 5) {
  stop("Usage: Rscript figS2.R <meta.tsv> <name> <depth_filter_raw> [stage] [outdir]\n  e.g. Rscript figS2.R ",
       file.path(DATA, "SM2/step1_integrate/SM2.full.metadata_updated.txt"), " SM2 200 PreClean ", OUTDIR)
}
meta_path <- args[1]
name_base <- args[2]
depth_filter_raw <- suppressWarnings(as.numeric(args[3]))
stage_keep <- if (length(args) >= 4) args[4] else NA_character_
# outdir holds the per-version metadata tables and a plots/ subdirectory
outdir <- if (length(args) >= 5) args[5] else OUTDIR

if (!is.finite(depth_filter_raw) || depth_filter_raw <= 0) {
  stop("depth_filter_raw must be positive (e.g. 100, 200, 1000). Got: ", args[3])
}

dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

tag_depth <- paste0(".minDepth", as.integer(depth_filter_raw))
tag_stage <- if (!is.na(stage_keep)) paste0(".stage", stage_keep) else ".allStages"

# -------------------------
# MetaQC function
# -------------------------
MetaQC <- function(meta,
                   name,
                   dirout,
                   depth_filter_raw,
                   tss_z_filter = -8,
                   acr_z_filter = -8,
                   org_thresh   = 0.1,
                   diff_z_keep  = -1) {

  # depth filter applied on raw unique
  if (!("unique" %in% names(meta))) meta$unique <- meta$total

  meta$unique <- suppressWarnings(as.numeric(meta$unique))

  # required columns for QC panels
  req <- c("total","unique","pTSS","FRiP","pOrg","tss_z","acr_z")
  miss <- setdiff(req, names(meta))
  if (length(miss) > 0) stop("[MetaQC] Missing columns: ", paste(miss, collapse=", "))

  if (!("dif" %in% names(meta))) meta$dif <- NA_real_
  if (!("qc_check" %in% names(meta))) meta$qc_check <- 1

  # keep finite unique
  meta <- meta[is.finite(meta$unique) & meta$unique > 0, , drop=FALSE]
  if (nrow(meta) == 0) stop("[MetaQC] No finite unique >0 rows for: ", name)

  # ---------- helpers ----------
  .thresh_from_z <- function(df, zcol, metric_col, zcut, floor_val) {
    below <- df[df[[zcol]] < zcut & is.finite(df[[metric_col]]), , drop=FALSE]
    if (nrow(below) == 0) return(floor_val)
    max(floor_val, max(below[[metric_col]], na.rm=TRUE))
  }

  .plot_density <- function(df, y, ylab, main, ylims, thresh, xlims=c(2.0, 6.5), h=c(0.2,0.05)) {
    df2 <- df[is.finite(df$unique) & df$unique > 0 & is.finite(df[[y]]), , drop=FALSE]
    if (nrow(df2) < 25) {
      plot.new(); title(main=paste0(main, " (too few points: ", nrow(df2), ")"))
      return(invisible(NULL))
    }
    den <- kde2d(log10(df2$unique), df2[[y]], n=300, h=h, lims=c(xlims, ylims))
    image(den, useRaster=TRUE, col=c("white", rev(magma(100))),
          xlab="Unique insertions (log10)", ylab=ylab, main=main)
    grid(lty=1, lwd=0.5, col="grey90")
    abline(h=thresh, col="red", lty=2, lwd=1)
    box()
  }

  plotDist <- function(x, main="") {
    x <- x[order(x$unique, decreasing = TRUE), , drop=FALSE]
    rank  <- log10(seq_len(nrow(x)))
    depth <- log10(x$unique + 1)

    keep_idx <- which(x$unique >= depth_filter_raw)
    if (length(keep_idx) == 0) stop("[MetaQC] No barcodes pass depth_filter_raw=", depth_filter_raw, " for ", name)

    cells_n <- max(keep_idx)
    knee    <- rank[cells_n]
    reads   <- as.integer(x$unique[cells_n])

    plot(rank[(cells_n+1):length(rank)], depth[(cells_n+1):length(depth)],
         type="l", lwd=2, col="grey75", main=main,
         xlim=range(rank), ylim=range(depth),
         xlab="Barcode rank (log10)", ylab="Unique insertions (log10)")
    lines(rank[1:cells_n], depth[1:cells_n], lwd=2, col="darkorchid4")
    grid()
    abline(v=knee, col="red", lty=2, lwd=1)
    abline(h=depth[cells_n], col="red", lty=2, lwd=1)
    text(x=min(rank) + 0.1, y=min(depth) + 0.2,
         labels=paste0("# kept=", cells_n, " | cutoff unique>=", reads),
         adj=c(0,0))
    x[seq_len(cells_n), , drop=FALSE]
  }

  # ---------- output filenames ----------
  tag <- paste0(tag_depth, tag_stage)
  pdf_out <- file.path(paste0(dirout,"/","plots"), paste0(name, tag, ".QC_FIGURES.pdf"))

  # ---------- start plotting ----------
  pdf(pdf_out, width=3, height=14)
  layout(matrix(1:5, nrow=5))

  # v1: depth
  meta.v1 <- plotDist(meta, main=paste0(name, " | depth >= ", depth_filter_raw))
  write.table(meta.v1, file=paste0(dirout, "/", name, tag, ".updated_metadata_v1.txt"),
              sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  # v2: pTSS
  meta.v1 <- meta.v1[order(meta.v1$tss_z, decreasing=TRUE), , drop=FALSE]
  tss_thresh <- .thresh_from_z(meta.v1, "tss_z", "pTSS", tss_z_filter, floor_val=0.2)
  meta.v2 <- meta.v1[meta.v1$pTSS >= tss_thresh, , drop=FALSE]
  .plot_density(meta.v1, "pTSS", "Fraction reads in TSS",
                paste0(name, " | pTSS (>= ", signif(tss_thresh,3), ")"),
                ylims=c(0,1), thresh=tss_thresh)
  legend("topright", legend=paste0("# cells=", nrow(meta.v2)), bty="n")
  write.table(meta.v2, file=paste0(dirout, "/", name, tag, ".updated_metadata_v2.txt"),
              sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  # v3: FRiP
  meta.v2 <- meta.v2[order(meta.v2$acr_z, decreasing=TRUE), , drop=FALSE]
  frip_thresh <- .thresh_from_z(meta.v2, "acr_z", "FRiP", acr_z_filter, floor_val=0.2)
  meta.v3 <- meta.v2[meta.v2$FRiP >= frip_thresh, , drop=FALSE]
  .plot_density(meta.v2, "FRiP", "Fraction insertions in ACRs (FRiP)",
                paste0(name, " | FRiP (>= ", signif(frip_thresh,3), ")"),
                ylims=c(0,1), thresh=frip_thresh)
  legend("topright", legend=paste0("# cells=", nrow(meta.v3)), bty="n")
  write.table(meta.v3, file=paste0(dirout, "/", name, tag, ".updated_metadata_v3.txt"),
              sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  # v4: pOrg (<=)
  df_org <- meta.v3[is.finite(meta.v3$pOrg), , drop=FALSE]
  meta.v4 <- df_org[df_org$pOrg <= org_thresh, , drop=FALSE]
  .plot_density(df_org, "pOrg", "Fraction organelle (pOrg)",
                paste0(name, " | pOrg (<= ", org_thresh, ")"),
                ylims=c(0,0.5), thresh=org_thresh, h=c(0.2,0.03))
  legend("topright", legend=paste0("# cells=", nrow(meta.v4)), bty="n")
  write.table(meta.v4, file=paste0(dirout, "/", name, tag, ".updated_metadata_v4.txt"),
              sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  # v5: dif (optional; skipped gracefully if mostly NA)
  meta.dif <- meta.v4[is.finite(meta.v4$dif), , drop=FALSE]
  if (nrow(meta.dif) >= 50) {
    z_dif <- as.numeric(scale(meta.dif$dif))
    dif_thresh <- min(meta.dif$dif[z_dif >= diff_z_keep], na.rm=TRUE)
    meta.v5 <- meta.dif[meta.dif$dif >= dif_thresh, , drop=FALSE]

    ymin <- min(-1, min(meta.dif$dif, na.rm=TRUE))
    ymax <- max( 1, max(meta.dif$dif, na.rm=TRUE))
    .plot_density(meta.dif, "dif", "dif (good - bad)",
                  paste0(name, " | dif (>= ", signif(dif_thresh,3), ")"),
                  ylims=c(ymin,ymax), thresh=dif_thresh, h=c(0.2,0.1))
    legend("topright", legend=paste0("# cells=", nrow(meta.v5)), bty="n")
  } else {
    meta.v5 <- meta.v4
    plot.new()
    title(main=paste0(name, " | dif (skipped: finite dif n=", nrow(meta.dif), ")"))
  }
  write.table(meta.v5, file=paste0(dirout, "/", name, tag, ".updated_metadata_v5.txt"),
              sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  # v6: qc_check (if present)
  meta.v6 <- meta.v5[meta.v5$qc_check == 1, , drop=FALSE]
  write.table(meta.v6, file=paste0(dirout, "/", name, tag, ".updated_metadata_v6.txt"),
              sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  dev.off()

  message("[MetaQC] ", name, " | v1=", nrow(meta.v1),
          " v2=", nrow(meta.v2),
          " v3=", nrow(meta.v3),
          " v4=", nrow(meta.v4),
          " v5=", nrow(meta.v5),
          " v6=", nrow(meta.v6),
          " | wrote ", pdf_out)

  # v1 (depth-filtered) and v6 (fully filtered) tables for downstream use
  list(v1=meta.v1, v6=meta.v6, pdf=pdf_out)
}

# -------------------------
# Load merged metadata
# -------------------------
a <- read_tsv(meta_path, show_col_types = FALSE, progress = FALSE)

# stage subset if requested
if (!is.na(stage_keep)) {
  a <- a %>% filter(stage == stage_keep)
  if (nrow(a) == 0) stop("No rows after stage filter stage='", stage_keep, "'")
}

# ensure unique
if (!("unique" %in% names(a))) a <- a %>% mutate(unique = total)

# -------------------------
# Define 3 subsets from the SAME input
# -------------------------
full_meta <- a
at_meta   <- a %>% filter(species == "At")
b73_meta  <- a %>% filter(species == "B73")

if (nrow(at_meta) == 0) warning("At subset is empty (species == 'At'). Check labels.")
if (nrow(b73_meta) == 0) warning("B73 subset is empty (species == 'B73'). Check labels.")

# -------------------------
# Run MetaQC on each subset (same function). The FRiP0.2 / FRiP0.4 tags are name-only: the FRiP
# threshold itself derives from acr_z inside MetaQC.
# -------------------------
qc_full <- MetaQC(full_meta, paste0(name_base, ".FRiP0.2.FULL"), outdir, depth_filter_raw)
qc_at   <- MetaQC(at_meta,   paste0(name_base, ".FRiP0.2.At"),   outdir, depth_filter_raw)
qc_b73  <- MetaQC(b73_meta,  paste0(name_base, ".FRiP0.4.B73"),  outdir, depth_filter_raw)

# -------------------------
# Build ONE long_all (compatible with violin + ECDF + KS)
# -------------------------
metric_specs <- tribble(
  ~metric,        ~metric_label,        ~transform,
  "total",        "total (log10+1)",     "log10p1",
  "tss",          "tss (log10+1)",       "log10p1",
  "acrs",         "acrs (log10+1)",      "log10p1",
  "log10nSites",  "log10nSites",         "identity",
  "FRiP",         "FRiP",                "identity",
  "pTSS",         "pTSS",                "identity",
  "pOrg",          "pOrg",                 "identity"
)

# Ensure dif exists in each v1 (otherwise NA)
ensure_cols <- function(df) {
  if (!("dif" %in% names(df))) df$dif <- NA_real_
  df
}

qc_full$v1 <- ensure_cols(qc_full$v1)
qc_at$v1   <- ensure_cols(qc_at$v1)
qc_b73$v1  <- ensure_cols(qc_b73$v1)

make_long <- function(df, group_label) {
  df %>%
    mutate(group = group_label) %>%
    select(group, any_of(metric_specs$metric)) %>%
    pivot_longer(cols = -group, names_to = "metric", values_to = "value") %>%
    left_join(metric_specs, by = "metric") %>%
    mutate(
      value = suppressWarnings(as.numeric(value)),
      value_plot = case_when(
        transform == "log10p1" ~ log10(value + 1),
        TRUE ~ value
      )
    ) %>%
    filter(is.finite(value_plot)) %>%
    select(group, metric, metric_label, transform, value, value_plot)
}

long_all <- bind_rows(
  make_long(qc_at$v1,   "At"),
  make_long(qc_b73$v1,  "B73"),
  make_long(qc_full$v1, "Full")
) %>%
  mutate(group = factor(group, levels = c("At", "B73", "Full")))


# -------------------------
# 1) Violin plot (panel D)
# -------------------------
out_violin <- file.path(outdir, "plots", paste0(name_base, tag_depth, tag_stage, ".VIOLIN.v1.pdf"))

subset_metrics <- c("tss", "FRiP", "total", "pOrg")

long_all_violion <- long_all
long_all_violion$value_plot[long_all$metric == 'pOrg'] <- pmin(long_all_violion$value_plot[long_all$metric == 'pOrg'], 0.1, na.rm = TRUE)


long_all_violion$metric_label <- factor(long_all_violion$metric_label,
                                        levels = c("total (log10+1)",
                                                   "tss (log10+1)",
                                                   "acrs (log10+1)",
                                                   "log10nSites",
                                                   "FRiP",
                                                   "pOrg"))

p_violin <- long_all_violion %>%
  filter(group %in% c("At", "B73")) %>%
  filter(metric %in% subset_metrics) %>%
  ggplot(aes(x=group, y=value_plot, fill = group, color = group)) +
  geom_violin(trim = TRUE, linewidth = 0.5, alpha = 0.2) +
  scale_fill_manual(values = c("At" = "#377eb8", "B73" = "#e41a1c")) +
  scale_color_manual(values = c("At" = "#377eb8", "B73" = "#e41a1c")) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.6, color = "grey20", notch = TRUE) +
  facet_wrap(~ metric_label, scales = "free", ncol = 1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(size = 9)) +
  labs(
    title = paste0(
      name_base, " | depth-filtered (v1) distributions | minDepth=", depth_filter_raw,
      if (!is.na(stage_keep)) paste0(" | stage=", stage_keep) else ""
    ),
    x = NULL, y = NULL
  )
if (SHOW_TEST_ANNOTATIONS) p_violin <- p_violin + stat_compare_means()

ggsave(out_violin, p_violin, width = 2.5, height = 14, useDingbats = FALSE)
message("Wrote: ", out_violin)

# -------------------------
# 2) KS tests (At vs B73) per metric_label; the table is written, the labels are drawn only if
#    SHOW_TEST_ANNOTATIONS is TRUE
# -------------------------
ks_tbl <- long_all %>%
  filter(group %in% c("At", "B73")) %>%
  group_by(metric_label, metric, transform) %>%
  summarise(
    n_At   = sum(group == "At"),
    n_B73  = sum(group == "B73"),
    ks_D = {
      x <- value_plot[group == "At"]; y <- value_plot[group == "B73"]
      if (length(x) >= 20 && length(y) >= 20) suppressWarnings(as.numeric(ks.test(x, y)$statistic)) else NA_real_
    },
    ks_p = {
      x <- value_plot[group == "At"]; y <- value_plot[group == "B73"]
      if (length(x) >= 20 && length(y) >= 20) suppressWarnings(ks.test(x, y)$p.value) else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    ks_label = paste0(
      "At vs B73: D=", sprintf("%.3f", ks_D), ", p=", format.pval(ks_p, digits = 2, eps = 1e-300), "\n"
    )
  )


out_ks <- file.path(outdir, "plots", paste0(name_base, tag_depth, tag_stage, ".KS_At_vs_B73.tsv"))
write_tsv(ks_tbl, out_ks)
message("Wrote: ", out_ks)

# -------------------------
# 3) ECDF plot (panel E), same long_all
# -------------------------

# Ensure the facet order is stable (only set levels that exist)
wanted_levels <- c(
  "total (log10+1)", "tss (log10+1)", "acrs (log10+1)",
  "log10nSites", "FRiP", "pOrg"
)

long_all <- long_all %>%
  mutate(metric_label = factor(metric_label, levels = intersect(wanted_levels, unique(metric_label))))

# Build the annotation table so that it matches exactly the plotted facets
ann_pos <- long_all %>%
  filter(metric %in% subset_metrics) %>%
  distinct(metric_label) %>%
  left_join(ks_tbl %>% select(metric_label, ks_label), by = "metric_label") %>%
  mutate(x = Inf, y = Inf)

ann_pos$metric_label <- factor(ann_pos$metric_label, levels = wanted_levels)

p_ecdf <- long_all %>%
  filter(group %in% c("At", "B73")) %>%
  filter(metric %in% subset_metrics) %>%
  ggplot(aes(x = value_plot, color = group)) +
  stat_ecdf(geom = "step", linewidth = 0.6) +
  scale_color_manual(values = c("At" = "#377eb8", "B73" = "#e41a1c")) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  ) +
  labs(
    title = paste0(
      name_base, " | ECDF (v1 depth-filtered) | minDepth=", depth_filter_raw,
      if (!is.na(stage_keep)) paste0(" | stage=", stage_keep) else ""
    ),
    x = NULL, y = "Cumulative fraction"
  ) +
  facet_wrap(~ metric_label, scales = "free_x", ncol = 1, drop = TRUE)

if (SHOW_TEST_ANNOTATIONS) {
  p_ecdf <- p_ecdf +
    geom_text(
      data = ann_pos,
      aes(x = x, y = y, label = ks_label),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.1,
      size = 2.8
    )
}

out_de <- file.path(outdir, "plots", paste0(name_base, tag_depth, tag_stage, ".D_E.pdf"))

pde <- p_violin | p_ecdf
ggsave(out_de, pde, width = 4, height = 8, useDingbats = FALSE)
message("Wrote: ", out_de)
