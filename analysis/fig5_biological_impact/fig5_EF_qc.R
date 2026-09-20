#!/usr/bin/env Rscript
# =============================================================================
# fig5_EF_qc.R  -  Fig 5 panel E: TSS enrichment (pTSS) per species, PreClean vs PostClean, on the
#   combined-genome (concatenated ZmATcombined) SM2v2 fixed cell set. Also rebuilds the SUPERSEDED
#   peak-count-normalised FRiP rarefaction panel from two precomputed caches (not a manuscript panel).
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/compare/SM2_{B73,At}.fixedSet.minDepth200.pre_post.txt
#         figures/main/fig5/frip_rarefaction.tsv + frip_full_points.tsv   (caches, rarefaction panel only)
# Output  figures/main/fig5/Fig5_P1_QC_TSS.{pdf,png} (panel E), Fig5_P1_QC_FRiPnorm.*, Fig5_P1_QC_block.*
# Run     Rscript analysis/fig5_biological_impact/fig5_EF_qc.R        (from the repo root)
# =============================================================================
#
# ONLY THE TSS PANEL IS A MANUSCRIPT PANEL (E). The FRiPnorm panel is SUPERSEDED and must not be
# published: it rarefies to the common peak count with a RANDOM draw of K peaks, which is not a
# matched comparison. The stage with FEWER peaks is the anchor and draws K of K (never
# subsampled), while the other throws away 65% of its peaks at random. There is no unbiased rule
# and the two sensible ones bracket ZERO (At at K = 65,379: random +0.161, top-K -0.068), so the
# At FRiP effect is NOT identifiable. Never quote "+0.162 at equal peak count". The replacement,
# which reports both rules side by side and says so on the panel, is fig5_EF_qc_fripfair.R
# (manuscript panel F).
#
# TSS is unaffected by any of this: it is scored against a FIXED TSS annotation, so it has no
# peak-set confound. That is exactly why the per-cell quality claim belongs to it.
# Depth (reads/barcode) is intentionally omitted (shown in other figures).
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SOC    <- file.path(DATA, "socrates")
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
SP_LEVELS     <- c("Maize (B73)", "Arabidopsis (At)")
sp_labels     <- c(B73 = "Maize (B73)", At = "Arabidopsis (At)")
stage_cols_pp <- c(PreClean = "#FF83FA", PostClean = "#43CD80")   # TSS panel (stage = PreClean/PostClean)
stage_cols    <- c(Pre = "#FF83FA", Post = "#43CD80")             # rarefaction (Pre/Post)

# =============================================================================
# PANEL E -- TSS enrichment (pTSS), per species, Pre vs Post
# =============================================================================
load1 <- function(f, sp) { x <- fread(file.path(SOC, "compare", f)); x[, species := sp]; x }
x <- rbind(load1("SM2_B73.fixedSet.minDepth200.pre_post.txt", "B73"),
           load1("SM2_At.fixedSet.minDepth200.pre_post.txt",  "At"))
x[, stage := factor(stage, levels = c("PreClean", "PostClean"))]
x[, species_lab := factor(sp_labels[species], levels = SP_LEVELS)]
xf <- x[present == TRUE & !is.na(pTSS)]
tss_med <- xf[, .(pTSS = median(pTSS)), by = .(species_lab, stage)]

pTSS <- ggplot(xf, aes(stage, pTSS, fill = stage)) +
  geom_violin(trim = FALSE, alpha = 0.75, width = 0.9, linewidth = 0.2, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.55, fill = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_text(data = tss_med, aes(label = sprintf("%.2f", pTSS)), vjust = -0.5, size = 2.8, fontface = "bold", colour = "grey15") +
  facet_wrap(~ species_lab) +
  scale_fill_manual(values = stage_cols_pp) + coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey70")) +
  labs(title = "TSS enrichment (pTSS) - maintained after cleaning",
       subtitle = "Fixed TSS reference: no peak-set confound. Dashed = QC gate (0.2).",
       x = NULL, y = "pTSS (reads near TSS)")

# Panel E is written before the cache guard below, so it is produced even when the two
# rarefaction caches are absent.
ggsave(file.path(OUTDIR, "Fig5_P1_QC_TSS.pdf"),      pTSS,  width = 5.6, height = 3.0)
ggsave(file.path(OUTDIR, "Fig5_P1_QC_TSS.png"),      pTSS,  width = 5.6, height = 3.0, dpi = 300)

# =============================================================================
# SUPERSEDED PANEL -- peak-count-normalized FRiP (rarefaction), per species
# =============================================================================
rf <- file.path(OUTDIR, "frip_rarefaction.tsv"); pf <- file.path(OUTDIR, "frip_full_points.tsv")
if (!file.exists(rf) || !file.exists(pf))
  stop("Missing ", rf, " and/or ", pf, ".\n",
       "  These two tables are PRECOMPUTED CACHES (the per-cell FRiP rarefaction curves and the\n",
       "  full-peak-set points of the superseded FRiPnorm panel). No script in this repository\n",
       "  produces them (the development script that wrote them is archived); they ship as cache\n",
       "  files with the companion data. Copy both into ", OUTDIR, " to build the FRiPnorm and block\n",
       "  figures. Panel E (Fig5_P1_QC_TSS) has already been written above.")
rare <- fread(rf); pts <- fread(pf)
rare[, `:=`(stage = factor(stage, levels = c("Pre", "Post")), species = factor(species, levels = SP_LEVELS))]
pts[,  `:=`(stage = factor(stage, levels = c("Pre", "Post")), species = factor(species, levels = SP_LEVELS))]
anch <- unique(pts[, .(species, Kanchor)])

pFRiP <- ggplot(rare, aes(K, med, colour = stage, fill = stage)) +
  geom_ribbon(aes(ymin = med - sd, ymax = med + sd), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = anch, aes(xintercept = Kanchor), linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point(data = pts, size = 2, shape = 21, colour = "grey20") +
  facet_wrap(~ species, scales = "fixed") +   # shared x AND y across both facets (identical scales for direct cross-panel comparison)
  scale_x_log10(labels = scales::comma) +
  scale_colour_manual(values = stage_cols) + scale_fill_manual(values = stage_cols) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = "grey70")) +
  labs(title = "Peak-count-normalized FRiP - maintained/gained at equal #peaks",
       subtitle = "Median FRiP vs #peaks (subsampled); dashed = equal peak count.",
       x = "Number of peaks (log10)", y = "Median per-cell FRiP", colour = "Stage", fill = "Stage")

# =============================================================================
# EXPORT -- the superseded panel + the combined block
# =============================================================================
ggsave(file.path(OUTDIR, "Fig5_P1_QC_FRiPnorm.pdf"), pFRiP, width = 5.6, height = 3.2)
ggsave(file.path(OUTDIR, "Fig5_P1_QC_FRiPnorm.png"), pFRiP, width = 5.6, height = 3.2, dpi = 300)

block <- pTSS / pFRiP + plot_layout(heights = c(1, 1.15)) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUTDIR, "Fig5_P1_QC_block.pdf"), block, width = 6.5, height = 6)
ggsave(file.path(OUTDIR, "Fig5_P1_QC_block.png"), block, width = 6.5, height = 6, dpi = 300)

cat("[done] wrote Fig5_P1_QC_{TSS,FRiPnorm,block}.{pdf,png} to ", OUTDIR, "/\n", sep = "")
