#!/usr/bin/env Rscript
# =============================================================================
# fig5_GH_umap.R  -  Fig 5 panel G: per-genome UMAPs across the three stages (PreClean / WD / ND)
#   on the plate-split objects, every stage coloured by each cell's PreClean cluster of ORIGIN.
#   Also writes an own-stage-cluster variant (not a manuscript panel).
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step2_cluster/
#           <stage>_<genome>.mQCv6.updated_metadata_v7.<cfg>.txt   (6 files, via the helper)
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P2_UMAP_{B73v5,TAIR10}_{origin,cluster}.{pdf,png}, part2_stage_summary.tsv
# Run     Rscript analysis/fig5_biological_impact/fig5_GH_umap.R        (from the repo root)
# =============================================================================
#
#   PreClean  SM2_*            raw
#   wd        Clean.SM2v2wd_*  AmbientMapper decontam WITH plate design (WD)
#   nd        Clean.SM2v2_*    AmbientMapper decontam design-free       (ND)
# wd and nd are INDEPENDENT treatments of the same raw input, NOT a chain.
# Each stage is a SEPARATE embedding (re-TF-IDF -> re-SVD -> re-UMAP). Facet side by side;
#   never overlay, and never compare cluster IDs by number across stages.
#
# THE HEADLINE COLOURING IS ORIGIN (each cell's PreClean cluster), because cellIDs are
# identical across stages -- that is an exact carry-over, not fuzzy matching, and it is the
# only thing that makes colour comparable across panels. The own-stage-cluster variant is
# also written as an individual part.
#
# WHAT THIS PANEL MAY AND MAY NOT CLAIM:
#   MAY: "populations are stable across modes; the Leiden boundaries move."
#        The origin colouring shows the same coherent blocks in all 3 embeddings.
#   MAY NOT: "nd destroys structure" -- not supported by these plots.
#   MAY NOT: any cell-type identity for a cluster; clusters stay as bare IDs here (cell-type
#        annotation lives on the consensus meta-cells, panel J).
#   At specifically: a CONTINUUM -- clusters are provisional GRADIENT BINS. A Leiden
#        seed-stability control has not been run, so At re-partitioning must not be read as a
#        cleaning effect.
# The metadata files are R-exported with row.names: the header has ONE FEWER field than the
#   data rows. fread warns and adds `V1`. Select BY NAME, never by position.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source("analysis/_helpers/fig5_part2_helpers.R")   # SOC, PLATE, CFG, STAGES, palettes, loaders

p2 <- load_part2_cached()

# --- provenance, printed on every run ------------------------------------------
summ <- p2[, .(cells = .N, clusters = uniqueN(LouvainClusters),
               median_reads = as.integer(median(total))),
           by = .(genome, stage)][order(genome, stage)]
summ[, net_pct_of_pre := round(100 * cells / cells[stage == "PreClean"], 1), by = genome]

cat("=== FIG 5 PANEL G | cells and clusters per genome x stage (frozen config, ", MQC, ") ===\n", sep = "")
print(summ)
cat("\n  NOTE `net_pct_of_pre` is a NET count and is NOT retention - every mode both\n",
    "  drops AND gains cells. Set-membership retention is in fig5_GH_flow.R.\n", sep = "")

cat("\n=== Cluster sizes ===\n")
for (g in names(CFG)) for (s in names(STAGES)) {
  x <- p2[genome == g & stage == s, .N, by = LouvainClusters][order(LouvainClusters)]
  cat(sprintf("  %-7s %-9s  n=%-6d  clusters: %s\n", g, s, sum(x$N),
              paste(sprintf("%s(%d)", x$LouvainClusters, x$N), collapse = " ")))
}

cat("\n=== CAVEATS carried on these panels ===\n",
    "  * No cell-type labels: clusters are bare IDs on purpose (annotation lives on the\n",
    "    consensus meta-cells, panel J).\n",
    "  * At (TAIR10) is a continuum; a Leiden seed-stability control has not been run.\n",
    "    Do not read At re-partitioning as a cleaning effect.\n",
    "  * Claim allowed: populations stable, boundaries move. NOT \"nd destroys structure\".\n", sep = "")

# --- facet strips carry n + cluster count so each panel is self-documenting ----
lab_dt <- summ[, .(genome, stage,
                   strip = sprintf("%s\nn = %s  |  %d clusters",
                                   stage, format(cells, big.mark = ","), clusters))]
p2 <- merge(p2, lab_dt, by = c("genome", "stage"), sort = FALSE)
p2[, strip := factor(strip, levels = lab_dt[order(genome, stage), unique(strip)])]

# Scale point size to cell count: At has ~14x fewer cells than B73 and reads as
# empty at a shared size. Cosmetic only -- no points are dropped or clipped.
pt_for <- function(d) {
  if (d[, .N, by = stage][, max(N)] < 3000) c(size = .8,  alpha = .75)
  else                                      c(size = .25, alpha = .55)
}

base_theme <- function() {
  theme_void(base_size = 10) +
    theme(legend.position = "right",
          strip.text = element_text(size = 8.5, face = "bold", lineheight = 1.15,
                                    margin = margin(b = 4)),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey35", lineheight = 1.2))
}

# --- headline: colour every stage by ORIGIN (PreClean cluster) -----------------
SUBTITLE_ORIGIN <- list(
  B73v5 = paste("Each cell keeps its PreClean cluster colour into wd and nd, so colour is comparable across panels.",
                "The same populations stay coherent blocks in all three embeddings:\nthe cell populations are stable, what moves is where Leiden draws the boundaries.",
                "Grey = present post-clean but absent pre-clean."),
  TAIR10 = paste("Arabidopsis is a CONTINUUM - these clusters are provisional gradient bins, not cell types.",
                 "A Leiden seed-stability control has not been run,\nso re-partitioning here cannot yet be attributed to cleaning.",
                 "Grey = present post-clean but absent pre-clean.")
)

umap_by_origin <- function(g, letter) {
  d <- copy(p2[genome == g]); pt <- pt_for(d)
  setorder(d, pre_lab)                      # draw "(new)" last so it stays visible
  ggplot(d, aes(umap1, umap2, colour = pre_lab)) +
    geom_point(size = pt[["size"]], alpha = pt[["alpha"]]) +
    facet_wrap(~strip, nrow = 1, scales = "free") +
    scale_colour_manual(values = c(PRE_PAL, `(new)` = NEW_COL),
                        name = "PreClean\ncluster", drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    base_theme() +
    labs(title = sprintf("%s. %s  -  coloured by cluster of ORIGIN", letter, CFG[[g]]$label),
         subtitle = SUBTITLE_ORIGIN[[g]])
}

# --- variant: colour by each stage's OWN cluster (individual part only) --------
umap_by_cluster <- function(g) {
  d <- copy(p2[genome == g]); pt <- pt_for(d)
  d[, cl := factor(LouvainClusters, levels = sort(unique(LouvainClusters)))]
  ggplot(d, aes(umap1, umap2, colour = cl)) +
    geom_point(size = pt[["size"]], alpha = pt[["alpha"]]) +
    facet_wrap(~strip, nrow = 1, scales = "free") +
    scale_colour_viridis_d(option = "turbo", name = "Cluster", end = .92) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    base_theme() +
    labs(title = CFG[[g]]$label,
         subtitle = "Coloured by each stage's OWN Leiden cluster. Separate embeddings - cluster IDs are NOT comparable across stages.")
}

# both genomes are manuscript panel G
pP2_umap_B73v5  <- umap_by_origin("B73v5",  "G")
pP2_umap_TAIR10 <- umap_by_origin("TAIR10", "G")
pP2_umapcl_B73v5  <- umap_by_cluster("B73v5")
pP2_umapcl_TAIR10 <- umap_by_cluster("TAIR10")

# --- save individual parts -----------------------------------------------------
save_part <- function(p, name, w = 10.5, h = 4.2) {
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
save_part(pP2_umap_B73v5,    "Fig5_P2_UMAP_B73v5_origin")
save_part(pP2_umap_TAIR10,   "Fig5_P2_UMAP_TAIR10_origin")
save_part(pP2_umapcl_B73v5,  "Fig5_P2_UMAP_B73v5_cluster")
save_part(pP2_umapcl_TAIR10, "Fig5_P2_UMAP_TAIR10_cluster")

fwrite(summ, file.path(OUTDIR, "part2_stage_summary.tsv"), sep = "\t")

cat("\n[done] wrote Fig5_P2_UMAP_{B73v5,TAIR10}_{origin,cluster}.{pdf,png} + part2_stage_summary.tsv to ",
    OUTDIR, "/\n", sep = "")
