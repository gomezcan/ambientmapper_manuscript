#!/usr/bin/env Rscript
###############################################################################
## 3_1_crossstage_plate.R  --  Part-2 Step-3 cross-stage compare, PLATE arm
##
## Plate-arm sibling of 3_1_cluster_crossstage_compare.R (indep, 2-stage Pre-vs-wd).
## The plate arm has THREE stages per genome -- Pre (reference), wd (design-guided),
## nd (design-free) -- and the decisive plate finding is QC-CONDITIONED: which Pre
## cluster does cleaning remove, and does it carry an ambient signature?
##
## For the At plate: Pre=5 clusters, and nd strips 90% of Pre cluster 1
## (lowest depth, highest organelle + doublet) while wd leaves it -> nd, design-free,
## specifically removes an ambient-signature population. This script generalizes that.
##
## Computes, per genome, with Pre as reference and each cleaned stage (wd, nd):
##   (A) Pre per-cluster QC signature  : median depth/nSites/FRiP/pTSS/pOrg/dif/doublet + n
##   (B) retention per Pre cluster     : kept / dropped / %kept in each stage
##       + Fisher enrichment (is this cluster removed MORE than the overall rate?)
##   (C) contingency Pre x stage (shared cells) + adjusted Rand index + mover fraction
##       (adj_rand / greedy-alignment / movers reused from the indep 3_1 script)
##   Plots (ggplot2, skipped if unavailable -> tables still written):
##     retention bars (Pre cluster x stage), Pre-UMAP colored by kept/dropped per stage,
##     per-stage contingency heatmap.
##
## Usage:
##   Rscript 3_1_crossstage_plate.R <genome> <outdir> <pre_meta_v7> <wd_meta_v7|NA> <nd_meta_v7|NA>
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  stop("Usage: Rscript 3_1_crossstage_plate.R <genome> <outdir> <pre_meta_v7> <wd_meta_v7|NA> <nd_meta_v7|NA>")
genome <- args[1]; outdir <- args[2]
pre_path <- args[3]
stage_paths <- list(wd = args[4], nd = args[5])
stage_paths <- stage_paths[!vapply(stage_paths, function(p) is.na(p) || p == "NA" || !nzchar(p), logical(1))]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plotdir <- file.path(outdir, "plots"); dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
op <- function(s) file.path(outdir, paste0(genome, s))
have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (!have_ggplot) message(" ! ggplot2 unavailable -> TSVs only, plots skipped (rerun with ggplot2 installed for PDFs)")

## ---- robust metadata read (write.table row.names=TRUE -> field1 is rownames) ----
read_meta <- function(p) {
  m <- read.table(p, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  m <- m[, colnames(m) != "...1", drop = FALSE]
  m$cellID  <- if ("cellID" %in% colnames(m)) as.character(m$cellID) else rownames(m)
  if (!"LouvainClusters" %in% colnames(m)) stop("LouvainClusters missing in ", p)
  m$cluster <- as.character(m$LouvainClusters)
  m <- m[!is.na(m$cluster) & nzchar(m$cluster), ]
  m
}
pre <- read_meta(pre_path)
message(sprintf(" - %s Pre: %d cells / %d clusters | cleaned stages: %s",
                genome, nrow(pre), length(unique(pre$cluster)),
                if (length(stage_paths)) paste(names(stage_paths), collapse = ",") else "(none)"))

QC <- intersect(c("total", "nSites", "log10nSites", "FRiP", "pTSS", "pOrg", "dif", "doubletscore"), colnames(pre))

## ---- (A) Pre per-cluster QC signature ---------------------------------------
a <- aggregate(pre[QC], list(pre_cluster = pre$cluster), function(x) round(median(x, na.rm = TRUE), 4))
a$n <- as.integer(table(pre$cluster)[as.character(a$pre_cluster)])
a <- a[order(a$pre_cluster), c("pre_cluster", "n", QC)]
write.table(a, op(".crossstage.pre_qc_by_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - (A) Pre QC-by-cluster:"); print(a, row.names = FALSE)

## ---- adjusted Rand (label-invariant) reused from indep 3_1 -------------------
adj_rand <- function(a, b) {
  ct <- table(a, b); cs <- function(x) sum(choose(x, 2)); n <- length(a)
  ii <- cs(as.vector(ct)); ai <- cs(rowSums(ct)); bj <- cs(colSums(ct))
  expected <- ai * bj / choose(n, 2); maxi <- (ai + bj) / 2; denom <- maxi - expected
  if (denom == 0) return(1); (ii - expected) / denom
}

pre_cl <- setNames(pre$cluster, pre$cellID)
overall <- list(); retention_all <- list(); ct_list <- list()

for (st in names(stage_paths)) {
  s <- read_meta(stage_paths[[st]]); s_cl <- setNames(s$cluster, s$cellID)
  shared  <- intersect(pre$cellID, s$cellID)
  dropped <- setdiff(pre$cellID, s$cellID)
  n_overall_kept <- length(shared); n_overall <- nrow(pre)

  ## (B) retention per Pre cluster + Fisher enrichment vs the rest
  pre$in_st <- pre$cellID %in% s$cellID
  rt <- as.data.frame.matrix(table(pre$cluster, ifelse(pre$in_st, "kept", "dropped")))
  if (!"kept"    %in% colnames(rt)) rt$kept <- 0L
  if (!"dropped" %in% colnames(rt)) rt$dropped <- 0L
  rt$pre_cluster <- rownames(rt); rt$n <- rt$kept + rt$dropped
  rt$pct_kept <- round(100 * rt$kept / rt$n, 1)
  tot_kept <- sum(rt$kept); tot_drop <- sum(rt$dropped)
  rt$fisher_p <- NA_real_; rt$odds_removed <- NA_real_
  for (i in seq_len(nrow(rt))) {
    m2 <- matrix(c(rt$dropped[i], rt$kept[i], tot_drop - rt$dropped[i], tot_kept - rt$kept[i]), nrow = 2)
    ft <- suppressWarnings(fisher.test(m2))
    rt$fisher_p[i] <- signif(ft$p.value, 3); rt$odds_removed[i] <- round(unname(ft$estimate), 3)
  }
  rt$stage <- st
  rt <- rt[order(-rt$pct_kept), c("stage", "pre_cluster", "n", "kept", "dropped", "pct_kept", "odds_removed", "fisher_p")]
  retention_all[[st]] <- rt

  ## (C) contingency + ARI + movers on shared cells
  sp <- factor(pre_cl[shared]); sq <- factor(s_cl[shared])
  ct <- table(Pre = sp, Stage = sq); ct_list[[st]] <- ct
  ari <- adj_rand(as.character(sp), as.character(sq))
  ctm <- as.matrix(ct); pre_sizes <- rowSums(ctm); post_sizes <- colSums(ctm)
  map_s2p <- setNames(rownames(ctm)[apply(ctm, 2, which.max)], colnames(ctm))
  mover_frac <- round(mean(map_s2p[as.character(sq)] != as.character(sp)), 4)

  most_removed <- rt$pre_cluster[which.min(rt$pct_kept)]
  overall[[st]] <- data.frame(
    genome = genome, stage = st, n_pre = n_overall, n_stage = nrow(s),
    n_shared = length(shared), n_dropped = length(dropped),
    frac_pre_dropped = round(length(dropped) / n_overall, 4),
    n_clusters_pre = length(unique(pre$cluster)), n_clusters_stage = length(unique(s$cluster)),
    ARI_shared = round(ari, 4), mover_frac = mover_frac,
    most_removed_precluster = most_removed,
    most_removed_pct_kept = rt$pct_kept[rt$pre_cluster == most_removed], stringsAsFactors = FALSE)
  write.table(as.data.frame.matrix(ct), op(paste0(".crossstage.", st, ".contingency.tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  message(sprintf(" - (%s) dropped %d (%.1f%% of Pre) | ARI=%.3f movers=%.1f%% | most-removed Pre cl %s (%.1f%% kept)",
                  st, length(dropped), 100 * overall[[st]]$frac_pre_dropped, ari, 100 * mover_frac,
                  most_removed, overall[[st]]$most_removed_pct_kept))
}

retention <- do.call(rbind, retention_all)
write.table(retention, op(".crossstage.retention_by_precluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary_tbl <- do.call(rbind, overall)
write.table(summary_tbl, op(".crossstage.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - (B) retention_by_precluster + summary written")

## ---- plots (ggplot2 only) ---------------------------------------------------
if (have_ggplot && length(stage_paths)) {
  suppressPackageStartupMessages(library(ggplot2))
  base_theme <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

  ## retention bars: Pre cluster x stage, %kept
  rb <- retention; rb$pre_cluster <- factor(rb$pre_cluster, levels = sort(unique(rb$pre_cluster)))
  p_ret <- ggplot(rb, aes(pre_cluster, pct_kept, fill = stage)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = pct_kept), position = position_dodge(0.8), vjust = -0.3, size = 2.6) +
    scale_fill_brewer(palette = "Set2") + ylim(0, 105) + base_theme +
    labs(title = paste0(genome, "  --  % of each Pre cluster kept after cleaning"),
         subtitle = "low bar = cleaning removes that Pre cluster (candidate ambient population)",
         x = "Pre cluster", y = "% cells kept")
  ggsave(file.path(plotdir, paste0(genome, ".crossstage.retention.pdf")), p_ret, width = 7, height = 4.5, device = cairo_pdf)

  ## Pre UMAP colored by kept/dropped per stage (spatial removal) -- needs umap coords
  if (all(c("umap1", "umap2") %in% colnames(pre))) {
    um_list <- lapply(names(stage_paths), function(st) {
      s <- read_meta(stage_paths[[st]])
      data.frame(umap1 = pre$umap1, umap2 = pre$umap2, pre_cluster = pre$cluster,
                 fate = ifelse(pre$cellID %in% s$cellID, "kept", "dropped"), stage = st, stringsAsFactors = FALSE)
    })
    um <- do.call(rbind, um_list)
    p_um <- ggplot(um, aes(umap1, umap2, color = fate)) +
      geom_point(size = 0.5, alpha = 0.7) +
      facet_wrap(~stage) +
      scale_color_manual(values = c(kept = "grey70", dropped = "#d7301f")) +
      base_theme + theme(axis.text = element_blank(), axis.ticks = element_blank()) +
      labs(title = paste0(genome, "  --  Pre UMAP: which cells each stage removes"),
           subtitle = "red = dropped by that cleaning; clustering of red = a targeted (ambient) population", x = NULL, y = NULL)
    ggsave(file.path(plotdir, paste0(genome, ".crossstage.removal_umap.pdf")), p_um,
           width = 3.6 * length(stage_paths) + 1, height = 4.2, device = cairo_pdf)
  }

  ## per-stage contingency heatmaps
  for (st in names(ct_list)) {
    ctdf <- as.data.frame(ct_list[[st]]); colnames(ctdf) <- c("Pre", "Stage", "n")
    ctdf$frac_of_pre <- ctdf$n / ave(ctdf$n, ctdf$Pre, FUN = sum)
    p_hm <- ggplot(ctdf, aes(Stage, Pre, fill = frac_of_pre)) +
      geom_tile(color = "grey85") + geom_text(aes(label = ifelse(n > 0, n, "")), size = 2.4) +
      scale_fill_viridis_c(name = "frac of\nPre cluster", limits = c(0, 1)) + base_theme +
      theme(panel.grid = element_blank()) +
      labs(title = sprintf("%s  Pre -> %s (shared cells)", genome, st), x = paste0(st, " cluster"), y = "Pre cluster")
    ggsave(file.path(plotdir, paste0(genome, ".crossstage.", st, ".heatmap.pdf")), p_hm, width = 6, height = 5, device = cairo_pdf)
  }
  message(" - plots written to ", plotdir)
}
message(" - DONE ", genome)
