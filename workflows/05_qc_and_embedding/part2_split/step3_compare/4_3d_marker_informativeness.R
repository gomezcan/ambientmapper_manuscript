#!/usr/bin/env Rscript
###############################################################################
## 4_3d_marker_informativeness.R  --  Part-2 Step-4 (markers): empirical marker evaluation.
##
## The payoff of the uncapped "test all" panel: from the 4_3 per-marker z-scores, score each candidate
## marker by how cell-type-SPECIFIC it actually is in the SM2v2 data, and compare the yield of the
## canonical (Marand/Ecker) markers vs the PlantscRNAdb-added markers -- per source and per cell type.
## This tells us which of the ~thousands of carried candidates the data supports, so the panel can be
## pruned empirically (rather than pre-capped by guesswork).
##
## Informativeness of a marker = max z across clusters (a cell-type-specific marker peaks high in one
## cluster; a flat/noise marker never does). n_clusters_ge_thr = how sharp (1 = clean single-cluster).
## Optional concordance: does the marker's peak cluster's annotation == the marker's own cell type?
##
## Usage:
##   Rscript 4_3d_marker_informativeness.R <marker_zscore.tsv> <canonical_bed> <out_prefix> [z_thr=1.5] [cluster_annotation.tsv]
##     marker_zscore.tsv : 4_3 output (long: stage, geneID, name, species, type, type_label, cluster, zscore)
##     canonical_bed     : ORIGINAL canonical panel (markers.maize.Marand2025.bed / markers.At.Ecker2025.bed)
##                         -> tags each marker source = canonical vs PlantscRNAdb
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: Rscript 4_3d_marker_informativeness.R <marker_zscore.tsv> <canonical_bed> <out_prefix> [z_thr] [cluster_annotation.tsv]")
zp <- args[1]; canon_bed <- args[2]; outpref <- args[3]
z_thr <- if (length(args) >= 4) as.numeric(args[4]) else 1.5
annp  <- if (length(args) >= 5) args[5] else NA_character_

z <- read.table(zp, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
z$zscore <- as.numeric(z$zscore)
canon <- read.table(canon_bed, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
canon_ids <- unique(canon$geneID)

## ---- per-marker aggregate over clusters ----
byg <- split(z, z$geneID)
res <- do.call(rbind, lapply(byg, function(d) {
  o <- order(-d$zscore)
  data.frame(geneID = d$geneID[1], name = d$name[1], type = d$type[1], type_label = d$type_label[1],
             stage = d$stage[1], max_z = round(d$zscore[o[1]], 3), peak_cluster = d$cluster[o[1]],
             n_clusters_ge_thr = sum(d$zscore >= z_thr, na.rm = TRUE), stringsAsFactors = FALSE)
}))
res$source      <- ifelse(res$geneID %in% canon_ids, "canonical", "PlantscRNAdb")
res$informative <- res$max_z >= z_thr

## ---- optional peak-concordance (marker's peak cluster annotated as the marker's own type) ----
if (!is.na(annp) && file.exists(annp)) {
  ann <- read.table(annp, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  top_by_cl <- setNames(ann$top_type, as.character(ann$cluster))
  res$peak_cluster_top_type <- top_by_cl[as.character(res$peak_cluster)]
  res$peak_concordant <- !is.na(res$peak_cluster_top_type) & res$peak_cluster_top_type == res$type_label
}
res <- res[order(res$type_label, -res$max_z), ]
write.table(res, paste0(outpref, ".marker_informativeness.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- by source (headline: are the PlantscRNAdb additions pulling their weight?) ----
agg_src <- do.call(rbind, lapply(split(res, res$source), function(d) data.frame(
  source = d$source[1], n_markers = nrow(d), n_informative = sum(d$informative),
  frac_informative = round(mean(d$informative), 3),
  mean_max_z = round(mean(d$max_z), 3), median_max_z = round(median(d$max_z), 3),
  stringsAsFactors = FALSE)))
write.table(agg_src, paste0(outpref, ".informativeness_by_source.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- by cell type x source ----
res$ts <- paste(res$type_label, res$source, sep = "|")
agg_ts <- do.call(rbind, lapply(split(res, res$ts), function(d) data.frame(
  type_label = d$type_label[1], source = d$source[1], n_markers = nrow(d),
  n_informative = sum(d$informative), frac_informative = round(mean(d$informative), 3),
  mean_max_z = round(mean(d$max_z), 3), stringsAsFactors = FALSE)))
agg_ts <- agg_ts[order(agg_ts$type_label, agg_ts$source), ]
write.table(agg_ts, paste0(outpref, ".informativeness_by_type.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- console headline (robust to a missing source) ----
fi <- setNames(agg_src$frac_informative, agg_src$source)
gf <- function(s) if (s %in% names(fi)) sprintf("%.0f%%", 100 * fi[s]) else "NA"
message(sprintf(" - %s | markers=%d | informative(max_z>=%.1f)=%d (%.0f%%) | canonical=%s vs PlantscRNAdb=%s",
                basename(outpref), nrow(res), z_thr, sum(res$informative),
                100 * mean(res$informative), gf("canonical"), gf("PlantscRNAdb")))
