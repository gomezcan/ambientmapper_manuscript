#!/usr/bin/env Rscript
###############################################################################
## 4_3k_cluster_annotation_from_rZ.R  --  INDEPENDENT cluster annotation from per-cell rZ.
##
## The clusters have NO assumed identity. For each cluster we ask: which marker cell-TYPE
## its highest-rZ markers belong to. Two outcomes, both informative:
##   * one type clearly wins            -> confident annotation = that type
##   * the top types tie at similar rZ  -> AMBIGUOUS: the cluster has no single identity
##                                          (a mix / continuum / organelle-ambient) -- itself a result.
##
## Method (per cluster cl, per metric):
##   markers that PEAK in cl  ->  for each type, best_rZ = max peak_rZ over its markers
##   rank types by best_rZ  ->  top_type = annotation ; ratio = 2nd_type_rZ / top_type_rZ
##   AMBIGUOUS if ratio >= amb_thr (default 0.80: the runner-up type is within 20% of the winner).
##
## Usage:
##   Rscript 4_3k_cluster_annotation_from_rZ.R <cluster_mean.tsv> <out_prefix> [metrics=euclid,geom] [amb_thr=0.80]
## Output: <out_prefix>.rZ_annotation.<metric>.tsv  +  console table.
###############################################################################

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 4_3k_cluster_annotation_from_rZ.R <cluster_mean.tsv> <out_prefix> [metrics] [amb_thr]")
cm_tsv <- args[1]; out_prefix <- args[2]
metrics <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("euclid", "geom")
amb_thr <- if (length(args) >= 4) as.numeric(args[4]) else 0.80
exclude <- if (length(args) >= 5) args[5] else ""     # regex on type_label to DROP (e.g. "dividing" = cell-cycle)

cm <- read.table(cm_tsv, header = TRUE, sep = "\t", quote = "", comment.char = "")
if (nzchar(exclude)) {
  n0 <- nrow(cm); cm <- cm[!grepl(exclude, cm$type_label, ignore.case = TRUE), ]
  message(sprintf(" - excluded %d markers matching type '%s'  (%d -> %d)", n0 - nrow(cm), exclude, n0, nrow(cm)))
}

for (metric in metrics) {
  pcol <- paste0("peak_cl_", metric); rcol <- paste0("peak_rZ_", metric)
  if (!all(c(pcol, rcol) %in% colnames(cm))) { message(" ! missing ", pcol, "/", rcol, " -- skip"); next }
  cls <- sort(unique(cm[[pcol]]))
  out <- do.call(rbind, lapply(cls, function(cl) {
    sub <- cm[cm[[pcol]] == cl, ]
    sub <- sub[order(-sub[[rcol]]), ]
    ## best rZ per TYPE among markers peaking in this cluster
    bytype <- tapply(sub[[rcol]], sub$type_label, max)
    bytype <- sort(bytype, decreasing = TRUE)
    top_type <- names(bytype)[1];  top_rZ <- as.numeric(bytype[1])
    snd_type <- if (length(bytype) >= 2) names(bytype)[2] else NA
    snd_rZ   <- if (length(bytype) >= 2) as.numeric(bytype[2]) else 0
    ratio    <- if (top_rZ > 0) snd_rZ / top_rZ else NA
    n_support <- sum(sub$type_label == top_type)            # how many markers back the winning type
    ## organelle tell: chloroplast/mito gene IDs among the cluster's top-5 markers
    top5 <- head(sub, 5)
    organelle <- sum(grepl("^ATC|^ATM|^ATCG|^ATMG", top5$geneID) | grepl("Pt|Mt", top5$geneID))
    data.frame(
      metric = metric, cluster = cl, n_markers = nrow(sub),
      top_type = top_type, top_rZ = round(top_rZ, 3), n_support = n_support,
      top_marker = sub$name[1], top_marker_type = sub$type_label[1],
      second_type = snd_type, second_rZ = round(snd_rZ, 3),
      ratio = round(ratio, 3),
      call = ifelse(!is.na(ratio) & ratio >= amb_thr, "AMBIGUOUS", "confident"),
      organelle_in_top5 = organelle,
      top3_types = paste(sprintf("%s(%.2f)", names(bytype)[1:min(3, length(bytype))],
                                 as.numeric(bytype[1:min(3, length(bytype))])), collapse = " ; ")
    )
  }))
  write.table(out, paste0(out_prefix, ".rZ_annotation.", metric, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n=====================  ", basename(out_prefix), " | metric = ", metric, "  =====================\n", sep = "")
  print(out[, c("cluster", "top_type", "top_rZ", "n_support", "second_type", "second_rZ", "ratio", "call", "organelle_in_top5")],
        row.names = FALSE)
}
cat("\n(ratio = 2nd_type_rZ / top_type_rZ ; >= ", amb_thr, " => AMBIGUOUS = competing types, no single identity)\n", sep = "")
