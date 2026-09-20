#!/usr/bin/env Rscript
###############################################################################
## 4_3c_annotation_diff.R  --  Part-2 Step-4 (markers), cross-stage annotation diff.
##
## Joins the per-cluster marker-annotation coherence (4_3 cluster_annotation.tsv) with the Step-3
## cross-stage shedding (3_1 dropped_by_precluster.tsv + alignment.tsv) to answer two questions:
##
##  (A) CONTAMINATION CONFIRMATION -- which Pre clusters are ambient/junk rather than real cell types?
##      Signature = DISPROPORTIONATE shedding (frac_dropped >> the genome's overall shed rate) AND
##      LOW marker coherence (weak top cell-type z / high annotation entropy). For maize this should
##      isolate Pre cl5 (91% shed) and cl12 (61% shed); for At (a continuum) we expect ~none.
##
##  (B) MARKER-COHERENCE DELTA -- does cleaning sharpen cell-type identity? Compare the distribution of
##      per-cluster top-type z + annotation entropy, Pre vs Post-wd (and Pre-excluding-contamination).
##
##  (C) MATCHED PAIRS -- for each Post cluster, its aligned Pre cluster (Step-3 max-overlap): did the
##      call change, did coherence rise (delta top-score) / ambiguity fall (delta entropy)?
##
## Dependency-light (base R; optional ggplot2 scatter under tryCatch).
##
## Usage:
##   Rscript 4_3c_annotation_diff.R <genome> <pre_annot.tsv> <post_annot.tsv> <crossstage_dir> <out_dir>
##     pre_annot  = 4_3_annotation_<PreName>/<PreName>.cluster_annotation.tsv
##     post_annot = 4_3_annotation_<PostName>/<PostName>.cluster_annotation.tsv
##     crossstage_dir holds <genome>.crossstage.{dropped_by_precluster,alignment,summary}.tsv (from 3_1)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  stop("Usage: Rscript 4_3c_annotation_diff.R <genome> <pre_annot.tsv> <post_annot.tsv> <crossstage_dir> <out_dir>")
genome <- args[1]; pre_p <- args[2]; post_p <- args[3]; xdir <- args[4]; outdir <- args[5]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
op <- function(s) file.path(outdir, paste0(genome, s))
rd <- function(p) {
  if (!file.exists(p)) stop("missing input: ", p)
  read.table(p, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
}
md_table <- function(df) {
  c(paste0("| ", paste(colnames(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

pre  <- rd(pre_p);  pre$cluster  <- as.character(pre$cluster)
post <- rd(post_p); post$cluster <- as.character(post$cluster)
drop <- rd(file.path(xdir, paste0(genome, ".crossstage.dropped_by_precluster.tsv")))
alg  <- rd(file.path(xdir, paste0(genome, ".crossstage.alignment.tsv")))
summ <- rd(file.path(xdir, paste0(genome, ".crossstage.summary.tsv")))
drop$pre_cluster <- as.character(drop$pre_cluster)
overall_shed <- as.numeric(summ$frac_pre_dropped[1])
shed_thr <- round(1.3 * overall_shed, 4)     # "disproportionate" = 1.3x the genome-wide shed rate

## ---- (A) Pre-cluster contamination table: coherence x shedding --------------------------
A <- merge(pre, drop, by.x = "cluster", by.y = "pre_cluster", all.x = TRUE)
med_score <- median(A$top_score, na.rm = TRUE)
med_ent   <- median(A$shannon_entropy, na.rm = TRUE)
A$disproportionate_shed <- !is.na(A$frac_dropped) & A$frac_dropped >= shed_thr
A$low_coherence         <- (A$top_score <= med_score) | (A$shannon_entropy >= med_ent)
A$contamination_flag    <- A$disproportionate_shed & A$low_coherence
A <- A[order(-A$frac_dropped), ]
Acols <- c("cluster", "n_cells", "frac_dropped", "top_type", "top_score", "margin",
           "shannon_entropy", "dominant_species", "disproportionate_shed", "low_coherence",
           "contamination_flag")
Acols <- Acols[Acols %in% colnames(A)]
write.table(A[, Acols], op(".step4.pre_contamination.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
flagged <- A$cluster[A$contamination_flag]

## ---- (B) marker-coherence delta Pre vs Post-wd -----------------------------------------
cohsum <- function(df, stg) data.frame(
  stage = stg, n_clusters = nrow(df),
  mean_top_score   = round(mean(df$top_score, na.rm = TRUE), 3),
  median_top_score = round(median(df$top_score, na.rm = TRUE), 3),
  mean_margin      = round(mean(df$margin, na.rm = TRUE), 3),
  mean_entropy     = round(mean(df$shannon_entropy, na.rm = TRUE), 3),
  median_entropy   = round(median(df$shannon_entropy, na.rm = TRUE), 3),
  stringsAsFactors = FALSE)
B <- rbind(
  cohsum(pre, "PreClean_all"),
  cohsum(pre[!(pre$cluster %in% flagged), , drop = FALSE], "PreClean_noContam"),
  cohsum(post, "PostClean_wd"))
write.table(B, op(".step4.coherence_delta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- (C) matched pairs Pre <-> Post (via Step-3 alignment) -----------------------------
alg$post_cluster <- as.character(alg$post_cluster); alg$aligned_pre <- as.character(alg$aligned_pre)
M <- merge(alg, post, by.x = "post_cluster", by.y = "cluster", all.x = TRUE)
M <- merge(M, pre, by.x = "aligned_pre", by.y = "cluster", all.x = TRUE, suffixes = c(".post", ".pre"))
M$delta_top_score <- round(M$top_score.post - M$top_score.pre, 3)
M$delta_entropy   <- round(M$shannon_entropy.post - M$shannon_entropy.pre, 3)
M$call_changed    <- M$top_type.post != M$top_type.pre
Mcols <- c("post_cluster", "aligned_pre", "overlap", "jaccard",
           "top_type.pre", "top_score.pre", "shannon_entropy.pre",
           "top_type.post", "top_score.post", "shannon_entropy.post",
           "delta_top_score", "delta_entropy", "call_changed")
Mcols <- Mcols[Mcols %in% colnames(M)]
M <- M[order(-M$overlap), ]
write.table(M[, Mcols], op(".step4.matched_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- results.md ------------------------------------------------------------------------
ft <- A[A$contamination_flag, c("cluster", "n_cells", "frac_dropped", "top_type", "top_score",
                                "shannon_entropy")]
md <- c(
  paste0("# 4_3c annotation diff — ", genome, " (Pre vs Post-wd)"), "",
  paste0("Overall Pre cells dropped by cleaning: **", round(100 * overall_shed, 1),
         "%**  (disproportionate-shed threshold = ", round(100 * shed_thr, 1),
         "%; low-coherence = top_score ≤ ", round(med_score, 2),
         " or entropy ≥ ", round(med_ent, 2), ")"), "",
  "## (A) Contamination confirmation",
  paste0("**Flagged Pre clusters (disproportionate shedding + low coherence): ",
         if (length(flagged)) paste(flagged, collapse = ", ") else "none",
         "**", if (!length(flagged)) "  (consistent with a continuum — no discrete contamination cluster)" else ""),
  "",
  if (nrow(ft)) md_table(ft) else "_none flagged_", "",
  "## (B) Marker-coherence delta (Pre vs Post-wd)",
  "Higher top_score + lower entropy = sharper cell-type identity.", "",
  md_table(B), "",
  "_Compare PostClean_wd against PreClean_noContam for the fair 'real-cluster' delta; against",
  "PreClean_all for the total effect (which includes removing the contamination clusters)._", "",
  "## Tables",
  paste0("- `", genome, ".step4.pre_contamination.tsv` — per Pre cluster: shedding x coherence + flags"),
  paste0("- `", genome, ".step4.coherence_delta.tsv` — Pre/Pre-noContam/Post coherence summary"),
  paste0("- `", genome, ".step4.matched_pairs.tsv` — Post->Pre aligned pairs, call change + coherence deltas"))
writeLines(md, op(".step4.annotation_diff.results.md"))

## ---- optional scatter: shedding vs coherence (contamination corner) ---------------------
tryCatch({
  suppressWarnings(suppressMessages(library(ggplot2)))
  A$flag <- ifelse(A$contamination_flag, "contamination", "retained")
  p <- ggplot(A, aes(frac_dropped, top_score, label = cluster)) +
    geom_vline(xintercept = shed_thr, linetype = 2, color = "grey60") +
    geom_hline(yintercept = med_score, linetype = 2, color = "grey60") +
    geom_point(aes(size = n_cells, color = flag)) +
    geom_text(vjust = -0.8, size = 3) +
    scale_color_manual(values = c(contamination = "firebrick", retained = "grey30")) +
    theme_bw(base_size = 11) +
    labs(title = paste0(genome, ": Pre-cluster shedding vs marker coherence"),
         x = "frac cells dropped by cleaning (Step-3)", y = "top cell-type z (4_3)")
  dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(outdir, "plots", paste0(genome, ".step4.shed_vs_coherence.pdf")),
         p, width = 7, height = 5.5, device = grDevices::cairo_pdf)
}, error = function(e) message(" ! scatter skipped: ", conditionMessage(e)))

message(sprintf(" - %s DONE | flagged contamination clusters: %s | coherence Pre_all=%.2f Pre_noContam=%.2f Post=%.2f",
                genome, if (length(flagged)) paste(flagged, collapse = ",") else "none",
                B$mean_top_score[B$stage == "PreClean_all"],
                B$mean_top_score[B$stage == "PreClean_noContam"],
                B$mean_top_score[B$stage == "PostClean_wd"]))
