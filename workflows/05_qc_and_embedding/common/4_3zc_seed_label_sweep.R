#!/usr/bin/env Rscript
## =================================================================================================
## 4_3zc_seed_label_sweep.R -- top-type labels across ALL 100 seed partitions (the _S100 sweep),
## both genomes, 3 stages: the full-distribution version of 4_3zb's 5-seed answer.
##
## WHY. The top_type of a meta-cell rides on the partition; 4_3zb showed the label churns between
## draws (At Pre 46% seed-pair agreement). This engine scores EVERY _S100 seed partition and
## emits, per CONSENSUS meta-cell, how stable its label actually is -- the quantity a heatmap can
## print next to the label, and the input for choosing which labels are worth testing against the
## embedding null (within-stage ranking only).
##
## HOW. Lean re-implementation of 4_3p's marker rZ (aggregate raw perkb -> scale to 1e4 common
## total -> Zi across meta-cells / Zj across ALL genes -> rZ_geom -> per-type mean over its
## markers, 'dividing' excluded, >=3 markers) that loads the perkb matrix ONCE per object and
## loops memberships. RECONCILIATION GATE (hard stop): scored on the FROZEN deliverable seed-1
## membership, the per-type matrix must match 4_3p's on-disk permetacell_type_rZ.geom.tsv at
## max|diff| == 0 on the 4-digit grid, and top_type must agree on every meta-cell.
##
## OUTPUTS (<s5>/rZ_annotation/consensus/):
##   <prefix>.seed_label_sweep.permc_labels.tsv.gz   seed x seed-metacell top labels (long)
##   <prefix>.seed_label_sweep.cells.tsv.gz          per cell: modal label, modal freq,
##                                                   consensus label, agree fraction
##   <prefix>.seed_label_sweep.metacell.tsv          per CONSENSUS meta-cell: label_support
##                                                   (frac of cell x seed = consensus label),
##                                                   modal seed label + freq, n_cells
##   <genome>.seed_label_sweep.summary.tsv           per stage: pairwise agreement, modal freq,
##                                                   consensus agreement (all over 100 seeds)
##
## Labels are DESCRIPTIVE (the gene null licenses almost none as annotation). Stability is a
##   property of the description. Cross-stage: report, never rank.
##
## Usage: Rscript 4_3zc_seed_label_sweep.R <s5_dir> <At|maize> [n_seeds=100]
## =================================================================================================

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Usage: Rscript 4_3zc_seed_label_sweep.R <s5_dir> <At|maize> [n_seeds]")
S5  <- args[1]; GEN <- args[2]
NS  <- if (length(args) >= 3) as.integer(args[3]) else 100L
stopifnot(GEN %in% c("At", "maize"))
msg <- function(...) message(sprintf(...))

PFX <- if (GEN == "At") {
  c(Pre = "SM2_At_TAIR10", wd = "Clean.SM2v2wd_At_TAIR10", nd = "Clean.SM2v2_At_TAIR10")
} else {
  c(Pre = "SM2_B73_B73v5", wd = "Clean.SM2v2wd_B73_B73v5", nd = "Clean.SM2v2_B73_B73v5")
}
RUN_FROZEN <- if (GEN == "At") "seacells_cps50_F14235_N14" else "seacells_cps50_F44658_N286"
RUN_S100   <- paste0(RUN_FROZEN, "_S100")
S3   <- file.path(dirname(S5), "step3_compare")
MKB  <- file.path(dirname(dirname(dirname(S5))), "_data/markers",
                  paste0("markers.", GEN, ".informative_top15.bed"))
if (!file.exists(MKB))                       # resolve _data relative to project base
  MKB <- file.path(sub("/SM2v2_plate.*$", "", normalizePath(S5)), "_data/markers",
                   paste0("markers.", GEN, ".informative_top15.bed"))
RZC  <- file.path(S5, "rZ_annotation")
CONS <- file.path(S5, "consensus")
OUTD <- file.path(RZC, "consensus")

markers <- read.table(MKB, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
markers$type_label <- paste0(ifelse(grepl("^Zm", markers$geneID), "Zm", "At"), ":", markers$type)

## per-type mean rZ + top type for ONE membership --------------------------------------------------
score_membership <- function(sm, cells_all, mem, mk_rows, type_of, types) {
  mc_lv <- unique(mem$mc)
  Ind <- sparseMatrix(i = match(mem$cellID, cells_all), j = match(mem$mc, mc_lv),
                      x = 1, dims = c(length(cells_all), length(mc_lv)))
  agg <- as.matrix(sm %*% Ind)                                   # genes x mc
  agg <- sweep(agg, 2, pmax(colSums(agg), .Machine$double.eps), "/") * 1e4
  cmean <- colMeans(agg)
  csd   <- apply(agg, 2, sd); csd[csd == 0 | is.na(csd)] <- 1
  smk   <- agg[mk_rows, , drop = FALSE]
  rmean <- rowMeans(smk)
  rsd   <- apply(smk, 1, sd); rsd[rsd == 0 | is.na(rsd)] <- 1
  Zi <- sweep(sweep(smk, 1, rmean, "-"), 1, rsd, "/")
  Zj <- sweep(sweep(smk, 2, cmean, "-"), 2, csd,  "/")
  Zi[Zi < 0] <- 0; Zj[Zj < 0] <- 0
  rZ <- sqrt(Zi * Zj)
  tm <- do.call(rbind, lapply(types, function(ty) colMeans(rZ[type_of == ty, , drop = FALSE])))
  rownames(tm) <- types; colnames(tm) <- mc_lv
  list(type_mean = tm, top = rownames(tm)[max.col(t(tm), ties.method = "first")], mc_lv = mc_lv)
}

sum_rows <- list()
for (st in names(PFX)) {
  prefix <- PFX[[st]]
  t0 <- Sys.time()
  sm <- readRDS(file.path(S3, paste0(prefix, ".plate.perkb.genes.sparse.rds")))

  frozen_c2s <- if (GEN == "At")
    file.path(S5, RUN_FROZEN, "seed1", paste0(prefix, ".seacells.cell_to_seacell.tsv")) else
    file.path(S5, RUN_FROZEN, paste0(prefix, ".seacells.cell_to_seacell.tsv"))
  read_mem <- function(f) { a <- read.delim(f); colnames(a)[1] <- "cellID"
                            data.frame(cellID = a$cellID, mc = as.character(a$SEACell)) }
  mem0  <- read_mem(frozen_c2s)
  cells <- intersect(colnames(sm), mem0$cellID)
  sm    <- sm[, cells, drop = FALSE]
  mk_rows <- intersect(rownames(sm), markers$geneID)
  type_of <- markers$type_label[match(mk_rows, markers$geneID)]
  keep_t  <- !grepl("dividing", type_of, ignore.case = TRUE)
  tt      <- table(type_of[keep_t])
  types   <- names(tt)[tt >= 3]
  mk_rows <- mk_rows[keep_t & type_of %in% types]; type_of <- markers$type_label[match(mk_rows, markers$geneID)]

  ## GATE: frozen seed-1 membership must reproduce 4_3p's on-disk type matrix ---------------------
  g0  <- score_membership(sm, cells, mem0[mem0$cellID %in% cells, ], mk_rows, type_of, types)
  ref_f <- if (GEN == "At") file.path(RZC, "seed1",  paste0(prefix, ".permetacell_type_rZ.geom.tsv"))
           else             file.path(RZC, "pooled", paste0(prefix, ".permetacell_type_rZ.geom.tsv"))
  ref <- read.delim(ref_f, check.names = FALSE)
  ref$metacell <- sub("^s[0-9]+:", "", ref$metacell)
  stopifnot(setequal(ref$metacell, g0$mc_lv), all(types %in% colnames(ref)))
  dd <- abs(t(as.matrix(ref[match(g0$mc_lv, ref$metacell), types])) - round(g0$type_mean, 4))
  ref_top <- ref$top_type[match(g0$mc_lv, ref$metacell)]
  same <- ref_top == g0$top
  ## argmax is order-dependent at ties: where labels disagree, require the two labels' values to
  ## be tied on the 4-digit grid -- then the matrix (the math) agrees and only the tiebreak differs.
  ties_ok <- all(vapply(which(!same), function(i) {
    v <- g0$type_mean[, i]
    ref_top[i] %in% names(v) && abs(v[ref_top[i]] - v[g0$top[i]]) <= 1e-4 + 1e-12
  }, TRUE))
  msg("%s %s | GATE vs 4_3p (frozen seed1): max|diff| = %.4g, top_type agree = %d/%d%s",
      GEN, st, max(dd), sum(same), length(same),
      if (all(same)) "" else if (ties_ok) " (disagreements are grid-ties: OK)" else " (REAL mismatch)")
  if (max(dd) > 0 || !ties_ok) stop("GATE FAIL: lean scorer does not reproduce 4_3p")

  ## sweep the _S100 partitions -------------------------------------------------------------------
  lab_idx <- matrix(NA_integer_, length(cells), NS)   # cells x seeds, integer type codes
  rownames(lab_idx) <- cells
  permc <- list()
  for (s in seq_len(NS)) {
    f <- file.path(S5, RUN_S100, paste0("seed", s), paste0(prefix, ".seacells.cell_to_seacell.tsv"))
    mm <- read_mem(f); mm <- mm[mm$cellID %in% cells, ]
    if (nrow(mm) != length(cells)) stop("seed ", s, ": cell set differs (", nrow(mm), ")")
    gg <- score_membership(sm, cells, mm, mk_rows, type_of, types)
    lab_idx[mm$cellID, s] <- match(gg$top, types)[match(mm$mc, gg$mc_lv)]
    permc[[s]] <- data.frame(stage = st, seed = s, metacell = gg$mc_lv, top_type = gg$top)
    if (s %% 25 == 0) msg("  %s %s: seed %d/%d (%.1f s elapsed)", GEN, st, s, NS,
                          as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  if (anyNA(lab_idx)) stop("unassigned cells in the label matrix")

  ## stability statistics -------------------------------------------------------------------------
  cnt <- t(apply(lab_idx, 1, tabulate, nbins = length(types)))    # cells x types
  pairw <- sum(cnt * (cnt - 1)) / (nrow(cnt) * NS * (NS - 1))     # mean per-cell seed-pair agreement
  modal_i <- max.col(cnt, ties.method = "first")
  modal_f <- cnt[cbind(seq_len(nrow(cnt)), modal_i)] / NS

  cc <- read.delim(file.path(CONS, paste0(prefix, ".consensus.cells.tsv")))
  hi <- read.delim(file.path(OUTD, paste0(prefix, ".consensus_rZ.heatmap_input.tsv")))
  cons_lab <- setNames(unique(hi[, c("metacell", "top_type")])$top_type,
                       unique(hi[, c("metacell", "top_type")])$metacell)
  cellmc   <- setNames(paste0("cMC-", cc$consensus_mc), cc$cellID)[cells]
  cellcons <- cons_lab[cellmc]
  cons_i   <- match(cellcons, types)
  agree_cons <- ifelse(is.na(cons_i), NA,
                       cnt[cbind(seq_len(nrow(cnt)), cons_i)] / NS)

  mc_tab <- do.call(rbind, lapply(split(seq_along(cells), cellmc), function(ix) {
    sub <- cnt[ix, , drop = FALSE]; tot <- colSums(sub)
    data.frame(n_cells = length(ix),
               consensus_label = cellcons[ix[1]],
               label_support = round(sum(sub[cbind(seq_along(ix), rep(cons_i[ix[1]], length(ix)))]) /
                                     (length(ix) * NS), 4),
               modal_seed_label = types[which.max(tot)],
               modal_seed_freq = round(max(tot) / (length(ix) * NS), 4))
  }))
  mc_tab <- data.frame(stage = st, metacell = rownames(mc_tab), mc_tab)
  mc_tab$labels_agree <- mc_tab$consensus_label == mc_tab$modal_seed_label

  gz1 <- gzfile(file.path(OUTD, paste0(prefix, ".seed_label_sweep.permc_labels.tsv.gz")), "w")
  write.table(do.call(rbind, permc), gz1, sep = "\t", quote = FALSE, row.names = FALSE); close(gz1)
  gz2 <- gzfile(file.path(OUTD, paste0(prefix, ".seed_label_sweep.cells.tsv.gz")), "w")
  write.table(data.frame(cellID = cells, consensus_mc = cellmc, consensus_label = cellcons,
                         modal_label = types[modal_i], modal_freq = round(modal_f, 4),
                         agree_with_consensus = round(agree_cons, 4)),
              gz2, sep = "\t", quote = FALSE, row.names = FALSE); close(gz2)
  write.table(mc_tab, file.path(OUTD, paste0(prefix, ".seed_label_sweep.metacell.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sum_rows[[st]] <- data.frame(genome = GEN, stage = st, n_seeds = NS, n_cells = length(cells),
                               pairwise_label_agreement = round(pairw, 4),
                               mean_modal_freq = round(mean(modal_f), 4),
                               mean_agree_with_consensus = round(mean(agree_cons, na.rm = TRUE), 4),
                               frac_mc_label_matches_modal = round(mean(mc_tab$labels_agree), 4),
                               median_mc_label_support = round(median(mc_tab$label_support), 4))
  msg("%s %-3s | pairwise %.3f | modal %.3f | vs-consensus %.3f | mc label==modal %d/%d | %.1f min",
      GEN, st, pairw, mean(modal_f), mean(agree_cons, na.rm = TRUE),
      sum(mc_tab$labels_agree), nrow(mc_tab),
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
out <- do.call(rbind, sum_rows)
write.table(out, file.path(OUTD, paste0(GEN, ".seed_label_sweep.summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("wrote %s", file.path(OUTD, paste0(GEN, ".seed_label_sweep.summary.tsv")))
