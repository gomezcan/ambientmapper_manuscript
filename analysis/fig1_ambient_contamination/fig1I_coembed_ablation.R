#!/usr/bin/env Rscript
# fig1I_coembed_ablation.R -- Figure 1, panel I: oracle feature ablation of the Fig 1F co-embedding
#   (the negative control for the cross-species co-projection). Re-embeds the Fig 1F barcode set with
#   each barcode given ONLY the bins of the genome its plate index says it came from, then re-scores
#   the species-mixing statistic. Writes the re-embedded object, per-nucleus metadata, the mixing
#   summary, the per-plate embeddability ladder and diagnostic UMAP renders (the panel is drawn from these).
# Inputs : data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_indep/coembed/SM2v2_coembed_Pre.coembed.meta_full.tsv
#          data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_indep/step0_qc/SM2_{B73v5,TAIR10}.raw.soc.rds
#          (optional) .../coembed/Pre_cluster/SM2v2_coembed_Pre.updated_metadata_v7.<cfg>.txt as the unablated comparator
# Run    : fig1I_coembed_ablation.sh (HPC: Socrates + Seurat + FNN). Usage line below; every argument is set by the wrapper.
#
# DESIGN NOTES
# The question. Fig 1F shows maize and Arabidopsis nuclei co-projecting in one embedding. They could
# co-project because their cell identities are similar across species (the library is maize
# above-ground organs against Arabidopsis whole seedlings, so above-ground identities genuinely
# overlap), not because of cross-genome (ambient / mis-assigned) signal.
# The test. Take the SAME barcodes Fig 1F is drawn from, give each one ONLY the bins of the genome its
# PLATE INDEX says it came from, and re-embed with the IDENTICAL frozen configuration. Species that
# SEPARATE after ablation mean the co-projection was carried by cross-genome signal, not by shared
# cell identity, depth structure, QC covariates or embedding geometry. Species that STILL MIX mean
# the identity alternative survives. The test is falsifiable.
# Why it is not circular. AmbientMapper appears nowhere in this pipeline: the inputs are the RAW
# (pre-QC, PreClean) per-genome count matrices, the mask is the physical plate index alone, and no
# decontamination output is read. The plate index is a trustworthy oracle: measured cross-plate
# barcode swapping is 0.04%.
# Why the ablation is done on the RAW objects, not inside the co-embed matrix. The co-embed build
# assembles the matrix as the UNION of the two per-genome QC-passed cell sets and zero-fills the
# genome a barcode did not pass QC on, so most Arabidopsis-plate barcodes carry an all-zero
# Arabidopsis block there as an artifact of the build. This script reads each barcode's own-genome
# bins from the raw step0_qc objects, which contain every barcode's real signal. Cross-genome bins
# are simply never included; that absence IS the ablation.
# What the raw data contain. On raw, unfiltered data an Arabidopsis-plate barcode from the Fig 1F
# set carries a median of 16 Arabidopsis reads in 16 bins, against 348 reads in 348 bins for a
# maize-plate barcode on maize. The ablated Arabidopsis population is small because those barcodes
# contain almost no Arabidopsis chromatin, not because a filter removed them. The ladder file
# reports this per plate.
# The one logical gap. The ablation removes ALL cross-genome signal at once: ambient chromatin from
# the shared pool AND conserved-sequence mis-mapping of a nucleus's own chromatin. On its own it
# cannot say which carried the co-projection; the discriminating evidence is the plate asymmetry in
# Fig 1 (per-read mis-mapping would be roughly symmetric, the observed off-target fraction is not).
# State the ablation and the asymmetry together.
# The embeddability gate, the ONLY cell filter, explicit and reported: keep a barcode iff it has
# >= min_bins non-zero own-genome bins (default 50, the pipeline's min.c convention). Everything it
# excludes is counted per plate in the ladder file. cleanData is called with min.c = min_bins, the
# same gate restated in the pipeline's idiom, never a stricter one; anything it still removes is
# reported per plate. The mixing statistic is computed on the EMBEDDED set from $UMAP, not on
# $Clusters: callClusters(m.clst) discards nuclei in small clusters, so clustering is run only to
# attach a label for the figure and never subsets anything.
# Held identical to the Fig 1F run (frozen for cross-object comparability, not re-optimised):
#   cleanData(min.c, min.t = 0.001, max.t = 0) -> tfidf(doL2 = TRUE)
#   -> reduceDims(SVD, n.pcs = pcs, cor.max = 0.6, num.var = 50% of tiles, scaleVar, doL2)
#   -> projectUMAP(k.near, m.dist) -> Socrates::callClusters(res, k.near, cl.method = 4, e.thresh = 3,
#      threshold = 3, m.clst);   pcs = 20, k_near = 30, min_dist = 0.3, res = 0.5, m.clst = 50 (min.c = 50)
#   These are ARGUMENTS here, not literals, so the caller states them (see COEMB_CFG in fig1.R).
#   min.t = 0.001 is the same parameter as Fig 1F but is applied over a different number of barcodes
#   because the gate changes the cell set; that is inherent to the design.
# Two different k values, never conflate them: k_near = 30 is the UMAP neighbourhood that BUILDS the
# embedding; mix_k = 15 is the MIXING-STATISTIC neighbourhood (the number quoted in the manuscript).
# The mixing statistic (identical to the grid-scan and Fig 5 definitions):
#   raw      = mean over nuclei of the fraction of its mix_k nearest neighbours carrying the OTHER
#              species label (label = plate of origin)
#   expected = 2p(1-p), p = minority plate fraction   (random-intermingling baseline)
#   obs/exp  = raw / expected; 1.0 = indistinguishable from random intermingling.
#   Report obs/exp, never raw: the gate changes which barcodes are embedded and hence p, and p moves
#   the null mechanically.
# "LouvainClusters" is Leiden: Socrates::callClusters(cl.method = 4) forwards to
# Seurat::FindClusters(algorithm = 4). The column name is the pipeline's historical name.
#
# Usage:
#   Rscript fig1I_coembed_ablation.R <meta_full_tsv> <raw_b73v5_rds> <raw_tair10_rds> \
#           <outdir> <prefix> <pcs> <k_near> <min_dist> <min_bins> <resolution> \
#           [seed=1] [m_clst=50] [mix_k=15] [own_min=200] [run_grid=0] [ref_v7_meta=NA]
#     meta_full_tsv  the Fig 1F barcode set + `Genome` oracle (SM2v2_coembed_Pre.coembed.meta_full.tsv).
#                    The co-embed .soc.rds is deliberately NOT an input (see the design notes).
#     raw_*_rds      step0_qc/<SM2_B73v5|SM2_TAIR10>.raw.soc.rds, pre-QC counts.
#     min_bins       embeddability gate, non-zero own-genome bins (default 50).
#     mix_k          neighbourhood for the mixing statistic (15).
#     own_min        own-genome READ floor for the second reported scope (200).
#     run_grid       1 -> also run the configuration grid for a Fig 5D-style band.
#     ref_v7_meta    optional unablated Fig 1F v7 metadata, scored through the SAME code path so both
#                    numbers in the panel share one definition.
#
# Outputs (all under <outdir>, created here):
#   <prefix>.ablated.SocObj_v7<tag>.rds              re-embedded Socrates object
#   <prefix>.ablated.updated_metadata_v7<tag>.txt    per-nucleus metadata + mixing
#   <prefix>.ablated.reduced_dimensions_v7<tag>.txt  SVD coordinates
#   <prefix>.ablated.mixing_summary<tag>.tsv         headline obs/exp numbers
#   <prefix>.ablated.embeddability_ladder<tag>.tsv   the per-plate attrition ladder
#   <prefix>.ablated.grid_metrics<tag>.tsv           only when run_grid = 1
#   plots/<prefix>.ablated.umap<tag>_{Genome,LouvainClusters}.png
# The metadata and reduced-dimension files follow the pipeline convention (row.names = TRUE), so the
# header carries ONE FEWER field than the data rows. Read them with a bare data.table::fread(f) or
# base read.table, never fread(f, header = TRUE), and select columns BY NAME.

suppressMessages(library(Socrates))
suppressMessages(library(igraph))
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(FNN))
suppressMessages(library(DelayedArray))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: Rscript fig1I_coembed_ablation.R <meta_full_tsv> <raw_b73v5_rds> <raw_tair10_rds> <outdir> <prefix> <pcs> <k_near> <min_dist> <min_bins> <resolution> [seed=1] [m_clst=50] [mix_k=15] [own_min=200] [run_grid=0] [ref_v7_meta=NA]")
}
meta_path  <- args[1]
raw_b73    <- args[2]
raw_at     <- args[3]
outdir     <- args[4]
out_prefix <- args[5]
pcs        <- as.integer(args[6])
k_near     <- as.integer(args[7])
min_dis    <- as.numeric(args[8])
min_bins   <- as.integer(args[9])                                    # embeddability gate
resolution <- as.numeric(args[10])
seed       <- if (length(args) >= 11) as.integer(args[11]) else 1L
m_clst     <- if (length(args) >= 12) as.integer(args[12]) else 50L  # labels only, never a filter
mix_k      <- if (length(args) >= 13) as.integer(args[13]) else 15L
own_min    <- if (length(args) >= 14) as.numeric(args[14]) else 200
run_grid   <- if (length(args) >= 15) as.integer(args[15]) else 0L
ref_meta_p <- if (length(args) >= 16 && !args[16] %in% c("", "NA", "none")) args[16] else NA_character_

## ---------------------------------------------------------------------------
## THE ORACLE MAP.  Plate of origin -> the ONE reference genome whose bins a barcode
## from that plate is allowed to keep. Bins of the other genome are never read.
## The genome tags are the feature-name prefixes the co-embed build (3_1_coembed_build.R) stamps on,
## reused here so the joint feature space is named identically to Fig 1F's.
## ---------------------------------------------------------------------------
PLATE_TO_GENOME <- c(At = "TAIR10", B73 = "B73v5")
GENOME_TAGS     <- c("B73v5", "TAIR10")
K_PC            <- 30L    # cross-check mixing k in SVD space (K_PC of the Fig 5 co-embed script)

## Plate of origin from the cellID suffix. Verbatim from 3_1_coembed_build.R.
plate_of <- function(k) ifelse(grepl("-SM2_At$",  k), "At",
                        ifelse(grepl("-SM2_B73$", k), "B73", NA_character_))

## Raw-object column key. Verbatim from 3_1_coembed_build.R. Turns
##   <barcode>-SM2_<plate>_<Zm_B73v5|AraTAIR10>_scifiATAC   (raw colname)
## into
##   <barcode>-SM2_<plate>                                  (== the Fig 1F cellID)
## WHY NOT THE BARE 26 bp BARCODE PREFIX. Both keys resolve the same barcodes here (8,602 of 8,852
## At, 16,670 of 16,672 B73), but the 26 bp prefix DISCARDS THE PLATE, so it would permit an At-plate
## barcode to join a B73-plate raw column if a prefix were ever shared. This key keeps the plate and
## is therefore collision-proof by construction, and it is an exact identity match on cellID.
strip_genome <- function(x) sub("_(Zm_B73v5|AraTAIR10)_scifiATAC$", "", x)

# ---- TFIDF (identical to the pipeline's clustering scripts) ------------------
tfidf <- function(obj, frequencies = TRUE, log_scale_tf = TRUE, scale_factor = 10000,
                  doL2 = FALSE, slotName = "residuals") {
  bmat <- obj$counts
  .safe_tfidf <- function(tf, idf, block_size = 2000e6) {
    tryCatch({ tf * idf }, error = function(e) {
      options(DelayedArray.block.size = block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      tf <- DelayedArray(tf); idf <- as.matrix(idf); tf * idf
    })
  }
  if (frequencies) tf <- t(t(bmat) / Matrix::colSums(bmat)) else tf <- bmat
  if (log_scale_tf) tf@x <- log1p(tf@x * (if (frequencies) scale_factor else 1))
  idf <- log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  tf_idf_counts <- .safe_tfidf(tf, idf)
  if (doL2) {
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x = 1/colNorm)
  }
  rownames(tf_idf_counts) <- rownames(bmat); colnames(tf_idf_counts) <- colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse = TRUE); obj$norm_method <- "tfidf"
  obj
}

# ---- the mixing statistic ---------------------------------------------------
## raw = mean fraction of a nucleus's k nearest neighbours carrying the OTHER plate
##       label;  expected = 2p(1-p) at minority fraction p;  obs/exp = raw / expected.
## Same construction as the co-embed grid scan (raw) + the per-genome grid scan (the 2p(1-p)
## normalisation) + the Fig 5 co-embed script.
mix_stat <- function(emb, lab, k) {
  emb <- as.matrix(emb)
  ok  <- stats::complete.cases(emb) & !is.na(lab)
  out <- list(n = sum(ok), n_a = NA_integer_, n_b = NA_integer_, p_minor = NA_real_,
              raw = NA_real_, expected = NA_real_, obs_exp = NA_real_,
              per_cell = rep(NA_real_, length(lab)))
  if (sum(ok) <= k || length(unique(lab[ok])) < 2L) return(out)
  e  <- emb[ok, , drop = FALSE]
  l  <- as.character(lab[ok])
  nn <- FNN::get.knn(e, k = k)$nn.index
  fr <- vapply(seq_len(nrow(nn)), function(i) mean(l[nn[i, ]] != l[i]), numeric(1))
  tb <- table(l)
  p  <- as.numeric(min(tb)) / length(l)
  ex <- 2 * p * (1 - p)
  out$n_a      <- as.integer(sum(l == "At"))
  out$n_b      <- as.integer(sum(l == "B73"))
  out$p_minor  <- p
  out$raw      <- mean(fr)
  out$expected <- ex
  out$obs_exp  <- if (is.finite(ex) && ex > 0) mean(fr) / ex else NA_real_
  out$per_cell[ok] <- fr
  out
}

mix_row <- function(object, space, k, scope, neighbour_set, m, note = "") {
  data.frame(object = object, space = space, k = k, scope = scope,
             neighbour_set = neighbour_set, n_cells = m$n, n_At = m$n_a, n_B73 = m$n_b,
             p_minority = m$p_minor, raw = m$raw, expected = m$expected,
             obs_exp = m$obs_exp, note = note, stringsAsFactors = FALSE)
}

# ---- compact UMAP plotting (diagnostic renders; the panel is assembled from these) ----
plot_umap <- function(df, colour, file, title, palette = NULL, legend = TRUE) {
  p <- ggplot(df, aes(x = umap1, y = umap2, colour = .data[[colour]])) +
    geom_point(size = 0.2, alpha = 0.5) +
    labs(title = title, x = "umap1", y = "umap2", colour = colour) +
    theme_bw(base_size = 10)
  if (!is.null(palette)) p <- p + scale_colour_manual(values = palette, na.value = "grey70")
  if (!legend) p <- p + theme(legend.position = "none")
  suppressMessages(ggsave(file, p, width = 6, height = 5, dpi = 300, device = "png"))
  invisible(p)
}

###############################################################################
## (1) THE Fig 1F BARCODE SET + THE PLATE ORACLE
###############################################################################
message("=== oracle feature ablation, rebuilt from the RAW per-genome objects ===")
message(" - Fig 1F barcode set : ", meta_path)
message(" - raw B73v5          : ", raw_b73)
message(" - raw TAIR10         : ", raw_at)

meta.data <- read.table(meta_path, header = TRUE, sep = "\t")
meta.data <- meta.data[, !colnames(meta.data) == "...1", drop = FALSE]
if (!"cellID" %in% colnames(meta.data)) stop("metadata has no cellID column: ", meta_path)
rownames(meta.data) <- meta.data$cellID
fig1f_cells <- as.character(meta.data$cellID)

tag <- paste0(".pcs_", pcs, ".k_near_", k_near, ".min_dis_", min_dis,
              ".minbins_", min_bins, ".res_", resolution,
              if (m_clst != 50L) paste0(".mclst_", m_clst) else "")
out_dir_plots <- file.path(outdir, "plots")
dir.create(outdir,        showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plots, showWarnings = FALSE, recursive = TRUE)

## Plate of origin from two independent sources, cross-checked. HARD FAIL on any
## disagreement or NA: if the oracle is not trustworthy there is no point ablating.
##   SOURCE 1 the `Genome` column, set by the co-embed build and preserved by the clustering step;
##   SOURCE 2 an independent re-derivation from the anchored cellID suffix.
if (!"Genome" %in% colnames(meta.data)) stop("metadata has no `Genome` column -- cannot resolve plate of origin")
plate_meta <- as.character(meta.data$Genome)
plate_id   <- plate_of(fig1f_cells)
if (any(is.na(plate_meta) | is.na(plate_id)))
  stop("plate of origin is NA for ", sum(is.na(plate_meta) | is.na(plate_id)),
       " barcodes. The design oracle must cover every barcode.")
if (any(plate_meta != plate_id))
  stop("`Genome` metadata and cellID plate suffix DISAGREE for ", sum(plate_meta != plate_id),
       " barcodes. Refusing to ablate on an inconsistent label.")
plate <- plate_meta
names(plate) <- fig1f_cells
if (!all(plate %in% names(PLATE_TO_GENOME)))
  stop("unexpected plate label(s): ", paste(setdiff(unique(plate), names(PLATE_TO_GENOME)), collapse = ", "))
message(" - plate oracle: metadata `Genome` and cellID suffix AGREE on all ",
        length(fig1f_cells), " barcodes  (",
        paste(names(table(plate)), table(plate), sep = "=", collapse = "  "), ")")

###############################################################################
## (2) THE ABLATION  --  the scientific core of this script
###############################################################################
## For each plate, pull ONLY that plate's barcodes and ONLY its own genome's bins out
## of the RAW (pre-QC) object. Cross-genome bins are never read, and that absence is
## the ablation -- there is no masking step, nothing to leak, and nothing that depends
## on the co-embed build. Feature names are then genome-prefixed exactly as the co-embed build
## does, so the joint space is named identically to Fig 1F's and the two blocks cannot
## collide on shared tile coordinates.
##
## The raw objects contain BOTH plates' barcodes (the TAIR10 raw object is
##   129,135 At + 160,414 B73 columns), so restricting to this plate's barcodes is
##   essential and is done by an exact cellID join, not by trusting the file name.
load_own_block <- function(rds, gtag, want_cells) {
  o <- readRDS(rds)
  m <- o$counts
  if (is.null(m)) stop("raw object has no $counts: ", rds)
  colnames(m) <- strip_genome(colnames(m))
  if (any(duplicated(colnames(m))))
    stop("raw object has duplicate barcode keys after strip_genome: ", rds)
  hit  <- want_cells[want_cells %in% colnames(m)]
  miss <- setdiff(want_cells, colnames(m))
  m <- m[, hit, drop = FALSE]
  rownames(m) <- paste0(gtag, ":", rownames(m))        # genome-prefixed, as the co-embed build does
  m <- m[Matrix::rowSums(m) > 0, , drop = FALSE]       # bins unused by these barcodes
  message(sprintf("   %-6s raw: matched %d of %d requested barcodes (%d unmatched) -> %d non-empty bins",
                  gtag, length(hit), length(want_cells), length(miss), nrow(m)))
  rm(o); invisible(gc(FALSE))
  list(m = m, matched = hit, missing = miss)
}

message(" - reading own-genome bins from the RAW objects (pre-QC, no filtering)")
B <- load_own_block(raw_b73, "B73v5",  fig1f_cells[plate == "B73"])
A <- load_own_block(raw_at,  "TAIR10", fig1f_cells[plate == "At"])

## Unmatched barcodes are REPORTED, never silently dropped. A Fig 1F barcode absent
## from its own genome's raw object cannot be given own-genome signal at all.
unmatched <- c(B$missing, A$missing)
if (length(unmatched) > 0) {
  message(" - NOTE: ", length(unmatched), " Fig 1F barcode(s) absent from their own genome's raw object:")
  for (pl in names(PLATE_TO_GENOME))
    message(sprintf("     %-3s plate: %d of %d", pl,
                    sum(plate[unmatched] == pl), sum(plate == pl)))
  message("   They cannot receive own-genome signal and are excluded before embedding.")
}

## Assemble the joint ablated matrix: [B73v5 bins ; TAIR10 bins] x (all matched cells),
## each barcode zero-filled on the genome it did not come from. Unlike the co-embed
## build, that zero block is now a DELIBERATE ABLATION rather than a QC artifact.
matched_cells <- c(B$matched, A$matched)
expand <- function(m, cells) {
  miss <- setdiff(cells, colnames(m))
  if (length(miss)) {
    z <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                      dims = c(nrow(m), length(miss)), dimnames = list(rownames(m), miss))
    m <- cbind(m, z)
  }
  m[, cells, drop = FALSE]
}
abl <- rbind(expand(B$m, matched_cells), expand(A$m, matched_cells))
feat_genome <- sub(":.*$", "", rownames(abl))
if (!all(feat_genome %in% GENOME_TAGS))
  stop("unexpected feature prefixes: ", paste(head(setdiff(feat_genome, GENOME_TAGS), 5), collapse = ", "))
plate_m <- plate[matched_cells]
rm(B, A); invisible(gc(FALSE))
message(" - joint ablated matrix: ", nrow(abl), " bins x ", ncol(abl), " barcodes")

## own-genome depth and complexity, on RAW data with no QC applied
own_reads <- Matrix::colSums(abl)
own_bins  <- Matrix::colSums(abl > 0)

## VERIFICATION of the ablation: no barcode may carry a single count on the other
## genome's block. Cheap, and it is the assertion a reviewer will look for.
leak_at  <- sum(abl[feat_genome == "B73v5",  plate_m == "At",  drop = FALSE])
leak_b73 <- sum(abl[feat_genome == "TAIR10", plate_m == "B73", drop = FALSE])
if (leak_at != 0 || leak_b73 != 0)
  stop("ABLATION FAILED: residual cross-genome counts (At on B73v5 = ", leak_at,
       ", B73 on TAIR10 = ", leak_b73, ")")
message(" - ablation verified: 0 cross-genome counts present by construction")
for (pl in names(PLATE_TO_GENOME)) {
  i <- plate_m == pl
  message(sprintf("   %-3s plate: n=%6d | median own reads %5.0f | median own bins %5.0f | max bins %6.0f",
                  pl, sum(i), median(own_reads[i]), median(own_bins[i]), max(own_bins[i])))
}

###############################################################################
## (3) THE EMBEDDABILITY GATE  --  the only cell filter, explicit and counted
###############################################################################
embeddable <- own_bins >= min_bins
message(" - embeddability gate: >= ", min_bins, " non-zero own-genome bins")
for (pl in names(PLATE_TO_GENOME)) {
  i <- plate_m == pl
  message(sprintf("   %-3s plate: %d of %d pass (%.1f%%), %d excluded",
                  pl, sum(i & embeddable), sum(i), 100 * sum(i & embeddable) / sum(i),
                  sum(i & !embeddable)))
}
if (sum(embeddable) <= mix_k)
  stop("only ", sum(embeddable), " barcodes pass the embeddability gate -- nothing to embed")

## the ladder, per plate, as a first-class result
ladder <- do.call(rbind, lapply(names(PLATE_TO_GENOME), function(pl) {
  in_f <- sum(plate == pl)
  i    <- plate_m == pl
  data.frame(plate = pl, kept_genome = unname(PLATE_TO_GENOME[[pl]]),
             step = c("in_Fig1F_set", "matched_in_raw", "ge_1_bin", "ge_20_bins",
                      paste0("ge_", min_bins, "_bins_EMBEDDED"),
                      paste0("ge_", own_min, "_reads")),
             n = c(in_f, sum(i), sum(i & own_bins >= 1), sum(i & own_bins >= 20),
                   sum(i & embeddable), sum(i & own_reads >= own_min)),
             pct_of_Fig1F = round(100 * c(in_f, sum(i), sum(i & own_bins >= 1),
                                          sum(i & own_bins >= 20), sum(i & embeddable),
                                          sum(i & own_reads >= own_min)) / in_f, 2),
             median_own_reads = median(own_reads[i]),
             median_own_bins  = median(own_bins[i]),
             stringsAsFactors = FALSE)
}))

abl_obj <- list(counts = abl[, embeddable, drop = FALSE],
                meta   = meta.data[matched_cells[embeddable], , drop = FALSE])
abl_obj$counts <- abl_obj$counts[Matrix::rowSums(abl_obj$counts) > 0, , drop = FALSE]
rm(abl); invisible(gc(FALSE))
message(" - into the frozen pipeline: ", nrow(abl_obj$counts), " bins x ",
        ncol(abl_obj$counts), " barcodes")

###############################################################################
## (4) RE-EMBED WITH THE FROZEN CONFIGURATION
###############################################################################
## min.c is set to min_bins -- the SAME gate already applied above, restated in the
## pipeline's idiom, NOT a second and stricter filter. Anything it still removes is
## reported per plate immediately below rather than absorbed silently.
message(" - cleanData min.c = ", min_bins, " (= the embeddability gate) | min.t = 0.001 | max.t = 0")
set.seed(seed)
soc.obj <- cleanData(abl_obj, min.c = min_bins, min.t = 0.001, max.t = 0, verbose = TRUE)
lost <- setdiff(colnames(abl_obj$counts), colnames(soc.obj$counts))
if (length(lost) > 0) {
  message(" - NOTE: cleanData removed ", length(lost), " barcode(s) BEYOND the gate ",
          "(min.t feature filtering can push a barcode below min.c). Reported, not hidden:")
  for (pl in names(PLATE_TO_GENOME))
    message(sprintf("     %-3s plate: %d", pl, sum(plate[lost] == pl)))
} else {
  message(" - cleanData removed no barcode beyond the gate")
}

soc.obj <- tfidf(soc.obj, doL2 = TRUE)
number.sites <- ceiling(nrow(soc.obj$counts) * 0.5)
soc.obj <- reduceDims(soc.obj, method = "SVD", n.pcs = pcs, cor.max = 0.6, num.var = number.sites,
                      verbose = TRUE, scaleVar = TRUE, doSTD = FALSE, doL1 = FALSE, doL2 = TRUE,
                      refit_residuals = FALSE)
if (ncol(soc.obj$PCA) < pcs)
  message(" - NOTE: reduceDims returned ", ncol(soc.obj$PCA), " of ", pcs,
          " requested components (cor.max = 0.6 dropped depth-correlated ones).")
soc.obj <- projectUMAP(soc.obj, verbose = TRUE, k.near = k_near, m.dist = min_dis)

## THE EMBEDDED SET. This, not $Clusters, is what the mixing statistic is computed
## on. Nothing below is allowed to shrink it.
emb_cells <- rownames(soc.obj$UMAP)
umap_xy   <- as.matrix(soc.obj$UMAP[emb_cells, c("umap1", "umap2"), drop = FALSE])
pca_xy    <- soc.obj$PCA[emb_cells, , drop = FALSE]
lab       <- unname(plate[emb_cells])
own_r     <- own_reads[emb_cells]
own_b     <- own_bins[emb_cells]
message(" - EMBEDDED: ", length(emb_cells), " nuclei (At ", sum(lab == "At"),
        ", B73 ", sum(lab == "B73"), ")")

## clusters: LABELS ONLY. m.clst discards small clusters, which on a minority plate of
## this size can remove most of it -- so callClusters must never gate the statistic,
## and a failure here must not kill the run.
clu <- NULL
ok_clu <- tryCatch({
  soc.obj <- Socrates::callClusters(soc.obj, res = resolution, k.near = k_near, verbose = TRUE,
                                    cleanCluster = FALSE, cl.method = 4, e.thresh = 3,
                                    threshold = 3, m.clst = m_clst)
  clu <- soc.obj$Clusters; TRUE
}, error = function(e) { message(" - callClusters SKIPPED (labels only): ", conditionMessage(e)); FALSE })
if (ok_clu) {
  message(" - ", length(unique(clu$LouvainClusters)), " clusters at res=", resolution,
          " m.clst=", m_clst, " | ", nrow(clu), " of ", length(emb_cells),
          " embedded nuclei carry a cluster label")
  message("   (the ", length(emb_cells) - nrow(clu), " unlabelled nuclei REMAIN in the ",
          "mixing statistic -- m.clst is a labelling choice, not a filter)")
}
cluster_lab <- rep(NA_character_, length(emb_cells))
if (ok_clu) cluster_lab[match(rownames(clu), emb_cells)] <- as.character(clu$LouvainClusters)

###############################################################################
## (5) THE MIXING STATISTIC, on the EMBEDDED set, twice
###############################################################################
in_scope <- own_r >= own_min
message(sprintf(" - mixing scopes: embedded = %d | own reads >= %g = %d (At %d, B73 %d)",
                length(emb_cells), own_min, sum(in_scope),
                sum(in_scope & lab == "At"), sum(in_scope & lab == "B73")))
if (min(sum(lab == "At"), sum(lab == "B73")) < 100)
  message(" - NOTE: one plate has <100 embedded nuclei; obs/exp rests on few nuclei and ",
          "must be labelled low-confidence.")

m_all_umap <- mix_stat(umap_xy, lab, mix_k)
m_sub_umap <- mix_stat(umap_xy[in_scope, , drop = FALSE], lab[in_scope], mix_k)
m_foc_umap <- m_all_umap
if (any(in_scope)) {
  fc <- m_all_umap$per_cell; fc[!in_scope] <- NA_real_
  ex <- m_all_umap$expected
  m_foc_umap <- list(n = sum(in_scope & !is.na(fc)),
                     n_a = sum(in_scope & lab == "At"), n_b = sum(in_scope & lab == "B73"),
                     p_minor = m_all_umap$p_minor, raw = mean(fc, na.rm = TRUE), expected = ex,
                     obs_exp = if (is.finite(ex) && ex > 0) mean(fc, na.rm = TRUE) / ex else NA_real_,
                     per_cell = fc)
}
m_all_pc <- mix_stat(pca_xy, lab, K_PC)
m_sub_pc <- mix_stat(pca_xy[in_scope, , drop = FALSE], lab[in_scope], K_PC)

scope_all <- paste0("embedded_ge_", min_bins, "_bins")
scope_sub <- paste0("embedded_and_ge_", own_min, "_reads")
summ <- rbind(
  mix_row("ablated", "umap", mix_k, scope_all, "within_scope", m_all_umap,
          "PRIMARY -- all embedded nuclei; NOT gated by callClusters"),
  mix_row("ablated", "umap", mix_k, scope_sub, "within_scope", m_sub_umap,
          "sparse-tail control: neighbourhoods among in-scope nuclei only"),
  mix_row("ablated", "umap", mix_k, scope_sub, "all_embedded", m_foc_umap,
          "sparse-tail control: focal nuclei in scope, neighbours from all embedded"),
  mix_row("ablated", "pc",   K_PC,  scope_all, "within_scope", m_all_pc,
          "cross-check in SVD space"),
  mix_row("ablated", "pc",   K_PC,  scope_sub, "within_scope", m_sub_pc,
          "cross-check in SVD space")
)

###############################################################################
## (5b) THE UNABLATED Fig 1F COMPARATOR, scored by the SAME code path
###############################################################################
if (!is.na(ref_meta_p)) {
  if (!file.exists(ref_meta_p)) stop("ref_v7_meta not found: ", ref_meta_p)
  message(" - reference (unablated) embedding: ", ref_meta_p)
  ref <- as.data.frame(data.table::fread(ref_meta_p))          # bare fread: row-name offset
  if (!all(c("umap1", "umap2", "Genome") %in% colnames(ref)))
    stop("ref_v7_meta lacks umap1/umap2/Genome")
  ref_id <- if ("cellID" %in% colnames(ref)) as.character(ref$cellID) else as.character(ref[[1]])
  ref_xy <- as.matrix(ref[, c("umap1", "umap2")])
  ref_lb <- as.character(ref$Genome)
  m_ref_all <- mix_stat(ref_xy, ref_lb, mix_k)
  shared    <- ref_id %in% emb_cells
  m_ref_sh  <- mix_stat(ref_xy[shared, , drop = FALSE], ref_lb[shared], mix_k)
  summ <- rbind(summ,
    mix_row("reference_unablated", "umap", mix_k, "all_Fig1F", "within_scope", m_ref_all,
            "Fig 1F embedding as published"),
    mix_row("reference_unablated", "umap", mix_k, "shared_with_ablated", "within_scope", m_ref_sh,
            "same nuclei as the ablated readout -- the like-for-like contrast"))
  message(sprintf(" - reference obs/exp (all Fig 1F) = %.3f | ablated obs/exp (embedded) = %.3f",
                  m_ref_all$obs_exp, m_all_umap$obs_exp))
}

cat("\n=== embeddability ladder (raw pre-QC own-genome signal, no filtering upstream) ===\n")
print(ladder, row.names = FALSE)
cat("\n=== mixing summary (report obs/exp, NEVER raw; 1.0 = random intermingling) ===\n")
print(format(summ, digits = 4), row.names = FALSE)
cat("\nNOTE: Report the ladder WITH the obs/exp number. The ablated Arabidopsis population\n",
    "  is small because those barcodes contain almost no Arabidopsis chromatin (see the\n",
    "  median own reads/bins above, measured on RAW pre-QC data), not because a filter\n",
    "  removed them and not because of the co-embed zero-fill.\n",
    "NOTE: This ablation removes ALL cross-genome signal, ambient chromatin and\n",
    "  conserved-sequence mis-mapping alike, so on its own it cannot say which of the\n",
    "  two carried the co-projection. Pair it with the Fig 1 plate asymmetry.\n", sep = "")

###############################################################################
## (6) OPTIONAL: the full configuration grid, for a Fig 5D-style band
###############################################################################
grid_df <- NULL
if (run_grid == 1L) {
  message("\n - run_grid = 1: repeating the co-embed grid on the ablated matrix")
  ## Structure from the co-embed grid scan: TF-IDF + SVD fixed ONCE at max(pcs_grid)+1, each
  ## combination subsets the PCs. Grid points are therefore NOT independent replicates and must
  ## never be tested against each other.
  pcs_grid      <- c(20, 25, 30, 40, 50)
  k_grid        <- c(15, 20, 30)
  min_dist_grid <- c(0.05, 0.15, 0.30)
  gsoc <- reduceDims(soc.obj, method = "SVD", n.pcs = max(pcs_grid) + 1, cor.max = 0.6,
                     num.var = number.sites, verbose = TRUE, scaleVar = TRUE, doSTD = FALSE,
                     doL1 = FALSE, doL2 = TRUE, refit_residuals = FALSE)
  common    <- intersect(rownames(gsoc$meta), rownames(gsoc$PCA))
  gsoc$meta <- gsoc$meta[common, , drop = FALSE]
  pca_full  <- gsoc$PCA[common, , drop = FALSE]
  avail_pcs <- ncol(pca_full)
  pcs_grid  <- pcs_grid[pcs_grid <= avail_pcs]
  if (length(pcs_grid) == 0) pcs_grid <- avail_pcs
  message("   available PCs after reduceDims: ", avail_pcs,
          " -> scanning pcs = ", paste(pcs_grid, collapse = ", "))
  if (mix_k > min(k_grid))
    message("   NOTE: mix_k (", mix_k, ") exceeds the smallest UMAP k_near (", min(k_grid),
            "); the grid scan sets diag_k = min(k_near). Reporting anyway, flagged.")

  grid <- expand.grid(pcs = pcs_grid, k_near = k_grid, min_dist = min_dist_grid,
                      stringsAsFactors = FALSE)
  glab <- unname(plate[common])
  gsc  <- own_reads[common] >= own_min
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    gi <- gsoc
    gi$PCA <- pca_full[, seq_len(grid$pcs[i]), drop = FALSE]
    gi <- projectUMAP(gi, m.dist = grid$min_dist[i], k.near = grid$k_near[i], metric = "cosine",
                      svd_slotName = "PCA", umap_slotName = "UMAP", verbose = FALSE, seed = seed)
    u  <- as.matrix(gi$UMAP[rownames(gi$meta), c("umap1", "umap2"), drop = FALSE])
    ma <- mix_stat(u, glab, mix_k)
    ms <- mix_stat(u[gsc, , drop = FALSE], glab[gsc], mix_k)
    rows[[i]] <- data.frame(pcs = grid$pcs[i], k_near = grid$k_near[i], min_dist = grid$min_dist[i],
                            mix_k = mix_k, n_cells = ma$n, p_minority = ma$p_minor,
                            raw_all = ma$raw, obs_exp_all = ma$obs_exp,
                            n_cells_scope = ms$n, raw_scope = ms$raw, obs_exp_scope = ms$obs_exp,
                            stringsAsFactors = FALSE)
    message(sprintf("   [%2d/%2d] pcs=%2d k=%2d min_dist=%.2f | obs/exp embedded=%.3f  own>=%g=%.3f",
                    i, nrow(grid), grid$pcs[i], grid$k_near[i], grid$min_dist[i],
                    ma$obs_exp, own_min, ms$obs_exp))
    rm(gi)
  }
  grid_df <- do.call(rbind, rows)
  message(sprintf("   grid obs/exp (embedded): %d configs, %.3f to %.3f",
                  nrow(grid_df), min(grid_df$obs_exp_all, na.rm = TRUE),
                  max(grid_df$obs_exp_all, na.rm = TRUE)))
  rm(gsoc, pca_full); invisible(gc(FALSE))
}

###############################################################################
## (7) WRITE EVERYTHING
###############################################################################
## per-nucleus metadata for EVERY EMBEDDED nucleus (not only the clustered ones)
per_nuc <- data.frame(
  cellID                  = emb_cells,
  Genome                  = lab,
  plate                   = lab,
  umap1                   = umap_xy[, 1],
  umap2                   = umap_xy[, 2],
  LouvainClusters         = cluster_lab,
  own_reads               = own_r,
  own_bins                = own_b,
  own_ge_min_reads        = in_scope,
  knn_other_frac_umap     = m_all_umap$per_cell,
  knn_other_frac_umap_scope = ifelse(in_scope, m_all_umap$per_cell, NA_real_),
  row.names               = emb_cells, stringsAsFactors = FALSE)
carry <- setdiff(colnames(meta.data), colnames(per_nuc))
if (length(carry)) per_nuc <- cbind(per_nuc, meta.data[emb_cells, carry, drop = FALSE])

f_rds  <- file.path(outdir, paste0(out_prefix, ".ablated.SocObj_v7",             tag, ".rds"))
f_meta <- file.path(outdir, paste0(out_prefix, ".ablated.updated_metadata_v7",   tag, ".txt"))
f_pca  <- file.path(outdir, paste0(out_prefix, ".ablated.reduced_dimensions_v7", tag, ".txt"))
f_mix  <- file.path(outdir, paste0(out_prefix, ".ablated.mixing_summary",        tag, ".tsv"))
f_lad  <- file.path(outdir, paste0(out_prefix, ".ablated.embeddability_ladder",  tag, ".tsv"))
f_grid <- file.path(outdir, paste0(out_prefix, ".ablated.grid_metrics",          tag, ".tsv"))

saveRDS(soc.obj, file = f_rds)
write.table(per_nuc, file = f_meta, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(pca_xy,  file = f_pca,  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(summ,    file = f_mix,  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(ladder,  file = f_lad,  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
if (!is.null(grid_df))
  write.table(grid_df, file = f_grid, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

pdf_df <- data.frame(umap1 = per_nuc$umap1, umap2 = per_nuc$umap2,
                     Genome = factor(lab, levels = c("At", "B73")),
                     LouvainClusters = factor(cluster_lab))
plot_umap(pdf_df, "Genome",
          file.path(out_dir_plots, paste0(out_prefix, ".ablated.umap", tag, "_Genome.png")),
          sprintf("oracle-ablated (raw-rebuilt) | obs/exp = %.3f (k=%d, n=%d)",
                  m_all_umap$obs_exp, mix_k, length(emb_cells)),
          # Fig 1 species palette: At blue, maize/B73 red (cols_species in analysis/_helpers/plotting.R).
          palette = c(At = "#377eb8", B73 = "#e41a1c"))
plot_umap(pdf_df, "LouvainClusters",
          file.path(out_dir_plots, paste0(out_prefix, ".ablated.umap", tag, "_LouvainClusters.png")),
          paste0("oracle-ablated (raw-rebuilt) | res=", resolution), legend = FALSE)

for (f in c(f_rds, f_meta, f_pca, f_mix, f_lad, if (!is.null(grid_df)) f_grid))
  if (!file.exists(f) || file.size(f) == 0) stop("FAILED to write: ", f)

message("\n - DONE. Wrote:")
for (f in c(f_rds, f_meta, f_pca, f_mix, f_lad, if (!is.null(grid_df)) f_grid))
  message("     ", f)
message(" - plots -> ", out_dir_plots)
