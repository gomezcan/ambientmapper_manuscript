#!/usr/bin/env Rscript
###############################################################################
## 4_3n_background_null_rZ.R  --  TEST A: covariate-corrected background-gene null for per-cell rZ.
##
## QUESTION. 4_3i gives each marker a cluster-mean rZ; 4_3k annotates a cluster by the type whose
## markers reach the highest rZ. Nothing in that chain says whether a peak rZ of 5.8 is remarkable.
## This script supplies the scale: push EVERY gene in the matrix through the identical rZ pipeline,
## then ask whether the panel's signal exceeds background (non-panel) genes at the SAME expression.
##
## TWO LEVELS, because they answer different questions:
##   (2) SET level  (cluster x type)  = HEADLINE. "Is cluster cl's TYPE call beyond chance?"
##       Observed = mean adjusted cluster-rZ of type T's markers in cl. Null = B draws of random
##       background sets of the SAME SIZE (a competitive gene-set test, and the unit the annotation
##       is actually made at). ~n_cl x n_type tests -> BH retains power.
##   (1) GENE level (per marker)      = detail. "Which individual markers carry the signal?"
##       Marker's adjusted spec vs the full background distribution, BH across markers.
##
## WHY THE COVARIATE CORRECTION (do not remove this).  rZ_geom = sqrt(max0(Zi)*max0(Zj)) and Zj is
## the cell's z ACROSS GENES -- i.e. largely an expression-level term, so rZ rises with expression
## for reasons that have nothing to do with cell-type identity. A naive null (random genes) is
## therefore dominated by low-expressed genes and would declare almost every marker significant.
##
##   EARLIER DESIGN, REJECTED: nearest-neighbour MATCHING (pool the n nearest background genes in
##   expression space). It fails silently in the sparse tails -- a synthetic decoy set placed at the
##   top of the expression range (where only ~1.6% of background genes live) had its "matched" pool
##   pulled down to lower-expressed genes and came out FALSELY SIGNIFICANT. Matching cannot match
##   what is not there.
##
##   ALSO REJECTED: standardising by a running median / MAD trend. rZ is ZERO-INFLATED by
##   construction -- sqrt(max0(Zi)*max0(Zj)) is exactly 0 whenever a gene is below its cell's
##   gene-wide mean, which holds in EVERY cell for lowly expressed genes (~60% of background genes
##   have peak rZ == 0), and the dependence on expression is a THRESHOLD, not a smooth trend. Both
##   location and scale are 0 across most of the covariate range, so the studentised statistic
##   degenerates (measured IQR 0.000).
##
##   ADOPTED: condition DISTRIBUTION-FREE. Express every gene as its MID-RANK QUANTILE u among
##   background genes in the same covariate neighbourhood. Ties (the zero block) share a mid-rank,
##   so u is uniform under the null BY CONSTRUCTION -- whatever the zero mass or the trend shape.
##   The background is then exchangeable everywhere, so the null draws uniformly from ALL background
##   genes, and the set-level p floor is 1/(1+B) rather than 1/(1+pool). Genes outside the background
##   expression range are flagged (`in_support`): no method can test what has no comparator.
##
## MODES (both run off one genome-wide pass):
##   mean     correct for expression level                     (PRIMARY)
##   mean_sd  correct for expression, then for cell-to-cell sd (CONSERVATIVE: isolates cluster-
##            STRUCTURED variance from plain noisiness; a call surviving this is the robust core)
##
## WHY THIS IS VALID ON A SMOOTHED MATRIX.  The unit of randomization is the GENE, not the cell.
## Every background gene is pushed through the exact same smoothed cell graph, the same clusters and
## the same cell counts as the markers it is compared against, so kNN-diffusion correlation is held
## fixed on both sides of the contrast and cannot inflate it. (A per-CELL test would be
## pseudoreplication -- that is the test we are deliberately NOT doing.) Scope is within-dataset:
## "beyond equally-expressed genes IN THIS OBJECT", not population inference -- one library, no
## biological replication. That is the correct scope for a cleaning-impact claim.
##
## REPORTED SENSITIVITIES (the primary call is mean + global BH; the others are printed alongside so
## a borderline result is visibly borderline rather than silently rounded to yes/no):
##   statistic   : mean u over all the type's markers   vs   mean of the top-k (default 5)
##   multiplicity: BH across all cluster x type pairs   vs   BH within each cluster
## Measured on SM2v2_plate: At is robust to both choices; maize is NOT (its best calls sit at
## p~1e-3 and flip on the multiplicity convention) -- which is itself the informative result.
##
## Usage:
##   Rscript 4_3n_background_null_rZ.R <out_dir> <prefix> <smoothed_rds> <meta_txt> <markers_bed> \
##           <stage> [cluster_col=LouvainClusters] [window=500] [mode=mean,mean_sd] \
##           [exclude=dividing] [metric=geom] [fdr=0.05] [B=10000] [min_markers=3] [seed=1] \
##           [cmp_tsv=] [topk=5]
##
## Outputs (<out_dir>/<prefix>...), genome-wide pass shared across modes:
##   .bgnull.genestats.tsv.gz       every gene: mean/sd expr, peak cluster, raw peak rZ, second, spec
##   .bgnull.<mode>.settest.tsv     HEADLINE: cluster x type, observed vs null, effect z, p/q
##   .bgnull.<mode>.annotation.tsv  per cluster: the type verdict (significant / ambiguous / none)
##   .bgnull.<mode>.markers.tsv     per marker: adjusted spec/peak vs background, p/q, in_support
##   .bgnull.<mode>.pdf             per-cluster type effect sizes + the trend that was removed
##   .bgnull.results.md             headline for every mode + reconciliation against 4_3i
###############################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("Usage: Rscript 4_3n_background_null_rZ.R <out_dir> <prefix> <smoothed_rds> <meta_txt> <markers_bed> <stage> [cluster_col] [window] [mode] [exclude] [metric] [fdr] [B] [min_markers] [seed] [cmp_tsv]")
out_dir      <- args[1]
prefix       <- args[2]
smoothed_rds <- args[3]
meta_txt     <- args[4]
markers_bed  <- args[5]
stage        <- args[6]
cluster_col  <- if (length(args) >=  7) args[7]              else "LouvainClusters"
window       <- if (length(args) >=  8) as.integer(args[8])  else 500L
modes        <- if (length(args) >=  9) strsplit(args[9], ",")[[1]] else c("mean", "mean_sd")
exclude      <- if (length(args) >= 10) args[10]             else "dividing"
metric       <- if (length(args) >= 11) args[11]             else "geom"
fdr_thr      <- if (length(args) >= 12) as.numeric(args[12]) else 0.05
B            <- if (length(args) >= 13) as.integer(args[13]) else 10000L
min_markers  <- if (length(args) >= 14) as.integer(args[14]) else 3L
seed         <- if (length(args) >= 15) as.integer(args[15]) else 1L
cmp_tsv      <- if (length(args) >= 16) args[16]             else ""
topk         <- if (length(args) >= 17) as.integer(args[17]) else 5L   # top-k sensitivity statistic

## A cluster x type row is UNTESTABLE if fewer than this fraction of its markers lie inside the
## background expression range -- outside it there is no valid comparator (see the set-test block).
min_support <- 0.5

if (!all(modes %in% c("mean", "mean_sd"))) stop("mode must be 'mean' and/or 'mean_sd'")
if (!metric %in% c("geom", "euclid"))      stop("metric must be 'geom' or 'euclid'")
set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
op <- function(suffix) file.path(out_dir, paste0(prefix, suffix))
message(" - 4_3n background null | ", prefix, " (", stage, ") | metric=", metric,
        " mode=", paste(modes, collapse = "+"), " window=", window, " B=", B, " fdr=", fdr_thr)

## ----------------------------------------------------------------- (1) markers ------------------
markers <- read.table(markers_bed, header = TRUE, sep = "\t", quote = "", comment.char = "")
markers <- markers[!duplicated(markers$geneID), ]
rownames(markers) <- markers$geneID
markers$species    <- ifelse(grepl("^Zm", markers$geneID), "Zm", "At")
markers$type_label <- paste0(markers$species, ":", markers$type)
message(" - panel markers (deduped): ", nrow(markers))

## ----------------------------------------------------------------- (2) metadata -----------------
b <- read.table(meta_txt, header = TRUE, sep = "\t", quote = "", comment.char = "")
if ("cellID" %in% colnames(b)) rownames(b) <- b$cellID
if (!cluster_col %in% colnames(b))
  stop("cluster column '", cluster_col, "' not in metadata; have: ", paste(colnames(b), collapse = ", "))
b$.cluster <- as.character(b[[cluster_col]])

## ------------------------------------------- (3) smoothed matrix (same intersection as 4_3i) ----
sm <- readRDS(smoothed_rds)                          # genes x cells, CP10k + kNN-diffused, perkb
cells <- intersect(colnames(sm), rownames(b))
if (length(cells) == 0) stop("no overlap between matrix colnames and metadata cellIDs")
sm <- sm[, cells, drop = FALSE]
b  <- b[cells, , drop = FALSE]
n_genes <- nrow(sm); n_cells <- ncol(sm)
message(" - matrix genes x cells: ", n_genes, " x ", n_cells)

## Zj background: per-cell mean/sd OVER ALL GENES (identical to 4_3i, incl. the gene-block trick)
cmean <- Matrix::colMeans(sm)
csumsq <- numeric(n_cells); blk <- 2000L
for (s in seq(1L, n_genes, by = blk)) {
  e <- min(s + blk - 1L, n_genes)
  csumsq <- csumsq + Matrix::colSums(sm[s:e, , drop = FALSE]^2)
}
cvar <- (csumsq - n_genes * cmean^2) / (n_genes - 1); cvar[cvar < 0] <- 0
csd  <- sqrt(cvar); csd[csd == 0 | is.na(csd)] <- 1
message(" - per-cell Zj background done (median cell sd = ", round(median(csd), 4), ")")

## ----------------------------------------------------------------- cluster weights --------------
cls <- unique(b$.cluster)
cls <- if (suppressWarnings(all(!is.na(as.numeric(cls))))) cls[order(as.numeric(cls))] else sort(cls)
n_cl <- length(cls)
Wc <- sparseMatrix(i = seq_len(n_cells), j = match(b$.cluster, cls), x = 1, dims = c(n_cells, n_cl))
Wc <- Wc %*% Diagonal(x = 1 / pmax(Matrix::colSums(Wc), 1))       # column cl averages over its cells
colnames(Wc) <- cls
message(" - clusters: ", paste(cls, collapse = ", "), " (n = ", paste(table(b$.cluster)[cls], collapse = ", "), ")")

## ------------------------------- (4) GENOME-WIDE pass: cluster-mean rZ for EVERY gene ------------
## Blocked so rZ is never materialized genome-wide: each block collapses straight to cluster means.
M      <- matrix(0, n_genes, n_cl, dimnames = list(rownames(sm), cls))
g_mean <- numeric(n_genes); g_sd <- numeric(n_genes)
gblk   <- 1000L
nb     <- length(seq(1L, n_genes, by = gblk))
message(" - genome-wide rZ pass in ", nb, " blocks of ", gblk, " genes ...")
bi <- 0L
for (s in seq(1L, n_genes, by = gblk)) {
  e   <- min(s + gblk - 1L, n_genes); bi <- bi + 1L
  X   <- as.matrix(sm[s:e, , drop = FALSE])
  rm_ <- rowMeans(X)
  rs_ <- sqrt((rowSums(X^2) - n_cells * rm_^2) / (n_cells - 1))
  g_mean[s:e] <- rm_; g_sd[s:e] <- rs_
  rs_[rs_ == 0 | is.na(rs_)] <- 1
  Zi <- sweep(sweep(X, 1, rm_, "-"), 1, rs_, "/")     # gene z across CELLS
  Zj <- sweep(sweep(X, 2, cmean, "-"), 2, csd, "/")   # cell z across GENES
  Zi[Zi < 0] <- 0; Zj[Zj < 0] <- 0
  rZ <- if (metric == "geom") sqrt(Zi * Zj) else sqrt(Zi^2 + Zj^2)
  M[s:e, ] <- as.matrix(rZ %*% Wc)
  rm(X, Zi, Zj, rZ)
  if (bi %% 10 == 0 || e == n_genes) { message("   block ", bi, "/", nb, " (gene ", e, ")"); gc(FALSE) }
}
message(" - genome-wide cluster-mean rZ done")

row_top2 <- function(Mx) {
  pk <- max.col(Mx, ties.method = "first")
  pv <- Mx[cbind(seq_len(nrow(Mx)), pk)]
  M2 <- Mx; M2[cbind(seq_len(nrow(Mx)), pk)] <- -Inf
  sk <- max.col(M2, ties.method = "first")
  sv <- M2[cbind(seq_len(nrow(M2)), sk)]
  sv[!is.finite(sv)] <- 0                             # single-cluster edge case (matches 4_3i)
  list(peak_cl = colnames(Mx)[pk], peak = pv, second = sv, spec = pv - sv)
}
S_raw <- row_top2(M)                                  # RAW scale: reconciles against 4_3i
genes <- rownames(sm)
is_mk <- genes %in% markers$geneID

gs <- data.frame(stage = stage, geneID = genes,
                 mean_expr = signif(g_mean, 5), sd_expr = signif(g_sd, 5),
                 peak_cluster = S_raw$peak_cl, peak_rZ = round(S_raw$peak, 4),
                 second_rZ = round(S_raw$second, 4), spec = round(S_raw$spec, 4),
                 is_marker = is_mk,
                 type_label = ifelse(is_mk, markers[genes, "type_label"], NA),
                 name = ifelse(is_mk, markers[genes, "name"], NA))
gzf <- gzfile(op(".bgnull.genestats.tsv.gz"), "w")
write.table(gs, gzf, sep = "\t", quote = FALSE, row.names = FALSE); close(gzf)
message(" - wrote genome-wide gene stats (", n_genes, " genes; ", sum(is_mk), " panel markers)")

## ---------------------------- background / marker index sets -------------------------------------
usable <- g_mean > 0 & g_sd > 0 & is.finite(g_mean) & is.finite(g_sd)
bg_i   <- which(!is_mk & usable)
mk_i   <- which(is_mk & usable)
if (length(bg_i) < 200) stop("background pool too small: ", length(bg_i))
if (!length(mk_i))      stop("no usable panel markers present in the matrix")
message(" - background pool: ", length(bg_i), " genes | testable markers: ", length(mk_i))

lm_ <- log10(g_mean + 1e-12); ls_ <- log10(g_sd + 1e-12)

## types eligible for the SET test (>= min_markers present, cell-cycle excluded) -------------------
mk_type <- markers[genes[mk_i], "type_label"]
drop_t  <- if (nzchar(exclude)) grepl(exclude, mk_type, ignore.case = TRUE) else rep(FALSE, length(mk_i))
if (nzchar(exclude)) message(" - excluding type regex '", exclude, "' (", sum(drop_t), " markers)")
tt        <- table(mk_type[!drop_t])
set_types <- names(tt)[tt >= min_markers]
message(" - set-testable types (>= ", min_markers, " markers): ", length(set_types))

## ------------------- covariate conditioning: LOCAL RANK (conditional quantile) --------------------
## The rZ statistic is ZERO-INFLATED by construction -- rZ = sqrt(max0(Zi)*max0(Zj)) is exactly 0
## whenever a gene sits below its cell's gene-wide mean, which holds in EVERY cell for lowly
## expressed genes. Empirically ~60% of background genes have peak rZ == 0, and the dependence on
## expression is a THRESHOLD (flat 0 across the lower deciles, then a steep climb), not a smooth
## trend. So location/scale standardisation is the wrong tool: median and MAD are both 0 over most
## of the covariate range (it collapses to a degenerate IQR of 0).
##
## Distribution-free replacement: express each gene as its MID-RANK QUANTILE among background genes
## in the same covariate neighbourhood. Ties (the zero block) share a mid-rank, so u stays calibrated;
## under the null u is uniform BY CONSTRUCTION whatever the shape of the trend or the zero mass.
##
## Neighbourhoods:
##   mode=mean     sliding window of `window` background genes by log-mean-expression RANK
##                 (rank-based, so a window always contains `window` real genes)
##   mode=mean_sd  2-D quantile bins: log-mean bins x log-sd bins within them
conditional_u <- function(Mx, grp, ref_list) {
  U <- matrix(NA_real_, nrow(Mx), ncol(Mx), dimnames = dimnames(Mx))
  for (k in seq_along(ref_list)) {
    gi <- which(grp == k); if (!length(gi)) next
    rf <- ref_list[[k]]; W <- length(rf); if (W < 20L) next
    for (ci in seq_len(ncol(Mx))) {
      wv  <- sort(Mx[rf, ci])
      m   <- Mx[gi, ci]
      nle <- findInterval(m, wv)                        # count of background <= m
      nlt <- findInterval(m, wv, left.open = TRUE)      # count of background <  m
      U[gi, ci] <- (nlt + nle) / (2 * W)                # mid-rank quantile, tie-safe
    }
  }
  U
}

## build (gene -> group) and (group -> background reference set) for a mode -------------------------
build_groups <- function(mode, W, n_anchor = 200L) {
  nbg <- length(bg_i)
  if (mode == "mean") {
    o    <- bg_i[order(lm_[bg_i])]                      # background sorted by expression
    xb   <- lm_[o]
    W    <- max(50L, min(W, nbg))
    half <- W %/% 2L
    anch <- unique(round(seq(half + 1L, nbg - half, length.out = min(n_anchor, max(1L, nbg - W)))))
    refs <- lapply(anch, function(a) o[(a - half):(a + half - 1L)])
    pos  <- findInterval(lm_, xb)                       # each gene's rank position in background
    ia   <- findInterval(pos, anch)                     # nearest anchor (anch is sorted)
    lo   <- pmax(ia, 1L); hi <- pmin(ia + 1L, length(anch))
    grp  <- ifelse(abs(pos - anch[lo]) <= abs(pos - anch[hi]), lo, hi)
    list(grp = grp, refs = refs, support = range(xb))
  } else {
    nq_m <- 20L; nq_s <- 5L
    qm   <- unique(quantile(lm_[bg_i], seq(0, 1, length.out = nq_m + 1)))
    bm   <- findInterval(lm_, qm, all.inside = TRUE)
    refs <- list(); grp <- integer(length(lm_)); k <- 0L
    for (i in seq_len(length(qm) - 1L)) {
      inb <- bg_i[bm[bg_i] == i]; alg <- which(bm == i)
      if (length(inb) < 50L) {                          # too thin to split on sd -> keep as one bin
        k <- k + 1L; refs[[k]] <- inb; grp[alg] <- k; next
      }
      qs <- unique(quantile(ls_[inb], seq(0, 1, length.out = nq_s + 1)))
      bs_bg <- findInterval(ls_[inb], qs, all.inside = TRUE)
      bs_al <- findInterval(ls_[alg], qs, all.inside = TRUE)
      for (j in seq_len(length(qs) - 1L)) {
        k <- k + 1L; refs[[k]] <- inb[bs_bg == j]; grp[alg[bs_al == j]] <- k
      }
    }
    list(grp = grp, refs = refs, support = range(lm_[bg_i]))
  }
}

## ============================ per-mode: correct, then test =======================================
run_mode <- function(mode) {
  message("\n===== mode: ", mode, " =====")
  G  <- build_groups(mode, window)
  UM <- conditional_u(M, G$grp, G$refs)                  # genes x clusters, conditional quantile
  SR <- cbind(spec = S_raw$spec, peak = S_raw$peak); rownames(SR) <- genes
  US <- conditional_u(SR, G$grp, G$refs)                 # genes x 2, for the gene-level test
  in_support <- lm_ >= G$support[1] & lm_ <= G$support[2]
  peak_cl_adj <- colnames(UM)[max.col(UM, ties.method = "first")]

  ## calibration check: background u must be ~uniform (mean 0.5, 5% above 0.95) whatever the trend
  message(sprintf(" - background conditional u: mean %.3f (target 0.500) | frac>0.95 = %.3f (target 0.050) | groups %d",
                  mean(UM[bg_i, ], na.rm = TRUE), mean(UM[bg_i, ] > 0.95, na.rm = TRUE), length(G$refs)))

  ## ================= (2) SET-LEVEL TEST: cluster x type  -- THE HEADLINE ========================
  ## The null replaces EACH marker by a random background gene from THAT MARKER'S OWN covariate
  ## neighbourhood -- it does NOT draw uniformly from the whole background.
  ##
  ## Why this matters (measured, not theoretical): markers of a type tend to share an expression
  ## range, so they share a reference window, so their u values are CORRELATED through that shared
  ## reference -- if a window happens to hold low background values in some cluster, every marker in
  ## it gets a high u at once. A null drawing uniformly from all background genes gives independent
  ## u's and is therefore far too narrow: on the synthetic fixture it produced null_sd 0.056 against
  ## an observed between-cluster scatter of ~0.18, and called a structure-free high-expression decoy
  ## set significant at z=4.4. Drawing each null gene from its marker's own window reproduces the
  ## set's covariate concentration, and with it the correct null variance.
  message(" - set test: ", length(set_types), " types x ", n_cl, " clusters, B=", B, " draws")
  nbg <- length(bg_i); set_rows <- list(); si <- 0L
  ## SENSITIVITY STATISTIC (top-k). The primary statistic is the MEAN u over all of a type's markers,
  ## which dilutes when only a subset of the panel applies to this tissue -- a real concern where the
  ## panel mixes sources (maize = Marand shoot markers augmented from other tissues). The top-k mean
  ## keeps only the k highest-u markers, so a signal carried by a subset survives. Computed off the
  ## SAME null draws, and the null uses the SAME statistic (else the comparison is meaningless).
  topk_mean <- function(v, k) {                       # column-wise mean of the k largest entries
    mm <- nrow(v); if (k >= mm) return(colMeans(v))
    s <- matrix(v[order(col(v), -v)], mm, ncol(v))    # each column sorted descending
    colMeans(s[seq_len(k), , drop = FALSE])
  }
  for (ty in set_types) {
    jj  <- which(mk_type == ty & !drop_t)
    m   <- length(jj)
    Uty    <- UM[mk_i[jj], , drop = FALSE]                       # m x n_cl (columns = clusters)
    obs    <- colMeans(Uty)
    obs_tk <- topk_mean(Uty, topk)
    refs_j <- lapply(jj, function(j) G$refs[[G$grp[mk_i[j]]]])   # each marker's own neighbourhood
    nulls    <- matrix(0, B, n_cl, dimnames = list(NULL, cls))
    nulls_tk <- matrix(0, B, n_cl, dimnames = list(NULL, cls))
    done <- 0L; blkB <- max(1L, min(B, as.integer(2e7 / max(m, 1))))
    while (done < B) {
      B0   <- min(blkB, B - done)
      pick <- matrix(0L, m, B0)
      for (r in seq_len(m)) pick[r, ] <- sample(refs_j[[r]], B0, replace = TRUE)
      for (ci in seq_len(n_cl)) {
        v <- UM[pick, ci]; dim(v) <- c(m, B0)
        nulls[(done + 1L):(done + B0), ci]    <- colMeans(v)
        nulls_tk[(done + 1L):(done + B0), ci] <- topk_mean(v, topk)
      }
      done <- done + B0
    }
    for (ci in seq_len(n_cl)) {
      nv <- nulls[, ci]; nvt <- nulls_tk[, ci]; si <- si + 1L
      set_rows[[si]] <- data.frame(
        stage = stage, mode = mode, metric = metric,
        cluster = cls[ci], type_label = ty, n_markers = m,
        frac_in_support = round(mean(in_support[mk_i[jj]]), 3),
        obs_mean_u = round(obs[ci], 4),
        obs_raw_rZ = round(mean(M[mk_i[jj], ci]), 4),
        null_mean = round(mean(nv), 4), null_sd = round(sd(nv), 4),
        null_q95  = round(as.numeric(quantile(nv, 0.95)), 4),
        effect_z  = round((obs[ci] - mean(nv)) / max(sd(nv), 1e-9), 3),
        p = (1 + sum(nv >= obs[ci])) / (1 + B),
        ## --- top-k sensitivity arm (same draws, same statistic on both sides) ---
        k_used = min(topk, m),
        obs_topk_u = round(obs_tk[ci], 4),
        null_topk_mean = round(mean(nvt), 4),
        effect_z_topk = round((obs_tk[ci] - mean(nvt)) / max(sd(nvt), 1e-9), 3),
        p_topk = (1 + sum(nvt >= obs_tk[ci])) / (1 + B))
    }
    message("   type ", ty, " (m=", m, ", k=", min(topk, m), ") done")
  }
  st <- do.call(rbind, set_rows)
  ## TWO multiplicity conventions, both reported -- they answer different questions and the choice
  ## matters exactly where signal is borderline (measured: At is robust to it, maize is not).
  ##   q            BH across ALL cluster x type pairs in the stage. Right for "WHICH pairs are
  ##                enriched?" (a discovery scan over the whole grid). CONSERVATIVE. Primary.
  ##   q_percluster BH within each cluster. Right for "does THIS cluster have an identity?", which
  ##                is what the annotation table actually claims, one verdict per cluster.
  st$q            <- p.adjust(st$p, method = "BH")
  st$q_percluster <- ave(st$p, st$cluster, FUN = function(p) p.adjust(p, method = "BH"))
  st$q_topk            <- p.adjust(st$p_topk, method = "BH")
  st$q_topk_percluster <- ave(st$p_topk, st$cluster, FUN = function(p) p.adjust(p, method = "BH"))
  ## A set whose genes lie outside the background expression range has no valid comparator: the
  ## neighbourhood has to reach down to lower-expressed genes and the set looks enriched for that
  ## reason alone. Measured on the synthetic fixture: an out-of-support decoy set was falsely
  ## significant in 3 of 6 replicates, while in-support decoys were 0/6. So this is enforced, not
  ## just reported -- such rows keep their numbers but cannot be called significant.
  st$untestable <- st$frac_in_support < min_support
  st$sig            <- st$q                 < fdr_thr & !st$untestable   # PRIMARY (global BH, mean)
  st$sig_percluster <- st$q_percluster      < fdr_thr & !st$untestable
  st$sig_topk       <- st$q_topk            < fdr_thr & !st$untestable
  st$sig_topk_percl <- st$q_topk_percluster < fdr_thr & !st$untestable
  if (any(st$untestable))
    message(" - [", mode, "] ", sum(st$untestable), " cluster x type rows marked UNTESTABLE (frac_in_support < ",
            min_support, ")")
  message(sprintf(" - [%s] significant: mean/globalBH %d | mean/perclusterBH %d | topk%d/globalBH %d | topk%d/perclusterBH %d  (of %d testable)",
                  mode, sum(st$sig), sum(st$sig_percluster), topk, sum(st$sig_topk), topk,
                  sum(st$sig_topk_percl), sum(!st$untestable)))
  st <- st[order(st$cluster, -st$effect_z), ]
  st$p <- signif(st$p, 4); st$q <- signif(st$q, 4)
  write.table(st, op(paste0(".bgnull.", mode, ".settest.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message(" - [", mode, "] set tests significant: ", sum(st$sig), "/", nrow(st),
          " (p floor = ", signif(1 / (1 + B), 3), ")")

  ## per-cluster annotation verdict ---------------------------------------------------------------
  ann <- do.call(rbind, lapply(cls, function(cl) {
    d  <- st[st$cluster == cl & !st$untestable, ]; d <- d[order(-d$effect_z), ]
    ds <- d[d$sig, ]
    e1 <- if (nrow(d)) d$effect_z[1] else NA
    e2 <- if (nrow(d) >= 2) d$effect_z[2] else NA
    data.frame(stage = stage, mode = mode, cluster = cl,
               n_types_tested = nrow(d), n_types_sig = nrow(ds),
               ## sensitivity: how the count moves under the other multiplicity / statistic choices
               n_sig_percluster = sum(d$sig_percluster), n_sig_topk = sum(d$sig_topk),
               n_sig_topk_percl = sum(d$sig_topk_percl),
               n_types_untestable = sum(st$cluster == cl & st$untestable),
               top_type = if (nrow(d)) d$type_label[1] else NA, top_effect_z = e1,
               top_q = if (nrow(d)) d$q[1] else NA,
               second_type = if (nrow(d) >= 2) d$type_label[2] else NA, second_effect_z = e2,
               ratio = if (!is.na(e2) && e1 > 0) round(e2 / e1, 3) else NA,
               verdict = if (!nrow(ds)) "no type beyond background" else
                         if (!is.na(e2) && e1 > 0 && (e2 / e1) >= 0.80) "AMBIGUOUS (competing types)" else
                         paste0("SIGNIFICANT: ", d$type_label[1]),
               all_sig_types = if (nrow(ds)) paste(sprintf("%s(z=%.1f,q=%.3g)", ds$type_label, ds$effect_z, ds$q),
                                                   collapse = " ; ") else "-")
  }))
  write.table(ann, op(paste0(".bgnull.", mode, ".annotation.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  ## ================= (1) GENE-LEVEL TEST: marker vs its own covariate neighbourhood ==============
  ## u is already the conditional quantile, so the upper-tail mid-p is simply 1 - u.
  ## Resolution is bounded by the neighbourhood size -> p floor ~ 1/(2*window).
  u_spec <- US[mk_i, "spec"]; u_peak <- US[mk_i, "peak"]
  p_spec <- 1 - u_spec; p_peak <- 1 - u_peak
  q_spec <- p.adjust(p_spec, method = "BH"); q_peak <- p.adjust(p_peak, method = "BH")

  mk <- data.frame(
    stage = stage, metric = metric, mode = mode,
    geneID = genes[mk_i], name = markers[genes[mk_i], "name"], type_label = mk_type,
    mean_expr = signif(g_mean[mk_i], 5), sd_expr = signif(g_sd[mk_i], 5),
    in_support = in_support[mk_i],
    peak_cluster_raw = S_raw$peak_cl[mk_i], peak_cluster_adj = peak_cl_adj[mk_i],
    peak_rZ_raw = round(S_raw$peak[mk_i], 4), spec_raw = round(S_raw$spec[mk_i], 4),
    u_spec = round(u_spec, 4), u_peak = round(u_peak, 4),
    p_spec = signif(p_spec, 4), q_spec = signif(q_spec, 4),
    p_peak = signif(p_peak, 4), q_peak = signif(q_peak, 4),
    sig_spec = q_spec < fdr_thr, sig_peak = q_peak < fdr_thr)
  mk <- mk[order(mk$peak_cluster_adj, mk$q_spec, -mk$u_spec), ]
  write.table(mk, op(paste0(".bgnull.", mode, ".markers.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message(" - [", mode, "] gene-level significant on spec: ", sum(mk$sig_spec), "/", nrow(mk),
          " (p floor = ", signif(1 / (1 + nbg), 3), ")")

  ## ---- figure ----------------------------------------------------------------------------------
  pdf(op(paste0(".bgnull.", mode, ".pdf")), width = 4 * min(n_cl, 3),
      height = 3.6 * (ceiling(n_cl / 3) + 1))
  par(mfrow = c(ceiling(n_cl / 3) + 1, min(n_cl, 3)), mar = c(4, 4, 3, 1))
  for (cl in cls) {
    d <- st[st$cluster == cl, ]; d <- d[order(-d$effect_z), ]
    if (!nrow(d)) { plot.new(); title(paste0("cluster ", cl, " (no types)")); next }
    dd <- d[seq_len(min(12L, nrow(d))), ]
    bp <- barplot(rev(dd$effect_z), horiz = TRUE, col = rev(ifelse(dd$sig, "firebrick", "grey70")),
                  border = NA, xlab = "effect (z vs size-matched background sets)", las = 1,
                  main = paste0("cluster ", cl, "  (", sum(d$sig), "/", nrow(d), " types sig)"))
    text(x = pmax(rev(dd$effect_z), 0), y = bp, labels = rev(sub("^[A-Za-z]+:", "", dd$type_label)),
         pos = 4, cex = 0.6, xpd = NA)
    abline(v = 0, col = "grey40")
  }
  ## why the conditioning is needed: raw rZ vs expression (the zero-inflated threshold)
  sub_bg <- sample(bg_i, min(4000L, length(bg_i)))
  plot(lm_[sub_bg], M[sub_bg, 1], pch = 16, cex = 0.25, col = "#00000022",
       xlab = "log10 mean expression", ylab = paste0("cluster-mean rZ (", cls[1], ")"),
       main = "raw rZ is expression-driven\n(zero-inflated, threshold-shaped)")
  points(lm_[mk_i], M[mk_i, 1], pch = 16, cex = 0.35, col = "#2166ac88")
  legend("topleft", c("background", "panel markers"), bty = "n", cex = 0.7,
         pch = 16, col = c("grey50", "#2166ac"))
  ## calibration: background conditional u must be flat (this is what makes the null valid)
  hist(UM[bg_i, ], breaks = 40, col = "grey85", border = "grey60", freq = FALSE,
       xlab = "conditional quantile u (background genes)", main = "null calibration\n(flat = valid)")
  abline(h = 1, col = "firebrick", lwd = 2, lty = 2)
  invisible(dev.off())

  cat("\n=========  ", prefix, " (", stage, ") | metric=", metric, " mode=", mode, "  =========\n", sep = "")
  print(ann[, c("cluster", "n_types_sig", "n_sig_percluster", "n_sig_topk", "n_sig_topk_percl",
                "top_type", "top_effect_z", "top_q", "verdict")], row.names = FALSE)
  list(st = st, ann = ann, mk = mk)
}

R <- lapply(modes, run_mode); names(R) <- modes

## ----------------------------------------------------------------- reconciliation + md -----------
## RAW peak rZ must reproduce 4_3i exactly (same matrix, same cells, same formula).
recon <- "not checked (no 4_3i cluster_mean TSV supplied)"
if (nzchar(cmp_tsv) && file.exists(cmp_tsv)) {
  cmp <- read.table(cmp_tsv, header = TRUE, sep = "\t", quote = "", comment.char = "")
  pc  <- paste0("peak_rZ_", metric); ref <- R[[1]]$mk
  if (pc %in% colnames(cmp)) {
    j <- match(cmp$geneID, ref$geneID); ok <- !is.na(j)
    if (sum(ok) > 2) {
      dd <- abs(cmp[[pc]][ok] - ref$peak_rZ_raw[j[ok]])
      recon <- sprintf("%d markers vs 4_3i: max|diff| = %.4g, r = %.6f%s", sum(ok), max(dd),
                       cor(cmp[[pc]][ok], ref$peak_rZ_raw[j[ok]]),
                       if (max(dd) < 1e-3) "  [PASS]" else "  [** MISMATCH -- investigate **]")
    }
  }
}

mode_block <- function(mn) {
  st <- R[[mn]]$st; ann <- R[[mn]]$ann; mk <- R[[mn]]$mk
  top <- head(st[st$sig, ][order(-st[st$sig, ]$effect_z), ], 15)
  c(paste0("## mode = ", mn, if (mn == "mean") "  (primary — expression corrected)"
           else "  (conservative — expression AND variability corrected)"),
    "",
    paste0("- **Set tests** (cluster x type) significant at q<", fdr_thr, ": **", sum(st$sig), " / ", nrow(st), "**"),
    paste0("- Clusters with a significant, non-ambiguous type call: **",
           sum(grepl("^SIGNIFICANT", ann$verdict)), " / ", nrow(ann), "**"),
    paste0("- Gene-level markers significant on adjusted `spec`: ", sum(mk$sig_spec), " / ", nrow(mk)),
    "",
    paste0("**Sensitivity** — significant set tests under each choice: mean+globalBH **", sum(st$sig),
           "**, mean+perclusterBH **", sum(st$sig_percluster), "**, top", topk, "+globalBH **",
           sum(st$sig_topk), "**, top", topk, "+perclusterBH **", sum(st$sig_topk_percl),
           "**. A call that moves between these is BORDERLINE, not established."),
    "",
    "### Per-cluster verdict",
    "| cluster | sig (mean/global) | sig (mean/per-cl) | sig (topk/global) | sig (topk/per-cl) | top type | effect z | q | verdict |",
    "|---|---|---|---|---|---|---|---|---|",
    apply(ann, 1, function(r) paste0("| ", r["cluster"], " | ", r["n_types_sig"], " | ", r["n_sig_percluster"],
                                     " | ", r["n_sig_topk"], " | ", r["n_sig_topk_percl"],
                                     " | ", r["top_type"], " | ", r["top_effect_z"], " | ", r["top_q"],
                                     " | ", r["verdict"], " |")),
    "",
    "### Strongest cluster x type enrichments",
    "| cluster | type | m | obs mean u | obs (raw rZ) | null mean u | effect z | q |",
    "|---|---|---|---|---|---|---|---|",
    if (!nrow(top)) "| (none) | | | | | | | |" else
      apply(top, 1, function(r) paste0("| ", r["cluster"], " | ", r["type_label"], " | ", r["n_markers"],
                                       " | ", r["obs_mean_u"], " | ", r["obs_raw_rZ"], " | ", r["null_mean"],
                                       " | ", r["effect_z"], " | ", r["q"], " |")),
    "")
}

md <- c(
  paste0("# 4_3n background-gene null (Test A) — ", prefix, " (", stage, ")"),
  "",
  paste0("Covariate-corrected background null for the per-cell reciprocal z of `4_3i`. Metric **", metric,
         "**, trend window **", window, "** genes, **", B, "** set-null draws, FDR ", fdr_thr, "."),
  "",
  "## Setup",
  paste0("- Smoothed perkb matrix: `", smoothed_rds, "`"),
  paste0("- Genes ", n_genes, " x cells ", n_cells, " | clusters: ", paste(cls, collapse = ", ")),
  paste0("- Panel: `", basename(markers_bed), "` — ", sum(is_mk), " of ", nrow(markers),
         " present; background ", length(bg_i), " non-panel genes"),
  paste0("- Types set-tested (>= ", min_markers, " markers, excluding `", exclude, "`): ", length(set_types)),
  paste0("- Reconciliation with 4_3i (raw peak rZ): ", recon),
  "",
  unlist(lapply(modes, mode_block)),
  "## Reading this",
  "- rZ rises with expression through its `Zj` term, and is ZERO-INFLATED (exactly 0 for genes below",
  "  their cell's gene-wide mean), so each gene is first converted to its MID-RANK QUANTILE `u` among",
  "  background genes at the same expression; `u` is uniform under the null by construction, whatever",
  "  the zero mass or trend shape. Two alternatives were tried and rejected: nearest-neighbour matching",
  "  (fails silently in sparse expression tails) and running median/MAD standardisation (degenerates,",
  "  since both location and scale are 0 over most of the range).",
  "- `null calibration` panel in the PDF must be FLAT — that is the assumption the p-values rest on.",
  "- The **set test** (cluster x type) is the headline: it is the level the annotation is made at, so BH",
  "  corrects over ~n_cluster x n_type tests rather than thousands of genes, and it retains power.",
  "  `effect_z` = (observed - null mean) / null sd, against random background sets of the SAME size.",
  "- The **gene test** is the detail layer — which individual markers carry the signal.",
  "- `in_support` flags genes outside the background expression range, where the trend is extrapolated flat.",
  "- The null is a gene permutation CONDITIONAL on the cell graph, so smoothing cannot inflate it;",
  "  a per-cell test on this matrix would be pseudoreplication and is deliberately not done.",
  "- Scope is within-object (one library): 'beyond equally-expressed genes here', not population inference.",
  "- `mode=mean` corrects expression level; `mode=mean_sd` also removes the variability trend, isolating",
  "  cluster-STRUCTURED variance — calls surviving it are the robust core."
)
writeLines(md, op(".bgnull.results.md"))
message("\n - wrote tables + pdfs + results md -> ", out_dir)
message(" - DONE")
