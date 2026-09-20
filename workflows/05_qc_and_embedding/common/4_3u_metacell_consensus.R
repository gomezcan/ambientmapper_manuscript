#!/usr/bin/env Rscript
## 4_3u_metacell_consensus.R -- CONSENSUS meta-cells from a multi-seed SEACells sweep.
## =================================================================================================
## A single SEACells run is one draw from a distribution: on At Pre two seeds agree at ARI 0.47.
## Shipping any one seed's partition as "the meta-cells" is therefore not defensible. This engine
## turns the 100-seed sweep (5_3 with SEEDS=100 OUT_TAG=_S100 MEMBERSHIP_ONLY=1) into ONE partition
## plus an explicit confidence attribute, so downstream work (Step 6 / DAR) has a unit it can name.
##
## THE OBJECT
##   F[i,j] = fraction of seeds placing cells i and j in the SAME meta-cell.  Built sparsely --
##   measured 1.8-2.3% non-zero on maize, 14-23% on At, ~19 s/stage for 14,427 cells.
##   Consensus partition = average-linkage hclust on (1 - F), cut at a PINNED k.
##
## THERE IS NO THRESHOLD ANYWHERE IN THIS ENGINE, AND THAT IS BY DESIGN.
##   F is a fraction; the cut is by k, not by height; loyalty and confidence are continuous. Pinning
##   k is precisely what removes the arbitrary number -- an emergent-k rule would need a cut height
##   or a PAC band, and that band would then be tunable against the answer. If a downstream consumer
##   wants "confident meta-cells only" it draws its own line on `loyalty`/`conf_meanF` and REPORTS
##   THE RETAINED FRACTION PER STAGE, which differs sharply (At Pre's continuum block has almost no
##   cells above loyalty 0.8; wd's largest keeps most of them). A confidence FILTER is a
##   stage-dependent treatment -- the `min.c=50` failure mode. Prefer weighting to filtering.
##
## WHY k IS PINNED AND NOT INFERRED.  `n_SEACells` is held at the deliverable's value (At 14,
##   maize 286) in EVERY stage, exactly as 5_3/5_7/5_8 do. Unmatched group counts MANUFACTURE stage
##   effects -- that is what killed the "maize nd lowest 17/17" claim in the count-matched controls.
##   Pinning k does NOT equalise meta-cell SIZE (At: 75.8 / 73.7 / 49.7 cells per meta-cell, because
##   nd lost 34% of its cells), and co-assignment falls ~0.136 per doubling of size, so `n_cells` is
##   emitted per meta-cell and MUST be conditioned on in any cross-stage per-meta-cell comparison.
##
## THE CONSENSUS IS NOT "THE AVERAGE SEED". Average linkage POOLS a continuum instead of carving
##   it: At Pre's largest consensus block is 369 cells against the largest seed block's 287. So it
##   is systematically COARSER where structure is absent and FINER where it is present. That is the
##   honest behaviour -- it says those 369 cells cannot be reliably separated -- but do not compare
##   consensus block sizes to seed block sizes as though they measured the same thing.
##
## MEASURED CHOICES (none of these are conventions)
##   * average linkage beat ward.D2 and complete on BOTH At stages, scored by ARI(consensus, seed):
##     Pre 0.555 / 0.516 / 0.538 ; wd 0.787 / 0.776 / 0.609. Linkage is re-scored every run.
##   * the cut does not chain: 0 singletons under all three linkages.
##   * ACCEPTANCE GATE -- mean ARI(consensus, seed) must EXCEED mean ARI(seed, seed). A consensus
##     that is no closer to the average seed than seeds are to each other is not representative and
##     nothing downstream should run on it. Measured At: 0.555 > 0.466 (Pre), 0.787 > 0.721 (wd).
##     HARD STOP on failure; set ALLOW_GATE_FAIL=1 only to inspect a known-bad case.
##   * split-half stability vs k is the DIAGNOSTIC that says how many units a stage really supports.
##     At: Pre declines above k~12 (0.85 -> 0.62 at k=20); wd and nd stay >=0.88 through k=24.
##     ⇒ k=14 sits AT the pre-clean object's resolution limit and INSIDE the cleaned ones'.
##     Read it on the cell-count-matched pair (Pre 1061 vs wd 1032). nd has 696 cells, so its
##     high-k stability is partly a small-block effect -- supporting evidence, not primary.
##
## OUTPUTS (<out_dir>/<prefix>.consensus.*)
##   cells.tsv      cellID, consensus_mc, loyalty, block_size, seed1_mc, LouvainClusters
##   metacells.tsv  per consensus meta-cell: n_cells, conf_meanF, stab, mean_loyalty, frac_F_ge0.9,
##                  dominant Leiden cluster + purity
##   summary.tsv    one row: gate numbers, F sparsity, high-confidence pair fraction
##   linkage.tsv    the linkage comparison actually scored on this object
##   ksweep.tsv     split-half stability vs k (if run)
##   seedpairs.tsv  pairwise ARI between seeds (the reproducibility distribution)
##   F.rds          sparse F + cell order (for heatmaps; ~2-3% dense)
##   heatmap.png    F ordered by the consensus dendrogram
##   results.md     human-readable
##
## `stab` = var(within-block F) / (mean*(1-mean)) in [0,1]. NOTE: It is a FUSION detector, NOT a quality
## score: the MEAN pair co-assignment of a block is forced by its per-run piece-size profile alone
## (identical for a boundary redrawn every run and one reproduced exactly), so reproducibility of
## the boundary lives only in the variance. High stab = the same split every run (fused cores);
## low stab = a different split every run (a continuum). Never rank stages on it.
##
## Usage:
##   Rscript 4_3u_metacell_consensus.R <out_dir> <prefix> <sweep_dir> <k> <stage> \
##           [n_seeds] [n_seedpairs] [ksweep_reps] [ksweep_maxn] [linkage]
##   n_seeds      0 = auto-detect seed*/ subdirs (default)
##   n_seedpairs  seed pairs sampled for the gate baseline; 0 = all C(n,2). default 500
##   ksweep_reps  split-half replicates; 0 = skip. default 25
##   ksweep_maxn  skip the sweep above this cell count. default 20000 (maize ~12 min; At ~1 min)
##   linkage      default "average" (re-scored against ward.D2/complete every run)
## =================================================================================================

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(TRUE)
if (length(args) < 5)
  stop("Usage: Rscript 4_3u_metacell_consensus.R <out_dir> <prefix> <sweep_dir> <k> <stage> ",
       "[n_seeds] [n_seedpairs] [ksweep_reps] [ksweep_maxn] [linkage]")
OUT    <- args[1]
PREFIX <- args[2]
SWEEP  <- args[3]
KPIN   <- as.integer(args[4])
STAGE  <- args[5]
NSEED_ARG <- if (length(args) >= 6) as.integer(args[6]) else 0L
NPAIR  <- if (length(args) >= 7) as.integer(args[7]) else 500L
KREPS  <- if (length(args) >= 8) as.integer(args[8]) else 25L
KMAXN  <- if (length(args) >= 9) as.integer(args[9]) else 20000L
LINK   <- if (length(args) >= 10 && nzchar(args[10])) args[10] else "average"

if (!is.finite(KPIN) || KPIN < 2) stop("k must be an integer >= 2, got: ", args[4])
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

msg <- function(...) message(sprintf(...))
msg("=== 4_3u consensus | %s (%s) | k pinned = %d | linkage = %s", PREFIX, STAGE, KPIN, LINK)

## ---- 1. load every seed's membership ------------------------------------------------------------
seed_dirs <- sort(Sys.glob(file.path(SWEEP, "seed*")))
seed_num  <- as.integer(sub(".*seed", "", basename(seed_dirs)))
seed_dirs <- seed_dirs[order(seed_num)]; seed_num <- sort(seed_num)
if (!length(seed_dirs)) stop("no seed*/ subdirectories under ", SWEEP)
if (NSEED_ARG > 0) { seed_dirs <- head(seed_dirs, NSEED_ARG); seed_num <- head(seed_num, NSEED_ARG) }
mem_file <- function(d) file.path(d, paste0(PREFIX, ".seacells.cell_to_seacell.tsv"))
miss <- seed_dirs[!file.exists(vapply(seed_dirs, mem_file, ""))]
if (length(miss)) stop("missing membership files in ", length(miss), " seed dir(s), e.g. ", miss[1])
NS <- length(seed_dirs)

a <- read.delim(mem_file(seed_dirs[1]), stringsAsFactors = FALSE)
cells <- a$index; n <- length(cells)
seed1_mc <- a$SEACell
louvain  <- if ("LouvainClusters" %in% names(a)) as.character(a$LouvainClusters) else rep(NA, n)

Lc <- matrix(NA_integer_, n, NS)          # integer label codes, per seed (labels are seed-local)
for (s in seq_len(NS)) {
  d <- read.delim(mem_file(seed_dirs[s]), stringsAsFactors = FALSE)
  i <- match(cells, d$index)
  if (anyNA(i)) stop("seed ", seed_num[s], " does not cover every cell of seed ", seed_num[1])
  Lc[, s] <- as.integer(factor(d$SEACell[i]))
}
ng <- apply(Lc, 2, max)
msg(" - %d seeds x %d cells | groups per seed: %s", NS, n,
    if (length(unique(ng)) == 1) as.character(ng[1]) else paste0(min(ng), "-", max(ng)))
if (length(unique(ng)) > 1)
  msg("   NOTE: group count VARIES across seeds -- cross-seed statistics assume it is pinned")

## ---- 2. sparse F --------------------------------------------------------------------------------
accumF <- function(cols) {
  Fs <- NULL
  for (s in cols) {
    A <- sparseMatrix(i = seq_len(n), j = Lc[, s], x = 1, dims = c(n, max(Lc[, s])))
    S <- tcrossprod(A)
    Fs <- if (is.null(Fs)) S else Fs + S
  }
  Fs <- Fs / length(cols); diag(Fs) <- 0; drop0(Fs)
}
t0 <- proc.time()[["elapsed"]]
Fm <- accumF(seq_len(NS))
v  <- Fm@x                                                  # each unordered pair appears twice
npair_all <- choose(n, 2)
pct_nz    <- 100 * (length(v) / 2) / npair_all
pct_all_hi <- 100 * (sum(v >= 0.9) / 2) / npair_all         # high-confidence pairs / ALL pairs
msg(" - F built in %.1fs | %.2f%% of pairs ever co-assigned | F>=0.9 = %.4f%% of ALL pairs",
    proc.time()[["elapsed"]] - t0, pct_nz, pct_all_hi)

## ---- 3. consensus partition, and the linkage comparison that justifies the choice -----------------
ari <- function(x, y) {
  t <- table(x, y); N <- sum(t)
  ii <- sum(choose(t, 2)); ai <- sum(choose(rowSums(t), 2)); bi <- sum(choose(colSums(t), 2))
  e <- ai * bi / choose(N, 2)
  (ii - e) / ((ai + bi) / 2 - e)
}
Fd <- as.matrix(Fm); D <- as.dist(1 - Fd)
lnk <- data.frame()
trees <- list()
for (lk in unique(c(LINK, "average", "ward.D2", "complete"))) {
  h <- hclust(D, lk); cl <- cutree(h, k = KPIN); trees[[lk]] <- list(h = h, cl = cl)
  wf <- vapply(sort(unique(cl)), function(g) {
    i <- which(cl == g); if (length(i) < 2) return(NA_real_)
    sub <- Fd[i, i]; mean(sub[upper.tri(sub)]) }, numeric(1))
  lnk <- rbind(lnk, data.frame(linkage = lk, n_singletons = sum(table(cl) == 1),
                               largest_block = max(table(cl)),
                               mean_within_F = round(mean(wf, na.rm = TRUE), 4),
                               ari_cons_seed = round(mean(vapply(seq_len(NS),
                                   function(s) ari(cl, Lc[, s]), numeric(1))), 4)))
}
cl <- trees[[LINK]]$cl; hc <- trees[[LINK]]$h
msg(" - linkage comparison (ARI to seeds): %s",
    paste(sprintf("%s %.3f", lnk$linkage, lnk$ari_cons_seed), collapse = " | "))
best <- lnk$linkage[which.max(lnk$ari_cons_seed)]
if (best != LINK)
  msg("   NOTE: %s scores higher (%.3f) than the requested %s (%.3f) on THIS object -- reported, not "
      , best, max(lnk$ari_cons_seed), LINK, lnk$ari_cons_seed[lnk$linkage == LINK])
if (best != LINK) msg("     switched automatically. Change the arg deliberately if you want it.")

## ---- 4. ACCEPTANCE GATE --------------------------------------------------------------------------
ari_cs <- vapply(seq_len(NS), function(s) ari(cl, Lc[, s]), numeric(1))
allp <- t(combn(NS, 2))
pick <- if (NPAIR > 0 && nrow(allp) > NPAIR) allp[sample(nrow(allp), NPAIR), , drop = FALSE] else allp
ari_ss <- apply(pick, 1, function(p) ari(Lc[, p[1]], Lc[, p[2]]))
gate <- mean(ari_cs) > mean(ari_ss)
msg(" - GATE: mean ARI(consensus, seed) = %.4f   vs   mean ARI(seed, seed) = %.4f  [%d pairs]  => %s",
    mean(ari_cs), mean(ari_ss), nrow(pick), if (gate) "PASS" else "FAIL")
if (!gate && !nzchar(Sys.getenv("ALLOW_GATE_FAIL")))
  stop("ACCEPTANCE GATE FAILED: the consensus is no closer to the average seed than seeds are to ",
       "each other, so it is not representative. Nothing downstream should use it. ",
       "Set ALLOW_GATE_FAIL=1 to write the outputs anyway for inspection.")

## ---- 5. per-cell and per-meta-cell attributes ----------------------------------------------------
loyalty <- numeric(n); bsize <- integer(n)
mcs <- sort(unique(cl))
per_mc <- data.frame()
for (g in mcs) {
  i <- which(cl == g); nb <- length(i); bsize[i] <- nb
  if (nb < 2) { loyalty[i] <- NA_real_
    per_mc <- rbind(per_mc, data.frame(consensus_mc = g, n_cells = nb, conf_meanF = NA_real_,
      stab = NA_real_, loyalty_sd = NA_real_, loyalty_min = NA_real_, frac_F_ge0.9 = NA_real_,
      dominant_cluster = louvain[i][1], purity = 1))
    next }
  sub <- Fd[i, i]; diag(sub) <- 0
  loyalty[i] <- rowSums(sub) / (nb - 1)
  pv <- sub[upper.tri(sub)]; m <- mean(pv)
  tb <- table(louvain[i])
  ## mean(loyalty) over a block is IDENTICALLY conf_meanF -- both reduce to the sum of the block's
  ## off-diagonal F over nb*(nb-1). It is not a second piece of evidence and is deliberately NOT
  ## emitted. What the per-cell layer adds that the pair layer cannot is the SPREAD: a block at
  ## conf 0.6 with loyalty_sd ~0 is uniformly middling, while the same conf with a large sd and a
  ## low loyalty_min is a solid core plus a halo. Those are different objects (cf. the pair-level
  ## bimodality that separated At Pre's fused contaminant core from its real one).
  per_mc <- rbind(per_mc, data.frame(
    consensus_mc = g, n_cells = nb, conf_meanF = round(m, 4),
    stab = round(if (m > 0 && m < 1) var(pv) / (m * (1 - m)) else NA_real_, 4),
    loyalty_sd = round(sd(loyalty[i]), 4), loyalty_min = round(min(loyalty[i]), 4),
    frac_F_ge0.9 = round(mean(pv >= 0.9), 4),
    dominant_cluster = names(tb)[which.max(tb)], purity = round(max(tb) / nb, 4)))
}
per_mc <- per_mc[order(-per_mc$n_cells), ]

cells_df <- data.frame(cellID = cells, consensus_mc = cl, loyalty = round(loyalty, 4),
                       block_size = bsize, seed1_mc = seed1_mc, LouvainClusters = louvain)

## ---- 6. split-half stability vs k (the DIAGNOSTIC behind "how many units does this stage support")
ks_tab <- NULL
if (KREPS > 0 && n <= KMAXN) {
  ## THE GRID MUST SCALE WITH k, NOT BE A FIXED SET OF SMALL NUMBERS. The question is whether the
  ## PINNED k exceeds what the stage supports, so the grid has to bracket it. A hardcoded 4,6,8,10,12
  ## answers that for At (k=14) and is meaningless for maize (k=286), where it would report on
  ## top-of-tree splits nobody asked about while costing the same as the informative points.
  KS <- sort(unique(round(KPIN * c(0.3, 0.5, 0.7, 0.85, 1, 1.15, 1.3, 1.5, 1.75))))
  KS <- KS[KS >= 2 & KS < n]
  ## Each replicate rebuilds F twice and runs two hclusts; on maize that is ~40 s, so 25 replicates
  ## x 3 stages would be ~50 min for an estimate that is already tight (286 blocks average out).
  if (n > 5000 && KREPS > 10) {
    msg("   note: n = %d > 5000, reducing split-half replicates %d -> 10 (the estimate is far ",
        n, KREPS); msg("   tighter at this block count; raise ksweep_reps deliberately to override)")
    KREPS <- 10L
  }
  msg(" - split-half stability: %d replicates over k = %s", KREPS, paste(KS, collapse = ","))
  M <- matrix(NA_real_, length(KS), KREPS)
  for (r in seq_len(KREPS)) {
    h <- sample(NS, floor(NS / 2))
    hA <- hclust(as.dist(1 - as.matrix(accumF(h))), LINK)
    hB <- hclust(as.dist(1 - as.matrix(accumF(setdiff(seq_len(NS), h)))), LINK)
    M[, r] <- vapply(KS, function(k) ari(cutree(hA, k), cutree(hB, k)), numeric(1))
  }
  ks_tab <- data.frame(k = KS, mean = round(rowMeans(M), 4), sd = round(apply(M, 1, sd), 4))
  msg("   %s", paste(sprintf("k=%d %.3f", ks_tab$k, ks_tab$mean), collapse = "  "))
} else if (KREPS > 0) {
  msg(" - split-half stability SKIPPED: n = %d > ksweep_maxn = %d", n, KMAXN)
}

## ---- 7. heatmap ----------------------------------------------------------------------------------
pal <- colorRampPalette(c("white", "#f7e9c6", "#e8a33d", "#a8341f", "#2b0b06"))(64)
if (n <= 3000) { idx <- seq_len(n); sub_note <- sprintf("all %d cells", n)
} else { mc_s <- sample(mcs, min(30, length(mcs)))       # WHOLE blocks, so structure survives
         idx <- which(cl %in% mc_s); sub_note <- sprintf("%d sampled meta-cells, %d cells",
                                                         length(mc_s), length(idx)) }
## NEVER let the graphics device lose the run. A compute node without a working cairo/X11 device
## makes png() throw, and everything above it -- F, the gate, the k-sweep, the whole reason the job
## was submitted -- would be discarded on a cosmetic failure. The heatmap is regenerable from F.rds.
tryCatch({
  Ms <- Fd[idx, idx]
  ord <- hclust(as.dist(1 - Ms), LINK)$order
  png(file.path(OUT, paste0(PREFIX, ".consensus.heatmap.png")), width = 900, height = 900, res = 110)
  par(mar = c(2, 2, 3, 1))
  image(Ms[ord, ord], col = pal, zlim = c(0, 1), axes = FALSE, useRaster = TRUE,
        main = sprintf("%s  %s  (%s)", PREFIX, STAGE, sub_note))
  box(); dev.off()
}, error = function(e) {
  msg("   NOTE: heatmap SKIPPED (%s: %s) -- regenerate later from %s.consensus.F.rds",
      class(e)[1], conditionMessage(e), PREFIX)
  try(dev.off(), silent = TRUE)
})

## ---- 8. write ------------------------------------------------------------------------------------
w <- function(x, s) write.table(x, file.path(OUT, paste0(PREFIX, ".consensus.", s)),
                                sep = "\t", quote = FALSE, row.names = FALSE)
w(cells_df, "cells.tsv"); w(per_mc, "metacells.tsv"); w(lnk, "linkage.tsv")
w(data.frame(seed_a = pick[, 1], seed_b = pick[, 2], ari = round(ari_ss, 5)), "seedpairs.tsv")
if (!is.null(ks_tab)) w(ks_tab, "ksweep.tsv")
summ <- data.frame(prefix = PREFIX, stage = STAGE, n_cells = n, n_seeds = NS, k_pinned = KPIN,
  linkage = LINK, pct_pairs_nonzero = round(pct_nz, 3), pct_all_pairs_F_ge0.9 = round(pct_all_hi, 5),
  cells_per_metacell = round(n / KPIN, 2), largest_block = max(per_mc$n_cells),
  mean_within_F = round(mean(per_mc$conf_meanF, na.rm = TRUE), 4),
  ari_cons_seed = round(mean(ari_cs), 4), ari_seed_seed = round(mean(ari_ss), 4),
  gate_pass = gate, n_singletons = sum(per_mc$n_cells == 1))
w(summ, "summary.tsv")
saveRDS(list(F = Fm, cells = cells, consensus = cl, seed1 = seed1_mc, k = KPIN, stage = STAGE),
        file.path(OUT, paste0(PREFIX, ".consensus.F.rds")))

md <- c(sprintf("# Consensus meta-cells -- %s (%s)", PREFIX, STAGE), "",
  sprintf("- %d cells, %d seeds, k **pinned at %d** (%.1f cells per meta-cell), linkage `%s`",
          n, NS, KPIN, n / KPIN, LINK),
  sprintf("- F: %.2f%% of pairs ever co-assigned; **%.4f%% of ALL pairs reach F>=0.9**",
          pct_nz, pct_all_hi),
  sprintf("- largest consensus block **%d cells** (%.1f%%); mean within-block F **%.3f**",
          max(per_mc$n_cells), 100 * max(per_mc$n_cells) / n, mean(per_mc$conf_meanF, na.rm = TRUE)),
  "",
  sprintf("## Acceptance gate: %s", if (gate) "**PASS**" else "**FAIL**"),
  sprintf("mean ARI(consensus, seed) = **%.4f** vs mean ARI(seed, seed) = **%.4f** over %d pairs.",
          mean(ari_cs), mean(ari_ss), nrow(pick)),
  "A consensus no closer to the average seed than seeds are to each other is not representative.",
  "", "## Linkage (re-scored on this object)", "",
  paste0("| ", paste(names(lnk), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(lnk)), collapse = "|"), "|"),
  apply(lnk, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
if (!is.null(ks_tab)) md <- c(md, "## Split-half stability vs k", "",
  "How many units does this stage actually support? Stability falls once k exceeds it.",
  "NOTE: Compare stages only at matched cell count -- smaller blocks reproduce more easily.", "",
  paste0("| k | ", paste(ks_tab$k, collapse = " | "), " |"),
  paste0("|---|", paste(rep("---", nrow(ks_tab)), collapse = "|"), "|"),
  paste0("| mean | ", paste(sprintf("%.3f", ks_tab$mean), collapse = " | "), " |"),
  paste0("| sd | ", paste(sprintf("%.3f", ks_tab$sd), collapse = " | "), " |"), "")
md <- c(md, "## Reading the outputs", "",
  "- `loyalty` (per cell) and `conf_meanF` (per meta-cell) are the confidence attributes.",
  "- NOTE: **Weight, do not filter.** A loyalty cut removes far more cells from Pre than from wd,",
  "  so it is a stage-dependent treatment. If you filter, report the retained fraction per stage.",
  "- NOTE: `stab` is a FUSION detector, not a quality score -- never rank stages on it.",
  "- NOTE: `n_cells` must be conditioned on in any cross-stage per-meta-cell comparison:",
  "  pinning k does not equalise size, and co-assignment falls ~0.136 per doubling.")
writeLines(md, file.path(OUT, paste0(PREFIX, ".consensus.results.md")))
msg(" - wrote -> %s/%s.consensus.{cells,metacells,summary,linkage,seedpairs%s}.tsv + F.rds + "
    , OUT, PREFIX, if (!is.null(ks_tab)) ",ksweep" else "")
msg("   heatmap.png + results.md")
