#!/usr/bin/env Rscript
## =================================================================================================
## 4_3y_seedpair_ami_ci.R -- AMI + bootstrap CIs for the 100-seed reproducibility headline.
##
## WHY THIS EXISTS. 4_3u emits the seed-pair ARI distribution (seedpairs.tsv) and its means, but
## it computes neither (a) AMI nor (b) bootstrap CIs on the headline (mean ARI over all 4,950
## seed pairs per stage, and the wd-Pre / wd-nd stage differences). This script computes both and
## writes the AMI pair distribution to disk.
##
## WHAT IT DOES (one invocation per GENOME -- stage differences need all 3 stages in one process):
##   * ARI is READ from 4_3u's <prefix>.consensus.seedpairs.tsv, never recomputed, so it cannot
##     drift from the acceptance gate. RECONCILIATION GATE (hard stop): mean of seedpairs.tsv must
##     equal summary.tsv's ari_seed_seed within rounding (5e-5).
##   * AMI is recomputed from the _S100 membership files (aligned by cell index; hard stop if the
##     cell sets differ across seeds). AMI = (MI - EMI) / ((H1+H2)/2 - EMI) -- arithmetic-mean
##     normalization, sklearn's default. EMI is the exact hypergeometric expectation (Vinh et al.
##     2010), precomputed once per stage as a term matrix over the observed group SIZES -- sizes
##     recur heavily (~50-cell meta-cells), so the whole 4,950-pair sweep costs seconds, not hours.
##   * Percentile bootstrap CIs (default B=10000): RESAMPLE RUNS (seeds) WITH REPLACEMENT, never
##     pairs -- pairs sharing a seed are dependent, so pair-level resampling is anticonservative.
##     Statistic = mean of the pair matrix over resampled positions whose underlying seeds differ
##     (a seed drawn twice contributes no self-pair). Stage differences resample each stage
##     independently: stages are separate SEACells runs, there is nothing to pair on.
##   * FIXTURE GATE first (~1 s; SKIP_FIXTURE=1 to skip): AMI(P,P)=1, invariance under label
##     permutation, AMI of independent random partitions ~ 0. Hard stop on failure.
##
## READING RULES: quote ARI as the headline and AMI as the robustness check (AMI sits above ARI
## and compresses the stage effect). The word is *reproducibility*, never *recovery* or *accuracy*.
##
## Usage:
##   Rscript 4_3y_seedpair_ami_ci.R <consensus_dir> <sweep_dir> <out_dir> <At|maize> [B]
## Outputs (<out_dir>/):
##   <genome>.seedpair_ami.tsv   stage, seed_a, seed_b, ami          (3 x 4,950 rows)
##   <genome>.ari_ami_ci.tsv     per-stage means + CIs, stage differences + CIs
##   <genome>.ari_ami_ci.md      human-readable, incl. reconciliation vs the recorded numbers
## =================================================================================================

args <- commandArgs(TRUE)
if (length(args) < 4)
  stop("Usage: Rscript 4_3y_seedpair_ami_ci.R <consensus_dir> <sweep_dir> <out_dir> <At|maize> [B]")
CONS  <- args[1]
SWEEP <- args[2]
OUT   <- args[3]
GEN   <- args[4]
B     <- if (length(args) >= 5) as.integer(args[5]) else 10000L
if (!GEN %in% c("At", "maize")) stop("genome must be 'At' or 'maize', got: ", GEN)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) message(sprintf(...))
msg("=== 4_3y seed-pair AMI + bootstrap CIs | genome = %s | B = %d", GEN, B)

PFX <- if (GEN == "At") {
  c(Pre = "SM2_At_TAIR10", wd = "Clean.SM2v2wd_At_TAIR10", nd = "Clean.SM2v2_At_TAIR10")
} else {
  c(Pre = "SM2_B73_B73v5", wd = "Clean.SM2v2wd_B73_B73v5", nd = "Clean.SM2v2_B73_B73v5")
}

## Previously recorded wd-Pre ARI estimates (first analysis pass) -- used for a REPORTED
## reconciliation only, never as a gate: the original bootstrap's B and RNG are unknown, so
## agreement is judged to Monte-Carlo tolerance, not equality.
REF_WD_PRE <- if (GEN == "At") {
  c(est = 0.253, lo = 0.239, hi = 0.266)
} else {
  c(est = -0.035, lo = -0.041, hi = -0.029)
}

## ---- AMI machinery ------------------------------------------------------------------------------

## Exact expected MI term for one (row-size, col-size) pair under the permutation model.
emi_term <- function(va, vb, n) {
  lo <- max(1L, va + vb - n); hi <- min(va, vb)
  if (hi < lo) return(0)
  nij <- lo:hi
  lp  <- lchoose(va, nij) + lchoose(n - va, vb - nij) - lchoose(n, vb)
  sum((nij / n) * log(n * nij / (as.numeric(va) * vb)) * exp(lp))
}

## Term matrix over all observed sizes (EMI = sum over group pairs of Tm[size_a, size_b]).
make_T <- function(sizes, n) {
  Tm <- matrix(0, length(sizes), length(sizes))
  for (x in seq_along(sizes)) for (y in seq_along(sizes))
    Tm[x, y] <- emi_term(sizes[x], sizes[y], n)
  Tm
}

## AMI of two integer label vectors (1..ki / 1..kj), with a precomputed size->index term matrix.
ami_pair <- function(ci, cj, n, Sz, Tm) {
  ki <- max(ci); kj <- max(cj)
  M  <- matrix(tabulate((ci - 1L) * kj + cj, nbins = ki * kj), nrow = ki, byrow = TRUE)
  a  <- rowSums(M); b <- colSums(M)
  pos <- which(M > 0)
  v   <- M[pos]
  ri  <- ((pos - 1L) %% ki) + 1L
  cx  <- ((pos - 1L) %/% ki) + 1L
  MI  <- sum(v / n * log(n * v / (a[ri] * b[cx])))
  an  <- a[a > 0]; bn <- b[b > 0]
  Ha  <- -sum((an / n) * log(an / n))
  Hb  <- -sum((bn / n) * log(bn / n))
  EMI <- sum(Tm[match(an, Sz), match(bn, Sz), drop = FALSE])
  den <- (Ha + Hb) / 2 - EMI
  if (abs(den) < 1e-12) return(0)
  (MI - EMI) / den
}

ami_standalone <- function(ci, cj) {          # fixture path: builds its own term matrix
  n  <- length(ci)
  Sz <- sort(unique(c(tabulate(ci), tabulate(cj)))); Sz <- Sz[Sz > 0]
  ami_pair(ci, cj, n, Sz, make_T(Sz, n))
}

## ---- fixture gate -------------------------------------------------------------------------------
if (!nzchar(Sys.getenv("SKIP_FIXTURE"))) {
  set.seed(99)
  n <- 2000L; k <- 14L
  p <- sample.int(k, n, TRUE)
  self <- ami_standalone(p, p)
  perm <- sample.int(k); relab <- ami_standalone(p, perm[p])
  q <- sample.int(k, n, TRUE); indep <- ami_standalone(p, q)
  msg("fixture: AMI(self) = %.6f | AMI(relabelled) = %.6f | AMI(independent) = %.6f",
      self, relab, indep)
  if (abs(self - 1) > 1e-9)  stop("FIXTURE FAIL: AMI of a partition with itself != 1")
  if (abs(relab - 1) > 1e-9) stop("FIXTURE FAIL: AMI not invariant under label permutation")
  if (abs(indep) > 0.05)     stop("FIXTURE FAIL: AMI of independent partitions not ~0: ", indep)
  msg("fixture gate PASSED")
}

## ---- load one stage -----------------------------------------------------------------------------
read_members <- function(prefix) {
  sd  <- Sys.glob(file.path(SWEEP, "seed*"))
  num <- as.integer(sub(".*seed", "", basename(sd)))
  sd  <- sd[order(num)]
  if (!length(sd)) stop("no seed*/ dirs under ", SWEEP)
  cells0 <- NULL; L <- vector("list", length(sd))
  for (s in seq_along(sd)) {
    f <- file.path(sd[s], paste0(prefix, ".seacells.cell_to_seacell.tsv"))
    if (!file.exists(f)) stop("missing membership file: ", f)
    a <- read.delim(f, stringsAsFactors = FALSE)
    if (!all(c("index", "SEACell") %in% names(a)))
      stop("membership file lacks index/SEACell columns: ", f)
    if (is.null(cells0)) { cells0 <- a$index; o <- seq_along(cells0) }
    else {
      o <- match(cells0, a$index)
      if (anyNA(o) || length(a$index) != length(cells0))
        stop("cell set differs across seeds at ", f, " -- the same-cells assumption is violated")
    }
    L[[s]] <- as.integer(factor(a$SEACell[o]))
  }
  list(mat = do.call(cbind, L), cells = cells0)
}

## ---- per-stage: ARI matrix (from seedpairs.tsv) + AMI matrix (recomputed) -----------------------
stages   <- names(PFX)
mats     <- list()
ami_long <- list()
mean_tab <- list()

for (st in stages) {
  prefix <- PFX[[st]]
  t0 <- Sys.time()

  sp <- read.delim(file.path(CONS, paste0(prefix, ".consensus.seedpairs.tsv")))
  if (!all(c("seed_a", "seed_b", "ari") %in% names(sp)))
    stop("unexpected seedpairs.tsv columns for ", prefix)
  NS <- max(sp$seed_a, sp$seed_b)
  if (nrow(sp) != choose(NS, 2))
    stop(prefix, ": seedpairs.tsv has ", nrow(sp), " rows, expected all C(", NS, ",2) -- ",
         "was 4_3u run with n_seedpairs=0?")

  smr <- read.delim(file.path(CONS, paste0(prefix, ".consensus.summary.tsv")))
  if (abs(mean(sp$ari) - smr$ari_seed_seed) >= 5e-5)
    stop(prefix, ": RECONCILIATION FAIL -- mean(seedpairs.tsv) ", sprintf("%.6f", mean(sp$ari)),
         " vs summary.tsv ari_seed_seed ", smr$ari_seed_seed)

  mem <- read_members(prefix)
  if (ncol(mem$mat) != NS)
    stop(prefix, ": ", ncol(mem$mat), " membership files vs ", NS, " seeds in seedpairs.tsv")
  n <- nrow(mem$mat)

  sizes <- sort(unique(unlist(lapply(seq_len(NS), function(s) {
    tt <- tabulate(mem$mat[, s]); tt[tt > 0] }))))
  Tm <- make_T(sizes, n)

  cmb <- combn(NS, 2)
  ami <- numeric(ncol(cmb))
  for (p in seq_len(ncol(cmb))) {
    ami[p] <- ami_pair(mem$mat[, cmb[1, p]], mem$mat[, cmb[2, p]], n, sizes, Tm)
    if (p %% 1000 == 0) msg("  %s: %d / %d pairs", st, p, ncol(cmb))
  }

  Ma <- matrix(NA_real_, NS, NS)
  Ma[cbind(sp$seed_a, sp$seed_b)] <- sp$ari
  Ma[cbind(sp$seed_b, sp$seed_a)] <- sp$ari
  Mm <- matrix(NA_real_, NS, NS)
  Mm[cbind(cmb[1, ], cmb[2, ])] <- ami
  Mm[cbind(cmb[2, ], cmb[1, ])] <- ami

  mats[[st]]     <- list(ARI = Ma, AMI = Mm, NS = NS)
  ami_long[[st]] <- data.frame(stage = st, seed_a = cmb[1, ], seed_b = cmb[2, ], ami = round(ami, 5))
  mean_tab[[st]] <- c(ari = mean(sp$ari), ami = mean(ami))
  msg("stage %-3s (%s): n_cells=%d  mean ARI=%.4f  mean AMI=%.4f  [%.1f s]",
      st, prefix, n, mean(sp$ari), mean(ami), as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

## ---- bootstrap over runs ------------------------------------------------------------------------
set.seed(2026)
bs <- list()
for (st in stages) {
  NS <- mats[[st]]$NS
  ut <- upper.tri(matrix(0, NS, NS))
  bsA <- numeric(B); bsM <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(NS, NS, replace = TRUE)
    MA <- mats[[st]]$ARI[idx, idx]; bsA[b] <- mean(MA[ut], na.rm = TRUE)
    MM <- mats[[st]]$AMI[idx, idx]; bsM[b] <- mean(MM[ut], na.rm = TRUE)
  }
  bs[[st]] <- list(ARI = bsA, AMI = bsM)
}

ci <- function(x) quantile(x, c(0.025, 0.975), names = FALSE)

rows <- list()
for (st in stages) for (m in c("ARI", "AMI")) {
  est <- if (m == "ARI") mean_tab[[st]]["ari"] else mean_tab[[st]]["ami"]
  q <- ci(bs[[st]][[m]])
  rows[[length(rows) + 1]] <- data.frame(genome = GEN, kind = "stage", name = st, metric = m,
                                         estimate = round(est, 4), ci_lo = round(q[1], 4),
                                         ci_hi = round(q[2], 4), B = B, n_pairs = choose(mats[[st]]$NS, 2))
}
diffs <- list(c("wd", "Pre"), c("nd", "Pre"), c("wd", "nd"))
for (d in diffs) for (m in c("ARI", "AMI")) {
  est <- (if (m == "ARI") mean_tab[[d[1]]]["ari"] - mean_tab[[d[2]]]["ari"]
          else            mean_tab[[d[1]]]["ami"] - mean_tab[[d[2]]]["ami"])
  q <- ci(bs[[d[1]]][[m]] - bs[[d[2]]][[m]])
  rows[[length(rows) + 1]] <- data.frame(genome = GEN, kind = "diff",
                                         name = paste0(d[1], "-", d[2]), metric = m,
                                         estimate = round(est, 4), ci_lo = round(q[1], 4),
                                         ci_hi = round(q[2], 4), B = B, n_pairs = NA)
}
tab <- do.call(rbind, rows)

## ---- write --------------------------------------------------------------------------------------
f_ami <- file.path(OUT, paste0(GEN, ".seedpair_ami.tsv"))
f_ci  <- file.path(OUT, paste0(GEN, ".ari_ami_ci.tsv"))
f_md  <- file.path(OUT, paste0(GEN, ".ari_ami_ci.md"))
write.table(do.call(rbind, ami_long), f_ami, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab, f_ci, sep = "\t", quote = FALSE, row.names = FALSE)

wdpre <- tab[tab$kind == "diff" & tab$name == "wd-Pre" & tab$metric == "ARI", ]
recon <- sprintf(paste0("wd-Pre ARI: %+.4f [%+.4f, %+.4f]  vs recorded %+.3f [%+.3f, %+.3f]  -> %s"),
                 wdpre$estimate, wdpre$ci_lo, wdpre$ci_hi,
                 REF_WD_PRE["est"], REF_WD_PRE["lo"], REF_WD_PRE["hi"],
                 if (abs(wdpre$ci_lo - REF_WD_PRE["lo"]) < 0.004 &&
                     abs(wdpre$ci_hi - REF_WD_PRE["hi"]) < 0.004) "CONSISTENT (within MC tolerance)"
                 else "WARNING: INCONSISTENT -- bootstrap scheme differs from the recorded one")

md <- c(sprintf("# %s -- seed-pair ARI/AMI with bootstrap CIs (4_3y, B=%d)", GEN, B),
        "",
        "ARI read from 4_3u seedpairs.tsv (reconciled vs summary.tsv, gate passed).",
        "AMI recomputed from the _S100 memberships (fixture gate passed).",
        "Bootstrap resamples RUNS with replacement; percentile CIs; stages independent.",
        "",
        "| name | metric | estimate | 2.5% | 97.5% |",
        "|---|---|---|---|---|",
        sprintf("| %s | %s | %+.4f | %+.4f | %+.4f |", tab$name, tab$metric,
                tab$estimate, tab$ci_lo, tab$ci_hi),
        "",
        paste0("**Reconciliation vs the previously recorded headline:** ", recon),
        "",
        "NOTE: The word is *reproducibility*, never *recovery* or *accuracy*.")
writeLines(md, f_md)

msg("---- %s", recon)
msg("wrote: %s", f_ami)
msg("wrote: %s", f_ci)
msg("wrote: %s", f_md)
