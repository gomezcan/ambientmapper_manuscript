#!/usr/bin/env Rscript
###############################################################################
## 3_1b_coembed_meta_enrich.R  --  STAGE 3.1b: give the co-embed meta the full v6 schema
##
## WHY THIS STEP EXISTS (found the hard way)
## 3_1 emits a co-embed metadata with only its own columns:
##     cellID Genome species stage in_B73v5 in_TAIR10 both log10nSites pOrg joint_depth
## That is enough for the 3_2 grid (which only correlates log10nSites / pOrg), but NOT for
## downstream consumers, and the failure mode is nasty -- late and non-obvious:
##
##   1. `Socrates::callClusters` dies in its "filtering clusters with low cell/read counts" step
##      with `invalid type (NULL) for variable 'sro.meta$nSites'`. It needs the raw `nSites`
##      column (the co-embed carries only `log10nSites`). This happens AFTER cleanData, SVD, UMAP
##      and the graph build have all succeeded -- ~2 minutes of work thrown away at the last step.
##   2. the Fig 1 script (analysis/fig1_ambient_contamination/fig1.R) gates panels F, G and H on
##      `total >= 500`, `pTSS >= 0.2`,
##      `FRiP >= 0.2` -- none of which exist either.
##
## Rather than patch each consumer for each missing column, restore the schema at the source:
## `2_3_cluster.R` runs this exact clustering engine successfully on the per-genome v6 metadata,
## so handing the co-embed that same column set puts it on the known-good path, and the Fig-1
## gate columns come along for free (callClusters merges meta into $Clusters, which is what the
## figure reads).
##
## QC-SOURCE RULE -- prefer B73v5, else TAIR10.
## Identical to 3_1's `pick()` for log10nSites/pOrg, so for any given cell every QC column comes
## from the SAME genome (mixing sources per column would be incoherent), and the grafted `nSites`
## is exactly consistent with the existing `log10nSites` (verified: v6 log10nSites == log10(nSites)).
## Note this is "prefer the MAJOR genome", not "the cell's own plate genome", and that is
## deliberate: an At-plate cell sitting in the maize blob is there *because* it is full of maize
## ambient reads. Gating it on its own (shallow) TAIR10 depth would drop exactly the cells that
## Fig 1F exists to show. Gating on the genome that actually drives its embedding position keeps them.
##
## Writes a NEW file -- the 3_1 output is never modified.
##
## Usage:
##   Rscript 3_1b_coembed_meta_enrich.R <coembed_meta_tsv> <v6_B73v5> <v6_TAIR10> <out_tsv>
###############################################################################

## Base R only -- no data.table. Deliberate: the v6 files carry the R row-name column shift
## (N header names, N+1 data fields) and base `read.table` has DOCUMENTED handling for exactly
## that case ("if the header line has one fewer field than the number of columns, the first column
## is used for the row names"), whereas fread's handling of the mismatch does not surface `cellID`
## (a real failure mode here). This also keeps the script runnable in a bare R install.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 3_1b_coembed_meta_enrich.R <coembed_meta_tsv> <v6_B73v5> <v6_TAIR10> <out_tsv>")
}
ce_path   <- args[1]
v6_B_path <- args[2]
v6_A_path <- args[3]
out_path  <- args[4]

for (p in c(ce_path, v6_B_path, v6_A_path)) {
  if (!file.exists(p)) stop("MISSING input: ", p)
}
if (normalizePath(ce_path, mustWork = TRUE) == suppressWarnings(normalizePath(out_path, mustWork = FALSE))) {
  stop("refusing to overwrite the 3_1 output in place -- choose a different <out_tsv>")
}

## ---- reader that tolerates the R rowname-column shift -----------------------
## v6 files: write.table(row.names = TRUE) -> N header names, N+1 data fields. read.table absorbs
## the extra leading field into row.names, so the named columns line up and `cellID` is a real
## column. 3_1's co-embed meta: written row.names = FALSE -> read normally. Both land the same way.
read_meta <- function(p) {
  d <- read.table(p, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cellID" %in% names(d)) {
    stop("no 'cellID' column found in ", p, " (got: ", paste(utils::head(names(d), 6), collapse = " "), " ...)")
  }
  d
}

## co-embed cellID is the genome-stripped key <barcode>-SM2_<plate> (3_1's strip_genome);
## per-genome v6 cellIDs still carry the mapping suffix. Strip to join.
strip_genome <- function(x) sub("_(Zm_B73v5|AraTAIR10)_scifiATAC$", "", x)

ce <- read_meta(ce_path)
message(" - co-embed meta: ", nrow(ce), " cells x ", ncol(ce), " cols (", paste(names(ce), collapse = " "), ")")

## ---- columns to graft: whatever v6 has that the co-embed lacks --------------
v6B_all <- read_meta(v6_B_path)
v6A_all <- read_meta(v6_A_path)
graft <- setdiff(intersect(names(v6B_all), names(v6A_all)), c("cellID", names(ce)))
if (!length(graft)) stop("nothing to graft -- co-embed meta already carries the v6 schema")
message(" - grafting ", length(graft), " columns: ", paste(graft, collapse = " "))

prep <- function(d, lab) {
  d <- d[, c("cellID", graft), drop = FALSE]
  d$cellID <- strip_genome(d$cellID)
  if (anyDuplicated(d$cellID)) stop("duplicate stripped cellIDs in ", lab, " -- key convention broken")
  message("   * ", lab, ": ", nrow(d), " cells")
  d
}
v6B <- prep(v6B_all, "B73v5")
v6A <- prep(v6A_all, "TAIR10")

## ---- reproduce 3_1's pick() EXACTLY -----------------------------------------
## 3_1 merged the step2 per-genome QC with all.x=TRUE, so its `qc$X.B` is NA precisely when the
## cell was absent from the step2 B73v5 OBJECT -- i.e. `pick(b,a)` == `ifelse(in_B73v5, b, a)`.
## That is NOT the same as "present in the v6 file": 20 Pre cells sit in v6-B73v5 but were dropped
## from the step2 object, and selecting on v6 availability instead of `in_B73v5` picks the wrong
## genome for them (measured: it breaks the nSites/log10nSites identity by up to 1.29 in log10).
## So select on the membership flags the co-embed itself recorded.
iB <- match(ce$cellID, v6B$cellID)
iA <- match(ce$cellID, v6A$cellID)
if ("in_B73v5" %in% names(ce)) {
  sel_B <- as.logical(ce$in_B73v5)
  bad <- (sel_B & is.na(iB)) | (!sel_B & is.na(iA))
  if (any(bad)) stop(sum(bad), " cells have in_B73v5 membership with no matching v6 row -- ",
                     "the v6 files do not correspond to the step2 objects the co-embed was built from")
} else {
  sel_B <- !is.na(iB)
  message(" - !! no in_B73v5 column; falling back to v6 availability. This can pick a different",
          " genome than 3_1 did -- the consistency check below is the guard.")
}
out <- ce
for (cc in graft) out[[cc]] <- ifelse(sel_B, v6B[[cc]][iB], v6A[[cc]][iA])

## ---- coverage: an all-NA join would silently empty every downstream gate ----
src <- ifelse(sel_B, "B73v5", ifelse(!is.na(iA), "TAIR10", "NONE"))
tb  <- table(src)
message(" - QC source per cell: ", paste(names(tb), tb, sep = "=", collapse = "  "))
n_bad <- sum(src == "NONE")
if (n_bad == nrow(ce)) stop("join matched nothing -- check the cellID key convention")
if (n_bad > 0) {
  message(" - !! ", n_bad, " cells matched NEITHER v6 file; their QC is NA and they will drop out",
          " of any QC gate. First few: ", paste(utils::head(ce$cellID[src == "NONE"], 3), collapse = ", "))
}

## ---- the two things that must hold before this is worth shipping -----------
## (a) nSites must be present and consistent with the log10nSites 3_1 already set.
if (all(c("nSites", "log10nSites") %in% names(out))) {
  ok <- !is.na(out$nSites) & !is.na(out$log10nSites) & out$nSites > 0
  d  <- max(abs(log10(out$nSites[ok]) - out$log10nSites[ok]))
  message(" - consistency check nSites vs log10nSites: max|log10(nSites) - log10nSites| = ", signif(d, 3),
          "  (n=", sum(ok), ")")
  ## This is a HARD guard, not a warning. v6 log10nSites is exactly log10(nSites), so if the grafted
  ## nSites came from the same genome 3_1 read log10nSites from, this is float-zero (~5e-15 measured).
  ## Anything larger means the two columns describe different genomes for the same cell -- which is
  ## precisely the bug that selecting on v6 availability instead of in_B73v5 introduces.
  if (d > 1e-6) {
    stop("nSites and log10nSites disagree by ", signif(d, 3), " in log10 -- the grafted QC came from a ",
         "different genome than 3_1's. Do NOT ship this metadata; fix the selector first.")
  }
}
## (b) preview the Fig-1 F/G/H gate so a silently-empty panel is caught here, not in the figure.
if (all(c("total", "pTSS", "FRiP") %in% names(out))) {
  gate <- !is.na(out$total) & out$total >= 500 &
          !is.na(out$pTSS)  & out$pTSS  >= 0.2 &
          !is.na(out$FRiP)  & out$FRiP  >= 0.2
  message(" - Fig-1 F/G/H gate (total>=500, pTSS>=0.2, FRiP>=0.2): ",
          sum(gate), "/", nrow(out), " cells retained")
  if ("Genome" %in% names(out)) {
    gt <- table(out$Genome[gate])
    message("   by plate-of-origin: ", paste(names(gt), gt, sep = "=", collapse = "  "))
  }
  if (sum(gate) == 0) stop("the Fig-1 gate would retain ZERO cells -- do not ship this")
}

write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message(" - DONE. Wrote ", nrow(out), " x ", ncol(out), " -> ", out_path)
