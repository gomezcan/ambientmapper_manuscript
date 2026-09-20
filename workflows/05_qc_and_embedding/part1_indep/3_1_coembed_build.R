#!/usr/bin/env Rscript
## ============================================================================
## 3_1_coembed_build.R  --  MULTI-REFERENCE co-embedding by feature concatenation
##   (part1_indep: the individual-mapping construction of Fig 1F/G/H + Fig S3)
## ----------------------------------------------------------------------------
## Builds a Fig-1-style cross-species co-projection from the TWO independent
## per-genome objects, WITHOUT combined-genome mapping. Every barcode in the indep
## arm was mapped to BOTH references, so we "extend the row of the matrix": stack the
## two per-genome 500bp-tile matrices (union of cells, zero-fill the missing genome
## block) into ONE joint object. A maize cell loads on the B73v5 block, an At cell on
## the TAIR10 block, contamination on both — the same wide matrix Fig 1's combined
## ZmATcombined object has, but assembled from the individual mappings the method
## (AmbientMapper, Fig 3) actually uses. This is the PRIMARY Fig 1F construction;
## the combined ZmATcombined object becomes the Fig-S1 "robust to mapping strategy" check.
##
## Mixing LABEL = plate-of-origin (design ground truth) -> $meta$Genome. So the joint
## object's genome_mixing (computed by 3_2_coembed_gridscan.R) is species/plate mixing,
## directly comparable to the combined Fig 1F / Fig S3 — NOT the per-genome plate-library
## mixing of the split-arm scan.
##
## Feature names are genome-PREFIXED ("B73v5:1_...", "TAIR10:1_...") so identical tile
## coordinates on the two genomes' chr1 do not collide.
##
## INPUT = the step2 per-genome objects (cleaned tile $counts + $Clusters = the QC-passed
## v7 cells + QC covariates). Using the cleaned tiles is EQUIVALENT to raw tiles here:
## min.t counts a B73v5 tile only over B73v5-block cells (At-only cells are zero), so
## per-genome and joint feature-filtering agree — and it is ~10x lighter. min.c is
## (re)applied on the JOINT by 3_2.
##
## Usage:
##   Rscript 3_1_coembed_build.R <b73v5_step2_rds> <tair10_step2_rds> <out_prefix> [stage]
## Emits <out_prefix>.coembed.soc.rds and <out_prefix>.coembed.meta.tsv
## (the .rds + .tsv pair that 3_2_coembed_gridscan.R consumes).
##
## Only the Matrix package is needed HERE; 3_2/3_3 need Socrates.
## ============================================================================
suppressWarnings(suppressMessages(library(Matrix)))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("Usage: Rscript 3_1_coembed_build.R <b73v5_step2_rds> <tair10_step2_rds> <out_prefix> [stage]")
b73_rds <- a[1]; at_rds <- a[2]; out_prefix <- a[3]
stage   <- if (length(a) >= 4) a[4] else "NA"

## --- barcode key helpers (strip the reference-genome tag from the cellID) ----
## cellID = <barcode>-SM2_<plate>_<Zm_B73v5|AraTAIR10>_scifiATAC
strip_genome <- function(x) sub("_(Zm_B73v5|AraTAIR10)_scifiATAC$", "", x)      # -> <barcode>-SM2_<plate>
plate_of     <- function(k) ifelse(grepl("-SM2_At$",  k), "At",
                            ifelse(grepl("-SM2_B73$", k), "B73", NA_character_))

## --- load one per-genome object: cleaned tile matrix restricted to v7 cells ---
load_block <- function(rds, tag) {
  o <- readRDS(rds)
  m <- o$counts                                         # tiles x cells (dgCMatrix)
  ## restrict to the Fig-5 v7 clustered (QC-passed) cells for cross-arm consistency
  keep <- if (!is.null(o$Clusters)) o$Clusters$cellID else colnames(m)
  m <- m[, colnames(m) %in% keep, drop = FALSE]
  ## QC covariates keyed by genome-stripped barcode (for the qc_leak axis in 1_2)
  md <- if (!is.null(o$Clusters)) o$Clusters else o$meta
  porg <- if ("pOrg" %in% names(md)) md$pOrg else md$ptmt      # real pOrg if present, else ptmt proxy
  qc <- data.frame(key = strip_genome(as.character(md$cellID)),
                   log10nSites = suppressWarnings(as.numeric(md$log10nSites)),
                   pOrg        = suppressWarnings(as.numeric(porg)),
                   stringsAsFactors = FALSE)
  qc <- qc[!duplicated(qc$key), ]
  colnames(m) <- strip_genome(colnames(m))              # genome-stripped barcode key
  m <- m[, !duplicated(colnames(m)), drop = FALSE]
  rownames(m) <- paste0(tag, ":", rownames(m))          # genome-prefixed features
  message(sprintf("  %s: %d tiles x %d cells (restricted to v7)", tag, nrow(m), ncol(m)))
  list(m = m, qc = qc)
}

message("Reading B73v5: ", b73_rds)
B <- load_block(b73_rds, "B73v5")
message("Reading TAIR10: ", at_rds)
A <- load_block(at_rds,  "TAIR10")

## --- union cell set; reindex each block to it (zero-fill missing genome block) -
cells <- union(colnames(B$m), colnames(A$m))
expand <- function(m, cells) {
  miss <- setdiff(cells, colnames(m))
  if (length(miss)) {
    z <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                      dims = c(nrow(m), length(miss)),
                      dimnames = list(rownames(m), miss))
    m <- cbind(m, z)
  }
  m[, cells, drop = FALSE]
}
joint <- rbind(expand(B$m, cells), expand(A$m, cells))   # [B73v5 tiles ; TAIR10 tiles] x union cells
stopifnot(identical(colnames(joint), cells))

## --- joint meta: plate=Genome ground-truth label + QC (prefer own/major genome) --
inB <- cells %in% colnames(B$m); inA <- cells %in% colnames(A$m)
qc  <- merge(data.frame(key = cells, stringsAsFactors = FALSE),
             B$qc, by = "key", all.x = TRUE, sort = FALSE)          # B73v5 QC first
qc  <- merge(qc, A$qc, by = "key", all.x = TRUE, sort = FALSE, suffixes = c(".B", ".A"))
qc  <- qc[match(cells, qc$key), ]
pick <- function(b, a) ifelse(!is.na(b), b, a)                      # prefer B73v5, else At
meta <- data.frame(
  cellID      = cells,
  Genome      = plate_of(cells),                                   # design ground truth (At / B73)
  species     = plate_of(cells),
  stage       = stage,
  in_B73v5    = inB, in_TAIR10 = inA, both = inB & inA,            # barnyard membership
  log10nSites = pick(qc$log10nSites.B, qc$log10nSites.A),
  pOrg        = pick(qc$pOrg.B, qc$pOrg.A),
  joint_depth = Matrix::colSums(joint),
  row.names   = cells, stringsAsFactors = FALSE)

## --- save as a soc-object list ($counts + $meta); 1_2 recomputes the rest --------
soc <- list(counts = joint, meta = meta)
saveRDS(soc, paste0(out_prefix, ".coembed.soc.rds"))
write.table(meta, paste0(out_prefix, ".coembed.meta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf(
  "\nco-embed [%s]: features=%d (B73v5 %d + TAIR10 %d) | cells=%d",
  stage, nrow(joint), nrow(B$m), nrow(A$m), length(cells)))
message(sprintf(
  "  membership: B73v5-only=%d  TAIR10-only=%d  BOTH(barnyard)=%d",
  sum(inB & !inA), sum(!inB & inA), sum(inB & inA)))
message(sprintf(
  "  plate(Genome) label: At=%d  B73=%d  (NA=%d)",
  sum(meta$Genome == "At", na.rm = TRUE), sum(meta$Genome == "B73", na.rm = TRUE),
  sum(is.na(meta$Genome))))
message("Wrote: ", out_prefix, ".coembed.soc.rds  +  .coembed.meta.tsv")
