## 1_1c_merge_and_qc.R — merge per-chunk soc.obj files, then run isCellv2 once.
##
## Reads <chunkdir>/<pool>.chunk0..(N-1).soc.rds (from 1_1b_per_chunk.R),
## merges them with Socrates::mergeSocratesRDS, runs isCellv2 on the combined
## soc.obj (global z-scores across all cells), and writes <out>.raw.soc.rds —
## drop-in replacement for the unchunked 1_1 output. Downstream 1_2 / 1_3
## consume it unchanged.
##
## Safety: before merging, asserts that chunks have pairwise-disjoint barcode
## sets (the hash partition in 1_1a guarantees this; the check is cheap
## insurance against a partition bug or stray duplicate barcode).
##
## Usage:
##   Rscript 1_1c_merge_and_qc.R <pool> <chunkdir> <out_prefix> [N=5]

suppressMessages(library(Socrates))
suppressMessages(library(Matrix))

## isCellv2() — copied verbatim from 1_1_QC_scifiATAC_data.R. Computes
## pTSS / FRiP / pOrg and their z-scores, then qc_check. Z-scores require
## global stats across all cells, which is why this step runs once on the
## merged soc.obj rather than per chunk.
isCellv2 <- function(obj, num.test = 20000, num.tn5 = NULL, num.ref = 1000,
                     background.cutoff = 100, min.pTSS = 0.2, min.FRiP = 0.2,
                     min.pTSS.z = -2, min.FRiP.z = -2, verbose = T)
{
  .RowVar <- function(x) {
    spm <- t(x)
    if (!methods::is(spm, "dgCMatrix")) {
      stop("Error: Input is not a 'dgCMatrix' sparse matrix.")
    }
    ans <- sapply(base::seq.int(spm@Dim[2]), function(j) {
      if (spm@p[j + 1] == spm@p[j]) {
        return(0)
      }
      mean <- base::sum(spm@x[(spm@p[j] + 1):spm@p[j + 1]]) / spm@Dim[1]
      sum((spm@x[(spm@p[j] + 1):spm@p[j + 1]] - mean)^2) +
        mean^2 * (spm@Dim[1] - (spm@p[j + 1] - spm@p[j]))
    }) / (spm@Dim[1] - 1)
    names(ans) <- spm@Dimnames[[2]]
    ans
  }

  if (verbose) message("Step 1: Converting count matrix to sparseMatrix format")
  tryCatch({
    if (is.data.frame(obj$counts)) {
      sparse_count_matrix <- obj$counts
      if (verbose) message(" - converting triplet format to sparseMatrix")
      sparse_count_matrix$V1 <- factor(sparse_count_matrix$V1)
      sparse_count_matrix$V2 <- factor(sparse_count_matrix$V2)
      sparse_count_matrix <- Matrix::sparseMatrix(
        i = as.numeric(sparse_count_matrix$V1),
        j = as.numeric(sparse_count_matrix$V2),
        x = as.numeric(sparse_count_matrix$V3),
        dimnames = list(levels(sparse_count_matrix$V1), levels(sparse_count_matrix$V2))
      )
    } else if (methods::is(obj$counts, "dgTMatrix")) {
      if (verbose) message(" - converting 'dgTMatrix' to 'dgCMatrix'")
      sparse_count_matrix <- as(obj$counts, "dgCMatrix")
    } else if (methods::is(obj$counts, "dgCMatrix")) {
      if (verbose) message(" - count matrix is already a 'dgCMatrix'")
      sparse_count_matrix <- obj$counts
    } else if (is.matrix(obj$counts)) {
      if (verbose) message(" - converting dense matrix to 'dgCMatrix'")
      sparse_count_matrix <- Matrix::Matrix(obj$counts, sparse = TRUE)
    } else {
      stop("Error: Unsupported format for count matrix.")
    }
  }, error = function(e) stop("Error during sparse matrix conversion: ", e$message))

  if (verbose) message("Step 2: Aligning metadata with count matrix")
  shared <- intersect(rownames(obj$meta), colnames(sparse_count_matrix))
  if (length(shared) == 0) stop("Error: No shared identifiers between metadata and count matrix.")
  sparse_count_matrix <- sparse_count_matrix[, shared]
  obj$meta <- obj$meta[shared, ]
  obj$meta <- obj$meta[order(obj$meta$nSites, decreasing = TRUE), ]

  if (verbose) message("Step 3: Calculating quality control metrics")
  obj$meta$pTSS <- obj$meta$tss  / obj$meta$total
  obj$meta$FRiP <- obj$meta$acrs / obj$meta$total
  obj$meta$pOrg <- obj$meta$ptmt / obj$meta$total

  if (verbose) message("Step 4: Setting quality control filters")
  obj$meta$tss_z   <- (obj$meta$pTSS - mean(obj$meta$pTSS)) / sd(obj$meta$pTSS)
  obj$meta$acr_z   <- (obj$meta$FRiP - mean(obj$meta$FRiP)) / sd(obj$meta$FRiP)
  obj$meta$sites_z <- (log10(obj$meta$nSites) - mean(log10(obj$meta$nSites))) / sd(log10(obj$meta$nSites))
  obj$meta$tss_z[is.na(obj$meta$tss_z)]     <- -10
  obj$meta$acr_z[is.na(obj$meta$acr_z)]     <- -10
  obj$meta$sites_z[is.na(obj$meta$sites_z)] <- -10
  obj$meta$qc_check <- ifelse(obj$meta$tss_z < min.pTSS.z | obj$meta$pTSS < min.pTSS, 0,
                              ifelse(obj$meta$acr_z < min.FRiP.z | obj$meta$FRiP < min.FRiP, 0, 1))
  return(obj)
}

# arguments
args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Rscript 1_1c_merge_and_qc.R <pool> <chunkdir> <out_prefix> [N=5]")
}
pool       <- as.character(args[1])
chunkdir   <- as.character(args[2])
out_prefix <- as.character(args[3])
N          <- if (length(args) >= 4) as.integer(args[4]) else 5L

# locate per-chunk soc.obj files (zero-indexed chunk0..chunk{N-1})
chunk_files <- file.path(chunkdir, sprintf("%s.chunk%d.soc.rds", pool, seq_len(N) - 1L))
missing <- chunk_files[!file.exists(chunk_files)]
if (length(missing) > 0) {
  stop("Missing chunk files:\n  ", paste(missing, collapse = "\n  "))
}

# Pairwise barcode disjointness check. mergeSocratesRDS factor()s barcodes
# and rebuilds the sparseMatrix from triplets; a duplicate barcode across
# chunks would silently sum its counts AND create duplicate meta rownames.
# Hash partition in 1_1a should guarantee disjointness; this is the assert.
message(" - checking pairwise barcode disjointness across ", N, " chunks ...")
chunk_bcs <- lapply(chunk_files, function(f) colnames(readRDS(f)$counts))
for (i in seq_along(chunk_bcs)) {
  for (j in seq_len(i - 1L)) {
    ov <- intersect(chunk_bcs[[i]], chunk_bcs[[j]])
    if (length(ov) > 0) {
      stop(sprintf("Chunks %d and %d share %d barcodes (first: %s) - partition is broken",
                   i - 1L, j - 1L, length(ov), ov[1]))
    }
  }
}
total_bcs <- sum(lengths(chunk_bcs))
message(" - OK: chunks are disjoint. Total barcodes across chunks: ", total_bcs)
rm(chunk_bcs); invisible(gc())

# Merge per-chunk soc.obj files.
# Quirk: mergeSocratesRDS's filenames= branch indexes all.rds[[x]] where x is
# the lapply value (the path) but the list names come from names(filenames).
# So names must equal values to avoid all.rds[[x]] being NULL.
message(" - merging with mergeSocratesRDS()")
fn <- setNames(chunk_files, chunk_files)
soc.obj <- mergeSocratesRDS(filenames = fn)

# Collapse per-chunk sampleID back to pool name.
if ("sampleID" %in% colnames(soc.obj$meta)) {
  soc.obj$meta$sampleID <- pool
}
message(" - merged soc.obj")
message("   cells    : ", nrow(soc.obj$meta))
message("   features : ", nrow(soc.obj$counts))

# Global z-score QC.
soc.obj <- isCellv2(soc.obj)

# Drop-in for 1_2 / 1_3.
saveRDS(soc.obj, file = paste0(out_prefix, ".raw.soc.rds"))
message(" - wrote ", out_prefix, ".raw.soc.rds")
