###################################################################################################
## 4_0_build_gene_body_matrix.R
##
## Build a gene-body accessibility matrix (gene x cell, raw Tn5 insertion counts) from one or more
## Tn5 BED files + the combined-genome GTF. This is the "gene activity" matrix consumed by the
## cluster-annotation steps; it mirrors the maize_282 reference all_pools.raw.genes.sparse.rds
## (rows = geneID, cols = cellID, values = Tn5 insertions overlapping the gene window).
##
## Socrates' generateMatrix() only builds peak- or fixed-bp-window matrices, so gene-body counts are
## built directly here, reusing the same findOverlaps(Tn5, features) logic the QC step (1_1) uses for
## ACRs -- only the feature set changes (gene windows instead of MACS peaks).
##
## Usage:
##   Rscript 4_0_build_gene_body_matrix.R <out_prefix> <gtf> <up_bp> <down_bp> <bed1> [bed2 ...]
##     out_prefix : matrix written to <out_prefix>.genes.sparse.rds
##     gtf        : GTF (gene_id dbxref) OR GFF3 -- routed by file extension. Combined Part-1 uses
##                  ZmATcombined.gtf; per-genome Part-2 uses B73v5.gff3 (NAM5) / TAIR10.gff3 (Ensembl r60).
##     up_bp      : bp added on the 5' (promoter) side of each gene body, strand-aware
##     down_bp    : bp added on the 3' side of each gene body, strand-aware
##     bed*       : Tn5 BED(s) (chr, start, end, cellID, strand); cells from all BEDs are pooled (cbind)
##
## Env:
##   GENE_LENGTH_NORM = raw (default) | perkb
##     raw   : raw Tn5 insertion counts per gene window (legacy; length-biased -- large genes like
##             kn1/knox accumulate more insertions at equal accessibility density).
##     perkb : divide each gene by its counted-window length -> insertions per kb. Removes the length
##             bias that lets a few large marker genes dominate the per-cell 4_3g enrichment classifier.
##     The output FILENAME is unchanged (driven by <out_prefix>); the caller encodes the mode in it.
##
## Cells are species-tagged in the BED V4 cellID (...-SM2_At_... / ...-SM2_B73_...), which matches the
## merged Socrates metadata cellID exactly -> the matrix aligns to the clustered metadata with no
## remapping. Genes are species-separable by seqname prefix (Zm_* = maize B73v5, At_* = Arabidopsis).
###################################################################################################

suppressMessages({
  library(data.table)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 4_0_build_gene_body_matrix.R <out_prefix> <gtf> <up_bp> <down_bp> <bed1> [bed2 ...]")
}
out_prefix <- args[1]
gtf        <- args[2]
up_bp      <- as.integer(args[3])
down_bp    <- as.integer(args[4])
beds       <- args[-(1:4)]
norm_mode  <- tolower(Sys.getenv("GENE_LENGTH_NORM", "raw"))     # raw (default) | perkb
if (!norm_mode %in% c("raw", "perkb"))
  stop("GENE_LENGTH_NORM must be 'raw' or 'perkb' (got '", norm_mode, "')")

message(" - out_prefix = ", out_prefix)
message(" - gtf        = ", gtf)
message(" - window     = +", up_bp, " bp (5', promoter) / +", down_bp, " bp (3'), strand-aware")
message(" - beds       = ", paste(basename(beds), collapse = ", "))
message(" - length norm= ", norm_mode, if (norm_mode == "perkb") "  (insertions per kb)" else "  (raw counts)")

## (1) gene-body windows from the annotation -------------------------------------
## Route by extension: combined Part-1 run uses a .gtf (gene_id dbxref); per-genome Part-2 uses
## .gff3 (B73v5 -> NAM5, TAIR10 -> Ensembl r60). The maize .gtf variants carry GFF3 keys in GTF
## syntax with NO gene_id, so per-genome MUST use the .gff3.
fmt <- if (grepl("\\.gff3?(\\.gz)?$", gtf, ignore.case = TRUE)) "gff3" else "gtf"
message(" - building TxDb from ", toupper(fmt), " (may take a few minutes) ...")
txdb <- if (fmt == "gff3") {
  suppressWarnings(makeTxDbFromGFF(gtf, format = "gff3"))
} else {
  suppressWarnings(makeTxDbFromGFF(gtf, format = "gtf", dbxrefTag = "gene_id"))
}
g <- genes(txdb)                 # one strand-aware range per gene_id (gene body)
names(g) <- sub("^gene:", "", names(g))   # Ensembl GFF3 prefixes gene IDs 'gene:AT1G...' -> strip to match panel/matrix
gene_levels <- names(g)
message("   genes = ", length(g))

## strand-aware extension: up_bp added on the 5' side, down_bp on the 3' side
is.minus  <- as.character(strand(g)) == "-"
left.ext  <- ifelse(is.minus, down_bp, up_bp)   # bp added on the lower-coordinate side
right.ext <- ifelse(is.minus, up_bp, down_bp)   # bp added on the higher-coordinate side
gext <- GRanges(
  seqnames = seqnames(g),
  ranges   = IRanges(start = pmax(1L, start(g) - left.ext), end = end(g) + right.ext),
  gene_id  = gene_levels
)

## (2) per-BED gene x cell counts, then cbind ------------------------------------
## Process one BED at a time so peak memory is bounded by the largest single BED (B73 ~47M rows),
## not the sum. Cells are disjoint across BEDs (species-tagged), so a plain cbind is correct.
mats <- vector("list", length(beds))
for (k in seq_along(beds)) {
  bf <- beds[k]
  message(" - [", k, "/", length(beds), "] reading ", basename(bf), " ...")
  bed <- fread(cmd = paste("gzip -dc", shQuote(bf)), header = FALSE, sep = "\t",
               select = 1:4, col.names = c("chr", "start", "end", "cell"))
  message("   insertions = ", nrow(bed))

  gr   <- GRanges(bed$chr, IRanges(bed$start + 1L, bed$end))   # BED 0-based -> 1-based 1bp Tn5 site
  hits <- findOverlaps(gr, gext, ignore.strand = TRUE)         # a site in N overlapping windows counts N times
  message("   genic insertions = ", length(hits),
          " (", round(100 * length(hits) / nrow(bed), 1), "%)")

  cf <- factor(bed$cell[queryHits(hits)])
  ## sparseMatrix() sums x over duplicated (i,j) pairs -> per-(gene,cell) insertion counts
  sub <- sparseMatrix(
    i = subjectHits(hits),
    j = as.integer(cf),
    x = 1,
    dims     = c(length(gext), nlevels(cf)),
    dimnames = list(gene_levels, levels(cf))
  )
  message("   ", basename(bf), ": ", ncol(sub), " cells with >=1 genic insertion")
  mats[[k]] <- sub
  rm(bed, gr, hits, cf, sub); gc()
}

## (3) combine ------------------------------------------------------------------
m <- if (length(mats) == 1) mats[[1]] else do.call(cbind, mats)   # identical gene rows across BEDs
m <- m[Matrix::rowSums(m) > 0, , drop = FALSE]
message(" - final matrix: ", nrow(m), " genes x ", ncol(m), " cells (nnz = ", length(m@x), ")")

## (3b) optional gene-length normalization ---------------------------------------
## Under a null of uniform accessibility, expected Tn5 insertions in a gene window scale with its
## length, so dividing by window length removes the bias that lets a few LARGE marker genes
## (kn1 ~8 kb, knox8 ~9 kb) inflate the per-cell 4_3g enrichment z in every cluster. Near no-op for
## the cluster-level 4_3 (its per-gene z across clusters already cancels any constant per-gene
## factor); the fix targets the absolute-vs-random-background per-cell classifier.
if (norm_mode == "perkb") {
  rn <- rownames(m); cn <- colnames(m)
  wlen <- setNames(width(gext), gext$gene_id)          # full counted-window length per gene (bp)
  wkb  <- wlen[rn] / 1000
  if (any(is.na(wkb) | wkb <= 0))
    stop("perkb: ", sum(is.na(wkb) | wkb <= 0), " genes have missing/non-positive window length")
  m <- Matrix::Diagonal(x = 1 / wkb) %*% m             # scale each gene row -> insertions per kb
  m <- as(m, "CsparseMatrix"); dimnames(m) <- list(rn, cn)
  message(" - perkb: normalized by window length (median ", round(median(wlen) / 1000, 2),
          " kb); values are now insertions/kb")
}

## (4) save ----------------------------------------------------------------------
saveRDS(m, file = paste0(out_prefix, ".genes.sparse.rds"))
message(" - wrote ", out_prefix, ".genes.sparse.rds  (norm=", norm_mode, ")")
