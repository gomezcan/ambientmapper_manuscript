###################################################################################################
## 5_1b_export_for_seacells.R  --  Step 5.1b: portable SEACells bundle
##
## Writes the ACR matrix + embedding + cell metadata in a format any Python environment can load,
## WITHOUT requiring anndata/scanpy on the machine that builds it. The R side (this project) and
## the Python side (SEACells) are therefore decoupled -- the local mambaforge R is bare and has no
## anndata, and forcing an h5ad here would make the export depend on the analysis environment.
##
## ORIENTATION: the .mtx is written CELLS x FEATURES, i.e. AnnData convention (obs = cells,
## var = ACRs). No transpose is needed on the Python side.
##
## SEACells needs exactly three things (Persad et al. 2023):
##   (1) raw count matrix          -> <p>.mtx.gz          [cells x ACRs, raw Tn5 insertion counts]
##   (2) low-dim representation    -> <p>.svd.tsv         [-> obsm['X_svd']]
##   (3) n_SEACells                -> caller's choice (5_3_run_seacells_plate.sh: CELLS_PER_SEACELL / N_SEACELLS)
## plus, for evaluation:
##       cell metadata             -> <p>.obs.tsv         [LouvainClusters + QC; purity needs a label]
##       feature names             -> <p>.var.tsv
##
## Purity must be fed CLUSTER ID, not cell type -- the annotation is the open question this layer
## exists to test, so using it as the evaluation label would be circular.
##
## Usage:
##   Rscript 5_1b_export_for_seacells.R <acrs.sparse.rds> <reduced_dims.txt> <metadata.txt> <out_prefix>
###################################################################################################

suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: Rscript 5_1b_export_for_seacells.R <acrs.sparse.rds> <reduced_dims.txt> <metadata.txt> <out_prefix>")
mat_f <- args[1]; rd_f <- args[2]; md_f <- args[3]; outp <- args[4]

m <- readRDS(mat_f)                                   # ACRs x cells
X <- as.matrix(read.table(rd_f, header = TRUE, sep = "\t", row.names = 1,
                          check.names = FALSE, comment.char = ""))
mode(X) <- "numeric"
MD <- read.table(md_f, header = TRUE, sep = "\t", row.names = 1,
                 check.names = FALSE, comment.char = "", quote = "")

## align all three on the intersection, in ONE canonical order
cells <- Reduce(intersect, list(colnames(m), rownames(X), rownames(MD)))
if (length(cells) < 50) stop("only ", length(cells), " cells shared across matrix/embedding/metadata")
if (length(cells) != ncol(m))
  message(" ! dropping ", ncol(m) - length(cells), " matrix cells absent from embedding/metadata")
m <- m[, cells, drop = FALSE]; X <- X[cells, , drop = FALSE]; MD <- MD[cells, , drop = FALSE]

message(" - ", basename(mat_f))
message("   ", nrow(m), " ACRs x ", length(cells), " cells | ", ncol(X), " SVD dims")

## (1) counts, CELLS x FEATURES (AnnData orientation)
mm <- as(t(m), "CsparseMatrix")
tmp <- paste0(outp, ".mtx")
invisible(writeMM(mm, tmp))
system2("gzip", c("-f", shQuote(tmp)))
message(" - wrote ", basename(tmp), ".gz  (cells x ACRs, ", format(length(mm@x), big.mark = ","), " nnz)")

## (2) embedding -> obsm['X_svd']  (row order == obs order)
write.table(data.frame(cellID = cells, X, check.names = FALSE),
            paste0(outp, ".svd.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## (3) obs / var
keep <- intersect(c("LouvainClusters", "total", "nSites", "FRiP", "pTSS", "pOrg", "doubletscore"),
                  colnames(MD))
obs <- data.frame(cellID = cells, MD[, keep, drop = FALSE], check.names = FALSE)
write.table(obs, paste0(outp, ".obs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(acrID = rownames(m)), paste0(outp, ".var.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote .svd.tsv / .obs.tsv / .var.tsv   (obs cols: ", paste(keep, collapse = ", "), ")")

## depth per cell -- the rarefaction target. SEACells will NOT do this: it takes ONE AnnData and
## cannot know Pre/wd/nd are meant to be compared, so an unmatched stage gets deeper meta-cells free.
d <- Matrix::colSums(m)
message(sprintf(" - counts/cell: median %.0f  mean %.0f  (RAREFACTION TARGET -- match across stages)",
                median(d), mean(d)))
