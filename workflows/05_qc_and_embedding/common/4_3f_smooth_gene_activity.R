###################################################################################################
## 4_3f_smooth_gene_activity.R   --   STEP 2 of the marker pipeline (gene activity -> SMOOTHED)
##
## Raw scATAC gene activity is ~0/1 per cell (At Pre: ~2.6M nonzeros over 32k genes x 179k cells),
## far too sparse for per-cell marker scoring. This step diffuses each cell's gene activity over a
## kNN-Markov graph in PC space, producing a dense, denoised "smoothed gene activity" matrix that is
## a REUSABLE artifact -- consumed by 4_3g (per-cell annotation) and later by DAR / marker-UMAPs.
##
## Method (Markov-affinity diffusion, from the maize_282 reference smooth.data, mem-safe as in 4_3b):
##   CP10k per cell  ->  symmetric kNN adjacency in PC space  ->  row-normalized transition A
##   ->  X <- A %*% X iterated `step` hops (never form A^step).
## Smooths ALL expressed genes so the artifact is general (annotation samples its own background from it).
##
## Usage:
##   Rscript 4_3f_smooth_gene_activity.R <out_dir> <prefix> <matrix_rds> <rd_txt> <stage> [k=25] [step=3] [seed=1]
## Output:
##   <out_dir>/<prefix>.genes.smoothed.rds   (genes x clustered-cells, CP10k + diffused; dense)
##   <out_dir>/<prefix>.genes.smoothed.info.tsv  (n_genes, n_cells, k, step, PC dims)
###################################################################################################

options(stringsAsFactors = FALSE)
suppressMessages({ library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  stop("Usage: Rscript 4_3f_smooth_gene_activity.R <out_dir> <prefix> <matrix_rds> <rd_txt> <stage> [k=25] [step=3] [seed=1]")
out_dir <- args[1]; prefix <- args[2]; matrix_rds <- args[3]; rd_txt <- args[4]; stage <- args[5]
k    <- if (length(args) >= 6) as.integer(args[6]) else 25L
step <- if (length(args) >= 7) as.integer(args[7]) else 3L
seed <- if (length(args) >= 8) as.integer(args[8]) else 1L

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
op <- function(s) file.path(out_dir, paste0(prefix, s))
set.seed(seed)
message(" - 4_3f smooth | prefix=", prefix, " stage=", stage, " | k=", k, " step=", step)

## kNN indices incl. self (col 1); RANN/FNN if present, else chunked exact base-R (from 4_3b) --------
get_knn_idx <- function(X, k) {
  if (requireNamespace("RANN", quietly = TRUE)) return(RANN::nn2(X, k = k)$nn.idx)
  if (requireNamespace("FNN",  quietly = TRUE)) return(cbind(seq_len(nrow(X)), FNN::get.knn(X, k = k - 1)$nn.index))
  n <- nrow(X); sq <- rowSums(X^2); idx <- matrix(0L, n, k); chunk <- 1024L
  for (s in seq(1L, n, by = chunk)) {
    e <- min(s + chunk - 1L, n)
    d2 <- outer(sq[s:e], sq, "+") - 2 * (X[s:e, , drop = FALSE] %*% t(X))
    for (r in seq_len(nrow(d2))) idx[s + r - 1L, ] <- order(d2[r, ])[seq_len(k)]
  }
  idx
}

## ---- inputs + align to the CLUSTERED cells (those in the reduced-dims manifold) -----------------
pcs <- read.table(rd_txt, header = TRUE, sep = "\t", quote = "", comment.char = "", row.names = 1)
m   <- readRDS(matrix_rds)
cells <- intersect(colnames(m), rownames(pcs))
if (length(cells) == 0) stop("no common cells between matrix and reduced_dims (check naming)")
pcs <- as.matrix(pcs[cells, , drop = FALSE])
m   <- m[Matrix::rowSums(m[, cells, drop = FALSE]) > 0, cells, drop = FALSE]   # drop all-zero genes
message(" - genes x cells (expressed, clustered): ", nrow(m), " x ", ncol(m), " | PC dims: ", ncol(pcs))

## ---- per-cell CP10k (all genes) ----------------------------------------------------------------
libcell <- Matrix::colSums(m); libcell[libcell == 0] <- 1
norm <- m %*% Matrix::Diagonal(x = 1 / libcell) * 1e4          # sparse genes x cells

## ---- kNN-Markov graph + iterative diffusion ----------------------------------------------------
message(" - building kNN graph (k=", k, ") ...")
nn <- get_knn_idx(pcs, k); n <- nrow(pcs)
A  <- sparseMatrix(i = rep(seq_len(n), each = ncol(nn)), j = as.vector(t(nn)), x = 1, dims = c(n, n))
A  <- A + Matrix::t(A)
A  <- Matrix::Diagonal(x = 1 / Matrix::rowSums(A)) %*% A       # row-normalized transition
message(" - diffusing ", nrow(norm), " genes x ", n, " cells, step=", step, " (dense ~",
        round(nrow(norm) * n * 8 / 1e9, 2), " GB) ...")
X <- as.matrix(Matrix::t(norm))                               # cells x genes (dense)
for (s in seq_len(step)) { X <- as.matrix(A %*% X); message("   hop ", s, "/", step, " done") }
sm <- t(X)                                                    # genes x cells
rownames(sm) <- rownames(norm); colnames(sm) <- cells

## ---- save --------------------------------------------------------------------------------------
saveRDS(sm, op(".genes.smoothed.rds"))
write.table(data.frame(prefix = prefix, stage = stage, n_genes = nrow(sm), n_cells = ncol(sm),
                       k = k, step = step, pc_dims = ncol(pcs)),
            op(".genes.smoothed.info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(" - wrote ", op(".genes.smoothed.rds"), " (", nrow(sm), " x ", ncol(sm), ")")
message(" - DONE ", prefix)
