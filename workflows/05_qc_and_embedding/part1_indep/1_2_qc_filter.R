###############################################################################
## 1_2_qc_filter.R -- Step 1_2: per-barcode QC statistics and reference-based cell call
## Reads the 1_1 raw Socrates object, computes pTSS / FRiP / pOrg and their z-scores over
## barcodes with total >= 100 reads, builds a background (bREF) and a good-cell (gREF)
## reference profile from the count matrix and calls a barcode a cell when its profile
## difference (dif = gREF - bREF) exceeds the background median. Output:
## <prefix>.updated_metadata.txt (all barcodes, with qc_check, dif and call columns).
## Usage: Rscript 1_2_qc_filter.R <prefix>.raw.soc.rds <prefix>
## Called by part1_indep/1_0_qc_run.sh (identical to common/1_2_filter_lowQC_cells_scifiATAC_data.R).
###############################################################################
# Load required libraries
suppressMessages(library(devtools))     # For package development utilities
suppressMessages(library(Seurat))       # For single-cell analysis
suppressMessages(library(Socrates))     # Specific library for single-cell ATAC-seq processing
suppressMessages(library(qlcMatrix))    # For sparse matrix operations

rm(list = ls())

# Custom TF-IDF function
tfidf <- function(obj,
                  frequencies = TRUE,
                  log_scale_tf = TRUE,
                  scale_factor = 10000,
                  doL2 = FALSE,
                  slotName = "residuals") {
  
  bmat <- obj$counts
  
  .safe_tfidf <- function(tf, idf, block_size = 2000e6) {
    tryCatch({
      tf * idf
    }, error = function(e) {
      options(DelayedArray.block.size = block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      tf <- DelayedArray(tf)
      idf <- as.matrix(idf)
      tf * idf
    })
  }
  
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    tf <- t(t(bmat) / Matrix::colSums(bmat))
  } else {
    tf <- bmat
  }
  
  # Log scale term frequency if requested
  if (log_scale_tf) {
    tf@x <- log1p(tf@x * scale_factor)
  }
  
  # Inverse document frequency
  idf <- log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # Compute TF-IDF
  tf_idf_counts <- .safe_tfidf(tf, idf)
  
  # Apply L2 normalization if requested
  if (doL2) {
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x = 1 / colNorm)
  }
  
  # Save results to object
  rownames(tf_idf_counts) <- rownames(bmat)
  colnames(tf_idf_counts) <- colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse = TRUE)
  obj$norm_method <- "tfidf"
  return(obj)
}


# function
RowVar <- function(x) {
  spm <- t(x)
  stopifnot(methods::is(spm, "dgCMatrix"))
  ans <- sapply(base::seq.int(spm@Dim[2]), function(j) {
    if (spm@p[j + 1] == spm@p[j]) {
      return(0)
    }
    mean <- base::sum(spm@x[(spm@p[j] + 1):spm@p[j +
                                                   1]])/spm@Dim[1]
    sum((spm@x[(spm@p[j] + 1):spm@p[j + 1]] - mean)^2) +
      mean^2 * (spm@Dim[1] - (spm@p[j + 1] - spm@p[j]))
  })/(spm@Dim[1] - 1)
  names(ans) <- spm@Dimnames[[2]]
  ans
}

# Load command-line arguments
args <- commandArgs(TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 1_2_qc_filter.R <input.rds> <prefix>")
}

input.dat <- as.character(args[1])
prefix <- as.character(args[2])

# Parameters

total_reads <- 100
pTSS <- 0.2
FRiP <- 0.2
max_sites <- 250000
max_cells <- 250000

# Load data
message(" - Loading data")
a <- readRDS(input.dat)

# Align metadata and counts
shared <- intersect(rownames(a$meta), colnames(a$counts))
a$counts <- a$counts[, shared]
a$meta <- a$meta[shared, ]
nrow(a$meta)

# Quality control, and Generate stats
message(" - Performing initial quality control")
a$meta <- a$meta[order(a$meta$nSites, decreasing = TRUE),]
a$meta$pTSS <- a$meta$tss/a$meta$total
a$meta$FRiP <- a$meta$acrs/a$meta$total
a$meta$pOrg <- a$meta$ptmt/a$meta$total

# update qc_check based on total_reads
table(a$meta$qc_check)
a$meta$qc_check <- ifelse(a$meta$total >= total_reads,1, 0)
table(a$meta$qc_check)

# set initial thresholds
message(" - setting filters")
a$meta$tss_z <- (a$meta$pTSS - mean(a$meta$pTSS[a$meta$qc_check==1]))/sd(a$meta$pTSS[a$meta$qc_check==1])
a$meta$acr_z <- (a$meta$FRiP - mean(a$meta$FRiP[a$meta$qc_check==1]))/sd(a$meta$FRiP[a$meta$qc_check==1])
a$meta$sites_z <- (log10(a$meta$nSites) - mean(log10(a$meta$nSites[a$meta$qc_check==1])))/sd(log10(a$meta$nSites[a$meta$qc_check==1]))

a$meta$tss_z[is.na(a$meta$tss_z)] <- -8
a$meta$acr_z[is.na(a$meta$acr_z)] <- -8
a$meta$sites_z[is.na(a$meta$sites_z)] <- -8 
a$meta$qc_check <- ifelse(a$meta$tss_z < -2 | a$meta$pTSS < 0.2, 0, 
                          ifelse(a$meta$acr_z < -2 | a$meta$FRiP < 0.2, 0,
                                 ifelse(a$meta$sites_z < -2, 0, a$meta$qc_check)))

table(a$meta$qc_check)

# Separate good and bad cells
message(" - Separating good and bad cells")
good.cells <- rownames(subset(a$meta, a$meta$qc_check==1))


gg <- a$counts[,colnames(a$counts) %in% good.cells]
gg <- gg[,Matrix::colSums(gg > 0) > 100]

a$meta$qc_check <- ifelse(rownames(a$meta) %in% colnames(gg), 1, 0)
bb <- a$counts[,! colnames(a$counts) %in% colnames(gg)]


# Filter for top sites
sites <- Matrix::rowMeans(gg > 0)
sites <- sites[order(sites, decreasing=T)]
num.sites <- min(c(max(a$meta$nSites), max_sites))

if(length(sites) < num.sites){
  num.sites <- length(sites)
}

gg <- gg[names(sites)[1:num.sites],]
gg <- gg[,Matrix::colSums(gg) > 0]

# Align bad cell matrix
bb <- bb[rownames(gg),]
bb <- bb[,Matrix::colSums(bb) > 0]

shared <- intersect(rownames(gg), rownames(bb))
gg <- gg[shared,]
bb <- bb[shared,]
gg <- gg[,Matrix::colSums(gg) > 0]
bb <- bb[,Matrix::colSums(bb) > 0]

# Do not use more than 250,000 'bad' cells
if(ncol(bb) > 250000){
  top <- Matrix::colSums(bb)
  top <- top[order(top, decreasing=T)]
  top <- names(top)[5001:250000] # skip cells at the boundary
}else{
  num.cells <- ncol(bb)
  top <- Matrix::colSums(bb)
  top <- top[order(top, decreasing=T)]
  if(num.cells > 100000){
    top <- names(top)[5001:100000]
  }else{
    top <- names(top)
  }
}


# clean up ref
bb <- bb[,top]
bb <- bb[Matrix::rowSums(bb)>0,]
shared.sites <- intersect(rownames(bb), rownames(gg))
bb <- bb[shared.sites,]
gg <- gg[shared.sites,]
bb <- bb[,Matrix::colSums(bb)>0]
gg <- gg[,Matrix::colSums(gg)>0]

# make references
num.good <- 1000
if(ncol(gg) < num.good){
  num.good <- ncol(gg)
}
top.gg <- Matrix::colSums(gg)
top.gg <- top.gg[order(top.gg, decreasing=T)]
top.gg.ids <- names(top.gg)[1:num.good]


# Normalize matrices
message(" - normalizing distributions and creating references")
sub.counts <- a$counts[,c(colnames(gg), colnames(bb))]
sub.counts <- sub.counts[rownames(gg),]
all.res <- tfidf(list(counts=sub.counts), doL2=T)$residuals
bb.norm <- all.res[,colnames(bb)]
gg.norm <- all.res[,colnames(gg)]

# Pick sites
res.ave <- Matrix::rowMeans(gg.norm)
res.res <- RowVar(gg.norm)
resis <- loess(res.res~res.ave)$residuals
names(resis) <- rownames(gg.norm)
#
top.sites <- names(resis[resis > 0])
bb.norm <- bb.norm[top.sites,]
gg.norm <- gg.norm[top.sites,]

# Make references
bad.ref <- Matrix::rowMeans(bb.norm)
good.ref <- Matrix::rowMeans(gg.norm) #gg.norm[,top.gg.ids]
bad.ref <- Matrix(matrix(c(bad.ref / (sqrt(sum(bad.ref^2)))), ncol=1), sparse=T)
good.ref <- Matrix(matrix(c(good.ref / (sqrt(sum(good.ref^2)))),ncol=1), sparse=T)

# Check each cell against ref
message(" - Estimating correlations")
b.ref <- corSparse(gg.norm, bad.ref)
g.ref <- corSparse(gg.norm, good.ref)
b.bad <- corSparse(bb.norm, bad.ref)
g.bad <- corSparse(bb.norm, good.ref)

# Combine results
refs <- data.frame(cbind(b.ref, g.ref))
bads <- data.frame(cbind(b.bad, g.bad))
colnames(refs) <-c("bREF", "gREF")
colnames(bads) <-c("bREF", "gREF")

refs$call <- ifelse(refs$bREF > refs$gREF, 0, 1)
bads$call <- ifelse(bads$bREF > bads$gREF, 0, 1)
rownames(refs) <- colnames(gg.norm)
rownames(bads) <- colnames(bb.norm)
all.refs <- rbind(refs, bads)

shared <- intersect(rownames(a$meta), rownames(all.refs))
meta <- a$meta[shared,]
all.refs <- all.refs[shared,]
all.refs <- cbind(meta, all.refs)

# Handle non-referenced cells
nonrefs <- a$meta[!rownames(a$meta) %in% rownames(all.refs),]
nonrefs$bREF <- NA
nonrefs$gREF <- NA
nonrefs$call <- NA
test <- rbind(all.refs, nonrefs)

# Get new call
test$dif <- test$gREF-test$bREF
top <- test$dif[test$qc_check==1]
bottom <- test$dif[test$qc_check==0]

cut.off <- median(bottom, na.rm=T)
test$call <- ifelse(test$dif > cut.off & test$qc_check == 1, 1, 0)
# Write output
write.table(
  test,
  file = paste0(prefix, ".updated_metadata.txt"),
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE,
  sep = "\t"
)

message(" - Completed step 2 successfully")
