###############################################################################
## 1_1_QC_scifiATAC_data.R -- Step 1_1: build the raw Socrates object for one library
## Loads a Tn5 BED + gene annotation + chr sizes (Socrates loadBEDandGenomeData), counts and
## removes organelle reads, calls ACRs with MACS2, builds per-barcode QC metadata, bins the
## genome into 500 bp tiles and runs isCellv2. Output: <prefix>.raw.soc.rds (plus the
## <prefix>.raw.before.soc.rds intermediate and <prefix>_macs2_temp/ peak files).
## Usage: Rscript 1_1_QC_scifiATAC_data.R <tn5.bed.gz> <prefix> <gff3|gtf> <chr_sizes> \
##          [genomesize=1.689e9] [org_scaffolds=At_Mt,At_Pt,Zm_Mt,Zm_Pt]
## Called by part2_split/step0_qc/1_QC_scifiATAC_SM2v2_plate.sh (plate arm); the chunked
## combined-genome path uses 1_1b_per_chunk.R + 1_1c_merge_and_qc.R instead.
###############################################################################

# load libraries
suppressMessages(library(Socrates))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

## main function's modufications
isCellv2 <- function(obj, num.test = 20000, num.tn5 = NULL, num.ref = 1000, 
                     background.cutoff = 100, min.pTSS = 0.2, min.FRiP = 0.2, 
                     min.pTSS.z = -2, min.FRiP.z = -2, verbose = T) 
{
  # Helper function to calculate row variance
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
  
  # Step 1: Converting count matrix to sparseMatrix
  if (verbose) {
    message("Step 1: Converting count matrix to sparseMatrix format")
  }
  
  tryCatch({
    if (is.data.frame(obj$counts)) {
      # Assuming the count matrix is provided as a data frame with columns V1, V2, V3 (triplet format)
      sparse_count_matrix <- obj$counts
      if (verbose) {
        message(" - converting triplet format to sparseMatrix")
      }
      sparse_count_matrix$V1 <- factor(sparse_count_matrix$V1)
      sparse_count_matrix$V2 <- factor(sparse_count_matrix$V2)
      sparse_count_matrix <- Matrix::sparseMatrix(
        i = as.numeric(sparse_count_matrix$V1), 
        j = as.numeric(sparse_count_matrix$V2), 
        x = as.numeric(sparse_count_matrix$V3), 
        dimnames = list(levels(sparse_count_matrix$V1), levels(sparse_count_matrix$V2))
      )
    } else if (methods::is(obj$counts, "dgTMatrix")) {
      # If obj$counts is already a sparse matrix in triplet format
      if (verbose) {
        message(" - converting 'dgTMatrix' to 'dgCMatrix'")
      }
      sparse_count_matrix <- as(obj$counts, "dgCMatrix")
    } else if (methods::is(obj$counts, "dgCMatrix")) {
      # If obj$counts is already a sparse matrix in compressed column format
      if (verbose) {
        message(" - count matrix is already a 'dgCMatrix'")
      }
      sparse_count_matrix <- obj$counts
    } else if (is.matrix(obj$counts)) {
      # If the count matrix is a dense matrix, convert it to a sparse matrix
      if (verbose) {
        message(" - converting dense matrix to 'dgCMatrix'")
      }
      sparse_count_matrix <- Matrix::Matrix(obj$counts, sparse = TRUE)
    } else {
      stop("Error: Unsupported format for count matrix. Expected a data frame, 'dgTMatrix', 'dgCMatrix', or dense matrix.")
    }
  }, error = function(e) {
    stop("Error during sparse matrix conversion: ", e$message)
  })
  
  # Step 2: Align metadata and counts
  if (verbose) {
    message("Step 2: Aligning metadata with count matrix")
  }
  shared <- intersect(rownames(obj$meta), colnames(sparse_count_matrix))
  if (length(shared) == 0) {
    stop("Error: No shared identifiers found between metadata and count matrix.")
  }
  sparse_count_matrix <- sparse_count_matrix[, shared]
  obj$meta <- obj$meta[shared, ]
  obj$meta <- obj$meta[order(obj$meta$nSites, decreasing = TRUE), ]
  
  # Step 3: Calculating quality control metrics
  if (verbose) {
    message("Step 3: Calculating quality control metrics")
  }
  obj$meta$pTSS <- obj$meta$tss / obj$meta$total
  obj$meta$FRiP <- obj$meta$acrs / obj$meta$total
  obj$meta$pOrg <- obj$meta$ptmt / obj$meta$total
  
  # Step 4: Setting QC filters
  if (verbose) {
    message("Step 4: Setting quality control filters")
  }
  obj$meta$tss_z <- (obj$meta$pTSS - mean(obj$meta$pTSS)) / sd(obj$meta$pTSS)
  obj$meta$acr_z <- (obj$meta$FRiP - mean(obj$meta$FRiP)) / sd(obj$meta$FRiP)
  obj$meta$sites_z <- (log10(obj$meta$nSites) - mean(log10(obj$meta$nSites))) / sd(log10(obj$meta$nSites))
  obj$meta$tss_z[is.na(obj$meta$tss_z)] <- -10
  obj$meta$acr_z[is.na(obj$meta$acr_z)] <- -10
  obj$meta$sites_z[is.na(obj$meta$sites_z)] <- -10
  obj$meta$qc_check <- ifelse(obj$meta$tss_z < min.pTSS.z | 
                                obj$meta$pTSS < min.pTSS, 0, ifelse(obj$meta$acr_z < 
                                                                      min.FRiP.z | obj$meta$FRiP < min.FRiP, 0, 1))
  
  # Proceed with the rest of the function as before...
  
  return(obj)
}



callACRs_parallel <- function(obj, ChrSize, genomesize = 1.6e+09, shift = -75, 
                              extsize = 150, 
                              output = "bulk_peaks", 
                              tempdir = "./macs2_temp", verbose = TRUE, fdr = 0.05) 
{
  # Load required packages
  library(parallel)
  library(data.table)
  
  if (verbose) {
    message(" - running MACS2 on bulk BED file in parallel by chromosome ...")
  }
  
  # Read chromosome size file
  chr_size <- fread(ChrSize, header = FALSE)
  colnames(chr_size) <- c("chr", "size")  # Assuming file has two columns: chr_name, chr_size
  
  # Read BED file
  bed <- fread(obj$bedpath, header = FALSE)
  colnames(bed) <- c("chr", "start", "end", "bc", "strand")  # Assuming typical BED format
  
  # Create temporary directory
  mac2temp <- tempdir
  if (file.exists(mac2temp)) {
    unlink(mac2temp, recursive = TRUE)
  }
  dir.create(mac2temp)
  
  # Split BED file by chromosome
  bed_split <- split(bed, bed$chr)
  
  # Function to run MACS2 for each chromosome
  run_macs2 <- function(chr_data, chr, chr_size, genomesize, shift, extsize, output, mac2temp, fdr, verbose) {
    # Write individual chromosome BED file
    chr_bed_path <- file.path(mac2temp, paste0(chr, "_temp.bed.gz"))
    fwrite(chr_data, chr_bed_path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    # Update genomesize for the specific chromosome
    chr_total_size <- chr_size$size[chr_size$chr == chr]
    total_genome_size <- sum(chr_size$size)
    new_chr_genomesize <- genomesize * (chr_total_size / total_genome_size)
    
    # Construct MACS2 command
    cmdline <- paste0("macs2 callpeak -t ", chr_bed_path, " -f BED -g ", 
                      new_chr_genomesize, " --keep-dup all -n ", output, "_", chr, 
                      " --nomodel --shift ", shift, " --extsize ", extsize, 
                      " --outdir ", mac2temp, " --qvalue ", fdr)
    
    # Run MACS2 command
    if (verbose) {
      message(" - Running MACS2 for chromosome: ", chr)
    }
    suppressMessages(system(cmdline))
  }
  
  # Get the list of chromosomes to process
  chrs <- names(bed_split)
  
  # Run MACS2 for each chromosome in parallel
  mclapply(chrs, function(chr) {
    run_macs2(bed_split[[chr]], chr, chr_size, genomesize, shift, extsize, output, mac2temp, fdr, verbose)
  }, mc.cores = detectCores() - 1)  # Use all available cores minus one to keep the system responsive
  
  # Combine all peak files into a single file
  all_peaks <- do.call(rbind, lapply(chrs, function(chr) {
    peak_file <- file.path(mac2temp, paste0(output, "_", chr, "_peaks.narrowPeak"))
    if (file.exists(peak_file)) {
      peaks <- fread(peak_file, header = FALSE)
      return(peaks)
    } else {
      warning("Peak file for chromosome ", chr, " not found.")
      return(NULL)
    }
  }))
  
  # Combine all _peaks.xls files into a single file
  all_peaks_xls <- do.call(rbind, lapply(chrs, function(chr) {
    peaks_xls_file <- file.path(mac2temp, paste0(output, "_", chr, "_peaks.xls"))
    if (file.exists(peaks_xls_file)) {
      peaks_xls <- fread(peaks_xls_file, header = TRUE)
      return(peaks_xls)
    } else {
      warning("Peaks XLS file for chromosome ", chr, " not found.")
      return(NULL)
    }
  }))
  
  # Combine all summits.bed files into a single file
  all_summits <- do.call(rbind, lapply(chrs, function(chr) {
    summits_file <- file.path(mac2temp, paste0(output, "_", chr, "_summits.bed"))
    if (file.exists(summits_file)) {
      summits <- fread(summits_file, header = FALSE)
      return(summits)
    } else {
      warning("Summits file for chromosome ", chr, " not found.")
      return(NULL)
    }
  }))
  
  # Store combined peaks, peaks_xls, and summits in the object
  obj$acr <- all_peaks
  # all_peaks_xls and 
  #obj$acr_xls <- all_peaks_xls
  #obj$summits <- all_summits
  
  # Write the combined peaks to final output files
  final_output_path_peaks <- file.path(mac2temp, paste0(output, "_combined_peaks.narrowPeak"))
  fwrite(all_peaks, final_output_path_peaks, sep = "\t", quote = FALSE, col.names = FALSE)
  
  final_output_path_peaks_xls <- file.path(mac2temp, paste0(output, "_combined_peaks.xls"))
  fwrite(all_peaks_xls, final_output_path_peaks_xls, sep = "\t", quote = FALSE, col.names = TRUE)
  
  final_output_path_summits <- file.path(mac2temp, paste0(output, "_combined_summits.bed"))
  fwrite(all_summits, final_output_path_summits, sep = "\t", quote = FALSE, col.names = FALSE)
  
  #   # Cleanup temporary files if needed
  #   if (verbose) {
  #     message(" - Cleanup temporary files")
  #   }
  #   unlink(mac2temp, recursive = TRUE)
  
  # Return the updated object
  return(obj)
}


loadBEDandGenomeData2 <- function(
    bed,
    ann,
    sizes,
    attribute = "gene_id",
    verbose = TRUE,
    is.fragment = FALSE
) {
  # -------------------------
  # Helpers
  # -------------------------
  .preRunChecks <- function(bed, ann, sizes, verbose = TRUE) {
    if (verbose) message("Running pre-check on input files and executable paths ...")
    
    if (!file.exists(bed)) stop("BED file does not exist: ", bed)
    if (verbose) message("BED file path = ", bed, " ... ok")
    
    if (!file.exists(ann)) stop("Annotation file does not exist: ", ann)
    if (verbose) message("GFF/GTF file path = ", ann, " ... ok")
    
    if (!file.exists(sizes)) stop("Chromosome sizes file does not exist: ", sizes)
    if (verbose) message("Chromosome sizes file path = ", sizes, " ... ok")
  }
  
  .readBed <- function(bed) {
    if (grepl("\\.gz$", bed, ignore.case = TRUE)) {
      read.table(gzfile(bed), header = FALSE, sep = "\t", quote = "", comment.char = "")
    } else {
      read.table(bed, header = FALSE, sep = "\t", quote = "", comment.char = "")
    }
  }
  
  .sanitizeStrict9Cols <- function(ann_in, ann_out = NULL) {
    # Ensures each non-comment line has exactly 9 tab-separated fields.
    # If NF>9, glue fields 10..NF into field 9 separated by spaces.
    if (is.null(ann_out)) {
      ann_out <- paste0(ann_in, ".strict9.tmp.gtf")
    }
    
    con_in <- if (grepl("\\.gz$", ann_in, ignore.case = TRUE)) gzfile(ann_in, "rt") else file(ann_in, "rt")
    on.exit(try(close(con_in), silent = TRUE), add = TRUE)
    
    con_out <- file(ann_out, "wt")
    on.exit(try(close(con_out), silent = TRUE), add = TRUE)
    
    while (length(line <- readLines(con_in, n = 1, warn = FALSE)) > 0) {
      if (startsWith(line, "#") || nchar(line) == 0) {
        writeLines(line, con_out)
        next
      }
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 9) next
      if (length(fields) > 9) {
        fields[9] <- paste(c(fields[9], fields[10:length(fields)]), collapse = " ")
        fields <- fields[1:9]
      }
      writeLines(paste(fields, collapse = "\t"), con_out)
    }
    
    ann_out
  }
  
  .detectAnnType <- function(ann_path) {
    # Accept .gtf/.gtf.gz and .gff/.gff3/.gff.gz/.gff3.gz
    ann_lc <- tolower(ann_path)
    if (grepl("\\.gtf(\\.gz)?$", ann_lc)) return("gtf")
    if (grepl("\\.gff3?(\\.gz)?$", ann_lc)) return("gff3")
    # fallback: assume gtf-like
    "gtf"
  }
  
  # -------------------------
  # Main
  # -------------------------
  .preRunChecks(bed, ann, sizes, verbose = verbose)
  
  bedpath <- bed
  annpath <- ann
  chrpath <- sizes
  
  if (verbose) message(" - loading data (this may take a while for big BED files) ...")
  obj <- .readBed(bedpath)
  
  # Expect at least 4 columns: chr, start, end, cell/barcode
  if (ncol(obj) < 4) {
    stop("BED must have >=4 columns (chr, start, end, cellID/barcode). Found: ", ncol(obj))
  }
  
  # Fragment-to-Tn5 insertion conversion (single bp sites at both ends)
  if (is.fragment) {
    if (verbose) message(" - converting fragment file to single-bp Tn5 insertion sites ...")
    
    # If fragments: V1 chr, V2 start, V3 end, V4 barcode/cell
    start.coordinates <- data.frame(
      V1 = obj$V1,
      V2 = obj$V2,
      V3 = obj$V2 + 1L,
      V4 = obj$V4,
      V5 = "+"
    )
    end.coordinates <- data.frame(
      V1 = obj$V1,
      V2 = obj$V3 - 1L,
      V3 = obj$V3,
      V4 = obj$V4,
      V5 = "-"
    )
    all.coordinates <- rbind(start.coordinates, end.coordinates)
    all.coordinates <- all.coordinates[order(all.coordinates$V1, all.coordinates$V2, decreasing = FALSE), ]
    obj <- all.coordinates[!duplicated(all.coordinates), , drop = FALSE]
  }
  
  # Sanitize annotation to strict 9 columns to avoid readGFF errors in some pipelines
  ann_strict <- .sanitizeStrict9Cols(annpath)
  on.exit({
    if (file.exists(ann_strict)) unlink(ann_strict)
  }, add = TRUE)
  
  anntype <- .detectAnnType(annpath)
  
  # Build TxDb
  if (verbose) message(" - building TxDb from annotation (format=", anntype, ", dbxrefTag=", attribute, ") ...")
  gff <- suppressWarnings(
    suppressMessages(
      GenomicFeatures::makeTxDbFromGFF(
        file = ann_strict,
        format = anntype,
        dbxrefTag = attribute
      )
    )
  )
  
  chrom <- read.table(chrpath, header = FALSE, sep = "\t", quote = "", comment.char = "")
  if (ncol(chrom) < 2) stop("Chromosome sizes file must have at least 2 columns: chr, size.")
  
  if (verbose) message(" - finished loading data")
  
  list(
    bed = obj,
    gff = gff,
    chr = chrom,
    bedpath = bedpath,
    annpath = annpath,
    chrpath = chrpath
  )
}


# arguments
args <- commandArgs(T)
if(length(args) < 4){stop("Rscript 1_1_QC_scifiATAC_data.R <bed> <prefix> <gff/gtf> <chr> [genomesize] [org_scaffolds_comma_sep]")}

# load data
bed <- as.character(args[1])
out <- as.character(args[2])
ann <- as.character(args[3])
chr <- as.character(args[4])

# optional per-genome params (arg5 = MACS genome size, arg6 = comma-sep organelle scaffolds).
# Defaults reproduce the ZmATcombined config, so existing combined-genome runs are unchanged.
genomesize    <- if(length(args) >= 5 && nzchar(args[5])) as.numeric(args[5]) else (1.6e9 + 0.089e9)
org_scaffolds <- if(length(args) >= 6 && nzchar(args[6])) strsplit(args[6], ",", fixed=TRUE)[[1]] else c("At_Mt","At_Pt","Zm_Mt","Zm_Pt")
message(" - split-genome params: genomesize=", format(genomesize, scientific=TRUE),
        " ; org_scaffolds=", paste(org_scaffolds, collapse=","))


# gtf -> Socrates built-in loader (proven on combined + B73v5); gff3 -> local gff3-capable loader (loadBEDandGenomeData2)
if(grepl("\\.gff3?(\\.gz)?$", ann, ignore.case=TRUE)){
  obj <- loadBEDandGenomeData2(bed, ann, chr, attribute="gene_id", verbose = TRUE)
} else {
  obj <- loadBEDandGenomeData(bed, ann, chr, attribute="gene_id",  verbose = TRUE)
}

# count organelle reads
obj <- countRemoveOrganelle(obj, org_scaffolds=org_scaffolds, remove_reads=T)

# call ACRs
obj <- callACRs_parallel(obj, chr, genomesize=genomesize, 
                         shift= -75, 
                         extsize=150,
                         fdr=0.1,
                         output=paste0(out,"_peaks"), 
                         tempdir=paste0(out, '_macs2_temp'), 
                         verbose=T)

# build metadata
obj <- buildMetaData(obj, tss.window=2000, verbose=TRUE)

# Save QC object
saveRDS(obj, file=paste0(out,".raw.before.soc.rds"))

# Generate sparse matrix
obj <- readRDS(paste0(out,".raw.before.soc.rds"))
sobj <- generateMatrix(obj, filtered=F,  windows=500, peaks=F, verbose=T)

# Add metadata for cells QC
# Convert to Socrates format for downstream analysis. 
soc.obj <- convertSparseData(sobj, verbose=T)
soc.obj <- isCellv2(soc.obj)

# save QC object
saveRDS(soc.obj, file=paste0(out,".raw.soc.rds"))