## 1_1b_per_chunk.R — chunked Socrates build for one BED chunk.
##
## Identical to the monolithic 1_1_QC_scifiATAC_data.R pipeline minus the final
## isCellv2 (which computes pTSS/FRiP z-scores across cells — global stats only
## available after merge, so isCellv2 runs once in 1_1c).
##
## Input  : one hash-partitioned BED chunk from 1_1a_chunk_bed.sh.
## Output : <prefix>.soc.rds  (post-convertSparseData soc.obj: meta + counts)
##          <prefix>.raw.before.soc.rds  (intermediate; cleaned up by 1_1c)
##
## Adaptations vs. 2_PopulationStress_maize reference:
##   - org_scaffolds = c("At_Mt","At_Pt","Zm_Mt","Zm_Pt") for combined Zm+At ref
##   - attribute = "gene_id" (matches existing 6_socrates Socrates run)
##
## Usage:
##   Rscript 1_1b_per_chunk.R <chunk_bed> <prefix> <gtf> <chr_sizes> <macs3_narrowPeak>

suppressMessages(library(Socrates))
suppressMessages(library(data.table))

args <- commandArgs(TRUE)
if (length(args) != 5) {
  stop("Rscript 1_1b_per_chunk.R <chunk_bed> <prefix> <gtf> <chr_sizes> <macs3_narrowPeak>")
}

bed       <- as.character(args[1])
out       <- as.character(args[2])
ann       <- as.character(args[3])
chr       <- as.character(args[4])
peak_path <- as.character(args[5])

message(" - chunk mode")
message("   BED       : ", bed)
message("   prefix    : ", out)
message("   peaks     : ", peak_path)

# load BED + annotation + chromosome sizes
obj <- loadBEDandGenomeData(bed, ann, chr, attribute = "gene_id")

# strip organelle reads (combined Zm + At reference)
obj <- countRemoveOrganelle(obj, org_scaffolds = c("At_Mt", "At_Pt", "Zm_Mt", "Zm_Pt"),
                            remove_reads = TRUE)

# load pre-computed MACS3 peaks (replaces inline callACRs_parallel)
message(" - loading MACS3 peaks from ", peak_path)
obj$acr <- fread(peak_path, header = FALSE)
message(" - loaded ", nrow(obj$acr), " peaks")

# per-cell metadata (chunk-local data.frame is ~1/N of the unchunked pool)
obj <- buildMetaData(obj, tss.window = 2000, verbose = TRUE)

# intermediate for sentinel/debug; cleaned up by 1_1c after .soc.rds is written.
saveRDS(obj, file = paste0(out, ".raw.before.soc.rds"))

# cell x window sparse matrix on the 500 bp grid (deterministic — every chunk
# uses the same window vocabulary so mergeSocratesRDS unions rows cleanly).
obj  <- readRDS(paste0(out, ".raw.before.soc.rds"))
sobj <- generateMatrix(obj, filtered = FALSE, windows = 500, peaks = FALSE, verbose = TRUE)

# convert to soc.obj — STOP HERE. isCellv2 deferred to merge stage.
soc.obj <- convertSparseData(sobj, verbose = TRUE)

saveRDS(soc.obj, file = paste0(out, ".soc.rds"))
message(" - chunk done.")
message("   cells (chunk)    : ", nrow(soc.obj$meta))
message("   features (chunk) : ", nrow(soc.obj$counts))
