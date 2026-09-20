#!/usr/bin/env Rscript
###############################################################################
## 2_0_0_Obj_integreation.R -- Step 2_0: merge the per-species Socrates objects of one stage
## Reads <sample_dir>/step0_qc/<prefix>*.raw.soc.rds and the matching *.updated_metadata.txt
## (At + B73 of one stage), tags species and stage (Clean. prefix -> PostClean, else PreClean)
## and merges them with Socrates::mergeSocratesRDS into one combined-genome object:
## <sample_dir>/step1_integrate/<prefix>.full.SocObj.rds + <prefix>.full.metadata_updated.txt
## Usage: Rscript 2_0_0_Obj_integreation.R <sample_dir> <prefix>
##        (SM2 SM2 for PreClean, SM2v2_clean Clean.SM2v2 for PostClean)
## Driver: 2_0_0_Obj_integreation_SM2v2.sh
###############################################################################

suppressPackageStartupMessages({
  library(Socrates)
  library(data.table)
  library(stringr)
})

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 2_0_0_Obj_integreation.R <sample_dir> <prefix>\n",
       "Example: Rscript 2_0_0_Obj_integreation.R SM2v2_clean Clean.SM2v2")
}

sample_dir <- args[1]  # e.g. "SM2" (directory containing RDS + metadata)
prefix     <- args[2]  # e.g. "Clean.SM2" or "SM2"

if (!dir.exists(sample_dir)) stop("sample_dir does not exist: ", sample_dir)

message(" - Reading objects to integrate from: ", sample_dir)
message(" - Prefix filter: ", prefix)

# -------------------------
# 1) Find Socrates objects
# -------------------------
# Expect filenames like:
#   <prefix>_At_ZmATcombined...raw.soc.rds
#   <prefix>_B73_ZmATcombined...raw.soc.rds
pattern_rds <- paste0("^", gsub("\\.", "\\\\.", prefix), ".*raw\\.soc\\.rds$")
rds_files <- list.files(path = file.path(sample_dir, "step0_qc"), pattern = pattern_rds, full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No .raw.soc.rds found with pattern: ", pattern_rds, " in ", sample_dir)
}

# Define object names = basename without suffix
obj_names <- basename(rds_files)
obj_names <- sub("\\.raw\\.soc\\.rds$", "", obj_names)

message(" - ", length(rds_files), " Socrates objects to merge:")
message("   ", paste(obj_names, collapse = ", "))

# -------------------------
# 2) Read and merge Socrates objects
# -------------------------
obj_DB <- setNames(lapply(rds_files, readRDS), obj_names)

obj_merge <- mergeSocratesRDS(obj.list = obj_DB)
message(" - Socrates merge done.")
print(obj_merge$counts[1:5,1:5])
message(" ... ")
message(" ... ")

# -------------------------
# 3) Read updated metadata files
# -------------------------
# Expect filenames like: <prefix>....updated_metadata.txt
pattern_meta <- paste0("^", gsub("\\.", "\\\\.", prefix), ".*updated_metadata\\.txt$")
meta_files <- list.files(path = file.path(sample_dir, "step0_qc"), pattern = pattern_meta, full.names = TRUE)

if (length(meta_files) == 0) {
  stop("No updated_metadata.txt found with pattern: ", pattern_meta, " in ", sample_dir)
}

read_meta <- function(f) {
  df <- read.table(f, header = TRUE, sep='\t')
  row.names(df) <- df$cellID
  df
}
meta_list <- lapply(meta_files, read_meta)
names(meta_list) <- sub("\\.updated_metadata\\.txt$", "", basename(meta_files))

# Add sample_id, species and stage, parsed from the object name (At vs B73; Clean. prefix = PostClean).
add_sample_id <- function(df, nm) {
  df$sample_id <- nm
  df$species <- fifelse(grepl("_At_", nm) | grepl("^SM.*_At", nm) | grepl("At", nm), "At",
                        fifelse(grepl("_B73_", nm) | grepl("B73", nm), "B73", NA_character_))
  df$stage <- fifelse(grepl("^Clean\\.", nm), "PostClean", "PreClean")
  df
}

meta_list <- Map(add_sample_id, meta_list, names(meta_list))

# bind with fill=TRUE so differing columns don’t break
meta_all <- rbindlist(meta_list, use.names = TRUE, fill = TRUE)
meta_all <- as.data.frame(meta_all)
row.names(meta_all) <- meta_all$cellID
print(meta_all[1:5,])

# Ensure rownames are cell IDs
if (is.null(rownames(meta_all)) || anyDuplicated(rownames(meta_all)) > 0) {
  # if duplicates exist, keep first occurrence
  meta_all <- meta_all[!duplicated(rownames(meta_all)), , drop = FALSE]
}

message(" - Metadata loaded: ", nrow(meta_all), " cells; ", ncol(meta_all), " columns")

# -------------------------
# 4) Align merged counts and metadata
# -------------------------
shared <- intersect(colnames(obj_merge$counts), rownames(meta_all))
message(" - Shared cells (counts ∩ meta): ", length(shared))

if (length(shared) == 0) {
  stop("No shared cell IDs between merged counts and metadata. Check cellID naming consistency.")
}

obj_merge$counts <- obj_merge$counts[, shared, drop = FALSE]
obj_merge$meta   <- meta_all[shared, , drop = FALSE]

# -------------------------
# 5) Save outputs
# -------------------------
dir.create(file.path(sample_dir, "step1_integrate"), showWarnings = FALSE, recursive = TRUE)
out_soc  <- file.path(sample_dir, "step1_integrate", paste0(prefix, ".full.SocObj.rds"))
out_meta <- file.path(sample_dir, "step1_integrate", paste0(prefix, ".full.metadata_updated.txt"))

saveRDS(obj_merge, file = out_soc)
fwrite(obj_merge$meta, file = out_meta, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

message(" - Saved merged Socrates object: ", out_soc)
message(" - Saved merged metadata:       ", out_meta)
