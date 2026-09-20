#!/usr/bin/env Rscript
###############################################################################
## 2_0_3_fixed_set_compare.R
##
## Fair pre-vs-post (cleaning) comparison on a FIXED barcode set, applying all
## quality filters but NO post-clean depth re-gating, with quality thresholds
## FROZEN from the pre-clean object.
##
## Why: AmbientMapper cleaning removes reads, so post-clean depth is mechanically
## <= pre-clean depth for every barcode. Re-applying an absolute depth floor
## post-clean therefore penalises a cell *for having been cleaned*. To isolate
## the effect of cleaning on per-cell quality we instead:
##   1. Call the cell set on RAW (pre-clean) depth   (total >= DEPTH).
##   2. Freeze the pTSS / FRiP / dif thresholds, and the qc_check z-baselines
##      (mean/sd of pTSS, FRiP), from the PRE-clean metadata.
##   3. Carry the identical barcode set into the post-clean metadata and apply
##      the frozen quality gates (everything except depth), so the only variable
##      that differs between stages is the data, not the thresholds.
##
## This is the complement to the standard 1_3 cascade: 1_3 re-calls cells per
## stage (the cell-COUNT story); this script holds the cell set fixed (the
## per-cell QUALITY story). Metadata-only — consumes the two 1_2 outputs
## (<pool>.updated_metadata.txt); no Socrates object is rebuilt.
##
## Output (long, one row per fixed-set cell x stage; stage in PreClean/PostClean):
##   <out_prefix>.fixedSet.minDepth<DEPTH>.pre_post.txt
##   <out_prefix>.fixedSet.thresholds.txt
##
## Usage:
##   Rscript 2_0_3_fixed_set_compare.R <pre_meta> <post_meta> <out_prefix> [DEPTH=200] [outdir=.]
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 2_0_3_fixed_set_compare.R <pre_meta> <post_meta> <out_prefix> [DEPTH=200] [outdir=.]")
}
pre_path   <- args[1]
post_path  <- args[2]
out_prefix <- args[3]
DEPTH      <- if (length(args) >= 4) as.numeric(args[4]) else 200
outdir     <- if (length(args) >= 5) args[5] else "."
if (!is.finite(DEPTH) || DEPTH <= 0) stop("DEPTH must be a positive number; got: ", args[4])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Filter constants — identical to 1_3_metaQC_scifiATAC_data.R so the frozen
# thresholds reproduce that cascade exactly on the pre-clean object.
tss_z_filter <- -8; acr_z_filter <- -8; org_thresh <- 0.2; diff_z_keep <- -1
req <- c("cellID","total","pTSS","FRiP","pOrg","dif","tss_z","acr_z","qc_check","nSites")

## ---- robust loader -------------------------------------------------------- #
## step 1_2 writes via R write.table(row.names=TRUE): the header carries one
## fewer field than each data row (leading unnamed rowname column). read.table
## absorbs that automatically but fread does not, so detect and adapt.
read_meta <- function(path) {
  h  <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]]
  d1 <- strsplit(readLines(path, n = 2L)[2L], "\t", fixed = TRUE)[[1]]
  nh <- length(h); nd <- length(d1)
  if (nd == nh + 1L) {
    dt <- fread(path, skip = 1L, header = FALSE, sep = "\t", showProgress = FALSE)
    setnames(dt, c("rowname", h))
  } else if (nd == nh) {
    dt <- fread(path, header = TRUE, sep = "\t", showProgress = FALSE)
  } else {
    stop("Unexpected column counts in ", path, ": header=", nh, " data=", nd)
  }
  miss <- setdiff(req, colnames(dt))
  if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))
  dt
}

message(" - loading pre  : ", pre_path)
pre  <- read_meta(pre_path)
message(" - loading post : ", post_path)
post <- read_meta(post_path)

## ---- frozen thresholds from PRE (replicates the 1_3 cascade) --------------- #
thresh_from_z <- function(df, zcol, metric, zcut, floor_val) {
  below <- df[df[[zcol]] < zcut]
  if (nrow(below) == 0L) return(floor_val)
  max(floor_val, max(below[[metric]], na.rm = TRUE))
}

v1 <- pre[is.finite(total) & total >= DEPTH]            # raw-depth cell call
if (nrow(v1) < 50L) stop("Too few pre barcodes with total >= ", DEPTH, ": ", nrow(v1))
tss_thr  <- thresh_from_z(v1, "tss_z", "pTSS", tss_z_filter, 0.2)
v2 <- v1[pTSS >= tss_thr]
frip_thr <- thresh_from_z(v2, "acr_z", "FRiP", acr_z_filter, 0.2)
v3 <- v2[FRiP >= frip_thr]
v4 <- v3[pOrg <= org_thresh]
difsub <- v4[is.finite(dif) & is.finite(total) & total > 0]
if (nrow(difsub) < 50L) stop("Too few pre cells with finite dif to set dif threshold: ", nrow(difsub))
z_dif   <- as.numeric(scale(difsub$dif))
dif_thr <- min(difsub$dif[z_dif >= diff_z_keep], na.rm = TRUE)
v5 <- difsub[dif >= dif_thr]
pre_v6 <- v5[qc_check == 1]

# frozen z-baselines for the qc_check z-component (pre object-wide mean/sd)
mu_tss  <- mean(pre$pTSS, na.rm = TRUE); sd_tss  <- sd(pre$pTSS, na.rm = TRUE)
mu_frip <- mean(pre$FRiP, na.rm = TRUE); sd_frip <- sd(pre$FRiP, na.rm = TRUE)

## ---- FIXED SET = barcodes called as cells on raw (pre) depth --------------- #
fixed_ids <- v1$cellID
setkey(pre, cellID); setkey(post, cellID)

# frozen quality gate applied to either stage (NO depth term)
apply_quality_frozen <- function(d) {
  tz <- (d$pTSS - mu_tss)  / sd_tss
  az <- (d$FRiP - mu_frip) / sd_frip
  qc_frozen <- as.integer(tz >= -2 & d$pTSS >= 0.2 & az >= -2 & d$FRiP >= 0.2)
  pass <- as.integer(d$pTSS >= tss_thr & d$FRiP >= frip_thr &
                     d$pOrg <= org_thresh & d$dif >= dif_thr & qc_frozen == 1L)
  list(tss_z_frozen = tz, acr_z_frozen = az, qc_frozen = qc_frozen,
       pass_quality_frozen = pass)
}

build_stage <- function(dt, ids, stage_label) {
  d <- dt[.(ids), nomatch = NA]                # keyed join; absent ids -> NA row
  q <- apply_quality_frozen(d)
  data.table(
    cellID  = ids,
    stage   = stage_label,
    present = as.integer(!is.na(d$total)),
    total   = d$total, nSites = d$nSites,
    pTSS = d$pTSS, FRiP = d$FRiP, pOrg = d$pOrg, dif = d$dif,
    tss_z = d$tss_z, acr_z = d$acr_z, qc_check = d$qc_check,
    tss_z_frozen = q$tss_z_frozen, acr_z_frozen = q$acr_z_frozen,
    qc_frozen = q$qc_frozen, pass_quality_frozen = q$pass_quality_frozen
  )
}

pre_rows  <- build_stage(pre,  fixed_ids, "PreClean")
post_rows <- build_stage(post, fixed_ids, "PostClean")

# in_pre_v6 == 1 marks the cells that pass the FULL pre-clean cascade (the
# canonical "good cells"); restrict downstream comparisons to these for the
# strict fixed-good-set view, or keep all raw-depth cells for the broad view.
v6_ids <- pre_v6$cellID
pre_rows[,  in_pre_v6 := as.integer(cellID %chin% v6_ids)]
post_rows[, in_pre_v6 := as.integer(cellID %chin% v6_ids)]

long <- rbindlist(list(pre_rows, post_rows), use.names = TRUE)
long[present == 0L, total := 0]                # lost post-clean -> 0 reads

## ---- write outputs --------------------------------------------------------- #
out_long <- file.path(outdir, paste0(out_prefix, ".fixedSet.minDepth", DEPTH, ".pre_post.txt"))
fwrite(long, out_long, sep = "\t")
message(" - wrote ", out_long, "  (", nrow(long), " rows = ",
        length(fixed_ids), " cells x 2 stages)")

thr <- data.table(
  param = c("depth_raw","tss_thr","frip_thr","org_thresh","dif_thr",
            "mu_pTSS_pre","sd_pTSS_pre","mu_FRiP_pre","sd_FRiP_pre",
            "n_fixed_set","n_pre_v6"),
  value = c(DEPTH, tss_thr, frip_thr, org_thresh, dif_thr,
            mu_tss, sd_tss, mu_frip, sd_frip, length(fixed_ids), nrow(pre_v6)))
out_thr <- file.path(outdir, paste0(out_prefix, ".fixedSet.thresholds.txt"))
fwrite(thr, out_thr, sep = "\t")
message(" - wrote ", out_thr)

## ---- summary --------------------------------------------------------------- #
np  <- sum(post_rows$present)
npq <- sum(post_rows$pass_quality_frozen, na.rm = TRUE)
med <- function(x) round(median(x, na.rm = TRUE), 4)
pres_pre  <- pre_rows[present == 1L]
pres_post <- post_rows[present == 1L]
message("\n==== FIXED-SET pre/post summary: ", out_prefix, " (raw depth >= ", DEPTH, ") ====")
message(" fixed set (raw-depth cells)      : ", length(fixed_ids))
message("   of which pass full pre QC (v6) : ", nrow(pre_v6))
message(" present post-clean              : ", np,  " (", round(100*np /length(fixed_ids),1), "%)")
message(" pass frozen quality post-clean  : ", npq, " (", round(100*npq/length(fixed_ids),1), "%)")
message(" median total  pre -> post       : ", med(pres_pre$total), " -> ", med(pres_post$total))
message(" median FRiP   pre -> post       : ", med(pres_pre$FRiP),  " -> ", med(pres_post$FRiP))
message(" median pTSS   pre -> post       : ", med(pres_pre$pTSS),  " -> ", med(pres_post$pTSS))
message(" median pOrg   pre -> post       : ", med(pres_pre$pOrg),  " -> ", med(pres_post$pOrg))
message("====================================================================\n")
