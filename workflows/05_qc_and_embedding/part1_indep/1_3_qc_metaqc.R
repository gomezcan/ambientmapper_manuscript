###############################################################################
## 1_3_qc_metaqc.R -- Step 1_3: meta-QC cascade on the 1_2 metadata (one library, one depth floor)
## Applies the depth floor, the pTSS and FRiP z-score filters, the organelle-fraction cap and
## the Diff filter to the 1_2 metadata and writes the surviving barcodes after each step as
## <prefix>.minDepth<D>.updated_metadata_v{1..6}.txt plus plots/<prefix>.minDepth<D>.QC_FIGURES.pdf.
## v6 (all filters + final qc_check) is the canonical cell set used downstream.
## Usage: Rscript 1_3_qc_metaqc.R <prefix>.updated_metadata.txt <prefix> <depth_filter>
##  - CLI arg: depth_filter (in raw unique insertions, not log10)
##  - depth_filter is embedded into ALL output filenames + plot PDF name
###############################################################################

suppressPackageStartupMessages({
  library(MASS)
  library(viridis)
})

rm(list=ls())

# -------------------------
# Args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 1_3_qc_metaqc.R <input> <name> <depth_filter>\n",
       "  <depth_filter> is minimum unique insertions (integer), e.g. 100 or 1000")
}
input <- as.character(args[1])
name  <- as.character(args[2])
depth_filter_raw <- suppressWarnings(as.numeric(args[3]))

if (!is.finite(depth_filter_raw) || depth_filter_raw <= 0) {
  stop("depth_filter must be a positive number of unique insertions (e.g. 100, 500, 1000). Got: ", args[3])
}

# depth_filter is applied on log10(unique + 1)
depth_filter <- log10(depth_filter_raw + 1)

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Filters (tune as needed)
# -------------------------
tss_z_filter  <- -8
acr_z_filter  <- -8
org_thresh    <- 0.2
diff_z_keep   <- -1

# filename tag
tag <- paste0(".minDepth", depth_filter_raw)

# -------------------------
# Rank-depth + cutoff by min depth
# -------------------------
plotDist <- function(x, main = "") {
  x <- x[order(x$unique, decreasing = TRUE), , drop = FALSE]
  rank  <- log10(seq_len(nrow(x)))
  depth <- log10(x$unique + 1)

  keep_idx <- which(depth >= depth_filter)
  if (length(keep_idx) == 0) {
    stop("No barcodes pass depth_filter=", depth_filter_raw,
         ". Check 'unique' column or lower the threshold.")
  }
  cells_n <- max(keep_idx)
  knee    <- rank[cells_n]
  reads   <- as.integer(x$unique[cells_n])

  plot(rank[(cells_n+1):length(rank)], depth[(cells_n+1):length(depth)],
       type="l", lwd=2, col="grey75", main=main,
       xlim=range(rank), ylim=range(depth),
       xlab="Barcode rank (log10)", ylab="Unique Tn5 insertions (log10)")
  lines(rank[1:cells_n], depth[1:cells_n], lwd=2, col="darkorchid4")
  grid()
  abline(v=knee, col="red", lty=2, lwd=1)
  abline(h=depth[cells_n], col="red", lty=2, lwd=1)
  text(x=min(rank) + 0.1, y=min(depth) + 0.2,
       labels=paste0("# cells=", cells_n, " | cutoff unique>=", reads),
       adj=c(0,0))

  x[seq_len(cells_n), , drop = FALSE]
}

# -------------------------
# Load metadata
# -------------------------
message(" - loading meta data for ", name)

a <- read.table(input, header = TRUE, sep = "\t", check.names = FALSE)

req <- c("total","pTSS","FRiP","pOrg","tss_z","acr_z","dif","qc_check")
missing <- setdiff(req, colnames(a))
if (length(missing) > 0) {
  stop("Missing required columns in metadata: ", paste(missing, collapse=", "))
}

# define unique (adjust if you have an actual unique column)
a$unique <- a$total

# -------------------------
# QC figure PDF
# -------------------------
pdf(file = file.path("plots", paste0(name, tag, ".QC_FIGURES.pdf")), width = 14, height = 3)
layout(matrix(1:5, nrow = 1))

# -------------------------
# v1: depth filter
# -------------------------
message(" - calling cells (depth) for ", name)
meta.v1 <- plotDist(a, main = paste0(name, " | depth >= ", depth_filter_raw))
message(" - number of cells after cell calling v1 = ", nrow(meta.v1))
write.table(meta.v1, file = paste0(name, tag, ".updated_metadata_v1.txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

.thresh_from_z <- function(df, zcol, metric_col, zcut, floor_val) {
  below <- df[df[[zcol]] < zcut, , drop=FALSE]
  if (nrow(below) == 0) return(floor_val)
  max(floor_val, max(below[[metric_col]], na.rm=TRUE))
}

.plot_density <- function(df, y, ylab, main, ylims, thresh, xlims=c(2.0, 6.5)) {
  df2 <- df[is.finite(df$unique) & df$unique > 0 & is.finite(df[[y]]), , drop=FALSE]
  if (nrow(df2) < 10) {
    plot.new(); title(main=paste0(main, " (too few finite points)"))
    return(invisible(NULL))
  }
  den <- kde2d(log10(df2$unique), df2[[y]], n=300, h=c(0.2, 0.05),
               lims=c(xlims, ylims))
  image(den, useRaster=TRUE, col=c("white", rev(magma(100))),
        xlab="Unique Tn5 insertions (log10)", ylab=ylab, main=main)
  grid(lty=1, lwd=0.5, col="grey90")
  abline(h=thresh, col="red", lty=2, lwd=1)
  box()
}

# -------------------------
# v2: pTSS filter
# -------------------------
meta.v1 <- meta.v1[order(meta.v1$tss_z, decreasing = TRUE), , drop=FALSE]
tss_thresh <- .thresh_from_z(meta.v1, "tss_z", "pTSS", tss_z_filter, floor_val=0.2)
meta.v2 <- meta.v1[meta.v1$pTSS >= tss_thresh, , drop=FALSE]

.plot_density(meta.v1, y="pTSS", ylab="Fraction reads in TSS",
              main=paste0(name, " | pTSS (keep >= ", signif(tss_thresh,3), ")"),
              ylims=c(0,1.0), thresh=tss_thresh)
legend("topright", legend=paste0("# cells=", nrow(meta.v2)), bty="n")

message(" - number of cells after pTSS filter v2 = ", nrow(meta.v2))
write.table(meta.v2, file = paste0(name, tag, ".updated_metadata_v2.txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# -------------------------
# v3: FRiP filter
# -------------------------
meta.v2 <- meta.v2[order(meta.v2$acr_z, decreasing = TRUE), , drop=FALSE]
frip_thresh <- .thresh_from_z(meta.v2, "acr_z", "FRiP", acr_z_filter, floor_val=0.2)
meta.v3 <- meta.v2[meta.v2$FRiP >= frip_thresh, , drop=FALSE]

.plot_density(meta.v2, y="FRiP", ylab="Fraction Tn5 insertions in ACRs (FRiP)",
              main=paste0(name, " | FRiP (keep >= ", signif(frip_thresh,3), ")"),
              ylims=c(0,1.0), thresh=frip_thresh)
legend("topright", legend=paste0("# cells=", nrow(meta.v3)), bty="n")

message(" - number of cells after FRiP filter v3 = ", nrow(meta.v3))
write.table(meta.v3, file = paste0(name, tag, ".updated_metadata_v3.txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# -------------------------
# v4: organelle filter (keep <= org_thresh)
# -------------------------
meta.v3 <- meta.v3[order(meta.v3$pOrg, decreasing = FALSE), , drop=FALSE]
meta.v4 <- meta.v3[meta.v3$pOrg <= org_thresh, , drop=FALSE]

df_org <- meta.v3[is.finite(meta.v3$unique) & meta.v3$unique > 0 & is.finite(meta.v3$pOrg), , drop=FALSE]
den <- kde2d(log10(df_org$unique), df_org$pOrg, n=300, h=c(0.2, 0.05),
             lims=c(c(2.0, 6.5), c(0, 0.5)))
image(den, useRaster=TRUE, col=c("white", rev(magma(100))),
      xlab="Unique Tn5 insertions (log10)", ylab="Fraction organelle (pOrg)",
      main=paste0(name, " | pOrg (keep <= ", org_thresh, ")"))
grid(lty=1, lwd=0.5, col="grey90")
abline(h=org_thresh, col="red", lty=2, lwd=1)
box()
legend("topright", legend=paste0("# cells=", nrow(meta.v4)), bty="n")

message(" - number of cells after pOrg filter v4 = ", nrow(meta.v4))
write.table(meta.v4, file = paste0(name, tag, ".updated_metadata_v4.txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# -------------------------
# v5: dif filter
# -------------------------
message("dif NA: ", sum(is.na(meta.v4$dif)),
        " | dif non-finite: ", sum(!is.finite(meta.v4$dif)),
        " | unique<=0: ", sum(!(meta.v4$unique > 0)))

meta.dif <- meta.v4[
  is.finite(meta.v4$dif) &
  is.finite(meta.v4$unique) &
  meta.v4$unique > 0,
  , drop = FALSE
]

if (nrow(meta.dif) < 50) {
  stop("Too few cells with finite dif/unique to run dif QC: n=", nrow(meta.dif))
}

z_dif <- as.numeric(scale(meta.dif$dif))
dif_thresh <- min(meta.dif$dif[z_dif >= diff_z_keep], na.rm = TRUE)

meta.v5 <- meta.dif[meta.dif$dif >= dif_thresh, , drop = FALSE]

ymin <- min(-1, min(meta.dif$dif, na.rm=TRUE))
ymax <- max( 1, max(meta.dif$dif, na.rm=TRUE))

den <- kde2d(log10(meta.dif$unique), meta.dif$dif,
             n = 300, h = c(0.2, 0.1),
             lims = c(c(2.0, 6.5), c(ymin, ymax)))

image(den, useRaster=TRUE, col=c("white", rev(magma(100))),
      xlab="Unique Tn5 insertions (log10)", ylab="dif (good - bad)",
      main=paste0(name, " | dif (keep >= ", signif(dif_thresh,3), ")"))
grid(lty=1, lwd=0.5, col="grey90")
abline(h=dif_thresh, col="red", lty=2, lwd=1)
box()
legend("topright", legend=paste0("# cells=", nrow(meta.v5)), bty="n")

message(" - number of cells after dif filter v5 = ", nrow(meta.v5))
write.table(meta.v5, file = paste0(name, tag, ".updated_metadata_v5.txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# -------------------------
# v6: final gate by qc_check
# -------------------------
meta.v6 <- meta.v5[meta.v5$qc_check == 1, , drop=FALSE]
write.table(meta.v6, file = paste0(name, tag, ".updated_metadata_v6.txt"),
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

dev.off()

message(" - number of cells after cell calling v1 = ", nrow(meta.v1))
message(" - number of cells after pTSS filter v2 = ", nrow(meta.v2))
message(" - number of cells after FRiP filter v3 = ", nrow(meta.v3))
message(" - number of cells after pOrg filter v4 = ", nrow(meta.v4))
message(" - number of cells after dif filter v5 = ", nrow(meta.v5))
message(" - number of cells after all filters v6 = ", nrow(meta.v6))
