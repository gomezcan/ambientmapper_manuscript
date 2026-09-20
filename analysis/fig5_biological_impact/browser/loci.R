#!/usr/bin/env Rscript
# =============================================================================
# browser/loci.R  -  STEP 4 of the Fig 5 panel L genome-browser chain: build the browser AUDITION
#   loci list, every marker gene that passes the panel-K selection gates (capped per type), with
#   GFF3 coordinates, written to figures/main/fig5/browser/plot_loci.tsv.
# Inputs  figures/main/fig5/Fig5_P3E_examples_contrast.tsv   (all genes, all metrics; fig5_K_examples.R)
#         figures/main/fig5/Fig5_P3E_examples_selection.tsv  (flag: current UMAP picks)
#         data/processed/scifiATAC_B73_Arabidopsis/socrates/_data/_GenomeInfo/{B73v5,TAIR10}.gff3
# Output  figures/main/fig5/browser/plot_loci.tsv   (REWRITTEN on every run)
# Run     Rscript analysis/fig5_biological_impact/browser/loci.R    (from the repo root, ~30 s)
# =============================================================================
#
# THE IDEA: render browser tracks for the WHOLE gate-passing candidate pool, then choose the
# final showcase markers/types for BOTH the panel-K UMAP figure and the browser exhibit from
# what the tracks show (an audition rather than "UMAP picks -> 3 loci").
# SELECTION INTEGRITY: the audition set is DETERMINISTIC -- exactly the genes passing the gates
#   of fig5_K_examples.R (positive wd contrast on rZ AND access, d_rZ > 0, >= MIN_ON_WD wd
#   meta-cells, wd seed support >= MIN_SUPPORT, loc_wd <= the per-genome quantile threshold),
#   top N_PER_TYPE per type by d_rZ, ALL qualifying types (no top-4 cut). Curating the final
#   figure from this stated pool keeps the "examples selected for the effect" caption honest;
#   hand-adding genes from outside it would not be.
# COORDINATES FROM THE GFF3, NEVER THE MARKER BED: vnd7's marker-BED entry sits 148 kb from its
#   GFF3 gene body, and the At marker BED is all 0-0 placeholders.
# After editing the final picks, mirror them in fig5_K_examples.R via PICK_TYPES / PICK_GENES
#   so the UMAP panel and the browser exhibit agree.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

# -------------------------
# CONFIG
# -------------------------
DATA     <- "data/processed/scifiATAC_B73_Arabidopsis"
FIG5     <- "figures/main/fig5"
BROWSER  <- file.path(FIG5, "browser")
CONTRAST <- file.path(FIG5, "Fig5_P3E_examples_contrast.tsv")
SELECT   <- file.path(FIG5, "Fig5_P3E_examples_selection.tsv")
OUT      <- file.path(BROWSER, "plot_loci.tsv")
GINFO    <- file.path(DATA, "socrates", "_data", "_GenomeInfo")
dir.create(BROWSER, showWarnings = FALSE, recursive = TRUE)

# gates -- MUST mirror fig5_K_examples.R exactly
MIN_ON_WD   <- 2L
MIN_SUPPORT <- 0.5
N_PER_TYPE  <- 3L            # audition cap per type (top by d_rZ; the UMAP
                             # panel's top-2 are always contained in the top-3)
# pad: ~1 kb around the gene body (wider windows, 8/4 kb, wasted panel space;
# tight windows read better and the PDFs are smaller)
PAD         <- c(B73v5 = 1000L, TAIR10 = 1000L)

# the gene annotations (GFF3) of the two references, placed under _GenomeInfo/
GFF <- c(B73v5  = file.path(GINFO, "B73v5.gff3"),
         TAIR10 = file.path(GINFO, "TAIR10.gff3"))

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", sub("^[A-Za-z]+:", "", x))
# (in_umap_panel flags refresh whenever this script reruns after the UMAP
#  figure's PICK_* change -- coordinates and the audition pool are unaffected)

# --- candidate pool: the panel-K gates, all qualifying types ------------------
W <- fread(CONTRAST)
if (!all(c("loc_thr", "wd_type_support") %in% names(W)))
  stop(CONTRAST, " predates the gated selection -- rerun fig5_K_examples.R")
cand <- W[c_rZ_wd > 0 & c_acc_wd > 0 & d_rZ > 0 & n_on_wd >= MIN_ON_WD &
          (is.na(wd_type_support) | wd_type_support >= MIN_SUPPORT) &
          !is.na(loc_wd) & loc_wd <= loc_thr]
if (nrow(cand) == 0L) stop("no gate-passing candidates -- wrong contrast TSV?")
setorder(cand, genome, type_label, -d_rZ)
cand <- cand[, head(.SD, N_PER_TYPE), by = .(genome, type_label)]
cand[, rank_in_type := seq_len(.N), by = .(genome, type_label)]

umap_picks <- character()
if (file.exists(SELECT)) umap_picks <- fread(SELECT)[, paste(genome, geneID)]
cand[, in_umap_panel := paste(genome, geneID) %in% umap_picks]

# --- gene coordinates from the GFF3s -----------------------------------------
read_genes <- function(g) {
  f <- GFF[[g]]
  if (!file.exists(f)) stop(g, ": no readable GFF3 at ", f,
                            " (place the reference gene annotation there)")
  message("  [gff] ", g, " <- ", f)
  gf <- fread(cmd = sprintf("awk -F'\t' '$3 == \"gene\"' %s", shQuote(f)),
              header = FALSE,
              col.names = c("chrom", "src", "feat", "start", "end",
                            "score", "strand", "frame", "attr"))
  gf[, geneID := sub("^gene:", "", sub("^.*ID=([^;]+).*$", "\\1", attr))]
  gf[, .(geneID, chrom, start, end, strand)]
}
coords <- rbindlist(lapply(unique(cand$genome), function(g)
  read_genes(g)[, genome := g]))

loci <- merge(cand, coords, by = c("genome", "geneID"), all.x = TRUE, sort = FALSE)
if (anyNA(loci$start))
  stop("gene(s) missing from the GFF3: ",
       paste(loci[is.na(start), paste0(genome, ":", geneID)], collapse = ", "))
if (anyDuplicated(loci[, .(genome, geneID)]) > 0)
  stop("duplicated gene after the GFF3 join")

# --- write plot_loci.tsv (locus names unique; groups match the bw tracks) -----
loci[, locus := safe(name)]
loci[duplicated(locus) | duplicated(locus, fromLast = TRUE),
     locus := paste0(locus, "_", safe(geneID))]
if (anyDuplicated(loci$locus) > 0) stop("locus names still not unique")
out <- loci[, .(locus, genome, group = safe(type_label), chrom, start, end,
                strand, pad = PAD[genome],
                d_rZ = round(d_rZ, 3), wd_support = wd_type_support,
                rank_in_type, in_umap_panel,
                note = sprintf("%s %s | c_rZ %s->%s | GFF3 coords", geneID,
                               type_label, round(c_rZ_PreClean, 2),
                               round(c_rZ_wd, 2)))]
setorder(out, genome, group, rank_in_type)
fwrite(out, OUT, sep = "\t")

cat("=== browser audition loci ===\n\n")
print(out[, .(loci = .N, in_umap_panel = sum(in_umap_panel)),
          by = .(genome, group)])
cat("\n", nrow(out), " loci -> ", OUT, "\n",
    "  gates identical to fig5_K_examples.R; top ", N_PER_TYPE,
    "/type by d_rZ; ALL qualifying types (no top-4 cut).\n",
    "  next: sbatch analysis/fig5_biological_impact/browser/plot.sh   (after the bigwigs)\n", sep = "")
