#!/usr/bin/env Rscript
# =============================================================================
# browser/groups.R  -  STEP 1 of the Fig 5 panel L genome-browser chain: cell -> group membership,
#   one file per arm, group = the CONSENSUS meta-cell's rZ type call.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/consensus/<prefix>.consensus.cells.tsv
#         .../SM2v2_plate/step5_metacell/rZ_annotation/consensus/<prefix>.permetacell_type_rZ.geom.tsv
#         figures/main/fig5/Fig5_P3E_examples_selection.tsv   (showcase flag, from fig5_K_examples.R)
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/browser/groups/<prefix>.cell2group.tsv (cellID -> group), groups/manifest.tsv
# Run     Rscript analysis/fig5_biological_impact/browser/groups.R    (from the repo root, ~5 s)
# =============================================================================
#
# Pseudobulk browser tracks per cell group, where a cell's group is the TYPE CALLED on its
# consensus meta-cell -- then Pre vs wd (vs nd) coverage is compared at the showcase loci of
# panel K. This script writes the memberships that browser/split.sh uses to produce the
# per-group BED.gz files.
# GROUPS: every type CALLED in that arm (argmax over the consensus meta-cells; partitions and
#   calls are PER ARM -- a type can exist in wd and not in Pre). The splitter also emits an
#   ALLcells aggregate per arm. The manifest flags the panel-K showcase types and thin groups
#   (< MIN_WARN cells -> browser noise).
# ARGMAX calls -- the track NAMES are descriptive annotation, not tested identity. Same
#   standing as the figure labels.
# Never compare a type's track across stages by eye without remembering the groups are
#   different cell sets chosen by different partitions.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

# -------------------------
# CONFIG
# -------------------------
DATA     <- "data/processed/scifiATAC_B73_Arabidopsis"
FIG5     <- "figures/main/fig5"
OUTDIR   <- file.path(FIG5, "browser", "groups")
SELTSV   <- file.path(FIG5, "Fig5_P3E_examples_selection.tsv")
MIN_WARN <- 50L
source("analysis/_helpers/fig5_part2_helpers.R")   # PLATE, CFG, STAGES
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MCD      <- file.path(PLATE, "step5_metacell")
RZDIR    <- file.path(MCD, "rZ_annotation", "consensus")
CONS     <- file.path(MCD, "consensus")

S_LEV <- names(STAGES)
prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", sub("^[A-Za-z]+:", "", x))

showcase <- if (file.exists(SELTSV)) fread(SELTSV)[, unique(type_label)] else {
  message("  [showcase] ", SELTSV, " absent -- no showcase flags")
  character()
}

manifest <- list()
for (g in names(CFG)) for (s in S_LEV) {
  pf   <- prefix_for(g, s)
  f_cc <- file.path(CONS,  sprintf("%s.consensus.cells.tsv", pf))
  f_rz <- file.path(RZDIR, sprintf("%s.permetacell_type_rZ.geom.tsv", pf))
  for (f in c(f_cc, f_rz)) if (!file.exists(f)) stop("missing input: ", f)

  cc <- fread(f_cc)[, .(cellID, metacell = paste0("cMC-", consensus_mc))]
  rz <- fread(f_rz)[, .(metacell, top_type)]
  if (!setequal(unique(cc$metacell), rz$metacell))
    stop(pf, ": consensus meta-cells and rZ calls disagree")
  m <- merge(cc, rz, by = "metacell")
  if (nrow(m) != nrow(cc)) stop(pf, ": call join lost cells")
  if (anyNA(m$top_type))   stop(pf, ": NA type call")
  m[, group := safe(top_type)]
  # two distinct types must never collapse to one safe name
  if (uniqueN(m[, .(top_type, group)]) != uniqueN(m$top_type))
    stop(pf, ": type -> safe-name collapse")

  fwrite(m[, .(cellID, group)],
         file.path(OUTDIR, paste0(pf, ".cell2group.tsv")),
         sep = "\t", col.names = FALSE)

  manifest[[pf]] <- m[, .(n_cells = .N), by = .(type = top_type, group)][
    , `:=`(prefix = pf, genome = g, stage = s,
           is_showcase = type %in% showcase,
           thin = n_cells < MIN_WARN)][]
}
manifest <- rbindlist(manifest)[
  order(genome, factor(stage, levels = S_LEV), -n_cells)]
setcolorder(manifest, c("prefix", "genome", "stage", "type", "group",
                        "n_cells", "is_showcase", "thin"))
fwrite(manifest, file.path(OUTDIR, "manifest.tsv"), sep = "\t")

cat("=== browser groups written ===\n\n")
print(manifest[, .(groups = .N, cells = sum(n_cells),
                   showcase_groups = sum(is_showcase),
                   thin_groups = sum(thin)), by = .(prefix)])
cat("\n  thin = < ", MIN_WARN, " cells: expect noisy tracks; the manifest has\n",
    "  per-group counts. Group names are ARGMAX type calls (descriptive).\n",
    "  Next: bash analysis/fig5_biological_impact/browser/split.sh   (then sbatch browser/bw.sh on the HPC)\n",
    sep = "")
