# =============================================================================
# fig5_part2_helpers.R  -  shared constants and loaders for Fig 5 panels G to L and Fig S8
#   (the plate-split per-genome Socrates objects under socrates/SM2v2_plate/).
# Source from the repo root: source("analysis/_helpers/fig5_part2_helpers.R")
# Sourced by analysis/fig5_biological_impact/fig5_{GH_umap,GH_flow,I_consensus,J_access_cache,
#   J_typeaccess,J_annotation_table,K_examples}.R, browser/groups.R and analysis/supplementary/figS8.R.
# Defines DATA (unless the caller set it first), SOC, PLATE, the stage/genome map (STAGES, MODES,
#   CFG), the origin palette, and read_meta() / load_part2() / load_part2_cached() /
#   cellset_summary(). Needs data.table and ggplot2.
# =============================================================================
#
# THE DESIGN: one object per (library -> its OWN reference), three stages each:
#     PreClean  SM2_*            raw
#     wd        Clean.SM2v2wd_*  AmbientMapper decontam WITH plate design (WD in the manuscript)
#     nd        Clean.SM2v2_*    AmbientMapper decontam design-free       (ND in the manuscript)
#
#   wd and nd are INDEPENDENT treatments of the SAME raw input, NOT a chain.
#   Never draw or describe Pre -> wd -> nd. Flows are Pre -> wd and Pre -> nd,
#   two separate two-stage comparisons sharing one source.
#
# The clustering config is held fixed, identical across the 3 stages within a genome, so any
# stage difference is cleaning/mode, not tuning. The frozen params exist ONLY under mQCv6,
# which is what pins the mQC version.
#
# At is a CONTINUUM: its clusters are provisional GRADIENT BINS, not cell types.
# Each stage is a SEPARATE embedding (re-TF-IDF -> re-SVD -> re-UMAP). Never overlay stages,
#   and never compare cluster IDs by number across stages.
# These metadata files are R-exported with row.names: the header has ONE FEWER field than the
#   data rows (field 1 is a duplicated cellID rowname). fread warns and adds `V1`. Always
#   select columns BY NAME, never by position.
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# -------------------------
# CONFIG (a sourcing script may set DATA before source(); the default is the repo layout)
# -------------------------
if (!exists("DATA")) DATA <- "data/processed/scifiATAC_B73_Arabidopsis"
SOC   <- file.path(DATA, "socrates")
PLATE <- file.path(SOC, "SM2v2_plate")

MQC    <- "mQCv6"
STAGES <- c(PreClean = "SM2", wd = "Clean.SM2v2wd", nd = "Clean.SM2v2")
MODES  <- c("wd", "nd")            # the two independent cleaning modes

CFG <- list(
  B73v5  = list(suffix = "B73_B73v5", lib = "B73", label = "Maize (B73) -> B73v5",
                cfg = "pcs_20.k_near_20.min_dis_0.05.minc_50.res_0.3"),
  TAIR10 = list(suffix = "At_TAIR10", lib = "At",  label = "Arabidopsis (At) -> TAIR10",
                cfg = "pcs_5.k_near_20.min_dis_0.3.minc_50.res_0.3.mclst_40")
)

# Colour by ORIGIN (PreClean cluster), shared across all stages and both plot
# types so a cluster keeps one colour everywhere. Pre has 5 clusters in both
# genomes; the extra slots cover a re-run that finds more.
PRE_PAL  <- setNames(viridisLite::turbo(7, end = .92)[1:7], as.character(1:7))
NEW_COL  <- "grey72"    # cells present post-clean but ABSENT pre-clean
DROP_COL <- "grey50"    # cells present pre-clean but dropped by cleaning

# --- load one (genome, stage) --------------------------------------------------
read_meta <- function(genome, stage) {
  k <- CFG[[genome]]
  f <- file.path(PLATE, "step2_cluster",
                 sprintf("%s_%s.%s.updated_metadata_v7.%s.txt",
                         STAGES[[stage]], k$suffix, MQC, k$cfg))
  if (!file.exists(f)) stop("missing metadata: ", f)
  d <- suppressWarnings(fread(f))          # warns: header 1 short (rowname col) -> V1
  d[, lib := fifelse(grepl("-SM2_At_", cellID), "At",
             fifelse(grepl("-SM2_B73_", cellID), "B73", NA_character_))]
  if (anyNA(d$lib)) stop("unparsed cellID suffix in ", f)
  # ASSERT the plate split actually took: a pooled object carries both libraries (the unsplit
  # SM2v2_indep objects were 34.5% At-well on B73v5). Do not remove.
  bad <- d[lib != k$lib, .N]
  if (bad > 0)
    stop(sprintf("PLATE SPLIT NOT APPLIED: %s/%s has %d cells from the %s library (expected 100%% %s). File: %s",
                 genome, stage, bad, setdiff(c("At", "B73"), k$lib), k$lib, f))
  d[, `:=`(genome = genome, stage = stage)]
  d[, .(cellID, umap1, umap2, LouvainClusters, total, nSites, pTSS, FRiP, genome, stage, lib)]
}

# --- load all 6 (genome x stage), with each cell's PreClean cluster attached ----
# `pre_cl` is the ORIGIN label: the cluster this cell belonged to before cleaning.
# NA => the cell is absent from PreClean (a "gained" cell). Carrying this label
# across stages is what makes colours comparable -- cluster IDs themselves are not.
load_part2 <- function() {
  all <- rbindlist(lapply(names(CFG), function(g)
    rbindlist(lapply(names(STAGES), function(s) read_meta(g, s)))))
  all[, stage := factor(stage, levels = names(STAGES))]
  key <- all[stage == "PreClean", .(genome, cellID, pre_cl = LouvainClusters)]
  all <- merge(all, key, by = c("genome", "cellID"), all.x = TRUE, sort = FALSE)
  all[, pre_lab := factor(fifelse(is.na(pre_cl), "(new)", as.character(pre_cl)),
                          levels = c(as.character(sort(unique(na.omit(all$pre_cl)))), "(new)"))]
  all[]
}

# --- cached loader -------------------------------------------------------------
# The UMAP and flow scripts each need the same 6 files; when both are sourced into
# one R session this loads them once. Returns a copy(): callers add columns (e.g.
# `strip`) and must not mutate the shared cache.
load_part2_cached <- function() {
  if (!exists(".P2_CACHE", envir = globalenv()))
    assign(".P2_CACHE", load_part2(), envir = globalenv())
  copy(get(".P2_CACHE", envir = globalenv()))
}

# --- cell-set accounting: shared / dropped / gained, per genome x mode ----------
# Set membership, NOT net counts: every mode both drops AND gains cells, and a
# net figure hides the gains (B73v5/wd nets -120 while gaining 222).
cellset_summary <- function(all) {
  rbindlist(lapply(names(CFG), function(g) rbindlist(lapply(MODES, function(m) {
    pre <- all[genome == g & stage == "PreClean", cellID]
    po  <- all[genome == g & stage == m,          cellID]
    data.table(genome = g, mode = m, pre = length(pre), post = length(po),
               shared = length(intersect(pre, po)),
               dropped = length(setdiff(pre, po)),
               gained  = length(setdiff(po, pre)))
  }))))[, `:=`(pct_pre_retained = round(100 * shared / pre, 1),
               pct_post_new     = round(100 * gained / post, 1))][]
}
