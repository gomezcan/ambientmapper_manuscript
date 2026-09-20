#!/usr/bin/env Rscript
# =============================================================================
# fig5_K_examples.R  -  Fig 5 panel K: showcase marker genes as meta-cell accessibility UMAPs,
#   PreClean | WD | ND, for cell types and genes with a Pre -> WD GAIN of on-type signal.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/{rZ_annotation/consensus,consensus}/
#         .../SM2v2_plate/step3_compare/<prefix>.plate.perkb.genes.sparse.rds, step2_cluster metadata (via read_meta)
#         .../socrates/_data/markers/markers.{maize,At}.informative_top15.bed
#         figures/main/fig5/Fig5_P3D_marker_access_pooled.tsv (fig5_J_access_cache.R), Fig5_P3B_annotation_metacells.tsv (fig5_J_annotation_table.R)
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P3E_examples_{B73v5,TAIR10}.{pdf,png} + _selection.tsv + _contrast.tsv (read by browser/)
# Run     Rscript analysis/fig5_biological_impact/fig5_K_examples.R   (repo root, ~3 min; run the two J scripts first)
# Part 3 internal lettering kept in titles and cross-references: E = this panel = manuscript K.
# =============================================================================
#
# EXAMPLES SELECTED FOR THE EFFECT -- ILLUSTRATIVE, NEVER EVIDENCE.
#   Genes are chosen by the largest Pre->wd gain in on-type contrast, so this
#   figure cannot also be the proof that cleaning helps; the quantitative
#   claims live in the aggregate panels (P3 B/D, S7, S8). The selection rule
#   is printed on the figure.
#
# SELECTION RULE (deterministic, no hand-picking):
#   For each marker gene and stage: on-type contrast =
#       mean(rZ over meta-cells CALLED the gene's type) - mean(rZ elsewhere)
#   from the consensus marker-rZ matrices; the same contrast is computed on the
#   per-gene 0-1 scaled pooled accessibility (cache of fig5_J_access_cache.R)
#   as a cross-check. Candidate genes must have their type CALLED in all three
#   stages and a POSITIVE wd contrast on both metrics; they are ranked by
#   delta_rZ = contrast(wd) - contrast(Pre). A SPATIAL gate additionally
#   requires the gene's wd signal to form a compact territory in the embedding
#   -- see the TOP_FRAC/MAX_GENE_LOC constants. Types are ranked by the mean
#   delta_rZ of their top-2 genes; the top N_TYPES types x top N_GENES genes are
#   drawn. Override with PICK_TYPES / PICK_GENES below (the shipped panel uses
#   the manual picks set there, all of them gate-passers).
#
# WHAT EACH UMAP SHOWS (META-CELL RENDERING: a per-cell rendering stayed
# unreadable, the dropout speckle defeats it; pooling into meta-cells is what
# defeats dropout, the same reason the B|D access half is legible):
#   META-CELLS ONLY -- the grey per-cell silhouette is not drawn (it was ~98% of
#   the PDF: ~25k maize cells x 9 facets as vector circles). The panel EXTENTS
#   are still set from the full cell set via an invisible geom_blank on the
#   per-facet range corners, so every dot keeps the exact position it would have
#   with the silhouette -- nothing is re-scaled or re-zoomed. Each stage is its
#   OWN embedding (never overlaid). Each DOT = one CONSENSUS meta-cell at the
#   MEDIAN position of its member cells, sized by its cell count, filled with the
#   gene's POOLED accessibility from the Access cache: summed over member cells,
#   meta-cell depth-normalised to a common total over ALL genes, then min-max
#   scaled 0-1 per gene ACROSS the genome's three stages -- the same quantity as
#   the B|D access half, at gene level. Shared-within-a-gene scale = the Pre vs
#   wd comparison is real; no scale is shared BETWEEN genes. BLACK OUTLINE =
#   meta-cell called the gene's type in that stage. The per-cell normalisation
#   question disappears here: the common-total scaling handles At's -75% read
#   removal at meta-cell level.
#
# OTHER TRAPS:
#   Types come from ARGMAX calls (the expression-matched gene null finds no
#     significant meta-cell x type pair in five of the six objects and 2 of 5,720
#     in maize ND) -- the type name on the header is descriptive annotation, not
#     tested identity.
#   wd and nd are INDEPENDENT treatments of the same raw input, not a chain.
#   The At marker BED carries DUMMY coordinates (chr "At", 0-0) -- the browser
#     step needs real TAIR10 loci from the annotation; the selection TSV flags
#     this. maize coordinates are real.
#   Contrast uses each stage's OWN calls (partitions differ per stage).
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2); library(patchwork)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source("analysis/_helpers/fig5_part2_helpers.R")   # SOC, PLATE, CFG, STAGES, read_meta()

STEM    <- "Fig5_P3E_examples"
MCD     <- file.path(PLATE, "step5_metacell")
RZDIR   <- file.path(MCD, "rZ_annotation", "consensus")
CONS    <- file.path(MCD, "consensus")
PERKB   <- file.path(PLATE, "step3_compare")
CACHE   <- file.path(OUTDIR, "Fig5_P3D_marker_access_pooled.tsv")
SUPPORT <- file.path(OUTDIR, "Fig5_P3B_annotation_metacells.tsv")
BEDS    <- c(B73v5  = file.path(SOC, "_data", "markers", "markers.maize.informative_top15.bed"),
             TAIR10 = file.path(SOC, "_data", "markers", "markers.At.informative_top15.bed"))
METRIC  <- "geom"
N_TYPES <- 4L      # example types per genome
N_GENES <- 2L      # genes per type
NORM_SC <- 1e4     # per-cell normalisation scale (used by the loc gate only)
DISP_Q  <- 0.95    # colour squish: fill = scaled / per-gene q95, values above
                   # squished to the top colour (a strict min-max let ONE outlier
                   # meta-cell crush a gene's whole panel to near-white, e.g.
                   # ccdp_1; same convention class as the rZ panels' p99 clip)
# Stability gates -- without them the first build picked a type carried by ONE
# meta-cell at 3% seed support (Zm:developing_pavement_cell), which is exactly
# what must not go to the genome browser as a showcase.
MIN_ON_WD   <- 2L    # >= this many meta-cells called the type in wd
MIN_SUPPORT <- 0.5   # mean wd label_support of the type's meta-cells
# Spatial-coherence gate (the first, ungated picks lit cells scattered all over
# the UMAP). A gene's wd TERRITORY = its top TOP_FRAC cells by per-cell-normalised
# signal; loc = median spread of those cells around their own (median) centre /
# median spread of ALL cells around the embedding centre. loc << 1 = a coherent
# patch, loc ~ 1 = scattered like background. Gated on wd ONLY -- Pre may
# legitimately be scattered (that tightening is part of the story; both stages
# are reported). Genes with < MIN_TOP positive wd cells fail the gate as
# unshowable.
# The gate is a PER-GENOME QUANTILE, not an absolute cut. Measured on the first
#   gated build: maize loc q10/50/90 = 0.50/0.93/1.12 but At = 0.67/0.95/1.04 --
#   At is a CONTINUUM, nothing there ever reaches maize-island compactness, and
#   an absolute 0.6 emptied the At candidate set. Keeping the most compact LOC_Q
#   share per genome adapts to each embedding's geometry.
TOP_FRAC <- 0.02
MIN_TOP  <- 15L      # 30 blanked most At markers (loc=NA): <30 positive cells
                     # in a 1k-cell object is common there, 15 still defines a
                     # territory and keeps At selectable
LOC_Q    <- 0.30     # keep genes in the most compact 30% of their genome

# Manual overrides (NULL = automatic selection by the rule above).
# FINAL PICKS: types and genes chosen by eye from the genome-browser audition
# (browser/loci.R renders every gate-passing locus; see browser/plot_loci.tsv
# and browser/plots/ under OUTDIR). All five are gate-passers, so the stated
# selection rule still bounds the pool; the browser tracks did the last-mile
# curation. Restore NULL to return to the automatic top-4 rule.
PICK_TYPES <- list(
  B73v5  = c("Zm:companion_cells", "Zm:leaf_primordia", "Zm:stomatal_precursor"),
  TAIR10 = c("At:bundle_sheath"))
PICK_GENES <- list(
  B73v5  = c("Zm00001eb115210",    # companion_cells,   c_rZ 2.05 -> 3.16
             "Zm00001eb220660",    # zyb9, leaf_primordia
             "Zm00001eb299050"),   # stomatal_precursor, c_rZ 0.36 -> 1.39
  TAIR10 = c("AT2G37170",          # PIP2B, bundle_sheath
             "AT5G41920"))         # SCL23, bundle_sheath

G_LEV <- vapply(CFG, `[[`, character(1), "label")
S_LEV <- names(STAGES)
FIXED <- c("metacell", "seed", "cluster", "n_cells", "top_type", "top_rZ", "margin")

# RdPu ramp for the meta-cell fill -- full 0-1, near-white at the gene's own
# minimum, dark purple at its maximum.
MC_PAL <- c("#FFF7F3", "#FDE0DD", "#FBB4B9", "#F768A1", "#DD3497", "#AE017E",
            "#49006A")
HDRCOL <- "grey25"
# Per-cell grey silhouette: OFF (it was the whole PDF-size problem). TRUE
# restores it exactly; with KEEP_CELL_EXTENTS the extents are identical either
# way, so the two renderings overlay dot-for-dot.
SHOW_CELLS <- FALSE
# TRUE  = panels keep the FULL cell set's ranges, so every meta-cell sits where
#         it did in the silhouetted version. Cost: the margin the grey cells
#         used to fill is now blank -- very visible on At, whose meta-cell
#         centroids occupy a small part of an embedding stretched by outliers.
# FALSE = panels range on the meta-cells only (auto), filling the space at the
#         price of a per-stage zoom that is no longer 1:1 with the silhouetted figure.
KEEP_CELL_EXTENTS <- TRUE

prefix_for <- function(g, s) sprintf("%s_%s", STAGES[[s]], CFG[[g]]$suffix)
if (!file.exists(CACHE))
  stop("Access cache missing -- run fig5_J_access_cache.R first: ", CACHE)

# =============================================================================
# 1. LOAD marker rZ matrices + calls; compute the on-type contrasts
# =============================================================================
mk_ids <- list(); rzmat <- list(); calls <- list()
for (g in names(CFG)) for (s in S_LEV) {
  pf <- prefix_for(g, s)
  ft <- file.path(RZDIR, sprintf("%s.permetacell_type_rZ.%s.tsv", pf, METRIC))
  fm <- file.path(RZDIR, sprintf("%s.permetacell_marker_rZ.%s.tsv", pf, METRIC))
  for (f in c(ft, fm)) if (!file.exists(f)) stop("missing upstream input: ", f)
  rz <- fread(ft)
  calls[[pf]] <- rz[, .(metacell, top_type)]
  m  <- fread(fm)
  mk_ids[[pf]] <- m[, .(geneID, name, type_label)][, `:=`(genome = g, stage = s)]
  M <- as.matrix(m[, -c("geneID", "name", "type_label")])
  rownames(M) <- m$geneID
  if (!setequal(colnames(M), rz$metacell)) stop(pf, ": marker/type meta-cells differ")
  rzmat[[pf]] <- M
}

# shared marker set per genome (same rule as the other Part-3 scripts)
gmap <- list()
for (g in names(CFG)) {
  sets   <- lapply(S_LEV, function(s) mk_ids[[prefix_for(g, s)]]$geneID)
  shared <- Reduce(intersect, sets)
  mk <- unique(rbindlist(mk_ids)[genome == g & geneID %in% shared,
                                 .(geneID, name, type_label)])
  if (anyDuplicated(mk$geneID) > 0) stop(g, ": a marker carries two type labels")
  gmap[[g]] <- mk[, genome := g]
}
gmap <- rbindlist(gmap)

# scaled pooled accessibility (per-gene 0-1 across stages) from the cache
pooled <- fread(CACHE)
rngs   <- pooled[, .(lo = min(access), hi = max(access)), by = .(genome, geneID)]
rngs[, flat := hi <= lo]
pooled <- merge(pooled, rngs, by = c("genome", "geneID"))
pooled[, scaled := fifelse(flat, 0, (access - lo) / (hi - lo))]

contrast_one <- function(g, s) {
  pf <- prefix_for(g, s)
  M  <- rzmat[[pf]]; cl <- calls[[pf]]
  A  <- dcast(pooled[genome == g & stage == s], geneID ~ metacell, value.var = "scaled")
  Am <- as.matrix(A[, -1]); rownames(Am) <- A$geneID
  gm <- gmap[genome == g]
  rbindlist(lapply(split(gm, by = "type_label"), function(gg) {
    on <- cl[top_type == gg$type_label[1], metacell]
    if (length(on) == 0L || length(on) == ncol(M)) return(NULL)
    off <- setdiff(cl$metacell, on)
    ids <- gg$geneID
    data.table(genome = g, stage = s, geneID = ids, type_label = gg$type_label[1],
               n_on = length(on),
               c_rZ  = rowMeans(M [ids, on,  drop = FALSE]) -
                       rowMeans(M [ids, off, drop = FALSE]),
               c_acc = rowMeans(Am[ids, on,  drop = FALSE]) -
                       rowMeans(Am[ids, off, drop = FALSE]))
  }))
}
ctr <- rbindlist(lapply(names(CFG), function(g)
  rbindlist(lapply(S_LEV, function(s) contrast_one(g, s)))))

W <- dcast(ctr, genome + geneID + type_label ~ stage,
           value.var = c("c_rZ", "c_acc", "n_on"))
W <- W[!is.na(c_rZ_PreClean) & !is.na(c_rZ_wd) & !is.na(c_rZ_nd)]   # called in all 3
W[, `:=`(d_rZ  = c_rZ_wd  - c_rZ_PreClean,
         d_acc = c_acc_wd - c_acc_PreClean)]
W <- merge(W, gmap[, .(genome, geneID, name)], by = c("genome", "geneID"))

# type-level wd seed support -- needed BEFORE selection (it is a gate)
if (file.exists(SUPPORT)) {
  sup <- fread(SUPPORT)[stage == "wd",
                        .(wd_type_support = round(mean(label_support), 3)),
                        by = .(genome, type_label = top_type)]
  W <- merge(W, sup, by = c("genome", "type_label"), all.x = TRUE, sort = FALSE)
} else {
  message("  [support] ", basename(SUPPORT), " absent -- the support gate is OFF",
          " (run fig5_J_annotation_table.R to enable it)")
  W[, wd_type_support := NA_real_]
}

# =============================================================================
# 2. PER-CELL VALUES for ALL shared markers (one perkb load per arm) and the
#    spatial-coherence metrics
# =============================================================================
VN <- list(); MD <- list()
for (g in names(CFG)) {
  ids <- gmap[genome == g, geneID]
  for (s in S_LEV) {
    pf <- prefix_for(g, s)
    frds <- file.path(PERKB, sprintf("%s.plate.perkb.genes.sparse.rds", pf))
    if (!file.exists(frds)) stop("missing perkb matrix: ", frds)
    t0 <- proc.time()[3]
    X  <- readRDS(frds)
    md <- read_meta(g, s)                       # step2_cluster cells + UMAP
    miss <- setdiff(md$cellID, colnames(X))
    if (length(miss) > 0)
      stop(pf, ": ", length(miss), " clustered cell(s) absent from the perkb matrix")
    bad <- setdiff(ids, rownames(X))
    if (length(bad) > 0)
      stop(pf, ": marker gene(s) absent from the perkb matrix: ",
           paste(bad, collapse = ", "))
    tot <- Matrix::colSums(X[, md$cellID, drop = FALSE])
    if (any(tot <= 0)) stop(pf, ": a clustered cell has zero gene-space signal")
    V <- as.matrix(X[ids, md$cellID, drop = FALSE])   # markers x cells
    rm(X); gc(verbose = FALSE)
    VN[[pf]] <- sweep(V, 2L, tot / NORM_SC, "/")      # per-cell normalised
    MD[[pf]] <- md
    message(sprintf("  [cells] %-28s %s cells x %d markers  (%.1f s)",
                    pf, format(nrow(md), big.mark = ","), length(ids),
                    proc.time()[3] - t0))
  }
}

# gene territory localisation (Pre + wd; the gate reads wd)
loc_arm <- function(g, s) {
  pf <- prefix_for(g, s)
  V <- VN[[pf]]; md <- MD[[pf]]
  gx <- md$umap1; gy <- md$umap2
  sp_all <- median(sqrt((gx - median(gx))^2 + (gy - median(gy))^2))
  K <- max(MIN_TOP, ceiling(TOP_FRAC * ncol(V)))
  rbindlist(lapply(rownames(V), function(id) {
    v <- V[id, ]; pos <- which(v > 0)
    if (length(pos) < MIN_TOP)
      return(data.table(genome = g, stage = s, geneID = id,
                        loc = NA_real_, n_pos = length(pos)))
    top <- pos[order(v[pos], decreasing = TRUE)[seq_len(min(K, length(pos)))]]
    cx <- median(gx[top]); cy <- median(gy[top])
    data.table(genome = g, stage = s, geneID = id,
               loc = median(sqrt((gx[top] - cx)^2 + (gy[top] - cy)^2)) / sp_all,
               n_pos = length(pos))
  }))
}
locs <- rbindlist(lapply(names(CFG), function(g)
  rbindlist(lapply(c("PreClean", "wd"), function(s) loc_arm(g, s)))))
locW <- dcast(locs, genome + geneID ~ stage, value.var = "loc")
setnames(locW, c("PreClean", "wd"), c("loc_PreClean", "loc_wd"))

# type territory compactness in the wd embedding (reported, not gated -- the
# gene-level gate is what shapes the picks)
compact <- rbindlist(lapply(names(CFG), function(g) {
  pf <- prefix_for(g, "wd")
  cc <- fread(file.path(CONS, sprintf("%s.consensus.cells.tsv", pf)))[
    , .(cellID, metacell = paste0("cMC-", consensus_mc))]
  cc <- merge(cc, calls[[pf]], by = "metacell")
  n0 <- nrow(cc)
  cc <- merge(cc, MD[[pf]][, .(cellID, umap1, umap2)], by = "cellID")
  if (nrow(cc) < 0.95 * n0)
    message(sprintf("  [compact] %s: %.1f%% of consensus cells missing UMAP coords",
                    pf, 100 * (1 - nrow(cc) / n0)))
  gx <- MD[[pf]]$umap1; gy <- MD[[pf]]$umap2
  sp_all <- median(sqrt((gx - median(gx))^2 + (gy - median(gy))^2))
  cc[, .(compact_wd = round(median(sqrt((umap1 - median(umap1))^2 +
                                        (umap2 - median(umap2))^2)) / sp_all, 3),
         n_type_cells = .N), by = .(type_label = top_type)][, genome := g]
}))

W <- merge(W, locW,    by = c("genome", "geneID"),    all.x = TRUE, sort = FALSE)
W <- merge(W, compact, by = c("genome", "type_label"), all.x = TRUE, sort = FALSE)
W[, loc_thr := quantile(loc_wd, LOC_Q, na.rm = TRUE), by = genome]
cat("-- wd territory localisation (loc: spread of top-signal cells / embedding spread) --\n")
print(W[, .(n_genes = .N,
            n_loc_NA = sum(is.na(loc_wd)),
            loc_thr = round(loc_thr[1], 3),
            gated_out = sum(is.na(loc_wd) | loc_wd > loc_thr),
            loc_wd_q10_50_90 = paste(round(quantile(loc_wd, c(.1, .5, .9),
                                                    na.rm = TRUE), 2),
                                     collapse = " / ")), by = genome])
cat("\n")

# =============================================================================
# 3. SELECT: top N_TYPES types by mean top-2 gene d_rZ; top N_GENES genes each
# =============================================================================
sel <- list()
for (g in names(CFG)) {
  # d_rZ > 0: a "gain" showcase must not contain a loss -- without this gate
  # sweet13a (d_rZ = -0.48) rode in as its type's second-best gene.
  cand <- W[genome == g & c_rZ_wd > 0 & c_acc_wd > 0 & d_rZ > 0 &
            n_on_wd >= MIN_ON_WD &
            (is.na(wd_type_support) | wd_type_support >= MIN_SUPPORT) &
            !is.na(loc_wd) & loc_wd <= loc_thr]
  if (nrow(cand) == 0L) stop(g, ": no candidate genes pass the positive-wd gate")
  trank <- cand[order(-d_rZ), .(score = mean(head(d_rZ, N_GENES)), n_genes = .N),
                by = type_label][order(-score)]
  tt <- PICK_TYPES[[g]]
  if (is.null(tt)) tt <- trank[n_genes >= N_GENES, head(type_label, N_TYPES)]
  if (length(tt) < N_TYPES)
    message(sprintf("  [select] %s: only %d type(s) qualify (asked %d)",
                    g, length(tt), N_TYPES))
  bad <- setdiff(tt, cand$type_label)
  if (length(bad) > 0) stop(g, ": picked type(s) not in candidates: ",
                            paste(bad, collapse = ", "))
  gg <- PICK_GENES[[g]]
  picked <- if (is.null(gg))
    cand[type_label %in% tt][order(type_label, -d_rZ), head(.SD, N_GENES),
                             by = type_label]
  else cand[geneID %in% gg]
  picked[, type_order := match(type_label, tt)]
  setorder(picked, type_order, -d_rZ)
  sel[[g]] <- picked
}
sel <- rbindlist(sel)

# marker BED coordinates for the browser step
# the BEDs carry MULTIPLE rows for some genes (maize: 292 rows / 221 genes) --
#   an un-deduplicated merge triplicated pan2 on the first build. One row per
#   gene: the widest span of its entries.
bed <- rbindlist(lapply(names(CFG), function(g) {
  b <- fread(BEDS[[g]])[, .(chr = chr[1], start = min(start), end = max(end),
                            n_bed_entries = .N), by = geneID]
  b[, genome := g]
}))
sel <- merge(sel, bed, by = c("genome", "geneID"), all.x = TRUE, sort = FALSE)
sel[, coords_dummy := genome == "TAIR10"]   # At BED is 0-0 placeholders
if (anyDuplicated(sel[, .(genome, geneID)]) > 0)
  stop("selection has duplicated genes after the joins")
setorder(sel, genome, type_order, -d_rZ)

# =============================================================================
# 4. meta-cell table for the picked genes -- centroid on the stage's embedding,
#    fill = the pooled per-gene 0-1 scaled accessibility (Access cache values)
# =============================================================================
mcfeat <- list()
for (g in names(CFG)) {
  ids <- sel[genome == g, geneID]
  if (length(ids) == 0L) next
  for (s in S_LEV) {
    pf <- prefix_for(g, s)
    cc <- fread(file.path(CONS, sprintf("%s.consensus.cells.tsv", pf)))[
      , .(cellID, metacell = paste0("cMC-", consensus_mc))]
    n0 <- nrow(cc)
    cc <- merge(cc, MD[[pf]][, .(cellID, umap1, umap2)], by = "cellID")
    if (nrow(cc) < 0.95 * n0)
      stop(pf, ": >5% of consensus cells have no UMAP coords -- wrong object?")
    pos <- cc[, .(mx = median(umap1), my = median(umap2), n_cells = .N),
              by = metacell]
    if (nrow(pos) != nrow(calls[[pf]]))
      stop(pf, ": meta-cell count mismatch between consensus and calls")
    pos <- merge(pos, calls[[pf]], by = "metacell")          # + top_type
    v  <- pooled[genome == g & stage == s & geneID %in% ids,
                 .(geneID, metacell, scaled)]
    dt <- merge(v, pos, by = "metacell", allow.cartesian = TRUE)
    if (nrow(dt) != length(ids) * nrow(pos))
      stop(pf, ": meta-cell x gene grid incomplete")
    mcfeat[[pf]] <- dt[, `:=`(genome = g, stage = s)]
  }
}
mcfeat <- rbindlist(mcfeat)

# display contrast: per-gene q95 squish over the gene's three stages together.
# Guard: a gene with > 95% zero meta-cells would get q95 = 0 -- fall back to its
# max so disp stays defined (and 0 if the gene is all-zero).
mcfeat[, q95 := {
  q <- quantile(scaled, DISP_Q)
  if (q <= 0) q <- max(scaled)
  if (q <= 0) q <- 1
  q
}, by = .(genome, geneID)]
mcfeat[, disp := pmin(scaled / q95, 1)]

# =============================================================================
# 5. FIGURES -- per genome: type blocks x (genes x 3 stages), screenshot-style
# =============================================================================
mk_block <- function(g, tt, first_block) {
  ss  <- sel[genome == g & type_label == tt]
  if (anyDuplicated(ss$geneID) > 0) stop(tt, ": duplicated gene in the block")
  d   <- mcfeat[genome == g & geneID %in% ss$geneID]
  d[, gene_f  := factor(geneID, levels = ss$geneID, labels = ss$name)]
  d[, stage_f := factor(stage, levels = S_LEV)]
  d[, on_type := top_type == tt]
  setorder(d, disp)                     # strongest meta-cells drawn last
  # Cell layer, replicated per gene facet. With SHOW_CELLS = FALSE only the
  # per-facet RANGE CORNERS are kept and fed to geom_blank: the panel extents
  # (and therefore every meta-cell's position) are identical to the silhouetted
  # rendering, at 4 rows per facet instead of ~25,000.
  bg <- rbindlist(lapply(S_LEV, function(s)
    MD[[prefix_for(g, s)]][, .(umap1, umap2, stage = s)]))
  if (!SHOW_CELLS)
    bg <- bg[, CJ(umap1 = range(umap1), umap2 = range(umap2)), by = stage]
  bgg <- rbindlist(lapply(seq_len(nrow(ss)), function(i)
    copy(bg)[, geneID := ss$geneID[i]]))
  bgg[, `:=`(gene_f  = factor(geneID, levels = ss$geneID, labels = ss$name),
             stage_f = factor(stage, levels = S_LEV))]
  ptbg   <- if (g == "B73v5") 0.10 else 0.30
  # maize: 286 dots tiled the whole embedding at the At size range -- kept small
  # so they stay individually readable where the meta-cells crowd (this was also
  # what let the grey silhouette show through, when it was drawn)
  sz_rng <- if (g == "B73v5") c(0.6, 2.4) else c(1.3, 4.5)
  ggplot() +
    (if (SHOW_CELLS)
       geom_point(data = bgg, aes(umap1, umap2), colour = "grey90",
                  size = ptbg, stroke = 0)
     else if (KEEP_CELL_EXTENTS) geom_blank(data = bgg, aes(umap1, umap2))
     else NULL) +
    geom_point(data = d, aes(mx, my, fill = disp, size = n_cells,
                             colour = on_type), shape = 21, stroke = 0.4) +
    scale_fill_gradientn(colours = MC_PAL, limits = c(0, 1),
                         name = "pooled marker\naccessibility\n(per-gene 0-1,\nscaled to q95;\n>q95 squished)") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey62"),
                        guide = "none") +
    scale_size(range = sz_rng, name = "cells in\nmeta-cell") +
    facet_grid(rows = vars(stage_f), cols = vars(gene_f),
               scales = "free", switch = "y") +
    theme_void(base_size = 8) +
    theme(strip.text.x = element_text(face = "italic", size = 7.2,
                                      margin = margin(1, 0, 2, 0)),
          strip.text.y.left = if (first_block)
            element_text(angle = 90, face = "bold", size = 7.2) else element_blank(),
          panel.spacing = unit(1.2, "pt"),
          legend.position = "right",
          legend.title = element_text(size = 5.4),
          legend.text = element_text(size = 5.2),
          legend.key.size = unit(0.34, "cm"),
          aspect.ratio = 1,
          plot.margin = margin(0, 4, 2, if (first_block) 2 else 4))
}
mk_header <- function(tt) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = HDRCOL) +
    annotate("text", x = 0.5, y = 0.5, label = sub("^[A-Za-z]+:", "", tt),
             colour = "white", fontface = "bold", size = 3.0) +
    coord_cartesian(expand = FALSE) + theme_void()
}

manual_note <- if (any(!vapply(PICK_TYPES, is.null, logical(1))) ||
                   any(!vapply(PICK_GENES, is.null, logical(1))))
  paste0("\n*** FINAL PICKS ARE MANUAL: types and genes chosen by eye from the 21-locus genome-browser audition ",
         "(all gate-passers; browser/plot_loci.tsv + browser/plots/) - the gates above still bound the pool. ***") else ""

sel_rule <- paste0(
  "EXAMPLES SELECTED FOR THE EFFECT - illustrative, not evidence (aggregate claims: panels B/D, S7, S8). Rule: on-type contrast = mean(marker rZ over meta-cells\n",
  "called the gene's type) - mean(elsewhere), per stage, on the consensus meta-cells; candidates need the type called in all 3 stages and positive wd contrast on\n",
  "rZ AND scaled accessibility, >= ", MIN_ON_WD, " wd meta-cells of the type, mean wd seed support >= ", MIN_SUPPORT, ", and a COMPACT wd territory: the spread of the gene's top ",
  100 * TOP_FRAC, "% signal\n",
  "cells (vs the whole embedding's spread) must sit in the most compact ", 100 * LOC_Q, "% of the genome's markers - a per-genome quantile, because At is a continuum and\n",
  "never reaches maize-island compactness (Pre may be scattered; that tightening is part of the story, both values reported). Only genes with delta = contrast(wd) -\n",
  "contrast(PreClean) > 0 qualify, ranked by delta;\n",
  "top ", N_TYPES, " types (by mean of their best ", N_GENES, " genes) x top ", N_GENES, " genes per type. Type names are ARGMAX calls (gene null: readable zero) - descriptive annotation.\n",
  "Each panel: that stage's OWN embedding (never overlaid); META-CELLS ONLY - the per-cell background is not drawn, though the axes still span the full cell set,\n",
  "so the dots sit exactly where they did with it. Each DOT = one CONSENSUS meta-cell at the median position of its member cells, sized\n",
  "by cell count, filled with the gene's POOLED accessibility (summed over member cells, meta-cell depth-normalised to a common total over all genes, min-max scaled\n",
  "0-1 per gene ACROSS the three stages - the same quantity as the B|D access half, at gene level; for COLOUR additionally scaled to the gene's q95 with values above\n",
  "squished, because a strict min-max lets one outlier meta-cell crush the rest of the panel). Pooling defeats the per-cell dropout speckle of the previous\n",
  "rendering; colour is relative to the gene's own max - comparable across a row, never between genes. BLACK OUTLINE = meta-cell called the gene's type in that\n",
  "stage (grey outline = other calls). Meta-cell partitions and positions differ per stage. wd and nd are independent treatments of the same raw input, not a chain.",
  manual_note)

fig_files <- character()
for (g in names(CFG)) {
  tts <- sel[genome == g, unique(type_label)]
  if (length(tts) == 0L) {
    message("  [figure] ", g, ": no qualifying types -- figure skipped")
    next
  }
  cols <- list()
  for (i in seq_along(tts)) {
    blk <- mk_block(g, tts[i], first_block = (i == 1L))
    cols[[i]] <- (mk_header(tts[i]) / blk) + plot_layout(heights = c(0.055, 1))
  }
  comp <- Reduce(`|`, cols) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("E. Gain of on-type marker signal after cleaning - %s", G_LEV[[g]]),
      subtitle = sel_rule,
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                    plot.subtitle = element_text(size = 6.2, colour = "grey35",
                                                 lineheight = 1.28)))
  for (ext in c("pdf", "png")) {
    f <- file.path(OUTDIR, sprintf("%s_%s.%s", STEM, g, ext))
    # default pdf device: cairo_pdf can fail silently without X11 (ggsave only warns)
    if (ext == "pdf") ggsave(f, comp, width = 14.5, height = 7.6, bg = "white")
    else              ggsave(f, comp, width = 14.5, height = 7.6, dpi = 300, bg = "white")
    fig_files <- c(fig_files, f)
  }
}

# =============================================================================
# REPORT + EXPORT
# =============================================================================
cat("=== FIG 5 PANEL K (Part 3 E) | showcase markers: Pre -> wd gain of on-type signal ===\n\n")
show_cols <- intersect(c("genome", "type_label", "name", "geneID",
                         "c_rZ_PreClean", "c_rZ_wd", "d_rZ", "d_acc",
                         "loc_PreClean", "loc_wd", "compact_wd",
                         "n_on_wd", "wd_type_support"), names(sel))
print(sel[, ..show_cols][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\n  contrast = on-type mean - off-type mean (rZ; _acc = scaled accessibility).\n",
    "  d_* = wd - PreClean. n_on_wd = meta-cells called the type in wd.\n", sep = "")

classic <- W[grepl("DCT|SSU|MDH|PDK|PEPC|NADP|CA[0-9]|ME[0-9]", name, ignore.case = TRUE)]
if (nrow(classic) > 0) {
  cat("\n-- classic C4/photosynthesis names in the marker panel (for the browser step) --\n")
  print(classic[order(-d_rZ), .(genome, name, geneID, type_label,
                                d_rZ = round(d_rZ, 3), c_rZ_wd = round(c_rZ_wd, 3))])
}
cat("\n  NOTE: At BED coordinates are 0-0 placeholders (coords_dummy=TRUE) -- the browser\n",
    "  chain pulls real TAIR10 loci from the GFF3 (browser/loci.R). maize coords are real.\n",
    sep = "")

written <- fig_files
f <- file.path(OUTDIR, paste0(STEM, "_selection.tsv"))
fwrite(sel[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)],
       f, sep = "\t")
written <- c(written, f)
f <- file.path(OUTDIR, paste0(STEM, "_contrast.tsv"))
fwrite(W[order(genome, -d_rZ),
         lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)],
       f, sep = "\t")
written <- c(written, f)

cat("\n--- output verification ---\n")
ok <- TRUE
for (f in written) {
  sz   <- if (file.exists(f)) file.size(f) else NA_integer_
  good <- !is.na(sz) && sz > (if (grepl("\\.tsv$", f)) 50 else 1000)
  ok   <- ok && good
  cat(sprintf("  %-4s %-40s %s\n", if (good) "OK" else "FAIL", basename(f),
              if (is.na(sz)) "missing" else format(sz, big.mark = ",")))
}
if (!ok) stop("one or more outputs failed to write")
cat("\n[done] ", STEM, "_{B73v5,TAIR10}.{pdf,png} + 2 TSVs -> ", OUTDIR, "/\n", sep = "")
