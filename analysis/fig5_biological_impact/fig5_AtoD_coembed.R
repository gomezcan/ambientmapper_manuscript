#!/usr/bin/env Rscript
# =============================================================================
# fig5_AtoD_coembed.R  -  Fig 5 panels A to D on the SM2v2 independent (multi-reference)
#   co-embedding, PreClean vs PostClean (WD): A cleaned co-embedding coloured by species,
#   B by cluster, C per-cluster minority-plate fraction, D species mixing obs/exp over the grid.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_indep/coembed/
#           {Pre,Post_wd}_cluster/SM2v2_coembed_<stage>.{updated_metadata,reduced_dimensions}_v7.<cfg>.txt
#           {Pre,Post_wd}_grid/UMAP_grid_scan.minc_50.metrics.tsv                   (panel D)
# Output  figures/main/fig5/Fig5_P1coembed_{A_umap_species,B_umap_cluster,C_percluster,D_mixing}.{pdf,png},
#         a cell-load panel (not a manuscript panel), the composite and 4 TSVs
# Run     Rscript analysis/fig5_biological_impact/fig5_AtoD_coembed.R      (from the repo root)
# =============================================================================
#
# THE CONSTRUCTION ("extend the row"): every barcode was mapped to BOTH references, so the two
# per-genome 500 bp-tile matrices are stacked (union of cells, zero-fill, genome-prefixed
# features) into one object. `Genome` = plate-of-origin design ground truth (At / B73).
# Built upstream by workflows/05_qc_and_embedding/part1_indep/3_1_coembed_build.R; QC schema
# grafted by 3_1b; clustered by 3_3 at a config shared by both stages. Fig 1F to H are the
# pre-clean view of this same object, so Fig 5 A to D live in the same space and measure the
# same quantity: SPECIES mixing, with plate-of-origin as the label.
#
# Panels A/B are POST-ONLY on purpose. The pre-clean UMAP of this object IS Fig 1F -- same
# co-embed, same config -- so a Pre|Post facet would re-print a published panel. The Pre->Post
# contrast is carried quantitatively by C and D, which both load Pre; A/B instead show what the
# cleaned space looks like.
#
# Pre and Post are SEPARATE embeddings (re-TF-IDF -> re-SVD -> re-UMAP). Facet them side by
# side; NEVER overlay them, and NEVER pair cluster IDs across stages -- cluster 3 in Pre has
# no relationship to cluster 3 in Post. Panel C therefore sorts clusters WITHIN each stage
# (it prints the real cluster ID, but the SORT is what carries the comparison, not the ID).
#
# REPORT obs/exp, NEVER raw mixing. Cleaning changes the minority fraction p by construction,
# so the raw baseline 2p(1-p) collapses mechanically and raw overstates the effect. This is
# the composition confound that invalidated the earlier combined-genome mixing panel.
#
# The metric depends on k AND on the cell set. Two conventions are computed and printed:
#     UMAP space, k=15  -- primary; ties the mixing number to panel A and to the 3_2 grid scan
#     PC space,   k=30  -- cross-check; the convention of the earlier per-reference version
#   Do not quote a number from one convention alongside a number from the other.
# ASCII arrows only in plot text -- the pdf device has no glyph for U+2192.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =============================================================================
# CONFIG
# =============================================================================
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
SOC    <- file.path(DATA, "socrates")
CE     <- file.path(SOC, "SM2v2_indep", "coembed")
OUTDIR <- "figures/main/fig5"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

STAGES       <- c(PreClean = "Pre", PostClean = "Post_wd")   # label -> on-disk stage token
PLOT_STAGES  <- names(STAGES)
MINORITY     <- "At"          # the minority plate-of-origin in this design
K_UMAP       <- 15L           # primary mixing k, in UMAP space (matches the 3_2 grid's diag_k)
K_PC         <- 30L           # cross-check mixing k, in PC space

SP_COLS    <- c(At = "#377eb8", B73 = "#e41a1c")             # house species colours (Fig 1/2/3)
STAGE_COLS <- c(PreClean = "#FF83FA", PostClean = "#43CD80") # house Pre/Post colours

CFG_PIN <- NULL   # e.g. "pcs_20.k_near_30.min_dis_0.3.minc_50.res_0.5"

# =============================================================================
# CONFIG RESOLUTION -- discovered from filenames so an upstream re-run needs no edit
# =============================================================================
stage_dir  <- function(st) file.path(CE, paste0(STAGES[[st]], "_cluster"))
obj_prefix <- function(st) sprintf("SM2v2_coembed_%s", STAGES[[st]])

resolve_cfg <- function(st) {
  d <- stage_dir(st)
  if (!dir.exists(d)) return(character(0))
  head <- paste0(obj_prefix(st), ".updated_metadata_v7.")
  fs   <- list.files(d, pattern = "\\.updated_metadata_v7\\..*\\.txt$")
  hit  <- fs[startsWith(fs, head)]
  unique(sub("\\.txt$", "", substring(hit, nchar(head) + 1L)))
}

cfgs <- lapply(PLOT_STAGES, resolve_cfg)
names(cfgs) <- PLOT_STAGES
missing <- PLOT_STAGES[vapply(cfgs, length, 1L) == 0L]
if (length(missing)) {
  stop(sprintf(paste0("missing co-embed clustering output for stage(s): %s\n",
                      "  expected under: %s\n",
                      "  Run upstream:  PCS=.. KNN=.. MD=.. RES=.. STAGE=%s ",
                      "sbatch workflows/05_qc_and_embedding/part1_indep/3_3_coembed_cluster.sh"),
               paste(missing, collapse = ", "),
               paste(vapply(missing, stage_dir, ""), collapse = ", "),
               STAGES[[missing[1]]]))
}

if (!is.null(CFG_PIN)) {
  if (!all(vapply(cfgs, function(x) CFG_PIN %in% x, TRUE)))
    stop("CFG_PIN '", CFG_PIN, "' is not present for every stage")
  CFG <- CFG_PIN
} else {
  shared <- Reduce(intersect, cfgs)
  if (length(shared) != 1L)
    stop(sprintf(paste0("cannot resolve a single clustering config shared by %s.\n",
                        "  Shared candidates: %s\n  Per stage: %s\n  Set CFG_PIN to choose."),
                 paste(PLOT_STAGES, collapse = " + "),
                 if (length(shared)) paste(shared, collapse = ", ") else "(none)",
                 paste(sprintf("%s={%s}", PLOT_STAGES,
                               vapply(cfgs, paste, "", collapse = "|")), collapse = "  ")))
  CFG <- shared
}
cat("=== Resolved co-embed clustering config (identical across stages, by construction) ===\n")
cat("  ", CFG, "\n\n", sep = "")

# =============================================================================
# LOADERS
# =============================================================================
## fread WITHOUT an explicit header= : these files come from write.table(row.names = TRUE), so
## the header is one field short and fread's auto-detect prepends "V1" for the rowname column.
## Passing header = TRUE defeats that fill logic and the named columns never appear.
ce_file <- function(st, what)
  file.path(stage_dir(st), sprintf("%s.%s_v7.%s.txt", obj_prefix(st), what, CFG))

read_meta <- function(st) {
  f <- ce_file(st, "updated_metadata")
  if (!file.exists(f)) stop("missing metadata: ", f)
  d <- suppressWarnings(fread(f))
  need <- c("cellID", "umap1", "umap2", "LouvainClusters", "total", "Genome")
  if (!all(need %in% names(d)))
    stop("missing columns in ", f, ": ", paste(setdiff(need, names(d)), collapse = ", "))
  # Genome is the plate-of-origin ground truth written by 3_1 -- use it directly rather than
  # parsing the cellID (the co-embed key is genome-stripped).
  if (!all(d$Genome %in% names(SP_COLS)))
    stop("unexpected Genome level(s) in ", f, ": ",
         paste(setdiff(unique(d$Genome), names(SP_COLS)), collapse = ", "))
  d[, stage := st]
  d[, .(cellID, umap1, umap2, LouvainClusters, total, Genome, stage)]
}

read_pcs <- function(st) {
  f <- ce_file(st, "reduced_dimensions")
  if (!file.exists(f)) stop("missing reduced dims: ", f)
  d <- suppressWarnings(fread(f))
  if (!"cellID" %in% names(d)) setnames(d, 1, "cellID")   # rowname column
  d
}

all <- rbindlist(lapply(PLOT_STAGES, read_meta))
all[, stage  := factor(stage,  levels = PLOT_STAGES)]
all[, Genome := factor(Genome, levels = names(SP_COLS))]

cat("=== Cells loaded ===\n")
print(dcast(all[, .N, by = .(stage, Genome)], stage ~ Genome, value.var = "N", fill = 0))
cat("\n")

# =============================================================================
# 1. COMPOSITION / LOAD -- how the two plate libraries survive cleaning
# =============================================================================
load_dt <- all[, .(n_cells = .N,
                   n_minority = sum(Genome == MINORITY),
                   n_major    = sum(Genome != MINORITY)), by = stage]
load_dt[, pct_minority := 100 * n_minority / n_cells]
load_dt[, stage := factor(stage, levels = PLOT_STAGES)]

cat("=== Plate-of-origin load ===\n")
print(load_dt[])
{
  pr <- load_dt[stage == "PreClean"]; po <- load_dt[stage == "PostClean"]
  cat(sprintf("\n  %s-plate cells : %d -> %d  (%+.1f%%)\n", MINORITY,
              pr$n_minority, po$n_minority,
              100 * (po$n_minority - pr$n_minority) / pr$n_minority))
  cat(sprintf("  major-plate cells: %d -> %d  (%+.1f%%)\n",
              pr$n_major, po$n_major, 100 * (po$n_major - pr$n_major) / pr$n_major))
  cat(sprintf("  minority share of object: %.1f%% -> %.1f%%\n",
              pr$pct_minority, po$pct_minority))
  cat("  NOTE: asymmetric removal is the point -- report the SHARE alongside the counts,\n",
      "  because the share is what moves the mixing null 2p(1-p).\n", sep = "")
}

# =============================================================================
# 2. SPECIES MIXING -- obs/exp, two conventions
# =============================================================================
# raw = mean over cells of the fraction of its k nearest neighbours carrying the OTHER
# plate-of-origin label. Under random labelling E[raw] = 2p(1-p), so obs/exp = raw / 2p(1-p)
# and obs/exp = 1.0 means proximity in the embedding carries NO species information.
#
# Exact kNN in base R (kept identical to the earlier per-reference implementation so the
# PC-space numbers remain comparable to it). Squared distance d(i,j) = |xi|^2 + |xj|^2 - 2 xi.xj;
# within a row |xi|^2 is constant, so ranking by (|xj|^2 - 2 xi.xj) is identical and that term is
# dropped. Blocked over rows to bound memory. Exact, not approximate.
knn_cross_frac <- function(X, lab, k, block = 512L) {
  n <- nrow(X); sq <- rowSums(X^2); out <- numeric(n)
  for (start in seq(1L, n, by = block)) {
    idx <- start:min(start + block - 1L, n)
    D <- -2 * tcrossprod(X[idx, , drop = FALSE], X)
    D <- sweep(D, 2L, sq, "+")
    for (i in seq_along(idx)) {
      d <- D[i, ]; d[idx[i]] <- Inf                    # exclude self
      nb <- which(d <= sort.int(d, partial = k)[k])
      if (length(nb) > k) nb <- nb[seq_len(k)]         # ties: keep the first k
      out[idx[i]] <- mean(lab[nb] != lab[idx[i]])
    }
  }
  mean(out)
}

mix_one <- function(st, space, k) {
  meta <- all[stage == st]
  if (space == "umap") {
    X   <- as.matrix(meta[, .(umap1, umap2)])
    lab <- as.integer(meta$Genome == MINORITY)
    ndim <- 2L
  } else {
    pcs <- read_pcs(st)
    m   <- merge(meta[, .(cellID, Genome)], pcs, by = "cellID", sort = FALSE)
    if (nrow(m) != nrow(meta))
      stop("cellID mismatch between metadata and reduced dims for stage ", st)
    pc_cols <- grep("^PC_", names(m), value = TRUE)
    if (!length(pc_cols)) stop("no PC_ columns for stage ", st)
    X    <- as.matrix(m[, ..pc_cols])
    lab  <- as.integer(m$Genome == MINORITY)
    ndim <- length(pc_cols)
  }
  raw <- knn_cross_frac(X, lab, k)
  p   <- mean(lab); expct <- 2 * p * (1 - p)
  data.table(stage = st, space = space, k = k, n_cells = nrow(X), n_dim = ndim,
             p_minority = min(p, 1 - p), raw = raw, expected = expct, obs_exp = raw / expct)
}

mix <- rbindlist(c(
  lapply(PLOT_STAGES, mix_one, space = "umap", k = K_UMAP),
  lapply(PLOT_STAGES, mix_one, space = "pc",   k = K_PC)
))
mix[, stage := factor(stage, levels = PLOT_STAGES)]

cat("\n=== Species mixing (plate-of-origin label), both conventions ===\n")
print(mix[, .(stage, space, k, n_cells, n_dim, p_minority = round(p_minority, 4),
              raw = round(raw, 4), expected = round(expected, 4),
              obs_exp = round(obs_exp, 3))])
for (sp in c("umap", "pc")) {
  a <- mix[space == sp & stage == "PreClean"]; b <- mix[space == sp & stage == "PostClean"]
  cat(sprintf("\n  %-4s space (k=%d):  raw %.3f -> %.3f (%.1fx)   |   obs/exp %.3f -> %.3f (%.2fx)\n",
              sp, a$k, a$raw, b$raw, a$raw / b$raw, a$obs_exp, b$obs_exp, a$obs_exp / b$obs_exp))
}
cat("\n  The obs/exp ratio is the claim; the raw ratio is raw x composition and overstates it.\n")
cat("  Robustness (3_2 grid, 36 configs, its own min.c=50 cell set -- a DIFFERENT cell set,\n")
cat("  so do not quote it in the same sentence as the numbers above): Pre obs/exp 0.845-0.877,\n")
cat("  Post 0.230-0.274, non-overlapping.\n")

mix_plot <- mix[space == "umap"]     # primary convention for the marked point on panel D

# ---- the same quantity across the whole embedding-parameter grid ------------------
## Panel D plots the DISTRIBUTION over the 36 grid combos, not a single bar: the claim is that
## Pre and Post do not overlap under ANY embedding parameters, and a two-bar chart cannot show
## that. The upstream 3_2 grid scan computed these; obs/exp is derived the documented way,
##   centered = raw - 2p(1-p)   =>   expected = raw - centered   =>   obs/exp = raw / expected
## which reproduces the upstream ranges exactly (Pre 0.845-0.877, Post 0.230-0.274).
##
## NO SIGNIFICANCE TEST ON THESE 36 POINTS. They are re-parameterizations of the SAME cells,
##   so a paired test over combos is pseudo-replication -- exactly the error that invalidated
##   the old combined-genome panel C (t-test over 45 combos, p < 2.2e-16 off sd = 0.007).
##   The separation of the ranges IS the evidence. Do not add a p-value here.
##
## DIFFERENT CELL SET from panels A/B/C. The grid ran pre-clustering, so it keeps cells the
##   cluster-size filter later drops (Pre 25,208 vs 24,914; Post 15,309 vs 15,187, ~1%). The two
##   are therefore drawn with DIFFERENT GLYPHS and both n's are stated -- never merge them into
##   one number, and never quote a grid value as if it were this figure's value.
GRID_DIR  <- c(PreClean = "Pre_grid", PostClean = "Post_wd_grid")
GRID_FILE <- "UMAP_grid_scan.minc_50.metrics.tsv"

read_grid <- function(st) {
  f <- file.path(CE, GRID_DIR[[st]], GRID_FILE)
  if (!file.exists(f)) return(NULL)
  d    <- fread(f)
  need <- c("genome_mixing", "genome_mixing_centered", "n_cells", "pcs", "k_near", "min_dist")
  if (!all(need %in% names(d)))
    stop("unexpected grid schema in ", f, ": missing ",
         paste(setdiff(need, names(d)), collapse = ", "))
  d[, expected := genome_mixing - genome_mixing_centered]
  d[, obs_exp  := genome_mixing / expected]
  d[, stage    := st]
  d[]
}
grid <- rbindlist(lapply(PLOT_STAGES, read_grid), fill = TRUE)
has_grid <- nrow(grid) > 0L && uniqueN(grid$stage) == length(PLOT_STAGES)

if (has_grid) {
  grid[, stage := factor(stage, levels = PLOT_STAGES)]
  cat("\n=== Mixing across the embedding grid (obs/exp, ", GRID_FILE, ") ===\n", sep = "")
  print(grid[, .(combos = .N, n_cells = unique(n_cells),
                 min = round(min(obs_exp), 3), median = round(median(obs_exp), 3),
                 max = round(max(obs_exp), 3)), by = stage])
  gp <- grid[stage == "PreClean", range(obs_exp)]; gq <- grid[stage == "PostClean", range(obs_exp)]
  cat(sprintf("  ranges %s: [%.3f, %.3f] vs [%.3f, %.3f]  -- gap %.3f\n",
              if (gp[1] > gq[2]) "DO NOT OVERLAP" else "OVERLAP (!)",
              gp[1], gp[2], gq[1], gq[2], gp[1] - gq[2]))
  cat("  This figure's own config, on its own (clustered) cell set: ",
      sprintf("Pre %.3f / Post %.3f\n", mix_plot[stage == "PreClean", obs_exp],
              mix_plot[stage == "PostClean", obs_exp]), sep = "")
  cat("  ^ different cell sets by ~1% -- plotted as different glyphs, never averaged together.\n")
} else {
  cat("\n!! grid scan not found for both stages -- panel D falls back to a single-value bar\n")
  cat("   expected: ", paste(file.path(CE, GRID_DIR, GRID_FILE), collapse = "\n             "), "\n")
}

# =============================================================================
# 3. PER-CLUSTER COMPOSITION -- the purity claim, within each stage
# =============================================================================
pc_dt <- all[, .(n = .N), by = .(stage, LouvainClusters, Genome)]
pc_dt <- dcast(pc_dt, stage + LouvainClusters ~ Genome, value.var = "n", fill = 0)
setnames(pc_dt, MINORITY, "n_minor", skip_absent = TRUE)
maj <- setdiff(names(SP_COLS), MINORITY)
pc_dt[, tot := n_minor + get(maj)]
pc_dt[, pct_minor := 100 * n_minor / tot]
pc_dt[, stage := factor(stage, levels = PLOT_STAGES)]
# Pre and Post are separate embeddings: cluster IDs are NOT comparable across stages, so rank
# clusters within their own stage rather than pairing IDs.
setorder(pc_dt, stage, -pct_minor)
pc_dt[, rank_in_stage := seq_len(.N), by = stage]

cat("\n=== Per-cluster minority-plate fraction (ranked within stage) ===\n")
print(pc_dt[, .(clusters = .N, min = round(min(pct_minor), 1),
                median = round(median(pct_minor), 1), max = round(max(pct_minor), 1),
                mixed_5_95 = sum(pct_minor > 5 & pct_minor < 95)), by = stage])

# =============================================================================
# PANELS
# =============================================================================
strip_lab <- merge(all[, .(cells = .N), by = stage],
                   load_dt[, .(stage, pct_minority)], by = "stage")
strip_lab[, lab := sprintf("%s  (n = %s, %.1f%% %s-plate)",
                           stage, comma(cells), pct_minority, MINORITY)]
lab_of <- setNames(strip_lab$lab, as.character(strip_lab$stage))
all[, stage_lab := factor(lab_of[as.character(stage)], levels = lab_of[PLOT_STAGES])]

# --- A/B. the CLEANED co-embedding, two colourings ----------------------------
## Post only, by design. The Pre co-embed UMAP IS Fig 1F -- this is the same object at the same
## config -- so facetting Pre|Post here would re-print a published panel. The Pre->Post contrast
## is not lost: it is carried quantitatively by panels D (mixing) and C (per-cluster), which
## both still load Pre. What A/B add instead is the payoff view -- what the cleaned space
## actually looks like, and how its clusters sit in it.
POST_ST <- "PostClean"
post    <- all[stage == POST_ST]
post_n  <- load_dt[stage == POST_ST]
top_at  <- pc_dt[stage == POST_ST][which.max(pct_minor)]   # the consolidated minority island
# The island does NOT hold every minority cell -- a real fraction still sits inside major
# clusters. Compute the split rather than implying "the At cells moved to the island".
at_in_island <- post[Genome == MINORITY & LouvainClusters == top_at$LouvainClusters, .N]
at_scattered <- post_n$n_minority - at_in_island

# draw the minority library last so it is never buried under the major mass
setorderv(post, "Genome", order = if (MINORITY == "At") 1L else -1L)
pA <- ggplot(post, aes(umap1, umap2, colour = Genome)) +
  geom_point(size = .22, alpha = .5) +
  scale_colour_manual(values = SP_COLS, name = "Species (plate of origin)") +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "A. Cleaned co-embedding, by species",
       subtitle = sprintf(paste0("PostClean (WD): n = %s cells, %.1f%% %s-plate.\n",
                                 "%.0f%% of %s cells now sit clear of the maize mass; %.0f%% still\n",
                                 "scatter inside it (cf. Fig 1F, same space, pre-clean)."),
                          comma(post_n$n_cells), post_n$pct_minority, MINORITY,
                          100 * at_in_island / post_n$n_minority, MINORITY,
                          100 * at_scattered  / post_n$n_minority),
       x = "UMAP 1", y = "UMAP 2")

# Cluster IDs are the object's own Leiden labels at the resolved config. They are meaningful
# only WITHIN this stage (Pre is a separate embedding), which is also why panel C never pairs
# them across stages.
post[, cl := factor(LouvainClusters, levels = sort(unique(LouvainClusters)))]
cl_cols <- viridisLite::turbo(nlevels(post$cl), begin = .05, end = .95)

## A Leiden cluster can project as TWO separated lobes -- Leiden partitions the kNN GRAPH,
## not the UMAP plane, so a single community need not be spatially contiguous. Here the At
## island (cluster 7) does exactly that. A single median label would then sit on one lobe and
## leave the other unnumbered, and for a well-separated pair the median can even land in the
## empty gap between them. So: try splitting each cluster in two and keep BOTH labels only when
## the lobes are far apart relative to their own scatter; otherwise fall back to one median.
## Medians, not means, so a cluster with a straggling tail still gets its label on the body.
## TWO conditions, both required, because either alone mislabels (measured on this object):
##   (1) gap > sep_mult x lobe radius -- "far apart". Alone it splits cluster 8 (ratio 6.8),
##       whose second "lobe" is 39 stragglers with no blob for the label to sit on.
##   (2) the smaller lobe is a real share of the cluster -- "both lobes substantial". Alone it
##       splits cluster 6 (339/352) which is one elongated band, merely cut in half.
## Measured ratios here: cl7 14.5 (157/256) . cl8 6.8 (39/270) . cl6 3.3 . all others <= 2.3,
## so only cluster 7 satisfies both -- which is the one a reader would otherwise ask about.
## Over-labelling is the worse failure: a stray duplicate number reads as a bug.
set.seed(1)   # kmeans(nstart=) uses the RNG -- pin it or label positions drift between runs
lobe_labels <- function(dt, min_n = 30L, min_frac = 0.25, sep_mult = 5) {
  rbindlist(lapply(levels(dt$cl), function(g) {
    m   <- as.matrix(dt[cl == g, .(umap1, umap2)])
    one <- data.table(cl = g, umap1 = median(m[, 1]), umap2 = median(m[, 2]))
    if (nrow(m) < 2L * min_n) return(one)
    km     <- kmeans(m, centers = 2L, nstart = 10L)
    gap    <- sqrt(sum((km$centers[1, ] - km$centers[2, ])^2))
    spread <- mean(sqrt(km$withinss / km$size))          # RMS radius of each lobe
    if (gap > sep_mult * spread &&
        min(km$size) >= max(min_n, min_frac * nrow(m)))
      data.table(cl = g, umap1 = km$centers[, 1], umap2 = km$centers[, 2]) else one
  }))
}
cent <- lobe_labels(post)
cent[, cl := factor(cl, levels = levels(post$cl))]

# how the minority actually distributes: island vs still-scattered
n_lobes <- cent[cl == as.character(top_at$LouvainClusters), .N]
cat(sprintf("\n=== %s cells in the Post embedding ===\n", MINORITY))
cat(sprintf("  in the island (cluster %s): %d of %d (%.1f%%), drawn as %d lobe(s)\n",
            top_at$LouvainClusters, at_in_island, post_n$n_minority,
            100 * at_in_island / post_n$n_minority, n_lobes))
cat(sprintf("  still scattered in other clusters: %d (%.1f%%)\n",
            at_scattered, 100 * at_scattered / post_n$n_minority))
cat("  NOTE: the scattered remainder is the co-embed analogue of the combined-genome\n",
    "  'scattered At' population - do not describe the island as containing all of them.\n", sep = "")

pB <- ggplot(post, aes(umap1, umap2, colour = cl)) +
  geom_point(size = .22, alpha = .5) +
  # white halo first, then the glyph (no ggrepel dependency)
  geom_text(data = cent, aes(umap1, umap2, label = cl), inherit.aes = FALSE,
            colour = "white", size = 3.6, fontface = "bold") +
  geom_text(data = cent, aes(umap1, umap2, label = cl), inherit.aes = FALSE,
            colour = "grey15", size = 3.0, fontface = "bold") +
  scale_colour_manual(values = setNames(cl_cols, levels(post$cl)), name = "Cluster") +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), nrow = 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "B. Cleaned co-embedding, by cluster",
       ## Build the subtitle line-by-line and keep each line short: this panel is HALF WIDTH in
       ## the composite and ggplot silently CLIPS an over-long subtitle line at the panel edge.
       ## Same trap as panel D.
       subtitle = paste(c(
         sprintf("%d Leiden clusters at the shared config.", nlevels(post$cl)),
         sprintf("Cluster %s is the %s island (%.1f%% %s-plate, n = %s);",
                 top_at$LouvainClusters, MINORITY, top_at$pct_minor,
                 MINORITY, comma(top_at$tot)),
         sprintf("every other cluster is <= %.1f%% %s-plate.",
                 pc_dt[stage == POST_ST][order(-pct_minor)][2]$pct_minor, MINORITY),
         if (n_lobes > 1)
           sprintf("It projects as %d separated lobes, both labelled.", n_lobes)),
         collapse = "\n"),
       x = "UMAP 1", y = "UMAP 2")

# --- cell load (not a manuscript panel; its counts are quoted in Results) ------
load_long <- rbind(
  load_dt[, .(stage, role = sprintf("%s plate (minority)", MINORITY), n = n_minority)],
  load_dt[, .(stage, role = sprintf("%s plate (major)", maj),        n = n_major)])
load_long[, stage := factor(stage, levels = PLOT_STAGES)]
pC <- ggplot(load_long, aes(stage, n, fill = role)) +
  geom_col(position = position_dodge(width = .75), width = .7, alpha = .9) +
  geom_text(aes(label = comma(n)), position = position_dodge(width = .75),
            vjust = -0.35, size = 2.9) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .14))) +
  scale_fill_manual(values = setNames(c(SP_COLS[[MINORITY]], SP_COLS[[maj]]),
                                      c(sprintf("%s plate (minority)", MINORITY),
                                        sprintf("%s plate (major)", maj))), name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank()) +
  labs(title = "Cell load by plate of origin",
       subtitle = "Removal is strongly asymmetric - it falls on the minority plate",
       x = NULL, y = "Cells")

# --- D. species mixing obs/exp, across the embedding grid ---------------------
## Distribution, not a bar: every grid combo is one embedding of the same cells, so the spread
## shows how little the result depends on the parameter choice, and the gap between the two
## clouds is the claim. Set SHOW_FIG_POINT <- FALSE to drop the this-figure marker and leave a
## pure grid panel (the two use slightly different cell sets -- see the note above).
SHOW_FIG_POINT <- TRUE

if (has_grid) {
  fold <- grid[stage == "PreClean",  median(obs_exp)] /
          grid[stage == "PostClean", median(obs_exp)]
  ## label sits a fixed distance ABOVE the cloud: anchoring it to max(obs_exp) with a vjust
  ## nudge put it inside the violin's KDE tail, overlapping the points it describes.
  gl <- grid[, .(lab = sprintf("%d combos\n%.3f-%.3f", .N, min(obs_exp), max(obs_exp)),
                 y = max(obs_exp) + 0.085), by = stage]
  ## LAYER ORDER MATTERS: ggplot types the x scale from the FIRST layer carrying x. An
  ## annotate() with numeric x placed first makes the scale continuous, and the factor data
  ## then errors with "Discrete value supplied to a continuous scale". geom_hline has no x, so
  ## it is safe anywhere; the annotate must come AFTER a discrete-x layer.
  pD <- ggplot(grid, aes(stage, obs_exp, fill = stage)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = .35, colour = "grey35") +
    ## the violin is the CONTAINER; the individual grid points are the message, so keep the
    ## violin narrow and pale and let the jitter dominate. (A KDE over 36 points spanning
    ## ~0.03 units carries little on its own -- the points and the printed range carry it.)
    geom_violin(width = .5, alpha = .4, colour = NA, scale = "width") +
    geom_jitter(width = .115, height = 0, size = 1.35, alpha = .75, colour = "grey20") +
    stat_summary(fun = median, geom = "crossbar", width = .3, linewidth = .35,
                 colour = "grey15") +   # no fatten= : deprecated in ggplot2 4.0
    ## right-aligned at the panel edge: left-aligned it collided with the PreClean range label,
    ## which sits just under the y=1 line because Pre's cloud is close to random.
    annotate("text", x = 2.52, y = 1.0, label = "random intermingling", hjust = 1, vjust = -0.5,
             size = 2.7, colour = "grey35") +
    geom_text(data = gl, aes(stage, y, label = lab), inherit.aes = FALSE,
              vjust = 0, size = 2.7, colour = "grey25", lineheight = .95) +
    scale_fill_manual(values = STAGE_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 1.22), expand = expansion(mult = c(0, .04))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank()) +
    labs(title = sprintf("D. Species mixing collapses (%.2fx)", fold),
         # keep each line short - this panel is half-width in the composite and long
         # subtitle lines are silently CLIPPED at the panel edge (see panel B).
         subtitle = paste(c(
           "obs/exp = cross-plate neighbours / 2p(1-p);",
           "1.0 = NO species information in the embedding.",
           sprintf("Each dot = one of %d embedding configs", grid[stage == "PreClean", .N]),
           "(pcs x k_near x min_dist); the ranges do not overlap.",
           if (SHOW_FIG_POINT) "Black diamond = this figure's config (see note)."),
           collapse = "\n"),
         x = NULL, y = "Mixing obs / exp")
  if (SHOW_FIG_POINT)
    pD <- pD + geom_point(data = mix_plot, aes(stage, obs_exp), inherit.aes = FALSE,
                          shape = 23, size = 2.9, fill = "black", colour = "black")
} else {
  ## fallback: the single-value bar, if the grid scan is not on disk
  pD <- ggplot(mix_plot, aes(stage, obs_exp, fill = stage)) +
    geom_col(width = .6, alpha = .9) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = .35, colour = "grey35") +
    annotate("text", x = 0.62, y = 1.0, label = "random intermingling", hjust = 0, vjust = -0.5,
             size = 2.7, colour = "grey35") +
    geom_text(aes(label = sprintf("%.3f", obs_exp)), vjust = -0.4, size = 3.1) +
    scale_fill_manual(values = STAGE_COLS, guide = "none") +
    scale_y_continuous(limits = c(0, 1.15), expand = expansion(mult = c(0, .04))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank()) +
    labs(title = sprintf("D. Species mixing collapses (%.2fx)",
                         mix_plot[stage == "PreClean", obs_exp] /
                         mix_plot[stage == "PostClean", obs_exp]),
         subtitle = sprintf(paste0("obs/exp = cross-plate neighbours / 2p(1-p);\n",
                                   "UMAP space, k = %d.\n",
                                   "1.0 = NO species information in the embedding.\n",
                                   "GRID SCAN NOT FOUND - single config only."), K_UMAP),
         x = NULL, y = "Mixing obs / exp")
}

# --- C. per-cluster minority fraction, ranked within stage --------------------
## Bars carry the object's REAL cluster ID, still sorted by fraction within the stage, so the
## Post bars can be read against panel B's labelled UMAP. The sort (not the ID) is what makes
## the two stages comparable -- IDs are per-embedding and are never paired across stages.
pc_dt[, stage_lab := factor(lab_of[as.character(stage)], levels = lab_of[PLOT_STAGES])]
pc_dt[, xkey := factor(paste(stage, LouvainClusters, sep = "|"),
                       levels = paste(stage, LouvainClusters, sep = "|"))]  # pc_dt is presorted
pE <- ggplot(pc_dt, aes(xkey, pct_minor, fill = stage)) +
  geom_col(width = .8, alpha = .9) +
  facet_wrap(~ stage_lab, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = function(v) sub("^[^|]*\\|", "", v)) +
  scale_fill_manual(values = STAGE_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .04))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank()) +
  labs(title = sprintf("C. Per-cluster %s-plate fraction", MINORITY),
       subtitle = paste("Clusters sorted by fraction within their own stage - Pre and Post are",
                        "separate embeddings, so cluster\nIDs cannot be paired across them.",
                        "Pre-clean no cluster is pure; post-clean the mixing concentrates",
                        "\ninto the single cluster labelled in panel B."),
       x = "Cluster ID (sorted within stage)",
       y = sprintf("%% %s-plate cells", MINORITY))

# =============================================================================
# ASSEMBLE
# =============================================================================
composite <- (pA | pB) / (pC | pD) / (pE) +
  plot_layout(heights = c(1.3, 1, 1)) +
  plot_annotation(
    title = "Fig 5 (A to D). Cleaning removes the cross-species signature in the multi-reference space",
    subtitle = paste("SM2v2 individual (multi-reference) co-embedding - the same space as Fig 1F/G/H.",
                     "A/B show the CLEANED object (Fig 1F is its\npre-clean counterpart, so it is not",
                     "repeated here); the load, mixing and per-cluster panels contrast PreClean vs",
                     "decontam-with-design (WD).",
                     "\nPlate of origin is the design ground truth; obs/exp is composition-corrected."),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

## Default pdf device with bg = "white", as in the other Fig 5 scripts. cairo_pdf is avoided:
## on a machine without X11 it fails to load and silently writes NOTHING (ggsave does not
## raise). That is also why plot text uses ASCII "->" and never U+2192 -- the default pdf
## device has no glyph for it.
written <- character(0)
save_both <- function(p, name, w, h, dpi = 300) {
  f_pdf <- file.path(OUTDIR, paste0(name, ".pdf"))
  f_png <- file.path(OUTDIR, paste0(name, ".png"))
  ggsave(f_pdf, p, width = w, height = h, bg = "white", limitsize = FALSE)
  ggsave(f_png, p, width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
  written <<- c(written, f_pdf, f_png)
}

save_both(composite, "Fig5.part1.coembed.composite", 12, 14)
save_both(pA, "Fig5_P1coembed_A_umap_species",  5.8, 5.2)
save_both(pB, "Fig5_P1coembed_B_umap_cluster",  5.8, 5.2)
save_both(pC, "Fig5_P1coembed_load",            5.5, 4.4)
save_both(pD, "Fig5_P1coembed_D_mixing",        5.5, 4.4)
save_both(pE, "Fig5_P1coembed_C_percluster",   10,   4.4)   # manuscript panel C

tsvs <- list(list(mix, "mixing"), list(load_dt, "load"), list(pc_dt, "percluster"))
if (has_grid) tsvs <- c(tsvs, list(list(grid, "mixing_grid")))
for (x in tsvs) {
  f <- file.path(OUTDIR, sprintf("Fig5_P1coembed_%s.tsv", x[[2]]))
  fwrite(x[[1]], f, sep = "\t"); written <- c(written, f)
}

## Verify rather than assert -- a device that fails to load writes nothing but does not error.
ok  <- file.exists(written) & file.size(written) > 0
cat("\n=== Wrote ===\n")
for (i in seq_along(written))
  cat(sprintf("  %-4s %s%s\n", if (ok[i]) "OK" else "FAIL", written[i],
              if (ok[i]) sprintf("  (%s)", format(structure(file.size(written[i]),
                                                            class = "object_size"),
                                                  units = "auto")) else ""))
if (!all(ok)) stop(sum(!ok), " output file(s) were not written -- see FAIL above")
