#!/usr/bin/env Rscript
# =============================================================================
# fig5_GH_flow.R  -  Fig 5 panel H: cluster-flow alluvials, PreClean cluster -> post-clean cluster,
#   per genome (Pre -> WD | Pre -> ND as independent facets), plus the cell-set accounting
#   (shared / dropped / gained; not a manuscript panel) on the plate-split objects.
# Inputs  data/processed/scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step2_cluster/
#           <stage>_<genome>.mQCv6.updated_metadata_v7.<cfg>.txt   (6 files, via the helper)
# Sources analysis/_helpers/fig5_part2_helpers.R
# Output  figures/main/fig5/Fig5_P2_flow_{B73v5,TAIR10}.{pdf,png}, Fig5_P2_retention.{pdf,png},
#         part2_{cellset_accounting,flow_table}.tsv
# Run     Rscript analysis/fig5_biological_impact/fig5_GH_flow.R        (from the repo root)
# =============================================================================
#
#   LEFT  = PreClean clusters (+ a "(new)" source for cells absent pre-clean)
#   RIGHT = post-clean clusters (+ a "(dropped)" sink for cells cleaning removed)
#   Ribbons are coloured by ORIGIN (PreClean cluster) -- the same palette fig5_GH_umap.R
#   uses, so a cluster keeps one colour across panels G and H.
#
# TWO INDEPENDENT FLOWS, NOT A CHAIN. wd and nd are alternative treatments of the SAME raw
#   input, so this draws Pre -> wd and Pre -> nd as separate facets sharing one source. There
#   is no wd -> nd transition and none is implied.
# POST-STAGE NODE BARS ARE NEUTRAL GREY ON PURPOSE. Colouring them from the origin palette
#   would give e.g. "nd cluster 5" the same colour as "PreClean cluster 5" and imply a
#   correspondence that does not exist -- post cluster IDs come from an independent Leiden
#   run. Only ribbons carry origin colour.
# RETENTION IS SET MEMBERSHIP, NOT NET COUNTS. Every mode both drops AND gains cells and a net
#   figure hides the gains (B73v5/wd nets -120 while gaining 222). The gains are not explained:
#   read removal cannot create cells; the suspected mechanism is the `min.t` feature-frequency
#   filter refitting on the cleaned matrix. Retention numbers are therefore descriptive.
# Drawn with hand-rolled ribbons (cosine-eased polygons) so that no alluvial package is needed.
#   Geometry only -- no statistics, no smoothing of the counts themselves.
# ASCII arrows only in plot text -- the pdf device has no glyph for U+2192.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

# -------------------------
# CONFIG
# -------------------------
DATA   <- "data/processed/scifiATAC_B73_Arabidopsis"
OUTDIR <- "figures/main/fig5"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
source("analysis/_helpers/fig5_part2_helpers.R")   # SOC, PLATE, CFG, STAGES, MODES, palettes, loaders

p2 <- load_part2_cached()
cs <- cellset_summary(p2)

cat("=== FIG 5 PANEL H | cell-set accounting (membership, NOT net counts) ===\n")
print(cs)
cat("\n=== Flags carried onto the accounting panel ===\n")
for (i in seq_len(nrow(cs))) with(cs[i], {
  if (gained > 0)
    cat(sprintf("  ! %s/%s gained %d cells (%.1f%% of its post set) that are ABSENT pre-clean.\n",
                genome, mode, gained, pct_post_new))
})
cat("    Read removal cannot create cells - suspect the min.t feature-frequency filter\n",
    "    refitting on the cleaned matrix (not verified upstream).\n",
    "    Corrected retention uses set membership: TAIR10/nd keeps 58.6% of Pre cells,\n",
    "    not the 65.6% a net count implies.\n", sep = "")

# --- flow table: PreClean cluster -> post cluster, per genome x mode -----------
flows_for <- function(g, m) {
  pre <- p2[genome == g & stage == "PreClean", .(cellID, from = as.character(LouvainClusters))]
  po  <- p2[genome == g & stage == m,          .(cellID, to   = as.character(LouvainClusters))]
  f <- merge(pre, po, by = "cellID", all = TRUE)
  f[is.na(from), from := "(new)"]
  f[is.na(to),   to   := "(dropped)"]
  f[, .N, by = .(from, to)]
}

# order nodes: numeric clusters ascending, then the special node last
node_levels <- function(x, special) {
  v <- setdiff(unique(x), special)
  c(as.character(sort(as.numeric(v))), if (special %in% x) special)
}

# stacked y-extent of each node on one axis
node_layout <- function(lv, sizes, gap) {
  y <- 0
  rbindlist(lapply(lv, function(n) {
    h <- sizes[[n]]; r <- data.table(node = n, y0 = y, y1 = y + h); y <<- y + h + gap; r
  }))
}

# one ribbon = cosine-eased polygon from (yl0,yl1) to (yr0,yr1)
ribbon_poly <- function(yl0, yl1, yr0, yr1, id, npts = 60) {
  t <- seq(0, 1, length.out = npts); s <- (1 - cos(pi * t)) / 2
  rbind(data.table(x = t,      y = yl1 + (yr1 - yl1) * s),
        data.table(x = rev(t), y = rev(yl0 + (yr0 - yl0) * s)))[, grp := id][]
}

# geometry for one (genome, mode); `mode` column lets the two modes share a facet
build_flow_geom <- function(g, m) {
  fl <- flows_for(g, m)
  fromlv <- node_levels(fl$from, "(new)")
  tolv   <- node_levels(fl$to,   "(dropped)")
  gap    <- sum(fl$N) * 0.02

  lsz <- fl[, .(n = sum(N)), by = from][, setNames(as.list(n), from)]
  rsz <- fl[, .(n = sum(N)), by = to  ][, setNames(as.list(n), to)]
  L <- node_layout(fromlv, lsz, gap); R <- node_layout(tolv, rsz, gap)

  # sub-stack flows within each node, ordered so ribbons cross as little as possible
  fl[, from := factor(from, levels = fromlv)][, to := factor(to, levels = tolv)]
  setorder(fl, from, to)
  fl[, `:=`(yl0 = L[match(from, node), y0] + cumsum(N) - N,
            yl1 = L[match(from, node), y0] + cumsum(N)), by = from]
  setorder(fl, to, from)
  fl[, `:=`(yr0 = R[match(to, node), y0] + cumsum(N) - N,
            yr1 = R[match(to, node), y0] + cumsum(N)), by = to]
  fl[, id := .I]

  polys <- rbindlist(lapply(seq_len(nrow(fl)), function(i)
    ribbon_poly(fl$yl0[i], fl$yl1[i], fl$yr0[i], fl$yr1[i], fl$id[i])))
  polys <- merge(polys, fl[, .(id, from, to)], by.x = "grp", by.y = "id", sort = FALSE)
  # colour by ORIGIN; flows out of "(new)" and into "(dropped)" get their own greys
  polys[, fillkey := fifelse(from == "(new)", "(new)",
                     fifelse(to == "(dropped)", "(dropped)", as.character(from)))]
  # grp must be unique once the two modes are rbind-ed into one facetted plot
  polys[, `:=`(mode = m, grp = paste(m, grp))]

  bars <- rbind(L[, .(node, y0, y1, side = "Pre")], R[, .(node, y0, y1, side = "Post")])
  bars[, `:=`(xmin = fifelse(side == "Pre", -0.045, 1.005),
              xmax = fifelse(side == "Pre", -0.005, 1.045))]
  bars[, fillkey := fifelse(side == "Pre", as.character(node),
                    fifelse(node == "(dropped)", "(dropped)", "(post)"))]
  bars[, lab := sprintf("%s (%s)", node, format(as.integer(y1 - y0), big.mark = ","))]
  bars[, mode := m]

  list(polys = polys, bars = bars)
}

FLOW_PAL <- c(PRE_PAL, `(new)` = NEW_COL, `(dropped)` = DROP_COL, `(post)` = "grey88")

SUBTITLE_FLOW <- list(
  B73v5 = paste("wd holds the partition almost 1:1; nd re-cuts the same cells into more clusters.",
                "Read with the origin-coloured UMAPs (panel G): the populations themselves stay put,\nso this is where boundaries land, not evidence that structure is lost."),
  TAIR10 = paste("Arabidopsis is a CONTINUUM - cut points are arbitrary and shuffle under any perturbation.",
                 "Without a Leiden seed-stability control (not run),\nthis panel cannot separate a cleaning effect from boundary instability. Provisional.")
)

flow_panel <- function(g, letter) {
  geo   <- lapply(MODES, function(m) build_flow_geom(g, m))
  polys <- rbindlist(lapply(geo, `[[`, "polys"))
  bars  <- rbindlist(lapply(geo, `[[`, "bars"))
  lv    <- sprintf("PreClean -> %s", MODES)      # ASCII: pdf device has no U+2192
  polys[, mode_lab := factor(sprintf("PreClean -> %s", mode), levels = lv)]
  bars[,  mode_lab := factor(sprintf("PreClean -> %s", mode), levels = lv)]

  ggplot() +
    geom_polygon(data = polys, aes(x, y, group = grp, fill = fillkey), alpha = .62) +
    geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = y0, ymax = y1, fill = fillkey),
              colour = "grey25", linewidth = .18) +
    geom_text(data = bars[side == "Pre"],  aes(x = -0.06, y = (y0 + y1) / 2, label = lab),
              hjust = 1, size = 2.4, colour = "grey15") +
    geom_text(data = bars[side == "Post"], aes(x = 1.06, y = (y0 + y1) / 2, label = lab),
              hjust = 0, size = 2.4, colour = "grey15") +
    facet_wrap(~mode_lab, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = FLOW_PAL, guide = "none") +
    scale_x_continuous(limits = c(-0.30, 1.30), breaks = c(0, 1),
                       labels = c("PreClean", "post-clean")) +
    theme_void(base_size = 10) +
    theme(axis.text.x = element_text(size = 8, face = "bold", colour = "grey20",
                                     margin = margin(t = 3)),
          strip.text = element_text(size = 9, face = "bold", margin = margin(b = 4)),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey35", lineheight = 1.2)) +
    labs(title = sprintf("%s. %s  -  cluster flow across cleaning", letter, CFG[[g]]$label),
         subtitle = SUBTITLE_FLOW[[g]])
}

# both genomes are manuscript panel H
pP2_flow_B73v5  <- flow_panel("B73v5",  "H")
pP2_flow_TAIR10 <- flow_panel("TAIR10", "H")

# --- cell-set accounting (not a manuscript panel) -------------------------------
# Pre  bar = shared + dropped   (what cleaning did to the pre-clean set)
# Post bar = shared + gained    (what the post-clean set is actually made of)
# Drawing both makes the gains impossible to hide behind a net count.
ret <- rbindlist(list(
  cs[, .(genome, mode, set = "PreClean",  part = "dropped", n = dropped)],
  cs[, .(genome, mode, set = "PreClean",  part = "shared",  n = shared)],
  cs[, .(genome, mode, set = "post-clean", part = "shared", n = shared)],
  cs[, .(genome, mode, set = "post-clean", part = "gained", n = gained)]
))
ret[, set  := factor(set,  levels = c("PreClean", "post-clean"))]
# position_stack draws the FIRST level on top, so putting `shared` last makes it
# the base of BOTH bars; dropped/gained then read as deltas off a common footing.
ret[, part := factor(part, levels = c("dropped", "gained", "shared"))]
ret[, mode := factor(mode, levels = MODES)]     # wd before nd, not alphabetical
ret[, genome_lab := factor(vapply(genome, function(g) CFG[[g]]$label, ""),
                           levels = vapply(names(CFG), function(g) CFG[[g]]$label, ""))]

lab <- cs[, .(genome, mode, ymax = pmax(pre, post),
              txt = sprintf("%.1f%% of PreClean kept\n%.1f%% of post set is new",
                            pct_pre_retained, pct_post_new))]
lab[, mode := factor(mode, levels = MODES)]
lab[, genome_lab := factor(vapply(genome, function(g) CFG[[g]]$label, ""),
                           levels = levels(ret$genome_lab))]

pP2_retention <- ggplot(ret, aes(set, n, fill = part)) +
  geom_col(width = .62, colour = "grey30", linewidth = .18) +
  geom_text(data = lab, aes(x = 1.5, y = ymax * 1.16, label = txt),
            inherit.aes = FALSE, size = 2.5, colour = "grey25", lineheight = 1.05) +
  facet_grid(genome_lab ~ mode, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(dropped = "#B2182B", shared = "grey78", gained = "#2166AC"),
                    name = NULL, breaks = c("shared", "dropped", "gained"),
                    labels = c(dropped = "dropped by cleaning", shared = "shared (kept)",
                               gained  = "gained (absent pre-clean)")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .22))) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 8, face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey35", lineheight = 1.2)) +
  labs(y = "cells",
       title = "Cell-set accounting - retention is set membership, not net counts",
       subtitle = paste("Every mode both DROPS and GAINS cells. Read removal cannot create cells:",
                        "the gains are not explained\n(suspected min.t feature-frequency refitting) and are not verified upstream.",
                        "Retention numbers are descriptive."))

# --- save individual parts -----------------------------------------------------
save_part <- function(p, name, w, h) {
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
save_part(pP2_flow_B73v5,  "Fig5_P2_flow_B73v5",  10.5, 5.0)
save_part(pP2_flow_TAIR10, "Fig5_P2_flow_TAIR10", 10.5, 5.0)
save_part(pP2_retention,   "Fig5_P2_retention",    7.5, 5.4)

fwrite(cs, file.path(OUTDIR, "part2_cellset_accounting.tsv"), sep = "\t")
ft <- rbindlist(lapply(names(CFG), function(g) rbindlist(lapply(MODES, function(m)
  flows_for(g, m)[, `:=`(genome = g, mode = m)]))))
fwrite(ft, file.path(OUTDIR, "part2_flow_table.tsv"), sep = "\t")

cat("\n[done] wrote Fig5_P2_flow_{B73v5,TAIR10}.*, Fig5_P2_retention.*,",
    "part2_{cellset_accounting,flow_table}.tsv to ", OUTDIR, "/\n", sep = " ")
