#!/usr/bin/env Rscript
# 06_48_purity_model.R
#
# Barcode-resolved, depth- & call-class-stratified purity analysis on the
# WASP-corrected BAMs (raw vs clean). Consumes 06_48_extract_barcode_block_counts.py
# output (per-(barcode,block) match counts). Anchor = AM call (top1=genome_1).
#
# Produces (all TSV unless noted):
#   *_barcode_purity.tsv.gz  per-barcode p_top1 / top1or2 / pairfrac + attrs
#   *_stratified.tsv         genotype x class x depth_bin x mode: read-weighted
#                            purity, barcode median/IQR, raw->clean delta
#   *_weak_doublet_diag.tsv  is weak_doublet a mislabeled singlet? (pairfrac,
#                            top1or2, p_top1 vs single_clean)
#   *_decomposition.tsv      within-barcode vs reweighting split of the raw->clean
#                            purity gain, per cell, block-bootstrap 95% CI
#   *_ambiguous_margin.tsv   ambiguous best/second/margin distribution
#   *_agg_betabinom.tsv.gz   (genotype x class x depth_bin x block x mode) match/total
#   *_model_coefs.tsv        glmmTMB beta-binomial coefficients (if glmmTMB present)
#   *_diagnostics.pdf, *_summary.md
#
# Class map: single_clean->singlet, doublet->doublet, weak_doublet->weak_doublet.
# singlet_like := singlet U weak_doublet (evidence-gated; used for the model only).
# Depth bins (on cells_calls n_reads, fixed for both modes): 200-500 ... >10000.

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(ggplot2); library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--raw-prefix",   type = "character"),
  make_option("--clean-prefix", type = "character"),
  make_option("--genotypes",    type = "character"),
  make_option("--sample",       type = "character", default = ""),
  make_option("--out-prefix",   type = "character"),
  make_option("--n-boot",       type = "integer", default = 200L),
  make_option("--display-floor",type = "integer", default = 20L),  # t_i floor for barcode-median
  make_option("--model-floor",  type = "integer", default = 5L)    # t_i floor for model/decomp
)))
for (k in c("raw-prefix", "clean-prefix", "genotypes", "out-prefix"))
  if (is.null(opt[[k]])) stop("missing --", k)
genos    <- strsplit(opt$genotypes, ",", fixed = TRUE)[[1]]
PREFIX   <- opt[["out-prefix"]]
NBOOT    <- opt[["n-boot"]]
DFLOOR   <- opt[["display-floor"]]
MFLOOR   <- opt[["model-floor"]]
set.seed(2025)
DBREAKS <- c(200, 500, 1000, 2000, 5000, 10000, Inf)
DLABS   <- c("200-500", "500-1000", "1000-2000", "2000-5000", "5000-10000", ">10000")
CLASSMAP <- c(single_clean = "singlet", doublet = "doublet", weak_doublet = "weak_doublet")

cat(sprintf("[06_48] sample=%s genotypes=%s n_boot=%d\n", opt$sample, opt$genotypes, NBOOT))

# ----------------------------------------------------------------- load
load_counts <- function(prefix, mode) {
  d <- fread(paste0(prefix, ".counts.tsv.gz"))
  d[, mode := mode]; d
}
counts <- rbind(load_counts(opt[["raw-prefix"]], "raw"),
                load_counts(opt[["clean-prefix"]], "clean"))
attrs <- fread(paste0(opt[["raw-prefix"]], ".barcode_attrs.tsv"))   # canonical (fixed both modes)
attrs[, class := CLASSMAP[call]]
attrs[, depth_bin := cut(raw_depth, breaks = DBREAKS, labels = DLABS, right = FALSE)]
cat(sprintf("[06_48] counts rows: %s ; barcodes in attrs: %s\n",
            format(nrow(counts), big.mark = ","), format(nrow(attrs), big.mark = ",")))

# ---------------------------------------------- per-barcode purity (aggregate blocks)
bc <- counts[, .(T1 = sum(t_top1), M1 = sum(m_top1),
                 Ptot = sum(p_tot), P1 = sum(p_top1), P2 = sum(p_top2),
                 Ttr = sum(t_truth), Mtr = sum(m_truth)), by = .(barcode, mode)]
bc[, p_top1  := fifelse(T1 > 0, M1 / T1, NA_real_)]
bc[, top1or2 := fifelse(Ptot > 0, (P1 + P2) / Ptot, NA_real_)]
bc[, pairfrac:= fifelse((P1 + P2) > 0, P1 / (P1 + P2), NA_real_)]
bc[, p_truth := fifelse(Ttr > 0, Mtr / Ttr, NA_real_)]
bc <- merge(bc, attrs[, .(barcode, call, class, top1, top2, truth, raw_depth, depth_bin)],
            by = "barcode")
bc <- bc[!is.na(class)]                                             # drop any stray
fwrite(bc, paste0(PREFIX, "_barcode_purity.tsv.gz"), sep = "\t")

# --------------------------------------------------------- stratified table
med_f <- function(x, w) { x <- x[w]; if (length(x)) median(x, na.rm = TRUE) else NA_real_ }
iqr_f <- function(x, w) { x <- x[w]; if (length(x) > 1) IQR(x, na.rm = TRUE) else NA_real_ }
strat_one <- bc[T1 > 0, .(
  n_bc        = .N,
  tot_info    = sum(T1),
  rw_purity   = sum(M1) / sum(T1),                          # read-weighted (pooled)
  bc_med_purity = med_f(p_top1, T1 >= DFLOOR),              # barcode-democratic
  bc_iqr_purity = iqr_f(p_top1, T1 >= DFLOOR),
  rw_top1or2  = fifelse(sum(Ptot) > 0, sum(P1 + P2) / sum(Ptot), NA_real_),
  rw_pairfrac = fifelse(sum(P1 + P2) > 0, sum(P1) / sum(P1 + P2), NA_real_)
), by = .(top1, class, depth_bin, mode)]
# pivot raw/clean + delta on read-weighted purity
sw <- dcast(strat_one, top1 + class + depth_bin ~ mode,
            value.var = c("n_bc", "tot_info", "rw_purity", "bc_med_purity", "rw_pairfrac"))
if (!"rw_purity_clean" %in% names(sw)) sw[, rw_purity_clean := NA_real_]
sw[, delta_rw := rw_purity_clean - rw_purity_raw]
sw[, delta_bc_med := bc_med_purity_clean - bc_med_purity_raw]
setorder(sw, top1, class, depth_bin)
fwrite(strat_one, paste0(PREFIX, "_stratified_long.tsv"), sep = "\t")
fwrite(sw,        paste0(PREFIX, "_stratified.tsv"),      sep = "\t")

# --------------------------------------------- weak_doublet-as-singlet diagnostic
# On CLEAN, barcodes with enough informative reads: does weak_doublet look like
# single_clean (top1 high) with an inert top2 (pairfrac ~1, top1or2 ~1)?
diag <- bc[mode == "clean" & T1 >= DFLOOR & class %in% c("singlet", "weak_doublet"),
           .(n = .N,
             med_p_top1  = median(p_top1, na.rm = TRUE),
             med_top1or2 = median(top1or2, na.rm = TRUE),
             med_pairfrac= median(pairfrac, na.rm = TRUE)), by = class]
fwrite(diag, paste0(PREFIX, "_weak_doublet_diag.tsv"), sep = "\t")

# ---------------------------------- within-barcode vs reweighting decomposition
# Point estimate per cell (top1 x class x depth_bin), plus block bootstrap CI.
# Convention (per spec): barcodes with T1_clean==0 -> w_clean=0, p_clean:=p_raw
# (pure reweighting/attrition; keeps the two-term split exact).
cellkey <- c("top1", "class", "depth_bin")
decomp_from <- function(bcm) {   # bcm: barcode-level with M1,T1 per mode (wide)
  w <- bcm
  w[is.na(T1_raw),   `:=`(T1_raw = 0,   M1_raw = 0)]
  w[is.na(T1_clean), `:=`(T1_clean = 0, M1_clean = 0)]
  w[, p_raw   := fifelse(T1_raw   > 0, M1_raw   / T1_raw,   NA_real_)]
  w[, p_clean := fifelse(T1_clean > 0, M1_clean / T1_clean, NA_real_)]
  w[is.na(p_clean), p_clean := p_raw]
  w[is.na(p_raw),   p_raw   := p_clean]
  w <- w[!(is.na(p_raw) & is.na(p_clean))]
  w[, wr := if (sum(T1_raw)   > 0) T1_raw   / sum(T1_raw)   else 0, by = cellkey]
  w[, wc := if (sum(T1_clean) > 0) T1_clean / sum(T1_clean) else 0, by = cellkey]
  w[, `:=`(wbar = (wr + wc) / 2, pbar = (p_raw + p_clean) / 2)]
  w[, .(P_raw    = sum(wr * p_raw),
        P_clean  = sum(wc * p_clean),
        d_within = sum(wbar * (p_clean - p_raw)),
        d_weight = sum((wc - wr) * pbar),
        n_bc     = .N), by = cellkey]
}
# point estimate
bc_pt <- dcast(bc[class != "ambiguous"], barcode + top1 + class + depth_bin ~ mode,
               value.var = c("M1", "T1"))
decomp <- decomp_from(copy(bc_pt))
decomp[, delta := P_clean - P_raw]        # == d_within + d_weight (check)

# block bootstrap: resample 1-Mb blocks (paired raw/clean), recompute per-barcode M/T
cnt_key <- merge(counts[, .(barcode, block, mode, t_top1, m_top1)],
                 attrs[, .(barcode, top1, class, depth_bin)], by = "barcode")
cnt_key <- cnt_key[class != "ambiguous"]
setkey(cnt_key, block)
blocks_all <- unique(cnt_key$block)
boot_stat <- function() {
  bs <- sample(blocks_all, length(blocks_all), replace = TRUE)
  mult <- data.table(block = names(table(bs)), w = as.integer(table(bs)))
  d <- cnt_key[mult, on = "block", nomatch = 0]
  agg <- d[, .(M = sum(m_top1 * w), T = sum(t_top1 * w)), by = .(barcode, top1, class, depth_bin, mode)]
  wide <- dcast(agg, barcode + top1 + class + depth_bin ~ mode, value.var = c("M", "T"))
  setnames(wide, c("M_raw","M_clean","T_raw","T_clean"),
           c("M1_raw","M1_clean","T1_raw","T1_clean"), skip_absent = TRUE)
  for (nm in c("M1_raw","M1_clean","T1_raw","T1_clean"))
    if (!nm %in% names(wide)) wide[, (nm) := 0]
  decomp_from(wide)[, .(top1, class, depth_bin, d_within, d_weight, delta = P_clean - P_raw)]
}
if (NBOOT > 0) {
  cat(sprintf("[06_48] decomposition block bootstrap: %d iters over %s blocks\n",
              NBOOT, format(length(blocks_all), big.mark = ",")))
  bt <- rbindlist(lapply(seq_len(NBOOT), function(i) { b <- boot_stat(); b[, iter := i]; b }))
  ci <- bt[, .(d_within_lo = quantile(d_within, .025, na.rm = TRUE),
               d_within_hi = quantile(d_within, .975, na.rm = TRUE),
               d_weight_lo = quantile(d_weight, .025, na.rm = TRUE),
               d_weight_hi = quantile(d_weight, .975, na.rm = TRUE),
               delta_lo    = quantile(delta,    .025, na.rm = TRUE),
               delta_hi    = quantile(delta,    .975, na.rm = TRUE),
               within_p_gt0= mean(d_within > 0, na.rm = TRUE)), by = cellkey]
  decomp <- merge(decomp, ci, by = cellkey, all.x = TRUE)
}
setorder(decomp, top1, class, depth_bin)
fwrite(decomp, paste0(PREFIX, "_decomposition.tsv"), sep = "\t")

# --------------------------------------------------- aggregated beta-binomial
agg <- merge(counts[, .(barcode, block, mode, t_top1, m_top1, p_tot, p_top1, p_top2)],
             attrs[, .(barcode, top1, class, depth_bin)], by = "barcode")
agg <- agg[class != "ambiguous"]
agg_bb <- agg[, .(match_top1 = sum(m_top1), total_top1 = sum(t_top1),
                  match_p12  = sum(p_top1 + p_top2), total_p12 = sum(p_tot)),
              by = .(top1, class, depth_bin, block, mode)]
fwrite(agg_bb, paste0(PREFIX, "_agg_betabinom.tsv.gz"), sep = "\t")

if (requireNamespace("glmmTMB", quietly = TRUE)) {
  suppressPackageStartupMessages(library(glmmTMB))
  fit_bb <- function(dat, mvar, tvar, tag) {
    dat <- dat[get(tvar) > 0]
    dat[, `:=`(match = get(mvar), nonmatch = get(tvar) - get(mvar))]
    dat[, bam_state := factor(mode, levels = c("raw", "clean"))]
    ok <- length(unique(dat$class)) >= 1 && nrow(dat) > 50
    if (!ok) return(NULL)
    form <- if (length(unique(dat$class)) > 1)
      cbind(match, nonmatch) ~ bam_state * class * depth_bin + top1 + (1 | block)
    else cbind(match, nonmatch) ~ bam_state * depth_bin + top1 + (1 | block)
    f <- tryCatch(glmmTMB(form, family = betabinomial(link = "logit"), data = dat),
                  error = function(e) { message("[06_48] model ", tag, " failed: ", conditionMessage(e)); NULL })
    if (is.null(f)) return(NULL)
    co <- as.data.table(coef(summary(f))$cond, keep.rownames = "term")
    co[, model := tag]; co[]
  }
  models <- rbindlist(list(
    fit_bb(agg_bb[class == "singlet"],                      "match_top1", "total_top1", "singlet_top1"),
    fit_bb(agg_bb[class %in% c("singlet", "weak_doublet")], "match_top1", "total_top1", "singletlike_top1"),
    fit_bb(agg_bb[class %in% c("doublet", "weak_doublet")], "match_top1", "total_top1", "doubletfam_top1"),
    fit_bb(agg_bb[class %in% c("doublet", "weak_doublet")], "match_p12",  "total_p12",  "doubletfam_top1or2")
  ), fill = TRUE)
  if (nrow(models)) fwrite(models, paste0(PREFIX, "_model_coefs.tsv"), sep = "\t")
  cat("[06_48] glmmTMB models fit:", paste(unique(models$model), collapse = ", "), "\n")
} else {
  cat("[06_48] glmmTMB NOT installed — skipped model; agg table written for external fit.\n")
  cat("[06_48]   install: conda install -n base -c conda-forge r-glmmtmb\n")
}

# ------------------------------------------------------- ambiguous margins
amb_path <- paste0(opt[["clean-prefix"]], ".ambiguous.tsv.gz")
if (file.exists(amb_path) && file.info(amb_path)$size > 50) {
  amb <- fread(amb_path)
  tcols <- paste0(genos, "_t"); mcols <- paste0(genos, "_m")
  a <- amb[, lapply(.SD, sum), by = barcode, .SDcols = c(tcols, mcols)]
  conc <- a[, setNames(lapply(genos, function(g)
    fifelse(get(paste0(g, "_t")) > 0, get(paste0(g, "_m")) / get(paste0(g, "_t")), NA_real_)), genos)]
  conc[, barcode := a$barcode]
  cm <- melt(conc, id.vars = "barcode", variable.name = "founder", value.name = "conc")
  cm <- cm[!is.na(conc)]
  marg <- cm[order(-conc), .(best = conc[1], second = conc[2],
                             best_ref = as.character(founder[1])), by = barcode]
  marg[, margin := best - second]
  marg <- merge(marg, attrs[, .(barcode, depth_bin, raw_depth)], by = "barcode")
  fwrite(marg, paste0(PREFIX, "_ambiguous_margin.tsv"), sep = "\t")
  amb_summary <- marg[, .(n = .N, med_best = median(best, na.rm = TRUE),
                          med_margin = median(margin, na.rm = TRUE)), by = depth_bin][order(depth_bin)]
} else {
  amb_summary <- data.table(depth_bin = character(), n = integer())
}

# ------------------------------------------------------------------- plots
plots <- list()
pd <- sw[!is.na(depth_bin)]
if (nrow(pd)) {
  pl <- melt(pd, id.vars = c("top1", "class", "depth_bin"),
             measure.vars = c("rw_purity_raw", "rw_purity_clean"),
             variable.name = "mode", value.name = "purity")
  pl[, mode := ifelse(mode == "rw_purity_raw", "raw", "clean")]
  plots[[1]] <- ggplot(pl, aes(depth_bin, purity, color = mode, group = mode)) +
    geom_line() + geom_point(size = 1) + facet_grid(class ~ top1) +
    labs(title = paste(opt$sample, "- read-weighted top1 purity by depth bin"),
         x = "raw_depth bin", y = "top1 purity") +
    theme_bw(8) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
}
if (nrow(decomp)) {
  dl <- melt(decomp, id.vars = cellkey, measure.vars = c("d_within", "d_weight"),
             variable.name = "component", value.name = "value")
  plots[[length(plots) + 1]] <- ggplot(dl[!is.na(depth_bin)],
      aes(depth_bin, value, fill = component)) +
    geom_col(position = "stack") + facet_grid(class ~ top1) +
    labs(title = "raw->clean purity gain: within-barcode vs reweighting",
         x = "raw_depth bin", y = "delta purity") +
    theme_bw(8) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
}
bcd <- bc[mode == "clean" & T1 >= DFLOOR & class %in% c("singlet", "doublet", "weak_doublet")]
if (nrow(bcd)) {
  plots[[length(plots) + 1]] <- ggplot(bcd, aes(p_top1, color = class)) +
    stat_ecdf() + labs(title = "barcode top1-purity ECDF (clean; weak_doublet vs singlet vs doublet)",
                       x = "top1 purity", y = "ECDF") + theme_bw(8)
}
if (length(plots)) {
  pdf(paste0(PREFIX, "_diagnostics.pdf"), width = 11, height = 4 * length(plots))
  print(Reduce(`/`, plots)); dev.off()
}

# ------------------------------------------------------------ summary md
sink(paste0(PREFIX, "_summary.md"))
cat(sprintf("# 06_48 barcode-resolved purity — %s (WASP-corrected)\n\n", opt$sample))
cat(sprintf("Anchor = AM call (top1=genome_1). Barcodes: %s. Depth = cells_calls n_reads (fixed both modes).\n\n",
            format(nrow(attrs), big.mark = ",")))
cat("## Barcodes per call class (n_reads>=200)\n\n")
cc <- attrs[, .N, by = class][order(-N)]
for (i in seq_len(nrow(cc))) cat(sprintf("- %s: %s\n", cc$class[i], format(cc$N[i], big.mark = ",")))
cat("\n## weak_doublet vs single_clean (clean, t_i>=", DFLOOR, ")\n\n", sep = "")
cat("| class | n | med p_top1 | med top1or2 | med pairfrac |\n|---|---|---|---|---|\n")
for (i in seq_len(nrow(diag)))
  cat(sprintf("| %s | %s | %.3f | %.3f | %.3f |\n", diag$class[i], format(diag$n[i], big.mark = ","),
              diag$med_p_top1[i], diag$med_top1or2[i], diag$med_pairfrac[i]))
cat("\n_Interpretation: weak_doublet pairfrac ~1 (+ top1or2 ~1, p_top1 ~ singlet) => genetically-clean singlets mislabeled; pairfrac ~0.5 => real doublets._\n")
cat("\n## Decomposition of raw->clean top1-purity gain (sum over cells)\n\n")
tot <- decomp[!is.na(delta), .(delta = sum(delta * n_bc) / sum(n_bc),
                               within = sum(d_within * n_bc) / sum(n_bc),
                               weight = sum(d_weight * n_bc) / sum(n_bc))]
cat(sprintf("- n-weighted mean delta=%+.4f  (within=%+.4f, reweight=%+.4f)\n",
            tot$delta, tot$within, tot$weight))
cat("\nFull tables: *_stratified.tsv, *_decomposition.tsv, *_weak_doublet_diag.tsv, *_ambiguous_margin.tsv\n")
sink()

cat(sprintf("[06_48] DONE — outputs at %s_*\n", PREFIX))
