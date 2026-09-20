#!/usr/bin/env Rscript
# 06_44_concordance.R
#
# 2d (count-based) ANALYSIS. Reference-anchored allele concordance + a GRM-style
# correlation, computed directly from per-(group x site) allele counts
# (06_44_extract_allele_counts.py) at KNOWN founder sites — no SNP discovery, no
# MAF / LD / missingness filters. Replaces the SNPRelate GRM (06_42) as the
# primary 2d readout; 06_42 is kept as a confirmatory supplement.
#
# Two estimates per (pseudo-bulk group G, reference founder F):
#   concordance_GF  = reads matching F's allele / panel-allele reads, over
#                     informative sites where F is homozygous. Interpretable %.
#   correlation_GF  = cor(alt_frac_G , dosage_F) over the same sites. The honest
#                     analog of the SNPRelate "Corr" GRM, on all covered sites.
# A correct genotyping signal => argmax_F (of both) == G, and CLEAN >= RAW on the
# self values. Uncertainty on self-value and on the raw->clean delta comes from a
# genomic-block bootstrap (LD-aware).
#
# Inputs:
#   --counts-raw / --counts-clean  group x site count TSV(.gz) from the engine
#   --genotype-panel  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
#                     -s <genotypes>  (no header; CHROM POS REF ALT then 1 GT/geno)
#   --genotypes       comma-separated, same order as bcftools -s
#   --out-prefix      output prefix
#   --min-reads       min panel-allele reads per site to use (default 2)
#   --n-boot          block-bootstrap iterations (default 200)
#   --block-mb        bootstrap block size in Mb (default 1)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts-raw",     type = "character"),
  make_option("--counts-clean",   type = "character"),
  make_option("--genotype-panel", type = "character"),
  make_option("--genotypes",      type = "character"),
  make_option("--out-prefix",     type = "character"),
  make_option("--min-reads",      type = "integer", default = 2L),
  make_option("--n-boot",         type = "integer", default = 200L),
  make_option("--block-mb",       type = "double",  default = 1.0)
)))
for (k in c("counts-raw", "counts-clean", "genotype-panel", "genotypes", "out-prefix"))
  if (is.null(opt[[k]])) stop("missing --", k)

genos    <- strsplit(opt$genotypes, ",", fixed = TRUE)[[1]]
PREFIX   <- opt[["out-prefix"]]
MIN_READS<- opt[["min-reads"]]
N_BOOT   <- opt[["n-boot"]]
BLOCK    <- opt[["block-mb"]] * 1e6
set.seed(2025)
cat(sprintf("[06_44] genotypes: %s\n", paste(genos, collapse = ", ")))

# ----------------------------------------------------------------- panel
cat("[06_44] reading genotype panel\n")
panel <- fread(opt$`genotype-panel`, header = FALSE)
if (ncol(panel) != 4L + length(genos))
  stop(sprintf("panel has %d cols, expected %d (4 + %d genotypes)",
               ncol(panel), 4L + length(genos), length(genos)))
setnames(panel, c("chrom", "pos", "ref", "alt", genos))
panel[, ref := toupper(ref)]
panel[, alt1 := toupper(tstrsplit(alt, ",", fixed = TRUE)[[1]])]

norm_gt <- function(x) {
  x <- gsub("\\|", "/", x)
  fcase(x %in% "0/0", "ref", x %in% "1/1", "alt",
        x %in% c("0/1", "1/0"), "het", default = NA_character_)
}
for (g in genos) {
  panel[, (paste0(g, "_call")) := norm_gt(get(g))]
  panel[, (paste0(g, "_dose")) := fcase(get(paste0(g, "_call")) == "ref", 0,
                                        get(paste0(g, "_call")) == "alt", 1,
                                        default = NA_real_)]
}
panel[, n_ref_hom := Reduce(`+`, lapply(genos, function(g) get(paste0(g, "_call")) == "ref" & !is.na(get(paste0(g, "_call")))))]
panel[, n_alt_hom := Reduce(`+`, lapply(genos, function(g) get(paste0(g, "_call")) == "alt" & !is.na(get(paste0(g, "_call")))))]
panel[, informative := n_ref_hom > 0L & n_alt_hom > 0L]
panel_inf <- panel[informative == TRUE,
                   c("chrom", "pos", "ref", "alt1",
                     paste0(genos, "_call"), paste0(genos, "_dose")), with = FALSE]
setkey(panel_inf, chrom, pos)
cat(sprintf("[06_44]   panel %s sites, %s informative\n",
            format(nrow(panel), big.mark = ","),
            format(nrow(panel_inf), big.mark = ",")))

# ------------------------------------------------- per-mode joined site table
load_join <- function(counts_file) {
  cnt <- fread(counts_file)                       # chrom pos ref alt group n_ref n_alt n_other
  cnt[, c("ref", "alt") := NULL]                  # panel ref/alt1 are authoritative
  m <- merge(cnt, panel_inf, by = c("chrom", "pos"))
  m[, total := n_ref + n_alt]
  m <- m[total >= MIN_READS]
  m[, alt_frac := n_alt / total]
  m[, block := paste0(chrom, "_", pos %/% BLOCK)]
  m[]
}

# group x founder concordance + correlation (point estimates)
matrices <- function(m) {
  conc <- corr <- matrix(NA_real_, length(genos), length(genos),
                         dimnames = list(genos, genos))
  for (G in genos) {
    sub <- m[group == G]
    if (!nrow(sub)) next
    for (F in genos) {
      fc <- sub[[paste0(F, "_call")]]
      fd <- sub[[paste0(F, "_dose")]]
      use <- fc %in% c("ref", "alt")
      if (sum(use) == 0) next
      mtch <- fifelse(fc == "ref", sub$n_ref, sub$n_alt)
      conc[G, F] <- sum(mtch[use]) / sum(sub$total[use])
      if (sum(use) > 2 && sd(sub$alt_frac[use]) > 0 && sd(fd[use]) > 0)
        corr[G, F] <- cor(sub$alt_frac[use], fd[use])
    }
  }
  list(conc = conc, corr = corr)
}

# per-(group,block) self sums for the bootstrap (self founder == group name)
self_blocks <- function(m) {
  out <- vector("list", length(genos))
  for (i in seq_along(genos)) {
    G <- genos[i]
    sub <- m[group == G]
    if (!nrow(sub)) next
    gc <- sub[[paste0(G, "_call")]]
    use <- gc %in% c("ref", "alt")
    sub <- sub[use]
    sub[, mtch := fifelse(gc[use] == "ref", n_ref, n_alt)]
    out[[i]] <- sub[, .(group = G, m = sum(mtch), t = sum(total)), by = block]
  }
  rbindlist(out)
}

cat("[06_44] loading + joining RAW\n");   m_raw   <- load_join(opt$`counts-raw`)
cat("[06_44] loading + joining CLEAN\n"); m_clean <- load_join(opt$`counts-clean`)

mat_raw <- matrices(m_raw); mat_clean <- matrices(m_clean)
sb_raw  <- self_blocks(m_raw); sb_clean <- self_blocks(m_clean)

# ----------------------------------------------------- coverage + matrices out
cov_tab <- rbindlist(list(
  m_raw[,   .(mode = "raw",   n_sites = .N, n_reads = sum(total)), by = group],
  m_clean[, .(mode = "clean", n_sites = .N, n_reads = sum(total)), by = group]))
fwrite(cov_tab, paste0(PREFIX, "_coverage.tsv"), sep = "\t")

mat_long <- function(mat, what, mode)
  data.table(mode = mode, metric = what,
             group = rep(rownames(mat), times = ncol(mat)),
             founder = rep(colnames(mat), each = nrow(mat)),
             value = as.numeric(mat))
mats <- rbindlist(list(
  mat_long(mat_raw$conc,   "concordance", "raw"),
  mat_long(mat_clean$conc, "concordance", "clean"),
  mat_long(mat_raw$corr,   "correlation", "raw"),
  mat_long(mat_clean$corr, "correlation", "clean")))
fwrite(mats, paste0(PREFIX, "_matrices_long.tsv"), sep = "\t")

# ------------------------------------------------------ best-ref + margin
best_for <- function(mat, mode) {
  rbindlist(lapply(rownames(mat), function(G) {
    v <- mat[G, ]; v_self <- v[G]; off <- v[setdiff(names(v), G)]
    j <- which.max(v)
    data.table(mode = mode, group = G, self = v_self,
               best_off = if (length(off)) max(off, na.rm = TRUE) else NA_real_,
               margin = v_self - (if (length(off)) max(off, na.rm = TRUE) else NA_real_),
               best_ref = names(v)[j], correct = names(v)[j] == G)
  }))
}
best <- rbindlist(list(best_for(mat_raw$conc, "raw"), best_for(mat_clean$conc, "clean")))
fwrite(best, paste0(PREFIX, "_best_ref.tsv"), sep = "\t")

# --------------------------------------------------- block bootstrap on self + delta
blocks_all <- union(unique(sb_raw$block), unique(sb_clean$block))
boot_self <- function(sb, samp) {       # samp: table(block)->mult ; returns named self-conc
  d <- sb[J(block = names(samp)), on = "block", nomatch = 0]
  d[, mult := as.integer(samp[block])]
  d[, .(self = sum(m * mult) / sum(t * mult)), by = group]
}
setkey(sb_raw, block); setkey(sb_clean, block)
boot <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  bs <- sample(blocks_all, length(blocks_all), replace = TRUE)
  samp <- table(bs)
  r <- boot_self(sb_raw, samp);   setnames(r, "self", "raw")
  c <- boot_self(sb_clean, samp); setnames(c, "self", "clean")
  z <- merge(r, c, by = "group", all = TRUE)
  z[, `:=`(delta = clean - raw, iter = b)]; z
}), fill = TRUE)

ci <- boot[, .(
  self_raw   = mean(raw,   na.rm = TRUE),
  raw_lo     = quantile(raw,   0.025, na.rm = TRUE),
  raw_hi     = quantile(raw,   0.975, na.rm = TRUE),
  self_clean = mean(clean, na.rm = TRUE),
  clean_lo   = quantile(clean, 0.025, na.rm = TRUE),
  clean_hi   = quantile(clean, 0.975, na.rm = TRUE),
  delta      = mean(delta, na.rm = TRUE),
  delta_lo   = quantile(delta, 0.025, na.rm = TRUE),
  delta_hi   = quantile(delta, 0.975, na.rm = TRUE),
  delta_p_gt0= mean(delta > 0, na.rm = TRUE)
), by = group]
# overlay point estimates from the full data
pt <- merge(best[mode == "raw",   .(group, self_raw_pt = self)],
            best[mode == "clean", .(group, self_clean_pt = self)], by = "group")
ci <- merge(ci, pt, by = "group", all.x = TRUE)
fwrite(ci, paste0(PREFIX, "_self_concordance_ci.tsv"), sep = "\t")

# ---------------------------------------------------------- ref-bias sidecar
read_bias <- function(counts_file, mode) {
  bp <- paste0(sub("\\.gz$", "", counts_file), ".bias.tsv")
  if (!file.exists(bp)) return(NULL)
  b <- fread(bp); b[, mode := mode]
  b[, `:=`(mean_nm = sum_nm / pmax(n_obs, 1), mean_mapq = sum_mapq / pmax(n_obs, 1))]
  b
}
bias <- rbindlist(list(read_bias(opt$`counts-raw`, "raw"),
                       read_bias(opt$`counts-clean`, "clean")), fill = TRUE)
if (nrow(bias)) fwrite(bias, paste0(PREFIX, "_refbias.tsv"), sep = "\t")

# ------------------------------------------------------------ markdown summary
sink(paste0(PREFIX, "_summary.md"))
cat("# 06_44 count-based 2d — reference-anchored allele concordance\n\n")
cat(sprintf("Genotypes: %s | informative panel sites: %s | min reads/site: %d | bootstrap: %d x %.0f-Mb blocks\n\n",
            paste(genos, collapse = ", "), format(nrow(panel_inf), big.mark = ","),
            MIN_READS, N_BOOT, BLOCK / 1e6))
cat("## Self-concordance (group vs its own founder), RAW -> CLEAN\n\n")
cat("| group | self_raw | self_clean | delta [95% CI] | P(delta>0) | correct(raw/clean) |\n")
cat("|-------|----------|------------|----------------|------------|--------------------|\n")
for (g in genos) {
  r <- ci[group == g]; if (!nrow(r)) next
  cr <- best[mode == "raw"   & group == g]$correct
  cc <- best[mode == "clean" & group == g]$correct
  dpt <- r$self_clean_pt - r$self_raw_pt          # point-estimate delta; CI from bootstrap
  cat(sprintf("| %s | %.3f | %.3f | %+.3f [%+.3f, %+.3f] | %.2f | %s/%s |\n",
              g, r$self_raw_pt, r$self_clean_pt, dpt, r$delta_lo, r$delta_hi,
              r$delta_p_gt0, cr, cc))
}
cat("\n## Coverage (sites x allele-obs per group)\n\n")
cat("| group | sites_raw | reads_raw | sites_clean | reads_clean |\n")
cat("|-------|-----------|-----------|-------------|-------------|\n")
for (g in genos) {
  cr <- cov_tab[mode == "raw" & group == g]; cc <- cov_tab[mode == "clean" & group == g]
  cat(sprintf("| %s | %s | %s | %s | %s |\n", g,
              format(if (nrow(cr)) cr$n_sites else 0, big.mark = ","),
              format(if (nrow(cr)) cr$n_reads else 0, big.mark = ","),
              format(if (nrow(cc)) cc$n_sites else 0, big.mark = ","),
              format(if (nrow(cc)) cc$n_reads else 0, big.mark = ",")))
}
cat(sprintf("\n[06_44] argmax-correct: raw %d/%d, clean %d/%d\n",
            sum(best[mode == "raw"]$correct), length(genos),
            sum(best[mode == "clean"]$correct), length(genos)))
sink()

# ------------------------------------------------------------------- plots
to_df <- function(mat, mode, metric) {
  d <- as.data.table(as.table(mat)); setnames(d, c("group", "founder", "value"))
  d[, `:=`(mode = mode, metric = metric)]; d
}
hm <- rbindlist(list(to_df(mat_raw$conc, "raw", "concordance"),
                     to_df(mat_clean$conc, "clean", "concordance")))
p1 <- ggplot(hm, aes(founder, group, fill = value)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", value)), size = 2.6) +
  facet_wrap(~ mode) + scale_fill_viridis_c(option = "mako", limits = c(0, 1)) +
  labs(title = "Allele concordance (group x founder)", x = "reference founder", y = "pseudo-bulk") +
  theme_bw(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))

cip <- melt(ci, id.vars = "group", measure.vars = c("self_raw_pt", "self_clean_pt"),
            variable.name = "mode", value.name = "self")
cip[, mode := ifelse(mode == "self_raw_pt", "raw", "clean")]
cib <- merge(cip, ci[, .(group, raw_lo, raw_hi, clean_lo, clean_hi)], by = "group")
cib[, `:=`(lo = ifelse(mode == "raw", raw_lo, clean_lo),
           hi = ifelse(mode == "raw", raw_hi, clean_hi))]
p2 <- ggplot(cib, aes(group, self, color = mode)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.4)) +
  labs(title = "Self-concordance with block-bootstrap 95% CI", y = "self-concordance", x = NULL) +
  theme_bw(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))

plots <- list(p1, p2)
if (nrow(bias)) {
  p3 <- ggplot(bias, aes(group, mean_nm, fill = interaction(obs_class, mode))) +
    geom_col(position = position_dodge()) +
    labs(title = "Reference-mapping bias: mean NM by observed allele class",
         subtitle = "alt > ref NM => ALT-carrying reads penalised by B73-only mapping",
         y = "mean NM / obs", x = NULL) +
    theme_bw(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  plots <- c(plots, list(p3))
}
pdf(paste0(PREFIX, "_diagnostics.pdf"), width = 10, height = 4 * length(plots))
print(Reduce(`/`, plots))
dev.off()

cat(sprintf("[06_44] DONE — argmax-correct raw %d/%d clean %d/%d ; outputs at %s_*\n",
            sum(best[mode == "raw"]$correct), length(genos),
            sum(best[mode == "clean"]$correct), length(genos), PREFIX))
