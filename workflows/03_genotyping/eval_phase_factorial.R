#!/usr/bin/env Rscript
# =============================================================================
# eval_phase_factorial.R — shared evaluator for the Phase 2/3/4 genotyping factorials
#
#   1. Phase tag dispatch (phase=phase2|phase3|phase4) selects the config
#      list, the metric function, and the input/output paths.
#   2. Two metric functions:
#        make_truth_metrics()              — Phase 2/4 (synthetic + full Root1
#                                            with operational truth)
#        make_root1_operational_metrics()  — Phase 3 (sub1k, no truth table)
#   3. Phase 4 includes a sub1k-vs-full consistency panel.
#
# Run from: ${PROJECT_ROOT}/5_AmbientDetection/ (the directory holding the synthetic/,
#           synthetic_disc/ and Root1_rep1/ genotyping_runs trees)
# Usage:    Rscript <repo>/workflows/03_genotyping/eval_phase_factorial.R <phase> [<date>]
#
#   phase ∈ {phase2, phase3, phase4}
#   date is optional and defaults to today's YYYY-MM-DD; use to re-eval an
#   older factorial run that was tagged with a different date.
#
# Outputs (the *_summary_metrics.tsv is the manuscript input):
#   Phase 2: synthetic/eval/phase2_<date>/phase2_summary_metrics.tsv    -> Fig. S5
#   Phase 3: Root1_rep1/eval/phase3_<date>/phase3_summary_metrics.tsv  -> Fig. S6
#   Phase 4: Root1_rep1/eval/phase4_<date>/phase4_summary_metrics.tsv  -> Fig. 4G
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(jsonlite)
})

# -------------------------
# CLI args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript eval_phase_factorial.R <phase2|phase3|phase4> [<date>]")
}
PHASE <- args[1]
if (!PHASE %in% c("phase2", "phase3", "phase4")) {
  stop("Unknown phase: ", PHASE, " (expected one of: phase2 phase3 phase4)")
}
EVAL_DATE <- if (length(args) >= 2) args[2] else format(Sys.Date(), "%Y-%m-%d")
FACTORIAL_TAG <- paste0("factorial_", PHASE, "_", EVAL_DATE)
MIN_READS <- 500  # exclude barcodes with fewer than this many reads

cat("\n=== eval_phase_factorial.R ===\n")
cat("  phase         =", PHASE, "\n")
cat("  eval date     =", EVAL_DATE, "\n")
cat("  factorial tag =", FACTORIAL_TAG, "\n")
cat("  min reads     =", MIN_READS, "\n")

# Shared call colours / labels (analysis/_helpers/fig4_helpers.R) are optional:
# nothing below depends on them, this script builds its own palette.
REPO_ROOT <- Sys.getenv("REPO_ROOT", unset = ".")
HELPERS <- file.path(REPO_ROOT, "analysis", "_helpers", "fig4_helpers.R")
if (file.exists(HELPERS)) source(HELPERS)

# -------------------------
# Phase config registries
# -------------------------
PHASE2_CONFIGS <- c(
  "C0"                   = "C0",
  "C1c_nofriend"         = "C1c_nofriend",
  "C2a_xmap"             = "C2a_xmap",
  "C2c_xmap_eta0"        = "C2c_xmap_eta0",
  "C3b_mq50"             = "C3b_mq50",
  "C4a_bic3"             = "C4a_bic3",
  "C4c_eta0"             = "C4c_eta0",
  "C4c_eta5"             = "C4c_eta5",
  "C4d_wamb005_wdon"     = "C4d_wamb005_wdon",
  "C4d_wamb05_wdon"      = "C4d_wamb05_wdon"
)

PHASE3_CONFIGS <- c(
  # Reference configs (Phase 1 single-knob, for in-figure comparison)
  "C0"                       = "C0",
  "C1c_nofriend"             = "C1c_nofriend",
  "C1g_naked"                = "C1g_naked",
  "C3b_mq50"                 = "C3b_mq50",
  "C4a_bic3"                 = "C4a_bic3",
  "C4d_wamb005_wdon"         = "C4d_wamb005_wdon",
  # Phase 3 stacked configs
  "S01_nofr_mq50"            = "S01_nofr_mq50",
  "S02_nofr_mq50_xaUL"       = "S02_nofr_mq50_xaUL",
  "S03_nofr_mq50_xmap"       = "S03_nofr_mq50_xmap",
  "S04_nofr_mq50_xmap_xaUL"  = "S04_nofr_mq50_xmap_xaUL",
  "S05_nofr_mq50_wamb"       = "S05_nofr_mq50_wamb",
  "S06_nofr_mq50_wamb_xaUL"  = "S06_nofr_mq50_wamb_xaUL",
  "S07_nofr_mq50_xmap_wamb"  = "S07_nofr_mq50_xmap_wamb",
  "S08_full_xa0"             = "S08_full_xa0",
  "S09_full_xaUL"            = "S09_full_xaUL",
  "S10_full_bic3"            = "S10_full_bic3",
  "S11_full_friendon"        = "S11_full_friendon",
  "S12_full_wamb01"          = "S12_full_wamb01",
  "S13_full_eta2"            = "S13_full_eta2"
)

# Phase 4 = C0 + the Phase 3 winner. The winner name is auto-detected from
# the dirs present under Root1_rep1/genotyping_runs/<factorial_tag>/.
PHASE4_BASE_CONFIGS <- c("C0" = "C0")

# -------------------------
# Phase-specific dispatch
# -------------------------
DATASETS_SYN <- c(
  "alpha_000",
  paste0("alpha_", sprintf("%03d", c(2, 5, 10, 20, 30, 40, 50)), "_Il14H"),
  paste0("alpha_", sprintf("%03d", c(2, 5, 10, 20, 30, 40, 50)), "_Ki11")
)

if (PHASE == "phase2") {
  TRACKS <- list("Track B" = "synthetic", "Track B-disc" = "synthetic_disc")
  CONFIGS <- PHASE2_CONFIGS
  EVALDIR <- file.path("synthetic", "eval", paste0("phase2_", EVAL_DATE))
  USE_TRUTH <- TRUE
  USE_BINS  <- FALSE

} else if (PHASE == "phase3") {
  TRACKS <- list("sub1k_A" = "Root1_rep1/sub1k_A",
                 "sub1k_B" = "Root1_rep1/sub1k_B")
  CONFIGS <- PHASE3_CONFIGS
  EVALDIR <- file.path("Root1_rep1", "eval", paste0("phase3_", EVAL_DATE))
  USE_TRUTH <- FALSE
  USE_BINS  <- TRUE

} else { # phase4
  TRACKS <- list("Root1_rep1_full" = "Root1_rep1")
  CONFIGS <- PHASE4_BASE_CONFIGS  # winner appended below after auto-detection
  EVALDIR <- file.path("Root1_rep1", "eval", paste0("phase4_", EVAL_DATE))
  USE_TRUTH <- FALSE
  USE_BINS  <- FALSE  # full Root1: no panel-style barcodes_with_bin
}

dir.create(EVALDIR, recursive = TRUE, showWarnings = FALSE)
cat("  output dir    =", EVALDIR, "\n\n")

# -------------------------
# Phase 4 winner auto-detect
# -------------------------
if (PHASE == "phase4") {
  base <- file.path("Root1_rep1", "genotyping_runs", FACTORIAL_TAG)
  if (dir.exists(base)) {
    found <- list.dirs(base, recursive = FALSE, full.names = FALSE)
    extras <- setdiff(found, names(CONFIGS))
    extras <- extras[extras != ""]
    if (length(extras) > 0) {
      for (e in extras) CONFIGS[[e]] <- e
      cat("  auto-detected phase4 winner config(s):", paste(extras, collapse = ", "), "\n")
    }
  }
}

# -------------------------
# Loaders
# -------------------------
load_truth_for_track <- function(base_dir) {
  barcoded <- file.path(base_dir, "barcoded")
  truth_list <- lapply(DATASETS_SYN, function(ds) {
    f <- file.path(barcoded, ds, "truth_table.tsv")
    if (!file.exists(f)) { warning("Missing truth: ", f); return(NULL) }
    dt <- fread(f); dt[, dataset := ds]; dt
  })
  truth <- rbindlist(truth_list[!sapply(truth_list, is.null)], fill = TRUE)
  setnames(truth,
           c("alpha", "genome_1", "genome_2", "rho"),
           c("true_alpha", "true_g1", "true_g2", "true_rho"),
           skip_absent = TRUE)
  truth
}

load_calls_synthetic <- function(base_dir) {
  rows <- list()
  n_missing <- 0
  for (ds in DATASETS_SYN) {
    for (m in names(CONFIGS)) {
      cfg_dir <- CONFIGS[[m]]
      f <- file.path(base_dir, ds, "genotyping_runs", FACTORIAL_TAG, cfg_dir,
                     paste0(ds, "_cells_calls.tsv.gz"))
      if (!file.exists(f)) {
        cat("  MISSING:", base_dir, "/", ds, "/", cfg_dir, "\n")
        n_missing <- n_missing + 1; next
      }
      dt <- fread(f, sep = "\t")
      dt[, `:=`(dataset = ds, method = m)]
      rows[[length(rows) + 1]] <- dt
    }
  }
  cat("  loaded", length(rows), "files (", n_missing, "missing)\n")
  if (length(rows) == 0) return(data.table())
  rbindlist(rows, fill = TRUE)
}

load_calls_root1_panel <- function(panel_dir) {
  # panel_dir is e.g. "Root1_rep1/sub1k_A" or "Root1_rep1"
  rows <- list()
  n_missing <- 0
  panel_name <- basename(panel_dir)
  for (m in names(CONFIGS)) {
    cfg_dir <- CONFIGS[[m]]
    f <- file.path(panel_dir, "genotyping_runs", FACTORIAL_TAG, cfg_dir,
                   paste0(panel_name, "_cells_calls.tsv.gz"))
    if (!file.exists(f)) {
      cat("  MISSING:", f, "\n")
      n_missing <- n_missing + 1; next
    }
    dt <- fread(f, sep = "\t")
    dt[, `:=`(dataset = panel_name, method = m)]
    rows[[length(rows) + 1]] <- dt
  }
  cat("  loaded", length(rows), "files for", panel_name, "(", n_missing, "missing)\n")
  if (length(rows) == 0) return(data.table())
  rbindlist(rows, fill = TRUE)
}

load_bins_for_panel <- function(panel_dir) {
  f <- file.path(panel_dir, "barcodes_with_bin.tsv")
  if (!file.exists(f)) {
    warning("Missing barcodes_with_bin.tsv for ", panel_dir)
    return(NULL)
  }
  dt <- fread(f, sep = "\t")
  # Expected columns: barcode, n_reads, bin
  setnames(dt, "bin", "depth_bin", skip_absent = TRUE)
  dt
}

# -------------------------
# Per-track load + merge
# -------------------------
per_track <- list()
for (track_name in names(TRACKS)) {
  base <- TRACKS[[track_name]]
  cat("Loading", track_name, "(", base, ")...\n")

  if (USE_TRUTH) {
    truth <- load_truth_for_track(base)
    calls <- load_calls_synthetic(base)
    if (nrow(calls) == 0) next

    # Rename cells_calls columns to avoid clashes with truth
    rename_cols <- intersect(c("alpha", "genome_1", "genome_2", "rho"), names(calls))
    rename_map <- c(alpha = "est_alpha", genome_1 = "called_g1",
                    genome_2 = "called_g2", rho = "est_rho")
    setnames(calls, rename_cols, rename_map[rename_cols])

    truth_keep <- c("barcode", "dataset", "type", "true_genome", "contaminant",
                    "true_alpha", "true_g1", "true_g2", "true_rho",
                    "n_cell_reads", "n_contam_reads", "n_total_reads",
                    "target_depth", "actual_depth")
    truth_keep <- intersect(truth_keep, names(truth))
    merged <- merge(calls, truth[, ..truth_keep],
                    by = c("barcode", "dataset"), all.x = TRUE)
    merged[, track := track_name]

    merged[, `:=`(
      is_singlet_call_strict = call %in% c("single_clean", "dirty_singlet"),
      is_singlet_call_relaxed = call %in% c("single_clean", "dirty_singlet") |
        (call == "weak_doublet" & called_g1 == "B73"),
      called_b73 = (!is.na(called_g1) & called_g1 == "B73"),
      is_doublet_call = (call == "doublet"),
      is_empty_call = (call == "empty")
    )]
    merged[, `:=`(
      is_correct_singlet_strict  = is_singlet_call_strict  & called_b73,
      is_correct_singlet_relaxed = is_singlet_call_relaxed & called_b73
    )]

    # Filter to barcodes with >= MIN_READS total reads
    if ("n_total_reads" %in% names(merged)) {
      n_before <- nrow(merged)
      merged <- merged[n_total_reads >= MIN_READS]
      cat("  filtered to n_total_reads >=", MIN_READS, ":",
          n_before, "->", nrow(merged), "rows\n")
    }

  } else {
    # Phase 3 / Phase 4: Root1 operational metrics
    calls <- load_calls_root1_panel(base)
    if (nrow(calls) == 0) next

    rename_cols <- intersect(c("alpha", "genome_1", "genome_2", "rho"), names(calls))
    rename_map <- c(alpha = "est_alpha", genome_1 = "called_g1",
                    genome_2 = "called_g2", rho = "est_rho")
    setnames(calls, rename_cols, rename_map[rename_cols])

    merged <- copy(calls)
    merged[, track := track_name]
    merged[, called_b73 := (!is.na(called_g1) & called_g1 == "B73")]
    merged[, is_singlet_call_strict := call %in% c("single_clean", "dirty_singlet")]
    merged[, is_singlet_call_relaxed := is_singlet_call_strict |
                                         (call == "weak_doublet" & called_b73)]
    merged[, is_correct_singlet_strict  := is_singlet_call_strict  & called_b73]
    merged[, is_correct_singlet_relaxed := is_singlet_call_relaxed & called_b73]
    merged[, is_doublet_call := (call == "doublet")]
    merged[, is_ambiguous_call := (call == "ambiguous")]
    merged[, is_wrong_genome := (!is_singlet_call_relaxed) & is_singlet_call_strict |
                                 ((!called_b73) & call %in% c("single_clean","dirty_singlet","weak_doublet"))]
    # Simpler wrong-genome definition: any singlet/weak_doublet call where g1 != B73
    merged[, is_wrong_genome := call %in% c("single_clean","dirty_singlet","weak_doublet") &
                                 (!called_b73)]

    # Attach depth bins for sub1k panels and filter on original barcode depth.
    # cells_calls barcode = "ACGT-sub1k_A"; barcodes_with_bin = "ACGT" (no suffix).
    # Strip the suffix to create a join key.
    if (USE_BINS) {
      bins <- load_bins_for_panel(base)
      if (!is.null(bins)) {
        merged[, barcode_short := sub("-[^-]+$", "", barcode)]
        merged <- merge(merged, bins[, .(barcode, total_nuclear_reads, depth_bin)],
                        by.x = "barcode_short", by.y = "barcode", all.x = TRUE)
        merged[, barcode_short := NULL]
        # Filter on original nuclear read depth (not post-filter n_reads)
        n_before <- nrow(merged)
        merged <- merged[!is.na(total_nuclear_reads) & total_nuclear_reads >= MIN_READS]
        cat("  filtered to total_nuclear_reads >=", MIN_READS, ":",
            n_before, "->", nrow(merged), "rows\n")
      }
    } else if ("n_reads" %in% names(merged)) {
      # Phase 4 full Root1 (no bins file): fall back to n_reads from cells_calls
      n_before <- nrow(merged)
      merged <- merged[n_reads >= MIN_READS]
      cat("  filtered to n_reads >=", MIN_READS, ":",
          n_before, "->", nrow(merged), "rows\n")
    }
  }

  per_track[[track_name]] <- merged
}

if (length(per_track) == 0) {
  stop("No data loaded — check that the genotyping runs exist under ", FACTORIAL_TAG)
}

merged_all <- rbindlist(per_track, fill = TRUE)
cat("\nMerged all tracks:", format(nrow(merged_all), big.mark = ","), "rows\n\n")

# =============================================================================
# Metric functions
# =============================================================================

# Phase 2/4: truth-based metrics
make_truth_metrics <- function(d) {
  s <- d[type == "singlet"]
  db <- d[type == "doublet"]
  em <- d[type == "empty"]
  list(
    n_singlets = nrow(s),
    n_singlet_calls_strict   = sum(s$is_singlet_call_strict),
    n_correct_singlet_strict = sum(s$is_correct_singlet_strict),
    n_singlet_calls_relaxed   = sum(s$is_singlet_call_relaxed),
    n_correct_singlet_relaxed = sum(s$is_correct_singlet_relaxed),
    n_doublets = nrow(db),
    n_doublet_correct = sum(db$is_doublet_call),
    n_empties = nrow(em),
    n_empty_correct = sum(em$is_empty_call)
  )
}

# Phase 3/4: Root1 operational metrics (no truth table; all real cells = true B73)
# Denominator = total cells in the input (or total cells in a depth bin)
make_root1_operational_metrics <- function(d) {
  total <- nrow(d)
  list(
    n_total                 = total,
    n_correct_strict        = sum(d$is_singlet_call_strict & d$called_b73),
    n_correct_relaxed       = sum(d$is_singlet_call_relaxed & d$called_b73),
    n_wrong_genome          = sum(d$is_wrong_genome),
    n_called_doublet        = sum(d$is_doublet_call),
    n_called_ambiguous      = sum(d$is_ambiguous_call),
    n_called_empty          = sum(d$call == "empty")
  )
}

# Aggregator
if (USE_TRUTH) {
  summary_metrics <- merged_all[, {
    m <- make_truth_metrics(.SD)
    sens_strict   <- if (m$n_singlets) m$n_correct_singlet_strict   / m$n_singlets   else NA_real_
    sens_relaxed  <- if (m$n_singlets) m$n_correct_singlet_relaxed  / m$n_singlets   else NA_real_
    prec_strict   <- if (m$n_singlet_calls_strict)  m$n_correct_singlet_strict  / m$n_singlet_calls_strict  else NA_real_
    prec_relaxed  <- if (m$n_singlet_calls_relaxed) m$n_correct_singlet_relaxed / m$n_singlet_calls_relaxed else NA_real_
    f1_strict     <- if (!is.na(sens_strict)  && !is.na(prec_strict)  && (sens_strict  + prec_strict)  > 0) 2*sens_strict*prec_strict/(sens_strict+prec_strict)   else NA_real_
    f1_relaxed    <- if (!is.na(sens_relaxed) && !is.na(prec_relaxed) && (sens_relaxed + prec_relaxed) > 0) 2*sens_relaxed*prec_relaxed/(sens_relaxed+prec_relaxed) else NA_real_
    dbl_rate      <- if (m$n_doublets) m$n_doublet_correct / m$n_doublets else NA_real_
    empty_rate    <- if (m$n_empties)  m$n_empty_correct   / m$n_empties  else NA_real_
    .(n_singlets   = m$n_singlets,
      n_doublets   = m$n_doublets,
      n_empties    = m$n_empties,
      sens_strict  = sens_strict,
      sens_relaxed = sens_relaxed,
      prec_strict  = prec_strict,
      prec_relaxed = prec_relaxed,
      f1_strict    = f1_strict,
      f1_relaxed   = f1_relaxed,
      doublet_detection_rate = dbl_rate,
      empty_detection_rate   = empty_rate)
  }, by = .(track, method, dataset)]

  summary_metrics[, contam_series := fcase(
    dataset == "alpha_000",        "none",
    grepl("Il14H$", dataset),       "Il14H",
    grepl("Ki11$",  dataset),       "Ki11",
    default = NA_character_
  )]
  summary_metrics[, true_alpha := as.numeric(sub("^alpha_(\\d{3}).*", "\\1", dataset)) / 100]

} else {
  # Operational metrics — Phase 3/4
  by_cols <- if (USE_BINS && "depth_bin" %in% names(merged_all)) {
    c("track", "method", "depth_bin")
  } else {
    c("track", "method")
  }
  summary_metrics <- merged_all[, {
    m <- make_root1_operational_metrics(.SD)
    sens_strict  <- if (m$n_total) m$n_correct_strict  / m$n_total else NA_real_
    sens_relaxed <- if (m$n_total) m$n_correct_relaxed / m$n_total else NA_real_
    n_singlet_calls <- sum(.SD$is_singlet_call_relaxed)
    prec_relaxed <- if (n_singlet_calls) m$n_correct_relaxed / n_singlet_calls else NA_real_
    f1_relaxed   <- if (!is.na(sens_relaxed) && !is.na(prec_relaxed) && (sens_relaxed+prec_relaxed)>0) 2*sens_relaxed*prec_relaxed/(sens_relaxed+prec_relaxed) else NA_real_
    .(n_total           = m$n_total,
      n_correct_strict  = m$n_correct_strict,
      n_correct_relaxed = m$n_correct_relaxed,
      n_wrong_genome    = m$n_wrong_genome,
      n_called_doublet  = m$n_called_doublet,
      n_called_ambig    = m$n_called_ambiguous,
      sens_strict       = sens_strict,
      sens_relaxed      = sens_relaxed,
      prec_relaxed      = prec_relaxed,
      f1_relaxed        = f1_relaxed,
      wrong_genome_rate = m$n_wrong_genome / m$n_total,
      doublet_rate      = m$n_called_doublet / m$n_total,
      ambig_rate        = m$n_called_ambiguous / m$n_total)
  }, by = by_cols]
}

OUT_TSV <- file.path(EVALDIR, paste0(PHASE, "_summary_metrics.tsv"))
fwrite(summary_metrics, OUT_TSV, sep = "\t")
cat("Wrote", OUT_TSV, "\n")

# -------------------------
# Console summary
# -------------------------
cat("\n=== Summary metrics ===\n")
if (USE_TRUTH) {
  print_cols <- c("track","method","dataset","contam_series","true_alpha",
                  "sens_strict","sens_relaxed","prec_strict","f1_strict","f1_relaxed")
} else {
  print_cols <- intersect(c("track","method","depth_bin","n_total",
                             "sens_relaxed","prec_relaxed","f1_relaxed",
                             "wrong_genome_rate","doublet_rate","ambig_rate"),
                          names(summary_metrics))
}
print(as.data.frame(summary_metrics[, ..print_cols]),
      row.names = FALSE, digits = 3)

# =============================================================================
# Figures
# =============================================================================
cat("\nGenerating figures...\n")
methods_seen <- unique(summary_metrics$method)
n_methods <- length(methods_seen)

# Generate a palette big enough for the configs
method_palette <- if (n_methods <= 10) {
  setNames(scales::hue_pal()(n_methods), methods_seen)
} else {
  setNames(scales::hue_pal(l = 60)(n_methods), methods_seen)
}

if (USE_TRUTH) {
  # Phase 2/4 — synthetic-style panels
  sm <- copy(summary_metrics)
  sm_nz <- sm[contam_series != "none"]
  sm_zero <- sm[contam_series == "none"]
  sm_il14h <- copy(sm_zero)[, contam_series := "Il14H"]
  sm_ki11  <- copy(sm_zero)[, contam_series := "Ki11"]
  sm_plot  <- rbindlist(list(sm_nz, sm_il14h, sm_ki11), fill = TRUE)
  sm_plot[, contam_series := factor(contam_series, levels = c("Il14H", "Ki11"))]

  pA <- ggplot(sm_plot,
               aes(x = true_alpha, y = sens_strict,
                   color = method, group = method)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
    facet_grid(track ~ contam_series) +
    scale_color_manual(values = method_palette) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       breaks = c(0,0.02,0.05,0.1,0.2,0.3,0.4,0.5)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    labs(x = "True contamination α", y = "Strict singlet sensitivity",
         color = NULL,
         title = paste0(PHASE, ": A — Strict sensitivity")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(color = guide_legend(ncol = 4))

  pB <- pA + aes(y = prec_strict) +
    labs(y = "Strict singlet precision", title = paste0(PHASE, ": B — Strict precision"))

  pC <- pA + aes(y = f1_strict) +
    labs(y = "Strict singlet F1", title = paste0(PHASE, ": C — Strict F1"))

  pD <- pA + aes(y = sens_relaxed) +
    labs(y = "Relaxed singlet sensitivity",
         title = paste0(PHASE, ": D — Relaxed sensitivity (incl. weak_doublet[B73])"))

  pE <- pA + aes(y = f1_relaxed) +
    labs(y = "Relaxed singlet F1",
         title = paste0(PHASE, ": E — Relaxed F1"))

  pF <- pA + aes(y = doublet_detection_rate) +
    labs(y = "True doublet detection rate",
         title = paste0(PHASE, ": F — Doublet preservation"))

  pdf_path <- file.path(EVALDIR, paste0(PHASE, "_compare.pdf"))
  pdf(pdf_path, width = 13, height = 9)
  print(pA); print(pB); print(pC); print(pD); print(pE); print(pF)
  dev.off()
  cat("Wrote", pdf_path, "\n")

} else {
  # Phase 3/4 — Root1 operational panels
  sm <- copy(summary_metrics)

  # Order methods by F1 descending for readability
  method_order <- sm[, .(mean_f1 = mean(f1_relaxed, na.rm = TRUE)), by = method][
    order(-mean_f1), method]
  sm[, method := factor(method, levels = method_order)]

  pA <- ggplot(sm, aes(x = method, y = f1_relaxed, fill = method)) +
    geom_col() +
    facet_wrap(~ track, scales = "free_x") +
    scale_fill_manual(values = method_palette, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    labs(x = NULL, y = "Relaxed F1",
         title = paste0(PHASE, ": A — Relaxed F1 by config")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  pB <- pA + aes(y = sens_relaxed) +
    labs(y = "Relaxed sensitivity",
         title = paste0(PHASE, ": B — Relaxed sensitivity"))

  pC <- pA + aes(y = prec_relaxed) +
    labs(y = "Relaxed precision",
         title = paste0(PHASE, ": C — Relaxed precision"))

  pD <- pA + aes(y = doublet_rate) +
    labs(y = "Fraction called doublet",
         title = paste0(PHASE, ": D — False-doublet rate (lower is better)"))

  pE <- pA + aes(y = wrong_genome_rate) +
    labs(y = "Fraction wrong-genome calls",
         title = paste0(PHASE, ": E — Wrong-genome rate (lower is better)"))

  pdf_path <- file.path(EVALDIR, paste0(PHASE, "_compare.pdf"))
  pdf(pdf_path, width = 13, height = 9)
  print(pA); print(pB); print(pC); print(pD); print(pE)

  # Phase 3 only: per-bin breakdown if depth_bin is present
  if (USE_BINS && "depth_bin" %in% names(sm) && any(!is.na(sm$depth_bin))) {
    pBin <- ggplot(sm[!is.na(depth_bin)],
                   aes(x = depth_bin, y = f1_relaxed,
                       color = method, group = method)) +
      geom_line() + geom_point(size = 2) +
      facet_wrap(~ track) +
      scale_color_manual(values = method_palette) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
      labs(x = "Depth bin", y = "Relaxed F1",
           title = paste0(PHASE, ": Per-bin F1 by config")) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 7)) +
      guides(color = guide_legend(ncol = 4))
    print(pBin)
  }

  dev.off()
  cat("Wrote", pdf_path, "\n")
}

# =============================================================================
# Phase 4 — sub1k vs full Root1 consistency check
# =============================================================================
if (PHASE == "phase4") {
  cat("\n=== Phase 4: sub1k vs full Root1 consistency check ===\n")
  cat("(Looks for matching configs in earlier Phase 1/3 outputs;\n")
  cat(" missing earlier outputs are skipped silently.)\n\n")

  consistency_rows <- list()
  for (cfg in names(CONFIGS)) {
    # Try Phase 1 sub1k C0 (factorial_2026-04-07) for the C0 config
    if (cfg == "C0") {
      for (panel in c("sub1k_A","sub1k_B")) {
        f <- file.path("Root1_rep1", panel, "genotyping_runs",
                        "factorial_2026-04-07", paste0("C0_alpha005"),
                        paste0(panel, "_cells_calls.tsv.gz"))
        if (file.exists(f)) {
          dt <- fread(f)
          dt[, `:=`(source = paste0("phase1_", panel), method = cfg)]
          consistency_rows[[length(consistency_rows)+1]] <- dt
        }
      }
    } else {
      # Look for Phase 3 sub1k stacked output for the winner
      for (panel in c("sub1k_A","sub1k_B")) {
        # Try both today's date and a few recent ones
        for (d in c(EVAL_DATE, format(Sys.Date()-1, "%Y-%m-%d"),
                    format(Sys.Date()-2, "%Y-%m-%d"))) {
          f <- file.path("Root1_rep1", panel, "genotyping_runs",
                          paste0("factorial_phase3_", d), cfg,
                          paste0(panel, "_cells_calls.tsv.gz"))
          if (file.exists(f)) {
            dt <- fread(f)
            dt[, `:=`(source = paste0("phase3_", panel), method = cfg)]
            consistency_rows[[length(consistency_rows)+1]] <- dt
            break
          }
        }
      }
    }
    # Phase 4 source is the merged_all already loaded above
    f4_rows <- merged_all[method == cfg]
    if (nrow(f4_rows) > 0) {
      x <- copy(f4_rows)
      x[, source := "phase4_full"]
      consistency_rows[[length(consistency_rows)+1]] <- x
    }
  }

  if (length(consistency_rows) > 0) {
    cons <- rbindlist(consistency_rows, fill = TRUE)

    # Filter to barcodes with >= MIN_READS reads
    if ("n_reads" %in% names(cons)) {
      n_before <- nrow(cons)
      cons <- cons[n_reads >= MIN_READS]
      cat("  consistency: filtered to n_reads >=", MIN_READS, ":",
          n_before, "->", nrow(cons), "rows\n")
    }

    # Recompute metrics on each (source, method) pair
    cons[, called_b73 := (!is.na(genome_1) & genome_1 == "B73") |
                          (!is.na(called_g1) & called_g1 == "B73")]
    cons[, is_singlet_call_relaxed := call %in% c("single_clean","dirty_singlet") |
                                       (call == "weak_doublet" & called_b73)]
    cons[, is_correct_singlet_relaxed := is_singlet_call_relaxed & called_b73]
    cons[, is_doublet_call := (call == "doublet")]

    cons_metrics <- cons[, .(
      n_total = .N,
      n_correct_relaxed = sum(is_correct_singlet_relaxed),
      n_doublet = sum(is_doublet_call),
      sens_relaxed = sum(is_correct_singlet_relaxed) / .N,
      doublet_rate = sum(is_doublet_call) / .N
    ), by = .(source, method)]

    fwrite(cons_metrics,
           file.path(EVALDIR, "phase4_consistency_metrics.tsv"), sep = "\t")
    cat("Wrote phase4_consistency_metrics.tsv\n")
    print(as.data.frame(cons_metrics[order(method, source)]),
          row.names = FALSE, digits = 3)

    # Side-by-side bar panel
    pCons <- ggplot(cons_metrics,
                    aes(x = source, y = sens_relaxed, fill = source)) +
      geom_col() +
      facet_wrap(~ method) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
      labs(x = NULL, y = "Relaxed sensitivity",
           title = "Phase 4: sub1k vs full Root1 consistency",
           fill = NULL) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

    pdf_path <- file.path(EVALDIR, "phase4_consistency.pdf")
    pdf(pdf_path, width = 11, height = 7)
    print(pCons)
    dev.off()
    cat("Wrote", pdf_path, "\n")
  } else {
    cat("(No earlier sub1k outputs found — skipping consistency panel)\n")
  }
}

cat("\n=== eval_phase_factorial.R complete ===\n")
cat("Output dir:", EVALDIR, "\n")
for (f in list.files(EVALDIR, full.names = FALSE)) cat("  ", f, "\n")
