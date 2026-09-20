# Figure 4: robustness and generalizability

Four scripts, one per group of panels; the manuscript figure is assembled by hand from their
outputs. All four source `analysis/_helpers/fig4_helpers.R` and **must be run from the repo
root**, since every path in them is relative to it. Inputs are read from `data/processed/`,
outputs are written to `figures/main/fig4/` (created on the fly).

| Script | Panels | Inputs (repo paths) | Run | Notes |
|---|---|---|---|---|
| `fig4_part1.R` | A, B, C (synthetic benchmark) | `data/processed/synthetic/{synthetic,synthetic_disc}/` (Track B, Track B-disc): `alpha_000/cell_map_ref_chunks/*_filtered.tsv.gz` (A); `<ds>/genotyping_runs/factorial_phase2_2026-04-09/C0/<ds>_cells_calls.tsv.gz` (B); `<ds>/decontam_without_design_alpha05_C0/<ds>_{pre,post}_barcode_genome_counts.tsv.gz` (C); 15 datasets per track, `alpha_000` and `alpha_{002,005,010,020,030,040,050}_{Il14H,Ki11}` | `Rscript analysis/fig4_robustness/fig4_part1.R` or `sbatch analysis/fig4_robustness/fig4_part1.sh` | A loads every filtered read chunk at alpha = 0 (16 GB is comfortable). Outputs `Fig4A_synth_density.pdf`, `Fig4B_synth_counts.pdf`, `Fig4C_synth_barnyard.pdf` and the composite `Fig4_part1_AtoC.{pdf,png}` |
| `fig4_part2.R` | D, E, F, G (maize root library, 26 NAM genomes) | `data/processed/marand2021_B73_root/Root1_rep1/`: `genotyping_runs/_archive/eval_xmap_v2_phi_B73.tsv` (D); `genotyping_runs/factorial_phase4_2026-04-28/{C0,S10_full_bic3}/Root1_rep1_cells_calls.tsv.gz` (E, F); `eval/phase4_2026-04-28/phase4_summary_metrics.tsv` (G, written by `workflows/03_genotyping/eval_phase_factorial.R phase4`) | `Rscript analysis/fig4_robustness/fig4_part2.R` or the `.sh` | E and F drop `low_reads` and keep `n_reads >= 250` (empty calls kept). phi in D is reported descriptively, as a measure of reference-panel redundancy, not as an applied correction: C0 runs with xmap off, `S10_full_bic3` (mq50, friend rescue OFF, xmap ON, w_amb 0.5, bic_margin 3) with xmap on. Outputs `Fig4_{D_phi,E_call_dist,F_genome1_by_call,G_wrong_genome}.pdf` and `Fig4_part2_DtoG.{pdf,png}` |
| `fig4_part3.R` | H, I, J, K (AmbientMapper vs Souporcell) | `data/processed/zhang2024/<s>/` for `B73Mo17_rep1`, `B73Mo17_rep2`, `multiGenotypes_rep1`: `genotyping_runs/4cfg_2026-05-01/C0/<s>_cells_calls.tsv.gz`; `souporcell/supervised/<s>.min500/{clusters.tsv, Genotype_ID_key.v2.txt}` (Souporcell v2.1 on the raw, pre-cleaning BAMs, `workflows/03b_variant_based_comparison/`) | `Rscript analysis/fig4_robustness/fig4_part3.R` or the `.sh` (24 GB) | B73/Mo17 replicates are pooled; the join key is the 26-character barcode prefix. `low_reads` is structurally absent after the Souporcell `.min500` filter. J and K use Souporcell singlets only, column-normalised. Outputs `Fig4_{H_B73Mo17_method,I_multi_method,J_B73Mo17_genotype,K_multi_genotype}.pdf` and `Fig4_part3_HtoK.{pdf,png}` |
| `fig4_part4.R` | L, M, N (allele purity before vs after cleaning) | `data/processed/zhang2024/<s>/diagnostics/06_48_barcode_purity/<s>_barcode_purity.tsv.gz` for the same three datasets (WASP-corrected barcode x 1-Mb-block purity on raw vs cleaned BAMs, from the 06_48 step of `workflows/03b_variant_based_comparison/`) | `Rscript analysis/fig4_robustness/fig4_part4.R` or the `.sh` | Purity is absolute (WASP-corrected). `weak_doublet` is folded into singlet, evidence in Table S4. Facets with fewer than 1,000 barcodes are dropped (the multi-genotype singlet facet). M and N are one patchwork object, `pM`, saved as `Fig4_MN_reads_removed.pdf` and lettered separately in the assembled figure. Also `Fig4_L_purity_abs.pdf` and `Fig4_part4_LtoN.{pdf,png}` |

## Helpers

`analysis/_helpers/fig4_helpers.R` provides `TRUE_GENOME`, `call_levels`, `call_colors`,
`genome_highlight_colors`, the genotyping-configuration encoding shared with Figs S5 and S6
(`CONFIG_ORDER`, `CONFIG_LABELS`, `CONFIG_GROUPS`, `CONFIG_COLORS`, `CONFIG_SHAPES`,
`CONFIG_PARAM_TABLE`), the synthetic-benchmark constants (`SYN_TRACKS`, `SYN_DATASETS`,
`SYN_ALPHA_LEVELS`, the palettes, `parse_synthetic_dataset()`) and the shared heatmap style
(`HEAT_PALETTE()`, `heat_theme`). `analysis/supplementary/figS5.R` (synthetic benchmark
validation, Phase 2) and `analysis/supplementary/figS6.R` (Root1 sub1k_B validation, Phase 3)
source the same file.

## Running on a cluster

The `.sh` files are SLURM wrappers. Submit them from the repo root after `mkdir -p _logs`; each
one changes to the repo root itself (`git rev-parse --show-toplevel`) before calling `Rscript`.
Activate the environment from `environment.yml` first; the commented placeholder in each wrapper
marks where that line belongs. R packages: tidyverse, patchwork, scales, data.table (hexbin is
optional in part 4 and falls back to `geom_bin2d`).
