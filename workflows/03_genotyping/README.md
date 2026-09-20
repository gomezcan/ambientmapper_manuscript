# 03 — Genotyping (AmbientMapper extract, filter, chunks, assign, genotyping)

Scripts that produced the per-barcode genotype calls (`*_cells_calls.tsv.gz`) behind the
manuscript, one subfolder per dataset, plus the shared configuration registry and the
factorial evaluator. Ported from the run-of-record scripts; parameter values are unchanged.
Only the scripts behind a result reported in the paper ship here; exploratory configurations
and archived runs do not.

## Conventions

- Every script defines `PROJECT_ROOT` (default: current directory) and works inside
  `${PROJECT_ROOT}/5_AmbientDetection/`, mirroring the original layout
  (`3_Mapping/`, `5_AmbientDetection/<dataset>/`, `5_AmbientDetection/configs/`).
  Export `PROJECT_ROOT` before `sbatch`.
- Scripts that source the shared helper or call a sibling script locate this repository
  through `REPO_ROOT="$(git rev-parse --show-toplevel)"`; submit from inside the checkout.
- The `--config` JSON files live in `config/` at the repo root (`SM2v2`, `Root1_rep1`,
  `B73Mo17_rep1`, `B73Mo17_rep2`, `multiGenotypes_rep1` `.ambientmapper.json`); only the two
  sub1k panel configs, which exist nowhere else, sit in `configs/` here. All carry the literal
  placeholder `${PROJECT_ROOT}` in every path. Copy them to
  `${PROJECT_ROOT}/5_AmbientDetection/configs/` (the scripts read `configs/<name>.json`
  relative to that directory) and substitute the absolute root once, e.g.
  `sed -i "s#\${PROJECT_ROOT}#$PROJECT_ROOT#g" *.json`.
- `GENOMES_ROOT` (directory holding `Zea/Zm_<NAM>_REFERENCE_NAM_*/`, MaizeGDB assemblies)
  and `BWA_INDEX_ROOT` (directory holding `Zea/NAN_Indexes/Index_Zm_<genome>_bwa`) must be
  exported for the synthetic benchmark scripts that read reference sequence or indexes.
- SLURM `--account` and `--mail-*` directives were removed; `# conda activate <env from
  environment.yml>` marks where the environment was activated. `AMBIENTMAPPER_REPO` may
  point at a local checkout of the tool so the run log records its commit.
- The dated tags inside output paths (`factorial_2026-04-07`, `factorial_phase2_2026-04-09`,
  `factorial_phase4_2026-04-28`, `4cfg_2026-05-01`) are directory names that the figure
  scripts resolve; they are parameters, not session notes.

## Shared files

| File | Role |
|---|---|
| `_genotyping_configs.sh` | Sourced registry of the 55 named genotyping configurations (C0, C1a to C1g, C2a to C2e, C3a to C3e, C4a/b/c/d, Cxmap_mq50, S01 to S13), the friend-rescue routing and the `ambientmapper genotyping` invocation. Winner-only mode is OFF and beta = 10 in every run that goes through it (comment above the invocation is the only record). |
| `eval_phase_factorial.R` | `Rscript eval_phase_factorial.R <phase2|phase3|phase4> [<date>]`, run from `${PROJECT_ROOT}/5_AmbientDetection/`. phase2 -> `synthetic/eval/phase2_<date>/phase2_summary_metrics.tsv` (Fig. S5); phase3 -> `Root1_rep1/eval/phase3_<date>/phase3_summary_metrics.tsv` (Fig. S6); phase4 -> `Root1_rep1/eval/phase4_<date>/phase4_summary_metrics.tsv` (Fig. 4G). Local R, no SLURM. |
| `configs/` | `Root1_rep1_sub1k_A.ambientmapper.json`, `Root1_rep1_sub1k_B.ambientmapper.json` (panel configs read by `02_09` to `02_11`). The five dataset JSONs are in `config/`; the SM2 sample list `SM2_AtB73.list.tsv` is in `data/metadata/barcode_lists/`. |

## Execution order per dataset

### `sm2/` — SM2, the combined-reference arm (Fig. 1B to E)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `01_01_pipeline_B73_At.sh configs/SM2_AtB73.list.tsv` | 4 h, 35 cpu, 40 G | `SM2/final/SM2_cells_calls.tsv.gz`, `SM2/cell_map_ref_chunks/*_filtered.tsv.gz` (input of `04_decontamination/sm2/01_02`) |

The sample list is `data/metadata/barcode_lists/SM2_AtB73.list.tsv`. The run of record named
the BAMs under `4_MappingCleaning/ambientmapper_input/`; that directory was later consolidated
into `3_Mapping/ambientmapper_input/` with the same file names, and the shipped list uses the
`3_Mapping/` location.

### `sm2v2/` — SM2v2, the independent-mapping arm (Fig. 2, 3, 5, S7, S8, Tables S1, S3)

Input BAMs: `3_Mapping/ambientmapper_input/SM2_{B73v5,TAIR10}_scifiATAC.mq10.BC.rmdup.mm.bam`
(produced by stage 02, `01_00_merge_SM2v2_inputs.sh`).

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `01_04_extract_SM2v2.sh` | 4 h, 16 cpu, 40 G | `SM2v2/qc/{B73,At}_QCMapping.txt` |
| 2 | `01_05_filter_SM2v2.sh` | 2 h, 16 cpu, 40 G | `SM2v2/filtered_QCFiles/filtered_*_QCMapping.txt` |
| 3 | `01_06_chunks_SM2v2.sh` | 1 h, 4 cpu, 16 G | `SM2v2/cell_map_ref_chunks/*_cell_map_ref_chunk_<N>.txt` |
| 4 | `01_07_prepare_SM2v2.sh` | 2 h, 8 cpu, 64 G | Parquet conversion of the filtered QCMapping files |
| 5 | `01_08a_assign_bootstrap_SM2v2.sh` | 4 h, 16 cpu, 64 G | `SM2v2/ExplorationReadLevel/global_{edges,ecdf}.npz` plus the first scored chunks |
| 6 | `01_08_assign_array_SM2v2.sh` | 6 h, 16 cpu, 64 G, array 0-15 | `SM2v2/cell_map_ref_chunks/*_filtered.tsv.gz` (assign at alpha 0.05, k 10) |
| 7 | `01_09_genotyping_SM2v2.sh` | 8 h, 8 cpu, 32 G | `SM2v2/genotyping_runs/SM2v2_C0/SM2v2_cells_calls.tsv.gz`, linked as `SM2v2/final/` for stage 04 |

`00_prepare_parquet.sh` (12 h, 8 cpu, 64 G, array 0-3) is the same Parquet conversion for
Root1_rep1, B73Mo17_rep1, B73Mo17_rep2 and multiGenotypes_rep1; it runs after their
`filter` step and before their `assign` step.

### `root1/` — Root1_rep1, 26 NAM genomes (Fig. 4D to G)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `02_01_pipeline_Root1.sh` | 20 h, 30 cpu, 160 G | extract, filter, chunks (the assign in this script used the legacy alpha 1e-6 bundle and was superseded) |
| 2 | `02_12a_pipeline_full_root1_alpha005.sh` | 48 h, 16 cpu, 128 G, resume-safe (about 120 h total) | assign at alpha 0.05 into `Root1_rep1/cell_map_ref_chunks/` and the friend-rescue knockout sibling `cell_map_ref_chunks_alpha005_friendwithout/` |
| 3 | `02_12b_repair_friendwithout.sh` (`RUN_REPAIR=1` after the audit) | 2 h, 2 cpu, 8 G | completes the knockout sibling if the relabel in step 2 stopped early |
| 4 | `02_12_genotyping_full_root1_phase4.sh` | 8 h, 16 cpu, 64 G, array 0-4 | `Root1_rep1/genotyping_runs/factorial_phase4_<date>/{C0,C3b_mq50,C4a_bic3,S04_nofr_mq50_xmap_xaUL,S10_full_bic3}/Root1_rep1_cells_calls.tsv.gz` (Fig. 4E, F) |
| 5 | `../eval_phase_factorial.R phase4 <date>` | local R | `Root1_rep1/eval/phase4_<date>/phase4_summary_metrics.tsv` (Fig. 4G) |

The phi table of Fig. 4D (`genotyping_runs/_archive/eval_xmap_v2_phi_B73.tsv`) comes from an
archived script on an archived run and is not part of this stage.

### `root1_sub1k/` — depth-balanced 1k-barcode panels (Fig. S6)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `02_07_subsample_root1_balanced.py` | local Python (pandas, numpy) | `Root1_rep1/sub1k_{A,B}/{barcodes,barcodes_by_genome,barcodes_with_bin}.tsv` from `3_Mapping/Root1_rep1/*_bc_counts.txt` |
| 2 | `02_08_filter_bams_subsample.sh` | 4 h, 2 cpu, 8 G, array 0-51 | per-genome panel BAMs under `sub1k_<panel>/bams/` |
| 3 | `02_09_pipeline_sub1k.sh` | 6 h, 16 cpu, 32 G, array 0-1 | extract to assign per panel at alpha 0.05 and alpha 1e-6, each with a `friendwith` and a `friendwithout` chunks dir |
| 4 | `02_10_genotyping_factorial_sub1k.sh` | 1 h, 8 cpu, 16 G, array 0-167 | Phase 1: 42 configs x 2 alphas x 2 panels under `genotyping_runs/factorial_2026-04-07/` |
| 5 | `02_11_genotyping_stacked_phase3_sub1k.sh` | 30 min, 8 cpu, 16 G, array 0-37 | Phase 3: 19 configs x 2 panels under `genotyping_runs/factorial_phase3_<date>/` |
| 6 | `../eval_phase_factorial.R phase3 <date>` | local R | `Root1_rep1/eval/phase3_<date>/phase3_summary_metrics.tsv` (Fig. S6, `TRACK = sub1k_B`) |

### `synthetic/` — synthetic benchmark, Track B and Track B-disc (Fig. 4A to C, S5)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `03_01_call_peaks.sh` | 4 h, 4 cpu, 32 G | `synthetic/peaks/B73_peaks_filtered.bed` (MACS3 on the Root1 B73 BAM, top 10,000 autosomal peaks >= 150 bp) |
| 2 | `03_02_find_orthologs.sh` (+ `.py`) | 2 h, 8 cpu, 16 G; needs `GENOMES_ROOT`, minimap2 | `synthetic/orthologs/ortholog_map.tsv`, `{B73,Il14H,Ki11}_peak_seqs.fa` (3,831 triplets) |
| 3 | `03_03_filter_disc_peaks.py` | local Python | `synthetic_disc/orthologs/` (843 peaks with >= 1 expected SNP per 75 bp read in both contaminants) |
| 4 | `03_04_simulate_reads.sh` | 2 h, 4 cpu, 8 G; needs ART | `synthetic/reads/{genome}_{1,2}.fq.gz` (75 bp PE, 100x, insert 166 +/- 30, seed 42) |
| 5 | `03_05_assign_barcodes.sh` (+ `.py`) | 1 h, 1 cpu, 16 G | `synthetic/barcoded/templates.tsv` and 15 datasets `alpha_000`, `alpha_{002..050}_{Il14H,Ki11}` with `truth_table.tsv` |
| 6 | `03_06_map_reads.sh` | 2 h, 8 cpu, 16 G, array 0-44; needs `BWA_INDEX_ROOT` | `synthetic/mapping/<dataset>/syn_{B73v5,Il14H,Ki11}.sorted.bam` |
| 7 | `03_07_extract_metrics.sh` (+ `.py`) | 2 h, 1 cpu, 8 G, array 0-14 | `synthetic/<dataset>/qc/*_QCMapping.txt` (barcodes parsed from read names) |
| 8 | `03_08_run_pipeline.sh` | 8 h, 16 cpu, 32 G, array 0-14 | per-dataset `config.json`, chunks and assign |
| 9 | `03_16_disc_pipeline.sh` | 12 h, 16 cpu, 32 G; needs ART, `BWA_INDEX_ROOT` | steps 4 to 8 for Track B-disc in one job (ART at 450x) |
| 10 | `03_09_rerun_assign.sh` | 6 h, 16 cpu, 32 G, array 0-29 | re-runs extract to assign on both tracks so `frag_loc` and the `rescued` class exist |
| 11 | `03_21a_synthetic_friend_relabel.sh` | 30 min, 2 cpu, 4 G, array 0-29 | `cell_map_ref_chunks_friendwithout/` siblings |
| 12 | `03_21_genotyping_synthetic_factorial.sh` | 2 h, 8 cpu, 16 G, array 0-299%60 | Phase 2: 10 configs x 15 datasets x 2 tracks under `genotyping_runs/factorial_phase2_<date>/` (Fig. 4B reads C0) |
| 13 | `../eval_phase_factorial.R phase2 <date>` | local R | `synthetic/eval/phase2_<date>/phase2_summary_metrics.tsv` (Fig. S5) |

### `zhang2024/` — B73Mo17_rep1, B73Mo17_rep2, multiGenotypes_rep1 (Fig. 4H to K)

Input BAMs: `3_Mapping/ambientmapper_input/<sample>_<genome>_scifiATAC.mq10.BC.rmdup.mm.bam`
(stage 02, scifi-demux step 2 merges).

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `04_01_pipeline_B73Mo17_Multi.sh` | 24 h, 16 cpu, 128 G, array 0-2 | extract, filter (min barcode freq 5), chunks (100 cells), assign at alpha 0.05 for all three libraries (run `00_prepare_parquet.sh` between filter and assign) |
| 2a | `04_04_scoring_array_B73Mo17_rep1.sh` | 1 d, 16 cpu, 128 G, array 0-3 | finishes the rep1 scoring with `--score-chunk-range` |
| 2b | `04_02_resub_multiGenotypes_scoring.sh` then `04_03_scoring_array_multiGenotypes.sh` | 7 d, 16 cpu, 128 G; 2 d, 16 cpu, 128 G, array 0-7 | finishes the multiGenotypes scoring (score batch 50, 6 workers) |
| 3a | `04_05a_genotyping_4cfg_B73Mo17.sh` | 36 h, 16 cpu, 64 G, array 0-7 | `<sample>/genotyping_runs/4cfg_2026-05-01/{C0,C2a_xmap,C3b_mq50,Cxmap_mq50}/` for rep1 and rep2 |
| 3b | `04_05b_genotyping_4cfg_multi.sh` | 48 h, 16 cpu, 128 G, array 0-3 | the same four configurations for multiGenotypes_rep1 |

Only the `C0` output of step 3 is displayed (Fig. 4H to K read `4cfg_2026-05-01/C0/<sample>_cells_calls.tsv.gz`);
the other three cells of the grid are held in reserve. `05_01_pipeline_MultiGenotype.sh`
(8 h, 35 cpu, 80 G) is the earlier single-shot `ambientmapper run` over multiGenotypes_rep1;
the sample-list TSV it took as `$1` is not preserved in `configs/`.

## Notes

- Friend (co-localization) rescue has no assign-step switch in AmbientMapper. Every
  `*friendwithout*` chunks directory is a post-hoc copy in which `assigned_class`
  `rescued` is rewritten to `ambiguous` (`02_09`, `02_12a`/`02_12b`, `03_21a`).
- The global cross-mapping profile (phi, `--no-xmap` omitted) is enabled in fifteen
  registry branches, including `S10_full_bic3` shown in Fig. 4G; it is off in `C0` and in
  the shipped SM2v2 run. Fig. 4D reports phi descriptively.
- The genotyping pass is not run-to-run deterministic at the level of `C_all` winner
  masses, so small cross-configuration contrasts on the synthetic panel are within
  run-to-run noise.
- Decontamination and BAM cleaning of these calls live in `workflows/04_decontamination/`;
  the SNP-based comparison (Souporcell, WASP, purity) in `workflows/03b_variant_based_comparison/`.
