# 04 — Decontamination (AmbientMapper decontam and clean-bams)

Scripts that turned the stage 03 genotype calls into per-read drop lists, barcode-level
contamination tables and cleaned BAMs, one subfolder per dataset. Ported from the
run-of-record scripts; parameter values are unchanged.

## Conventions

- Every script defines `PROJECT_ROOT` (default: current directory) and works inside
  `${PROJECT_ROOT}/5_AmbientDetection/` (the combined-genome scripts inside
  `${PROJECT_ROOT}/6_socrates/`), mirroring the original layout. Export `PROJECT_ROOT`
  before `sbatch`; scripts that call a sibling locate this checkout through
  `REPO_ROOT="$(git rev-parse --show-toplevel)"`.
- Two AmbientMapper modes per library: `decontam_with_design_*` (WD, plate design file
  guides the allowed genome per well, `--ambiguous-policy design_rescue`) and
  `decontam_without_design_*` (ND, genotype inferred from the calls,
  `--ambiguous-policy top1_rescue`). Shared parameters in every run: `--decontam-alpha 0.05`,
  `--doublet-policy top1`, `--indist-policy top1`, `--safe-keep-delta-as 3`,
  `--min-reads-post-clean 100`, `--min-allowed-frac-post-clean 0.90`, `--chunksize 1000000`.
- The plate design read by the WD passes (`configs/PlateDesign_<sample>_ATAC.txt`, i.e.
  `PlateDesign_SM2_ATAC.txt` and `PlateDesign_SM2v2_ATAC.txt`, the latter a symlink to the
  former in the run tree) is the two-block AmbientMapper design of record, shipped as
  `data/metadata/plate_designs/PlateDesign_scifiATAC_B73_Arabidopsis.txt` (B73 wells
  columns 1 to 4, Arabidopsis wells columns 9 to 12). Install it under
  `${PROJECT_ROOT}/5_AmbientDetection/configs/` under both names. Do not use
  `config/PlateDesign_SM2_ATAC.txt` for this: despite the identical file name it is the
  three-block demultiplexing design (with the MuDR wells), which `design_rescue` does not accept.
- Key outputs per run directory: `*_cells_calls.decontam.tsv.gz` (calls plus metrics),
  `*_barcode_policy.tsv.gz`, `*_pre/post_barcode_genome_counts.tsv.gz`,
  `*_pre/post_barcode_composition.tsv.gz`, `*_reads_to_drop.tsv.gz`, `*_decontam_params.json`.
- SLURM `--account` and `--mail-*` directives were removed; `# conda activate <env from
  environment.yml>` marks where the environment was activated.

## Two tool-level notes (apply to every run below)

- `--min-reads-post-clean 100` is a reporting flag, not an enforced gate: the post-clean
  metrics are computed after `reads_to_drop.tsv.gz` is written and reach only the summary
  tables, and `clean-bams` filters on `read_id` alone, so a barcode flagged
  `keep_postclean = False` keeps its surviving reads in the cleaned BAM (cleaning is
  read-level throughout).
- `--min-allowed-frac-post-clean 0.90` is inert: the post-clean allowed fraction is computed
  after the drop list is applied, so numerator and denominator are the same read set and
  the value is 1 by construction, and the threshold removed zero barcodes in every run.

## Execution order per dataset

### `sm2/` — SM2, the combined-reference arm (Fig. 1B to E)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `01_02_decontam_B73_At.sh SM2` | 2 h, 10 cpu, 40 G | `SM2/decontam_with_design_alpha05_v2/` (the WD tables read by Fig. 1B to E) and `SM2/decontam_without_design_alpha05_v2/` (no consumer) |
| 2 | `01_03_clean_bams_B73_At.sh SM2` | 3 h, 5 cpu, 40 G | `SM2/clean_bams_alpha05/*.Clean.bam` + `*.tn5.bed.gz` from the combined-genome (ZmATcombined) BAMs |

`01_03` names its inputs under `4_MappingCleaning/`, the directory of record at the time; it
was later consolidated into `3_Mapping/` (`3_Mapping/SM2_{At,B73}/` for the BAMs,
`3_Mapping/_archive/1_6_scifi_makeTn5bed.py` for the Tn5 BED script).

### `sm2v2/` — SM2v2, the independent-mapping arm (Fig. 2, 3, 5, S7, S8, Tables S1, S3)

Input: `SM2v2/final/SM2v2_cells_calls.tsv.gz` and `SM2v2/cell_map_ref_chunks/*_filtered.tsv.gz`
from `workflows/03_genotyping/sm2v2/01_09`.

| # | Script | Resources | Produces |
|---|---|---|---|
| 1a | `01_10a_decontam_SM2v2_with_design.sh` | 24 h, 10 cpu, 40 G | `SM2v2/decontam_with_design_alpha05_v2/` (WD: Fig. 2, Fig. 3, Table S1, Table S3) |
| 1b | `01_10b_decontam_SM2v2_without_design.sh` (in parallel with 1a) | 24 h, 10 cpu, 40 G | `SM2v2/decontam_without_design_alpha05_v2/` (ND: `*_decontam_params.json` in Table S1) |
| 2 | `01_11_clean_bams_SM2v2.sh` | 2 h, 8 cpu, 16 G, array 0-3 | `SM2v2/clean_bams_alpha05_C0_{wd,nd}/SM2_{B73v5,TAIR10}_*.Clean.bam` + `.tn5.bed.gz` (WD and ND arms of stage 05); uses the per-chromosome sidecar `workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py` |

The pre-clean tn5 BEDs and the WD/ND BED regeneration (`01_12`, `01_11c`) ship with stage 02.

### `combined_genome/` — SM2v2 drop list applied to the concatenated-reference BAMs (Fig. 5E, F; Fig. S2, S3)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `0_04_split_reads_to_drop_SM2v2.sh` | local, minutes | `6_socrates/_data/SM2v2_reads_to_drop.clean{At,B73}.tsv.gz` (WD drop list split by target BAM) |
| 2 | `0_05_clean_bams_SM2v2_combined.sh` | 3 h, 5 cpu, 40 G, array 0-1 | `6_socrates/_data/SM2v2_clean_bams_combined/Clean.SM2v2_{At,B73}_ZmATcombined_*.bam` |

The matching tn5 BED step (`0_06_make_tn5bed_SM2v2_combined.sh`) ships with stage 02.

### `synthetic/` — both synthetic tracks (Fig. 4C)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `03_22_decontam_synthetic_nodesign.sh` | 45 min, 8 cpu, 16 G, array 0-29 | `<track>/<dataset>/decontam_without_design_alpha05_C0/` for 15 datasets x 2 tracks, on the Phase 2 `C0` calls (`factorial_phase2_2026-04-09`); Fig. 4C reads the pre/post `barcode_genome_counts` |

### `zhang2024/` — B73Mo17_rep1, B73Mo17_rep2, multiGenotypes_rep1 (Fig. 4L to N, Table S1)

Input: `<sample>/genotyping_runs/4cfg_2026-05-01/C0/<sample>_cells_calls.tsv.gz` from
`workflows/03_genotyping/zhang2024/04_05a` and `04_05b`.

| # | Script | Resources | Produces |
|---|---|---|---|
| 1a | `04_06_decontam_B73Mo17_C0_nd.sh` | 48 h, 10 cpu, 40 G, array 0-1 | `B73Mo17_rep{1,2}/decontam_without_design_alpha05_C0/` |
| 1b | `05_03_decontam_MultiGenotype_C0_nd.sh` | 72 h, 10 cpu, 80 G | `multiGenotypes_rep1/decontam_without_design_alpha05_C0/` |
| 2a | `04_07_clean_bams_B73Mo17.sh` | 6 h, 5 cpu, 32 G, array 0-3 | `B73Mo17_rep{1,2}/clean_bams_alpha05_C0_nd/*.Clean.bam` + `.tn5.bed.gz` |
| 2b | `05_04_clean_bams_multi.sh` | 8 h, 5 cpu, 48 G, array 0-6 | `multiGenotypes_rep1/clean_bams_alpha05_C0_nd/*.Clean.bam` + `.tn5.bed.gz` |
| 3a | `04_07c_tn5bed_regen_B73Mo17_rep1.sh` | 4 h, 12 cpu, 24 G, array 0-1 | regenerated `.tn5.bed.gz` for rep1 after the sidecar `sort -u` fix |
| 3b | `05_04c_tn5bed_regen_multi.sh` | 6 h, 12 cpu, 32 G, array 0-6 | regenerated `.tn5.bed.gz` for multiGenotypes_rep1 |

The Clean BAMs of step 2 are the "clean" arm of the WASP and purity chain
(`workflows/03b_variant_based_comparison/`, Fig. 4L to N); the `*_decontam_params.json` of
step 1 feed Table S1. Steps 2 and 3 call `3_Mapping/_archive/1_6_scifi_makeTn5bed.py` and
`workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py`, both shipped with stage 02.
Only the ND mode was run on these libraries (no plate design was used). Root1_rep1 was not
decontaminated for any panel.
