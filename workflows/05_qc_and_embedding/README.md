# 05 QC and embedding (Socrates)

Single-cell QC, cell calling, embedding, clustering and downstream annotation of the
scifi-ATAC B73 + Arabidopsis library (SM2 / SM2v2) with Socrates, run before and after
AmbientMapper cleaning. Three tracks share the engines in `common/`:

| Track | Data directory | What it asks | Manuscript panels |
|---|---|---|---|
| A. Combined-genome SM2, pre and post | `6_socrates/SM2/` (PreClean), `6_socrates/SM2v2_clean/` (PostClean), `6_socrates/compare/` | per-cell QC before vs after cleaning on the concatenated ZmATcombined reference | Fig. 5E, 5F, Fig. S2, Fig. S3 |
| B. Independent co-embedding (`part1_indep/`) | `6_socrates/SM2v2_indep/` | cross-species co-projection built from the two single-reference mappings | Fig. 1F to 1I, Fig. 5A to 5D, Fig. S3 |
| C. Plate arm (`part2_split/` + `common/`) | `6_socrates/SM2v2_plate/` | each plate mapped to its expected genome, the normal user workflow, Pre vs WD vs ND | Fig. 5G to 5L, Fig. S7, Fig. S8 |

The directory mirrors the source tree `6_socrates/0_scripts/` one to one (`setup/`, `common/`,
`part1_combined/`, `part1_indep/`, `part2_split/step*`), so the drivers' relative calls
(`${SCRIPTS}/common/...`, `${SCRIPTS}/part1_indep/...`) resolve unchanged.

## Pointing the scripts at data and at each other

Every shell driver anchors on two variables, both overridable from the environment:

```bash
export PROJECT_ROOT=/path/to/project      # the directory that holds 6_socrates/ (and 3_Mapping/ for the two setup scripts)
export SCRIPTS=/path/to/repo/workflows/05_qc_and_embedding   # optional; default ${PROJECT_ROOT}/6_socrates/0_scripts
```

`BASE="${PROJECT_ROOT}/6_socrates"` and all inputs and outputs keep their original relative
layout under it (`_data/_BED_files/`, `_data/_PeakFiles/`, `_data/_GenomeInfo/`, `_data/markers/`,
`SM2/`, `SM2v2_clean/`, `SM2v2_indep/`, `SM2v2_plate/`, `compare/`, `_logs/`). The simplest setup is
a symlink `6_socrates/0_scripts -> workflows/05_qc_and_embedding`, after which the `sbatch
0_scripts/...` commands quoted in the script headers work verbatim when submitted from
`6_socrates/`. The R engines take every path as a command-line argument; the one Python builder
with a hard-wired root (`common/4_1c_merge_db_markers.py`) reads `PROJECT_ROOT` too.

The figure scripts in `analysis/` read the same objects from
`data/processed/scifiATAC_B73_Arabidopsis/socrates/` (subtrees `SM2/`, `SM2v2_clean/`,
`SM2v2_indep/`, `SM2v2_plate/`, `compare/`, `_data/`), so copy or symlink `6_socrates/` there
once this stage has run.

Inputs that come from other stages: the per-genome and combined-genome tn5 BEDs
(`02_mapping`, `04_decontamination/combined_genome/`), the AmbientMapper input BAMs used only
to derive chromosome sizes (`0_2_setup_indep.sh`), and the genome annotations, which live
outside the project tree and are passed through `AT_GFF3` / `B73V5_GFF3` (Ensembl Plants TAIR10
release 60 and the MaizeGDB Zm-B73-REFERENCE-NAM-5.0 GFF3).

## Environment

The R + Socrates stack is the repo `environment.yml` (`conda activate ambientmapper-manuscript`);
each driver carries that activation as a commented placeholder where the original activated a
site-specific env. R packages used here: Socrates (github.com/plantformatics/Socrates), Seurat,
Matrix, data.table, FNN, MASS, viridis, ggplot2, dplyr, cowplot, tidyverse, qlcMatrix, devtools.
Socrates' `loadBEDandGenomeData()` requires `macs2` on `PATH` even when peaks are precomputed
(the drivers carry a `# module load macs2` placeholder). `0_07_macs3_SM2v2.sh` needs `macs3`,
`0_2_setup_indep.sh` needs `samtools`, `0_1_bam_to_tn5bed.py` needs `pysam` and `pigz`.

SEACells (step 5 of the plate arm) runs in its own environment because it pins scanpy, anndata,
numpy and jax: create it once, interactively, with `bash setup/0_10_install_seacells_env.sh`
(Python 3.10, `pip install SEACells`, plus `ipywidgets`, which SEACells imports through
`tqdm.notebook` and which every headless job fails without). `5_3_run_seacells_plate.sh` expects
that env as `seacells`.

## Stage tokens in file names

| Token in object / file names | Stage | Meaning |
|---|---|---|
| `SM2` (`SM2_At`, `SM2_B73`, `SM2_B73v5`, `SM2_TAIR10`, `SM2_At_TAIR10`, `SM2_B73_B73v5`) | PreClean | raw reads, before AmbientMapper |
| `Clean.SM2v2wd` (`Clean.SM2v2wd_B73v5`, `Clean.SM2v2wd_At_TAIR10`, ...) | WD | AmbientMapper with the plate design file (design-guided) |
| `Clean.SM2v2` (`Clean.SM2v2_B73v5`, `Clean.SM2v2_At_TAIR10`, ...) | ND | AmbientMapper without a design file (design-free) |

The mapping above is the one the independent (`SM2v2_indep/`) and plate (`SM2v2_plate/`) tracks
use. In the combined-genome track the cleaned objects also carry the plain `Clean.SM2v2_At` /
`Clean.SM2v2_B73` prefix; they are built by the decontamination stage's
`combined_genome/0_04_split_reads_to_drop_SM2v2.sh`, whose input path names the with-design
`reads_to_drop` list, while the project notes describe that pair as design-free. Check the
decontamination stage README before reading WD or ND into those two names.

Other tokens: `ZmATcombined` = the concatenated maize + Arabidopsis reference (track A);
`B73v5` / `TAIR10` = single-reference mappings (tracks B, C); `minDepth200` = the depth floor of
the meta-QC cascade; `updated_metadata_v1..v6` = the meta-QC cascade outputs (v6 = final cell
set); `mQCv6` = clustered on the v6 cell set; `_v7` = clustered metadata / objects (adds
`umap1`, `umap2`, `LouvainClusters`); `pcs_<p>.k_near_<k>.min_dis_<d>[.minc_<c>][.res_<r>][.mclst_<m>]`
= the frozen embedding and clustering configuration. `LouvainClusters` is the Socrates column
name; the algorithm is Leiden (`callClusters(cl.method = 4)`, Seurat `FindClusters(algorithm = 4)`).

## Track A: combined-genome SM2, pre and post

Run from `6_socrates/`, in this order (resources from the `#SBATCH` lines):

| # | Script | Does | Resources |
|---|---|---|---|
| A0 | `setup/0_07_macs3_SM2v2.sh` (array 0-3) | MACS3 peaks per sample (SM2_At, SM2_B73, Clean.SM2v2_At, Clean.SM2v2_B73) into `_data/_PeakFiles/` | 4 h, 16 GB, 4 cpu |
| A1 | `common/submit_chunked_pipeline_SM2v2.sh` (`POOL=<sample> OUTDIR=SM2|SM2v2_clean [N=5] [DEPTH=200]`) | submits the four-stage `afterok` chain below, once per sample | login node |
| A1a | `common/1_1a_chunk_bed.sh` | hash-partitions the tn5 BED into N barcode-disjoint chunks | 8 h, 8 GB, 6 cpu |
| A1b | `common/1_1b_per_chunk.sh` -> `1_1b_per_chunk.R` (array 0..N-1) | Socrates object per chunk (organelle removal, precomputed peaks, 500 bp tiles) | 12 h, 130 GB, 5 cpu |
| A1c | `common/1_1c_merge_and_qc.sh` -> `1_1c_merge_and_qc.R` | `mergeSocratesRDS` + `isCellv2` -> `step0_qc/<POOL>.raw.soc.rds` | 6 h, 180 GB, 5 cpu |
| A1d | `common/1_2_1_3_run.sh` -> `1_2_filter_lowQC_cells_scifiATAC_data.R`, `1_3_metaQC_scifiATAC_data.R` | QC filter + meta-QC cascade -> `updated_metadata.txt`, `minDepth200.updated_metadata_v1..v6.txt`, QC PDF | 2 h, 20 GB, 2 cpu |
| A2 | `part1_combined/2_0_0_Obj_integreation_SM2v2.sh` -> `2_0_0_Obj_integreation.R` | merges At + B73 objects of one stage -> `step1_integrate/<prefix>.full.SocObj.rds` | 1 h, 50 GB, 1 cpu |
| A3 | `part1_combined/2_0_3_fixed_set_compare.sh` -> `2_0_3_fixed_set_compare.R` | fixed-barcode pre vs post QC comparison (thresholds frozen from Pre) -> `compare/*.fixedSet.minDepth200.pre_post.txt` | 30 min, 24 GB, 2 cpu |

The meta-QC update (Fig. S2) and the UMAP grid scan (Fig. S3) that follow A2 ship as analysis
scripts (`analysis/supplementary/figS2.R`, `figS3.R`, `figS3A_replot.R`). The combined-genome
clustering (`3_0_0_Normalization_clustering`) is not ported: the Fig. 1F to 1H artwork now comes
from track B. `common/1_1_QC_scifiATAC_data.R` is the unchunked form of A1b + A1c and is called
directly by the plate driver.

## Track B: independent co-embedding (`part1_indep/`)

Self-contained: the QC engines are copied in as `1_1_qc_build.R`, `1_2_qc_filter.R`,
`1_3_qc_metaqc.R` (identical to `common/1_1`, `1_2`, `1_3`). Order:

| # | Script | Does | Resources |
|---|---|---|---|
| B0 | `0_1_bam_to_tn5bed.py`, `bash 0_2_setup_indep.sh` | per-chrom parallel BAM -> tn5 BED; chr sizes from the BAM headers, annotation and BED symlinks | interactive |
| B1 | `1_0_qc_run.sh` (array 0-3; add 4,5 for WD) -> `1_1`, `1_2`, `1_3` | per-genome QC, Pre + ND (+ WD) x {B73v5, TAIR10} -> `SM2v2_indep/step0_qc/` | 3 h, 40 GB, 10 cpu |
| B2a | `2_1_gridscan_pergenome.sh` -> `2_1_gridscan_pergenome.R` (array 0-1) | pcs x k_near x min_dist grid on the Pre objects; score = kNN preservation - 0.5 x max QC correlation; library mixing reported as obs/exp, never scored | 2 h, 48 GB, 6 cpu |
| B2b | `2_2_resolution_scan.sh` -> `2_2_resolution_scan.R` (array 0-1) | pcs x resolution bootstrap-stability scan (marker-free) | 10 h, 120 GB, 10 cpu |
| B2c | `2_3_cluster.sh` -> `2_3_cluster.R` (array 0-1) | final per-genome clustering at the frozen config -> `step2_cluster/*.full.SocObj_v7.*.rds` | 6 h, 80 GB, 8 cpu |
| B3a | `3_1_coembed_build.sh` -> `3_1_coembed_build.R` (array 0-1) | stacks the two per-genome tile matrices into one joint object; `Genome` = plate of origin | 2 h, 64 GB, 4 cpu |
| B3b | `3_2_coembed_gridscan.sh` -> `3_2_coembed_gridscan.R` (array 0-1) | Fig. S3 grid on the co-embedding (kNN preservation, genome mixing, QC leakage) | 8 h, 96 GB, 8 cpu |
| B3c | `PCS=.. KNN=.. MD=.. RES=.. STAGE=Pre sbatch 3_3_coembed_cluster.sh` -> `3_1b_coembed_meta_enrich.R`, `3_3_coembed_cluster.R` | grafts the v6 QC schema onto the co-embed metadata, then clusters -> `coembed/Pre_cluster/*.updated_metadata_v7.*.txt` (Fig. 1F to 1H input) | 4 h, 80 GB, 8 cpu |

The Fig. 1I feature-ablation control (`3_4_coembed_ablation`) ships as
`analysis/fig1_ambient_contamination/fig1I_coembed_ablation.R`. The configuration search that
froze the B2 settings (`2_config_search/`) is reference material and is not ported.

## Track C: plate arm (`part2_split/` with `common/` engines)

| # | Script | Does | Resources |
|---|---|---|---|
| C0 | `step0_qc/0_09_split_beds_by_plate.sh` (array 0-5) | splits each per-genome BED by plate tag: At plate -> TAIR10, B73 plate -> B73v5; cross-plate barcodes discarded by design | 2 h, 8 GB, 4 cpu |
| C1 | `step0_qc/1_QC_scifiATAC_SM2v2_plate.sh` (array 0-3; add 4,5 for ND) -> `common/1_1`, `1_2`, `1_3` | per-plate QC with the isCell model refit on the plate's own barcodes -> `SM2v2_plate/step0_qc/` | 3 h, 40 GB, 10 cpu |
| C2a | `step1_cluster_opt/3_0_1_scan_plate.sh` -> `3_0_1_resolution_scan.R` (array 0-1) | pcs x resolution stability scan on the plate Pre objects | 6 h, 90 GB, 10 cpu |
| C2b | `step1_cluster_opt/3_0_1b_umap_panels.sh` -> `3_0_1b_umap_panels.R`; `3_0_1c_umap_kgrid.sh` -> `3_0_1c_umap_kgrid.R` | UMAP panels over pcs x res (and k_near x m.clst) coloured by cluster and by QC, the by-eye continuum check for At | 1.5 h, 90 GB, 10 cpu each |
| C3 | `step2_cluster/3_0_cluster_plate_v6.sh` -> `3_0_0b_cluster_pergenome.R` (array 0-5) | final clustering at the frozen plate configs (At pcs5/k20/md0.3/minc50/res0.3/m.clst40; B73 pcs20/k20/md0.05/minc50/res0.3) -> `step2_cluster/*.mQCv6.*_v7.*` | 4 h, 120 GB, 10 cpu |
| C3' | `step2_cluster/3_0_eval_plate_v5v6.sh` (array 0-11) | the v5 vs v6 comparison that settled the cell set (v6) | 4 h, 120 GB, 10 cpu |
| C4 | `bash step3_compare/3_1_crossstage_plate.sh` -> `3_1_crossstage_plate.R` | per Pre cluster: QC signature, retention in WD and ND, Fisher enrichment, ARI | local |
| C5a | `step3_compare/4_3_eval_plate.sh` -> `common/4_0_build_gene_body_matrix.R`, `common/4_3_marker_accessibility.R` (array 0-2 At, 3-5 B73) | gene-body matrix from the plate BED (raw or `GENE_LENGTH_NORM=perkb`), pseudobulk marker z-scores and cluster annotation | 4 h, 48 GB, 4 cpu |
| C5b | `step3_compare/4_3f_smooth_plate.sh` -> `common/4_3f_smooth_gene_activity.R` | kNN-Markov smoothing of gene activity in the v7 PC space | 2 h, 64 GB, 4 cpu |
| C5c | `step3_compare/4_3g_annotate_plate.sh` -> `common/4_3g_cell_annotation.R` | per-cell enrichment classifier on the smoothed activity | 1 h, 48 GB, 4 cpu |
| C5d | `step3_compare/4_3i_reciprocal_percell_plate.sh` -> `common/4_3i_reciprocal_zscore_percell.R` (perkb) | per-cell reciprocal z (Zi across cells, Zj across genes) | 2 h, 64 GB, 4 cpu |
| C5e | `step3_compare/4_3n_bgnull_plate.sh` -> `common/4_3n_background_null_rZ.R` | expression-conditioned background-gene null for the type calls | 4 h, 64 GB, 4 cpu |
| C5f | `step3_compare/4_3cd_postprocess.sh` -> `4_3c_annotation_diff.R`, `4_3d_marker_informativeness.R`; `4_3e_build_informative_panel.py` | cross-stage annotation diff, empirical marker informativeness, the informative_top15 panels | 30 min, 8 GB, 1 cpu |
| C6a | `step5_metacell/5_0_make_consensus_peaks.sh` -> `common/5_0_make_consensus_peaks.R` (array 0-1) | one frozen consensus ACR set per genome (Pre-primary; union emitted for the reverse arm) | 1 h, 16 GB, 2 cpu |
| C6b | `step5_metacell/5_1_build_acr_matrix_plate.sh` -> `common/5_1_build_acr_matrix.R`, `5_1b_export_for_seacells.R` (array 0-1) | ACR x cell matrices with the feature set frozen on Pre, plus the portable SEACells bundle | 4 h, 64 GB, 4 cpu |
| C6c | `step5_metacell/5_3_run_seacells_plate.sh` -> `common/5_3_run_seacells.py` (array 0-5; `SEEDS=100 OUT_TAG=_S100 MEMBERSHIP_ONLY=1` for the replicate sweep) | SEACells meta-cells; `RAREFY_TO` and `N_SEACELLS` must take one value across the three stages of a genome | 8 h, 96 GB, 8 cpu |
| C6d | `step5_metacell/5_5_metacell_rZ_plate.sh` -> `common/4_3p_metacell_rZ.R` (array 0-5) | marker reciprocal z at meta-cell level (raw perkb, never the smoothed matrix) | 2 h, 64 GB, 4 cpu |
| C6e | `step5_metacell/5_9_consensus_plate.sh` -> `common/4_3u_metacell_consensus.R` (array 0-5) | consensus partition from the 100-seed sweep, k pinned per genome, acceptance gate, no threshold anywhere | 2 h, 32 GB, 4 cpu (runs locally) |
| C6f | `bash step5_metacell/5_10_seedpair_ci_plate.sh At|maize` -> `common/4_3y_seedpair_ami_ci.R` | seed-pair ARI / AMI with run-resampling bootstrap CIs | local |
| C6g | `step5_metacell/5_12_consensus_rZ_plate.sh` -> `common/4_3p_metacell_rZ.R`, `4_3z_consensus_rZ_join.R` (array 0-5) | rZ on the consensus meta-cells joined with their confidence attributes -> `rZ_annotation/consensus/*.heatmap_input.tsv` (Fig. 5I) | 2 h, 32 GB, 4 cpu |
| C6h | `Rscript common/4_3zc_seed_label_sweep.R <s5_dir> At|maize` | label stability of every consensus meta-cell across the 100 seed partitions | local |

`common/4_3h_reciprocal_zscore.R` (cluster-level reciprocal z) and
`common/4_3k_cluster_annotation_from_rZ.R` (independent cluster annotation from the per-cell rZ
tables) have no driver; they were run locally on the `4_3_eval_plate` and `4_3i` outputs.

Marker resource preparation (`_data/markers/`): `common/4_1a_extract_at_markers.py` and
`4_1b_curate_at_markers.py` build the Arabidopsis panel from the Ecker et al. 2025 atlas
supplementary table, `4_1c_merge_db_markers.py` merges the PlantscRNAdb 4.0 high-confidence
markers into the curated At (Ecker) and maize (Marand) panels with a v4 to v5 gene lift, and
`4_3e_build_informative_panel.py` prunes that augmented panel to the `informative_top15` sets the
annotation steps score over.

## Reading the statistics

Comments kept verbatim in the engines define the quantities the manuscript quotes:
the observed / expected mixing (`2_1_gridscan_pergenome.R`, `3_2_coembed_gridscan.R`: raw mixing
is composition-confounded through the 2p(1-p) baseline, so obs/exp is reported and mixing is
never part of the selection score); the frozen configuration held fixed across stages for
comparability rather than re-optimised (the Fig. 5 co-embedding is not the score maximiser);
the statistical unit of the SEACells sweep (runs are the replication unit, seed pairs are not
independent, pooled seeds are pseudo-replicates, cross-stage comparisons are descriptive); and
the annotation policy (meta-cell labels are descriptive top types, not tested identities).

## Not ported

`setup/0_02`, `0_03` (B73Mo17 staging, no panel), `setup/0_04` to `0_06` (decontamination and
mapping stages), `part1_combined/2_0_1_Metadata_update_and_QC.v2.R`, `2_0_2_umap_grid_scan.{R,sh}`,
`2_0_2_umap_grid_scan_panelA_replot.R` and `part1_indep/3_4_coembed_ablation.{R,sh}` (analysis
scripts), the combined-genome clustering and its evaluations (`3_0_0_*`, `3_0_eval_*`, `4_2_*`),
the indep-arm drivers superseded by `part1_indep/` (`step0_qc/0_08`, `1_QC_*_indep`, `1_2_rerun_*`,
`step0b_metaqc/`, `step1_cluster_opt/2_0_2b`, `2_0_2c`, `2_0_3_*`, `3_0_2` to `3_0_5`,
`step2_cluster/3_0_eval_pergenome.sh`, `step3_compare/3_1_cluster_crossstage_compare.R`,
`3_1_eval_crossstage.sh`, `4_3_eval_pergenome.sh`), the SVD-null and gene-null layers
(`5_7*`, `5_8`, `5_11`, `5_13`, `5_14`, `common/4_3t`, `4_3za`), the meta-cell diagnostics
(`5_4a`, `5_4b`, `common/4_3j`, `4_3l`, `4_3m`, `4_3o`, `4_3q`, `4_3v`, `4_3w`, `4_3x`, `4_3zb`)
and the feature-filter sweeps (`common/5_1c`, `5_1s`, `5_2s*`).
