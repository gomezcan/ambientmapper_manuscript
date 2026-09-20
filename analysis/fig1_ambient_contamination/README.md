# Figure 1: ambient contamination in the interspecies barnyard

Panel A is a schematic drawn in Illustrator and has no script. All scripts run from the
repository root and source `analysis/_helpers/plotting.R`.

| Script | Panels | Inputs (repo paths, under `data/processed/scifiATAC_B73_Arabidopsis/`) | Run |
|---|---|---|---|
| `fig1.R` | B, C, D, E | `SM2/decontam_with_design_alpha05_v2/SM2_{barcode_policy, pre_barcode_genome_counts, pre_barcode_composition, cells_calls.decontam, barcode_postclean}.tsv.gz` | `Rscript analysis/fig1_ambient_contamination/fig1.R` |
| `fig1.R` | F, G, H | `socrates/SM2v2_indep/coembed/Pre_cluster/SM2v2_coembed_Pre.updated_metadata_v7.pcs_20.k_near_30.min_dis_0.3.minc_50.res_0.5.txt` | same call |
| `fig1I_coembed_ablation.R` + `.sh` | I | `socrates/SM2v2_indep/coembed/SM2v2_coembed_Pre.coembed.meta_full.tsv`, `socrates/SM2v2_indep/step0_qc/SM2_{B73v5,TAIR10}.raw.soc.rds`, optional unablated comparator = the F to H metadata above | `mkdir -p _logs && sbatch analysis/fig1_ambient_contamination/fig1I_coembed_ablation.sh` |

## Outputs

- `figures/main/fig1/Fig1{B,C,D,E,F,G,H}_*.{pdf,png}` plus `Fig1BtoH_preview.png` (composite preview; the manuscript figure is assembled by hand).
- `figures/main/fig1/coembed_ablation/`: the re-embedded Socrates object, per-nucleus metadata (`*.ablated.updated_metadata_v7.*.txt`), the mixing summary (`*.ablated.mixing_summary.*.tsv`), the per-plate embeddability ladder (`*.ablated.embeddability_ladder.*.tsv`) and diagnostic UMAP renders in `plots/`. Panel I was assembled from the ablated metadata and the `_Genome` render.

## Where the inputs come from

- Panels B to E read the design-aware AmbientMapper decontam run in `SM2/` (`workflows/04_decontamination/sm2/`, preceded by `workflows/03_genotyping/sm2/`).
- Panels F to I read Socrates objects built from the `SM2v2/` run by the independent (multi-reference) co-embedding chain in `workflows/05_qc_and_embedding/part1_indep/` (QC objects from stage 1, the co-embed barcode set from stage 3.1, the clustered metadata from stage 3.3). The concatenated-reference co-embedding is used only by Fig S1.

## Notes

- The co-embedding configuration `pcs_20.k_near_30.min_dis_0.3.minc_50.res_0.5` is held fixed for comparability with Fig S1, it is not the score maximiser. Fig 1I re-embeds with the identical configuration so that the ablation is the only difference between panels F and I.
- Clustering in G and H is Leiden (`Socrates::callClusters(cl.method = 4)` forwards to `Seurat::FindClusters(algorithm = 4)`). `LouvainClusters` is only the column name Socrates writes.
- Panel E plots the count-based pseudocounted dominance ratio `(k1 + 1) / (k2 + 1)` from winner reads, not AmbientMapper's `ratio_top1_top2`; the dashed line is the detection ceiling `ratio = reads + 1`. Panel D counts barcodes with `total_reads > 10`, panel C uses `total_reads > 0`.
- Panels F to H show nuclei with `total >= 500`, `pTSS >= 0.2`, `FRiP >= 0.2` in clusters of at least 100 cells.
- Fig 1I is HPC only (R with Socrates, Seurat and FNN; about 180 GB peak while the raw B73v5 object is loaded). AmbientMapper is not in its pipeline: the inputs are the raw pre-QC per-genome matrices and the only mask is the plate index. Its mixing statistic is obs/exp against the random-intermingling baseline `2p(1-p)`; report obs/exp together with the embeddability ladder, never the raw value.
