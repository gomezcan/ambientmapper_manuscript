# fig5_biological_impact

Fig. 5, the biological effect of cleaning on the maize and Arabidopsis library. Run from the
repository root. Inputs live under `data/processed/scifiATAC_B73_Arabidopsis/socrates/` (layout in
`data/processed/scifiATAC_B73_Arabidopsis/README.md`); outputs go to `figures/main/fig5/`. Scripts
G to L and Fig. S8 source `analysis/_helpers/fig5_part2_helpers.R`.

| Script | Panels | Inputs | Notes |
|---|---|---|---|
| `fig5_AtoD_coembed.R` | A, B, C, D | `SM2v2_indep/coembed/{Pre_cluster,Post_wd_cluster}/*` metadata and reduced dimensions; `{Pre_grid,Post_wd_grid}/UMAP_grid_scan.minc_50.metrics.tsv` | mixing is reported as observed over expected under random intermingling, never raw; the embedding configuration was held fixed for comparability and is not the score maximiser |
| `fig5_EF_qc.R` | E | `compare/SM2_{B73,At}.fixedSet.minDepth200.pre_post.txt`; caches `frip_rarefaction.tsv`, `frip_full_points.tsv` in the output directory | stops with a message if the two caches are absent |
| `fig5_EF_prefilter_beds.sh` then `fig5_EF_qc_fripfair.R` | F | combined-genome Tn5 BEDs `_data/_BED_files/{SM2,Clean.SM2v2}_{At,B73}_ZmATcombined*`, peaks in `_data/_PeakFiles/`, the paired cell lists from `compare/` | the shell step writes a scratch directory that the R step reads and deletes |
| `fig5_GH_umap.R`, `fig5_GH_flow.R` | G, H | `SM2v2_plate/step2_cluster/*.updated_metadata_v7.*.txt` (six files, three stages by two genomes) | stage tokens `SM2` = PreClean, `Clean.SM2v2wd` = WD, `Clean.SM2v2` = ND |
| `fig5_I_consensus.R` | I | `SM2v2_plate/step5_metacell/consensus/*.consensus.{summary,ksweep,cells,metacells}.tsv`, `*.consensus.F.rds` | |
| `fig5_J_access_cache.R` | cache for J and K | `SM2v2_plate/step3_compare/*.plate.perkb.genes.sparse.rds`, `consensus/*.consensus.cells.tsv`, `_data/markers/*` | writes `Fig5_P3D_marker_access_pooled.tsv` and `_genes.tsv`; run before `fig5_J_typeaccess.R` and `fig5_K_examples.R` |
| `fig5_J_annotation_table.R` | table for K | `step5_metacell/rZ_annotation/consensus/*` | writes `Fig5_P3B_annotation_metacells.tsv`; run before `fig5_K_examples.R` |
| `fig5_J_typeaccess.R` | J | rZ tables and the access cache | |
| `fig5_K_examples.R` | K | rZ tables, perkb matrices, cluster metadata, the two caches above | the example genes were selected for the effect and are illustrative, not evidence; the aggregate panels carry the claim |
| `browser/` | L | `consensus/*.consensus.cells.tsv`, rZ type tables, `Fig5_P3E_examples_{selection,contrast}.tsv` from `fig5_K_examples.R`, plate-split BEDs, `_GenomeInfo/{B73v5,TAIR10}.gff3` and chromosome sizes | HPC only |

## Run order

1. `fig5_AtoD_coembed.R`
2. `fig5_EF_qc.R`; `fig5_EF_prefilter_beds.sh <scratch>` then `fig5_EF_qc_fripfair.R <scratch>`
3. `fig5_GH_umap.R`, `fig5_GH_flow.R`
4. `fig5_I_consensus.R`
5. `fig5_J_access_cache.R`, `fig5_J_annotation_table.R`, then `fig5_J_typeaccess.R`
6. `fig5_K_examples.R`
7. `browser/env.sh` (creates the conda environment), `groups.R`, `split.sh`, `bw.sh`, `loci.R`, `plot.sh`

## HPC only

The browser chain (panel L) runs as a SLURM array and needs bedtools, wigToBigWig and
pyGenomeTracks 3.9 in its own environment, plus the GFF3 annotations of both genomes. Steps 1 to 6
run on a workstation.
