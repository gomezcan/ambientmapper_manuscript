# Supplementary figures

Scripts run from the repository root. Inputs are read from `data/processed/<dataset>/...` (each
`data/processed/*/README.md` lists where the files come from); outputs go to
`figures/supplementary/figS<N>/`. `...` below stands for `data/processed/scifiATAC_B73_Arabidopsis`.

| Script | Panels | Inputs (repo paths) | Run | Notes |
|---|---|---|---|---|
| `figS1.R` | S1 a, b | `.../bed/SM2_{At,B73}/SM2_{At,B73}_ZmATcombined_scifiATAC.mq10.tn5.bed.gz` | `Rscript analysis/supplementary/figS1.R` | Concatenated-reference (B73v5 + TAIR10) Tn5 BEDs, one per plate half, from `3_Mapping/SM2_{At,B73}/bed/`; not in the public deposit. Species of a read = chromosome prefix (`At`, `Zm`). |
| `figS2.R` | S2 A to E | `.../socrates/SM2/step1_integrate/SM2.full.metadata_updated.txt` | `Rscript analysis/supplementary/figS2.R <meta> SM2 200 PreClean figures/supplementary/figS2` | A to C = `SM2.FRiP0.2.FULL`, `.FRiP0.2.At`, `.FRiP0.4.B73` `*.QC_FIGURES.pdf` (the FRiP tags are name-only); D and E = `SM2.minDepth200.stagePreClean.D_E.pdf`. Test annotations (Wilcoxon, KS) are computed but not drawn. Also writes the `updated_metadata_v1..v6` tables; the Full v4 table is the input of `figS3.R`. |
| `figS3.R` + `figS3.sh` | S3 B, C (and the metrics behind A) | `.../socrates/SM2/step1_integrate/SM2.full.SocObj.rds`, `.../socrates/SM2/step2_metaqc/SM2.FRiP0.2.FULL.minDepth200.stagePreClean.updated_metadata_v4.txt` | `sbatch analysis/supplementary/figS3.sh` | HPC only (Socrates, about 30 GB). B and C are pages of `plots/UMAP_grid_scan.Genome.pdf`; writes `UMAP_grid_scan.{metrics,best}.tsv`. The wrapper also runs `figS3A_replot.R`. |
| `figS3A_replot.R` | S3 A | `.../socrates/SM2/step3_embedscan/UMAP_grid_scan.metrics.tsv` (or the `figS3.R` output) | `Rscript analysis/supplementary/figS3A_replot.R` | Metrics only. The highlighted "used" configuration (pcs 20, k 30, min_dist 0.3) was held fixed for comparability across objects and is not the score maximiser. |
| none | S4 | | | Schematic, one panel, assembled by hand; no script. |
| `figS5.R` + `figS5.sh` | S5 | see script header | see script header | Synthetic benchmark. |
| `figS6.R` + `figS6.sh` | S6 | see script header | see script header | |
| `figS7.R`, `figS7_fripnorm.R`, `figS7_prep_beds.sh` | S7 A to F | `scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step0_qc/*.minDepth200.updated_metadata_v1.txt` and `*_macs2_temp/*_peaks_combined_peaks.narrowPeak` (six each), plate-split Tn5 BEDs for panel D | run `figS7_prep_beds.sh <scratch>` then `figS7_fripnorm.R <scratch>` to build the panel D caches, then `figS7.R` |
| `figS8.R` | S8 A to C | `scifiATAC_B73_Arabidopsis/socrates/SM2v2_plate/step5_metacell/consensus/*.consensus.seedpairs.tsv`, `*.ari_ami_ci.tsv`, `*.seedpair_ami.tsv`; `step5_metacell/rZ_annotation/consensus/*.seed_label_sweep.*.tsv` | sources `analysis/_helpers/fig5_part2_helpers.R`; values are descriptive, runs are the replication unit |
