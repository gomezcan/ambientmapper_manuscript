# scifiATAC_B73_Arabidopsis

Files the analysis scripts read, by subfolder. Producing stage in brackets.

`SM2/decontam_with_design_alpha05_v2/` [04_decontamination/sm2]: `SM2_barcode_policy.tsv.gz`,
`SM2_pre_barcode_genome_counts.tsv.gz`, `SM2_pre_barcode_composition.tsv.gz`,
`SM2_cells_calls.decontam.tsv.gz`, `SM2_barcode_postclean.tsv.gz`. Read by `analysis/fig1_*/fig1.R` (Fig. 1B to E).

`SM2v2/decontam_with_design_alpha05_v2/` [04_decontamination/sm2v2]: `SM2v2_barcode_policy`,
`SM2v2_pre_barcode_genome_counts`, `SM2v2_post_barcode_genome_counts`, `SM2v2_pre_barcode_composition`,
`SM2v2_post_barcode_composition`, `SM2v2_cells_calls.decontam`, `SM2v2_barcode_postclean` (all `.tsv.gz`),
`SM2v2_decontam_params.json`. Read by Fig. 2, Fig. 3, Tables S1 and S3.
`SM2v2/decontam_without_design_alpha05_v2/`: `SM2v2_decontam_params.json` (Table S1).

`socrates/SM2v2_indep/` [05, independent co-embedding]: `coembed/{Pre_cluster,Post_wd_cluster}/*.updated_metadata_v7.*.txt`
and `*.reduced_dimensions_v7.*.txt`, `coembed/{Pre_grid,Post_wd_grid}/UMAP_grid_scan.minc_50.metrics.tsv`,
`coembed/SM2v2_coembed_Pre.coembed.meta_full.tsv`, `step0_qc/SM2_{B73v5,TAIR10}.raw.soc.rds`. Fig. 1F to I, Fig. 5A to D.

`socrates/SM2/` and `socrates/SM2v2_clean/` [05, combined genome]: `step1_integrate/SM2.full.{SocObj.rds,metadata_updated.txt}`,
`step2_metaqc/*.updated_metadata_v4.txt`, `step3_embedscan/UMAP_grid_scan.metrics.tsv`; `socrates/compare/SM2_{B73,At}.fixedSet.minDepth200.pre_post.txt`;
`socrates/_data/_PeakFiles/*/*_peaks.narrowPeak`; `socrates/_data/_BED_files/{SM2,Clean.SM2v2}_{At,B73}_ZmATcombined_*.tn5.bed.gz`. Fig. S2, S3, Fig. 5E, F.

`socrates/SM2v2_plate/` [05, plate arm]: `step2_cluster/*.updated_metadata_v7.*.txt`, `step0_qc/*.updated_metadata_v1.txt`
and `*_macs2_temp/*_peaks_combined_peaks.narrowPeak`, `step3_compare/*.plate.perkb.genes.sparse.rds`,
`step5_metacell/consensus/*`, `step5_metacell/rZ_annotation/consensus/*`; `socrates/_data/_BED_files/*_{At_TAIR10,B73_B73v5}_*.tn5.bed.gz`,
`socrates/_data/markers/*`, `socrates/_data/_GenomeInfo/{B73v5,TAIR10}.gff3` and chromosome sizes. Fig. 5G to L, S7, S8.

`bed/SM2_{At,B73}/SM2_{At,B73}_ZmATcombined_scifiATAC.mq10.tn5.bed.gz` [02_mapping, legacy route]: Fig. S1.
