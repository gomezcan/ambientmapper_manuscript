# Configuration files

Parameter records of the runs behind the manuscript. Placeholders: `${PROJECT_ROOT}` is the analysis
tree that holds `1_RawData/`, `2_CleanReads/`, `3_Mapping/` and `5_AmbientDetection/`; `${GENOMES_DIR}`
holds the BWA indexes (`Zea/NAN_Indexes/Index_Zm_<genome>_bwa`, `Arabidopsis/Index_AraTAIR10_bwa`).

## AmbientMapper

| File | Dataset | What it records |
|---|---|---|
| `ambientmapper_scifiATAC_SM2v2.yaml` | SM2v2 (maize B73 + Arabidopsis) | The shipped SM2v2 run: assign and genotyping flags read from the run log and from the shared wrapper (`workflows/03_genotyping/_genotyping_configs.sh`). Behind Fig 1F to I, 2, 3, 5, S4, Tables S1 and S3. |
| `ambientmapper_scifiATAC.yaml` | SM2 v1 | The earlier run read by Fig 1B to E only. |
| `ambientmapper_marand2021.yaml` | Root1_rep1 (26 NAM genomes) | Sample config of the zero-contamination control (Fig 4D to G, S6). |
| `SM2v2.ambientmapper.json`, `Root1_rep1.ambientmapper.json`, `B73Mo17_rep1.ambientmapper.json`, `B73Mo17_rep2.ambientmapper.json`, `multiGenotypes_rep1.ambientmapper.json` | as named | The `--config` JSON files passed to `ambientmapper` (sample, workdir, min_barcode_freq, chunk_size_cells, genome to BAM), absolute paths replaced by `${PROJECT_ROOT}`. |
| `genome_map_scifiATAC.tsv` | subsampled sandbox | Format example of the scifi-demux step 2 genome map; not a paper run. The paper run plans are `workflows/02_mapping/zhang2024/run_plan.*.tsv`. |

The genotyping factorial configurations (C0 to S13) are defined in `workflows/03_genotyping/_genotyping_configs.sh`.
Decontamination flags are recorded per run in the `*_decontam_params.json` next to each run's outputs
and in the `workflows/04_decontamination/` scripts.

## Plate designs, pools and the Tn5 layout

| File | Format | Dataset |
|---|---|---|
| `PlateDesign_SM2_ATAC.txt` | range (`{sample}\t{well_ranges}`), 2 groups: `B73` columns 1 to 4, `At` columns 9 to 12 | SM2 design of record: the plate design passed to the AmbientMapper WD decontamination pass and the Table S1 input. Carries no MuDR block (the Table S1 script rejects the 3-group files below). |
| `PlateDesign_SM2_ATAC.legacy_demux.txt` | range, 3 groups (`SM2_B73`, `SM2_MUDR`, `SM2_At`) | SM2 legacy demultiplexing design; the sample names become the `SM2_B73` and `SM2_At` library prefixes used downstream. MuDR wells (columns 5 to 8) are separated at demultiplexing and never mapped. |
| `PlateDesign_scifi_At_B73_rep1.txt` | range, 3 groups (`At_B73_B73`, `At_B73_MUDR`, `At_B73_At`) | The same SM2 plate in scifi-demux naming; see the provenance note in `workflows/01_preprocessing/README.md`. |
| `PlateDesign_scifi_B73Mo17_rep1.txt`, `PlateDesign_scifi_B73Mo17_rep2.txt`, `PlateDesign_scifi_multi_genotypes.txt` | 1:1 (`{well}\t{well_id}`, 96 rows) | Zhang et al. 2024 libraries, per-well demultiplexing (scifi-demux step 1). |
| `Well_to_Genotype_multiGenotypes_rep1.txt` | `{well}\t{genotype}` | Well to genotype assignment of the seven-genotype library (Table S1 input). |
| `Pools_scifi_*.txt` | one name per line | Well (96) or sample (3) names consumed by the mapping run plans. |
| `96well_Tn5_bc_layout.txt` | 8 x 12 table of `<Tn5 A>_<Tn5 B>` 5-mers | Tn5 well-barcode layout used by both demultiplexing routes. |

## Barcode lists (`../data/metadata/barcode_lists/`)

`sample, genome, bam, workdir` tables of the first AmbientMapper runs (`SM2_AtB73.list.tsv` = SM2 v1,
`Root1_Rep1.list.tsv` = Root1). The absolute HPC prefixes of the originals were replaced by
`${PROJECT_ROOT}/...`, and `4_MappingCleaning/` by `3_Mapping/`: the BAM-cleaning directory was merged
into `3_Mapping/` during the project, so `4_MappingCleaning/ambientmapper_input/` is now
`3_Mapping/ambientmapper_input/` and `4_MappingCleaning/Root1_rep1/` is `3_Mapping/Root1_rep1/`.
`AtB73_B73.list.tsv` is an early configuration pointing at pooled BAMs (`Pools_bams/`, a directory
that no longer exists); no panel reads it.

## SRA metadata

`scifi_Metadata_sra.clean.txt` (SampleID, Run, LibraryLayout) and `scifi_Metadata_sra.only.txt` (the
nine SRR ids) drive the download and renaming scripts in `workflows/01_preprocessing/zhang2024/` and the
Table S2 script; `../data/metadata/sra_metadata.txt` lists the SRA download URLs. `srr_root.txt` is the
run list of the maize root atlas as downloaded (four runs across two GEO samples); the paper uses the
single run SRR12331466 (Root1), and the Table S2 script reads the file with that restriction.

## Parameter glossary (AmbientMapper genotyping)

- `bic_margin`: BIC margin for the singlet vs doublet decision
- `empty.bic_margin`, `empty.top1_max`, `empty.ratio12_max`: empty-barcode gate
- `topk_genomes`: number of top genomes kept per barcode for read re-classification (`topk_reclass`)
- `single_mass_min`, `ratio_top1_top2_min`, `doublet_minor_min`: singlet and doublet call thresholds
- `max_alpha`, `alpha_grid`, `rho_grid`: contamination-fraction grid of the mixture fit
- `winner_only`, `beta`: winner-only scoring (off in every run) and Dirichlet pseudocount
- `winner_discount`: winner-mass discount of ambiguous reads (`winner_ratio` mode)
- `xmap`: cross-mapping profile phi; reported descriptively in Fig 4D, not applied in the shipped SM2v2 run
