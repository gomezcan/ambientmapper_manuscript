# 03b — Variant-based comparison (Souporcell, WASP, barcode purity)

The SNP-based reference against which AmbientMapper's calls are compared on the Zhang et al.
2024 libraries (B73Mo17_rep1, B73Mo17_rep2, multiGenotypes_rep1): supervised Souporcell
genotyping (Fig. 4H to K), WASP correction of reference-mapping bias, and the
barcode-resolved allele-purity model on raw versus cleaned BAMs (Fig. 4L to N, Table S4).
Ported from `5_Genotyping/0_scripts/` of the analysis project, plus the three
`06_44_*` scripts whose preparation outputs the WASP chain depends on.

## Conventions

- All scripts run from `${PROJECT_ROOT}/5_Genotyping/` and write to
  `5_Genotyping/<SAMPLE>/...`. They read AmbientMapper outputs from
  `${PROJECT_ROOT}/5_AmbientDetection/<SAMPLE>/` and BAMs from
  `${PROJECT_ROOT}/3_Mapping/ambientmapper_input/`. Export `PROJECT_ROOT` before `sbatch`.
- Sibling scripts are located through `REPO_ROOT="$(git rev-parse --show-toplevel)"`, so
  submit from inside this checkout.
- Environment variables replacing site paths: `GENOMES_ROOT` (holds
  `Zea/Zm_B73_REFERENCE_NAM_5.0/`), `BWA_INDEX_ROOT` (holds `Zea/NAN_Indexes/`) or `BWA_IDX`
  directly, `LOCALINSTALL` (default `$HOME/LocalInstall`; holds `souporcell_release.sif` and
  the `WASP` clone), `WASP_DIR`, `BWA_MODULE` (default `bwa/0.7.17-mil4ns7`, the bwa 0.7.17
  that mapped the BAMs).
- Souporcell is **v2.1** (`souporcell_release.sif`, supervised mode with `--known_genotypes`,
  `--skip_remap`, `--no_umi`, barcodes with >= 500 fragments). WASP is the van de Geijn et
  al. 2015 pipeline cloned from `github.com/bmvdgeijn/WASP` by `06_46 STEP=setup`.

## Configs shipped (`configs/`)

| File | Content |
|---|---|
| `Pools_DB.by_Sample.txt` | genotype -> sample (no header); gives Souporcell its `k` and known-genotype names |
| `Well_to_Genotype_multiGenotypes_rep1.txt` | well -> genotype -> genotype_rep for the 96-well multiGenotypes plate; plate-of-origin truth for `06_46`, `06_47` and `06_48` |

Not shipped (data, named here with their provenance):

- `configs/Final_Mo17_relative_to_B73.vcf` (+ `.vcf.gz`, `.tbi` from `06_00`): B73 and Mo17
  genotypes relative to B73 NAM 5.0, used for both B73Mo17 replicates.
- `configs/25NAM_full.vcf.gz` (+ `.tbi`): the 25 NAM founder genotypes relative to B73
  NAM 5.0, filtered by `06_01` to the seven pooled genotypes of multiGenotypes_rep1.
  Both VCFs were taken from the group's NAM-founder SNP set (`25NAM_snps`); the public
  source is to be named in the Data Availability statement.
- `Genotype_ID_key.v2.txt` is an output of `06_03` (written next to `clusters.tsv` under
  `<SAMPLE>/souporcell/supervised/<SAMPLE>.min500/`), not a config; Fig. 4H to K read both.

## Execution order (per SAMPLE)

| # | Script | Resources | Produces |
|---|---|---|---|
| 1 | `06_00_prep_mo17_vcf.sh` | 30 min, 1 cpu, 4 G | `configs/Final_Mo17_relative_to_B73.vcf.gz` + `.tbi` (read by steps 4 to 6) |
| 2 | `06_01_souporcell.sh <SAMPLE> [500]` | 48 h, 14 cpu, 160 G | `<SAMPLE>/souporcell/supervised/<SAMPLE>.min500/{clusters.tsv, cluster_genotypes.vcf, ...}`; the raw-BAM run (`round1_raw`) is the Fig. 4H to K input |
| 3 | `06_03_assign_genotype_by_pearson.R` | called by step 2 | `ref_clust_pearson_correlations.v2.tsv`, `Genotype_ID_key.v2.txt` |
| 4 | `06_44_run.sh` (`STEP=extract MODE=raw`, `STEP=extract MODE=clean`, then `STEP=analyze`; calls `06_44_extract_allele_counts.py` and `06_44_concordance.R`) | 36 h, 16 cpu, 120 G | its `prep()` builds the marker-site panel (`diagnostics/06_44_concordance/multi_7geno_sites.vcf.gz` for multiGenotypes_rep1, `genotype_panel.tsv`) and the barcode-to-genome table (`bc_to_genome1.tsv`, B73Mo17 replicates) that the WASP steps below read; the allele-concordance panel it computes itself is not in the manuscript |
| 5 | `06_46_wasp_pilot.sh` (`STEP=setup` on a login node; `STEP=wasp MODE=raw|clean`; `STEP=analyze`) | 48 h, 16 cpu, 72 G | per-chrom SNP files `diagnostics/06_46_wasp/wasp_snps/` (needed by step 5) and, for B73Mo17_rep2, the monolithic WASP BAMs `<SAMPLE>_{raw,clean}.wasp.bam` |
| 6 | `06_47_wasp_chunked_prep.sh` -> `06_47_wasp_chunked_array.sh` (array 1-10) -> `06_47_wasp_chunked_combine.sh`, each with `MODE=raw` and `MODE=clean`, then `06_47_wasp_chunked_analyze.sh` | 2 h/4 cpu/16 G; 12 h/16 cpu/96 G per chrom; 4 h/8 cpu/32 G; 2 h/4 cpu/32 G | chunked WASP BAMs `diagnostics/06_47_wasp_chunked/<SAMPLE>_{raw,clean}.wasp.bam` for B73Mo17_rep1 and multiGenotypes_rep1 (the monolithic run exceeded 72 G on multiGenotypes) |
| 7 | `06_48_run.sh` (`STEP=extract MODE=raw`, `STEP=extract MODE=clean`, then `STEP=analyze`) | 4 h, 16 cpu, 120 G | `diagnostics/06_48_barcode_purity/<SAMPLE>_barcode_purity.tsv.gz` (Fig. 4L to N) and `<SAMPLE>_weak_doublet_diag.tsv` (Table S4), via `06_48_extract_barcode_block_counts.py` and `06_48_purity_model.R` |

The "clean" arm of steps 4 to 7 reads the AmbientMapper-cleaned BAMs
`5_AmbientDetection/<SAMPLE>/clean_bams_alpha05_C0_nd/*.Clean.bam` produced by
`workflows/04_decontamination/zhang2024/{04_07,05_04}`; the "raw" arm reads the stage 02 BAMs.
`06_48` reads the barcode attributes (call class, top1, top2, depth) from
`5_AmbientDetection/<SAMPLE>/genotyping_runs/4cfg_2026-05-01/C0/<SAMPLE>_cells_calls.tsv.gz`.

## Note on the `06_44_*` scripts

`06_44_run.sh`, `06_44_extract_allele_counts.py` and `06_44_concordance.R` are exploratory: the
pooled allele-concordance panel they compute is not in the paper. They ship here anyway because the WASP chain depends on what
`06_44_run.sh prep()` builds: the marker-site panel (`multi_7geno_sites.vcf.gz`, the 25NAM VCF
restricted to the seven pooled genotypes, biallelic SNPs with minor allele frequency >= 0.01;
`genotype_panel.tsv`) and the barcode-to-genome table (`bc_to_genome1.tsv`, AmbientMapper
`genome_1` call per barcode with >= 200 reads). `06_46 STEP=wasp` and `06_47 combine` also
re-count alleles with the `06_44` extractor, and the `analyze` steps of `06_46`/`06_47` run
`06_44_concordance.R`; those outputs feed no manuscript panel.
