# 02 mapping: demultiplexed FASTQs to per-genome cleaned BAMs and Tn5 insertion BEDs

Every library was aligned with BWA-MEM 0.7.17 (`-M`) to each reference genome of its panel separately
(multi-reference mapping), then cleaned to MAPQ >= 10 proper pairs, barcode-corrected, deduplicated
with Picard MarkDuplicates (`BARCODE_TAG=BC`) and filtered for multi-mapping reads (MAPQ < 30 with more
than one alternative hit within 3 edits). The route differs by library and is fixed by the BAM `@PG`
headers: legacy Perl chain (Picard 2.18.29, samtools 1.18) for SM2 and Root1, scifi-demux step 2
(Picard 3.4.0, samtools 1.22) for the Zhang et al. 2024 libraries.

Placeholders: `PROJECT_ROOT` and `REPO_ROOT` as in `01_preprocessing/`, plus `GENOMES_DIR` for the BWA
indexes (`Zea/NAN_Indexes/Index_Zm_<genome>_bwa`, `Arabidopsis/Index_AraTAIR10_bwa`).

## zhang2024/ (scifi-demux step 2)

| Order | Script | Does | Resources |
|---|---|---|---|
| 1 | `0_generate_run_plans.sh` | writes `run_plan.<dataset>.tsv` (one row per well x genome) from `config/Pools_scifi_*.txt` and the genome lists in this directory | seconds |
| 2 | `2_01_step2_map.sh <run_plan> <dataset>` | `scifi-demux step2 run --mode hpc --mapq-min 10`: BWA-MEM, sort, BC tagging and MAPQ filter, dedup, fixBC, Tn5 BED, one well x genome per array task (1 to 192 for B73Mo17, 1 to 672 for multiGenotypes) | 4 h, 8 CPU, 20 GB per task |
| 3 | `3_01_merge_well_bams.sh <dataset> <genome_list>` | `samtools merge` of the per-well BAMs into `3_Mapping/ambientmapper_input/<dataset>_<genome>_scifiATAC.mq10.BC.rmdup.mm.bam` | 4 h, 8 CPU, 16 GB |
| 4 | `cleanup_per_well_bams.sh <dataset>` | removes the per-well BAM intermediates once the merged BAMs exist | minutes |

Genome lists: `Genome_list_scifi_B73_Mo17` (B73v5, Mo17) and `Genome_list_MultipleGenome` (B73v5, B97,
M162W, Mo18W, Ky21, Oh7B, Tzi8). The shipped run plans carry the literal `${GENOMES_DIR}` in `ref_path`.

## legacy_sm2_root1/ (BWA + Perl chain)

| Order | Script | Does | Applies to |
|---|---|---|---|
| 1 | `1_1_align.scifi.ATAC.NAM.sh <lib>` | `bwa mem -M` of `<lib>_R{1,3}.bc1.fastq.gz` against `Index_Zm_<genome>_bwa`, one genome per array task read from `GENOME_LIST` (`Genome_list_NAM_{1,2,3}.txt` were run in turn for the 26 NAM genomes of the Root1 panel) | Root1 |
| 1 | `1_1_align.scifi.sh <lib>`, `1_1_align.scifi.at.sh <lib>` | `bwa mem -M` against the B73v5 index and against the TAIR10 index (single-reference wrappers, archived as used) | SM2, see note |
| 2 | `1_scifi_processBAM.ATAC.sh <lib>` | the cleaning chain, one raw BAM per array task: sort; `1_1_scifi_modufy_BC_flag.pl` (barcode from the read name into `BC:Z:`) and `samtools view -q 10 -f 3`; `1_2_scifi_countBCs.BAM.pl` (barcodes seen on more than 5 reads); `1_3_scifi_correctBCs.10x.v2.pl` (10x segment corrected within two substitutions to the closest whitelist entry, Tn5 halves within one); `1_4_scifi_correctBAM.pl`; Picard MarkDuplicates; `1_5_scifi_fixBC.pl` (multi-mapping filter, library tag appended to the barcode, per-barcode counts); `1_6_scifi_makeTn5bed.py` (Tn5 cut sites, +4 forward / -5 reverse) | SM2, Root1 |

`1_3_scifi_correctBCs.10x.v2.pl` reads `tn5_bcs.txt` (shipped) and `737K-cratac-v1.txt` (the 10x
Genomics scATAC barcode whitelist, not shipped) from its working directory; the driver `cd`s there.

Note on the SM2 alignment: which BWA wrapper produced the SM2 raw BAMs is not recorded.
`1_1_align.scifi.sh` and `1_1_align.scifi.at.sh` are the two archived single-reference wrappers (they name
their outputs `<lib>_scifiATAC.raw.bam` and `<lib>_At_scifiATAC.raw.bam`), while the SM2 BAMs on disk follow
the `<lib>_<genome>_scifiATAC.raw.bam` pattern of `1_1_align.scifi.ATAC.NAM.sh` with the genome list
`Genome_list_SMs.txt` (Zm_B73v5, AraTAIR10, ZmATcombined). All three wrappers ship. The `ZmATcombined`
(concatenated B73v5 + TAIR10) mapping feeds Fig S1 to S3 and Fig 5E, F only; the per-reference BAMs are
the AmbientMapper input.

## tn5bed/ (Tn5 insertion BEDs for the Socrates QC step)

| Order | Script | Does |
|---|---|---|
| 1 | `01_00_merge_SM2v2_inputs.sh` | `samtools merge` of the per-plate-half SM2 BAMs (SM2_B73 + SM2_At) into one BAM per reference, `3_Mapping/ambientmapper_input/SM2_{B73v5,TAIR10}_scifiATAC.mq10.BC.rmdup.mm.bam`, the AmbientMapper input |
| 2 | `00_bam_to_tn5bed_parallel.py --bam --out` | per-chromosome parallel Tn5 BED extractor (cut = start + 4 on the forward strand, end - 5 on the reverse); records deduplicated on the full (chrom, cut, barcode, strand) key and concatenated in BAM header order |
| 3 | `01_12_tn5bed_preclean_SM2v2.sh` | PreClean per-genome BEDs from the merged AmbientMapper input BAMs (array 0 = B73v5, 1 = TAIR10) |
| 4 | `01_11c_tn5bed_regen_SM2v2.sh` | PostClean per-genome BEDs from the cleaned BAMs of both decontamination passes (array 0 to 3 = ND/WD x B73v5/TAIR10), overwriting any existing BED |
| 5 | `0_06_make_tn5bed_SM2v2_combined.sh` | PostClean BEDs of the cleaned concatenated-reference BAMs (`Clean.SM2v2_{At,B73}_ZmATcombined_*`, written by `workflows/04_decontamination/combined_genome/`) |

Steps 3 and 4 verify that records are not collapsed to distinct (chrom, start) positions (a correct
Tn5 BED has many barcodes sharing insertion sites). Resources: 12 CPU, 24 GB, 2 to 4 h per task.
