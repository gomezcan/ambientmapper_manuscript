# 01 preprocessing: raw reads to demultiplexed, barcode-corrected FASTQs

Two routes were used. Which library went through which is settled by the `@PG` header lines of the
final BAMs (Picard MarkDuplicates 2.18.29 and samtools 1.18 in the legacy chain, Picard 3.4.0 and
samtools 1.22 in scifi-demux; BWA-MEM 0.7.17 in both):

| Libraries | Route | Scripts |
|---|---|---|
| SM2 (maize B73 + Arabidopsis, in-house) and Root1 (Marand et al. 2021, SRR12331466) | legacy chain: UMI-tools + cutadapt + custom Python demultiplexing | `legacy_sm2_root1/` |
| B73Mo17_rep1, B73Mo17_rep2, multiGenotypes_rep1 (Zhang et al. 2024, SRR25320539 to SRR25320547) | scifi-demux step 1 | `zhang2024/` |

Unresolved: the project's own preprocessing notes also describe a scifi-demux route for SM2 (library
id `scifi_At_B73_rep1_1`, design `config/PlateDesign_scifi_At_B73_rep1.txt`, and an archived scifi-demux
run plan for `SM2_ATAC`). The BAM headers of the SM2 files the manuscript reads say legacy chain, so that is
what is documented here; the scifi-demux pass over SM2 is not part of the paper.

Both routes implement the same demultiplexing rule: the 16 bp 10x cell barcode is transferred from the
index read (R2) to the read names of R1 and R3; the Tn5 well barcodes (5 bp + 5 bp) are read at error
rate 0.2 against the mosaic end and corrected with at most one mismatch per half (N counts as a mismatch)
to the closest well of `config/96well_Tn5_bc_layout.txt`; reads are assigned to samples or wells by the
plate design.

Placeholders: `PROJECT_ROOT` (analysis tree with `1_RawData/`, `2_CleanReads/`) must be exported;
`REPO_ROOT` (this repository) defaults to `git rev-parse --show-toplevel`. SLURM headers keep the
resources but not the account or e-mail lines; `# conda activate <env from environment.yml>` marks where
the environment was activated.

## zhang2024/ (scifi-demux)

| Order | Script | Does | Resources |
|---|---|---|---|
| 1 | `1_prefecth_raw_reds.sh` | prefetch + fasterq-dump of the nine SRA runs listed in `config/scifi_Metadata_sra.only.txt`, one per array task | 12 h, 20 CPU, 10 GB |
| 2 | `rename_sample.sh` | rename `<SRR>_*.fastq.gz` to the library ids of `config/scifi_Metadata_sra.clean.txt` | seconds |
| 3 | `1_01_demux_step1.sh` | `scifi-demux step1 plan`: 100 chunks for one run (`LIB`, `DESIGN` environment variables) | 3 h, 25 CPU, 10 GB |
| 4 | `1_02_demux_step1.sh` | `scifi-demux step1 run --mode hpc`, array 1 to 100 (same `LIB`, `DESIGN`) | 3 h x 100 tasks, 20 CPU, 10 GB |
| 5 | `3_01_merge_fastq_wells.sh` | phase 1 `scifi-demux step1 merge` per run; phase 2 gzip concatenation of the three runs of each replicate into `2_CleanReads/combined/<dataset>_rep<N>/<well>_R{1,3}.bc1.bc2.fastq.gz` | 8 h, 4 CPU, 8 GB |

Steps 3 and 4 run once per sequencing run (nine times) with `LIB=scifi_<dataset>_rep<N>_<run>` and the
matching `config/PlateDesign_scifi_*.txt` (1:1 well format, 96 wells).

## legacy_sm2_root1/ (UMI-tools, cutadapt, custom demultiplexing)

| Order | Script | Does | Applies to | Resources |
|---|---|---|---|---|
| 1 | `1_1_UMItools.ATAC.R1_parallel.sh`, `1_1_UMItools.ATAC.R3_parallel.sh` | `umi_tools extract`: 16 bp barcode from R2 into the R1 and R3 read names (`<lib>_R{1,3}.bc1.fastq.gz`), 20 seqkit chunks in parallel | SM2, Root1 | 12 to 24 h, 30 CPU, 80 GB |
| 2 | `1_2_Cutadapt.ATAC.sh` | cutadapt: 5 bp Tn5 barcodes cut from R1 and R3 (error rate 0.2, mosaic end `AGATGTGTATAAGAGACAG`) and appended to the read name (`.bc1.bc2`) | SM2 only | 4 h, 4 CPU, 80 GB |
| 3 | `2_1_0_split_chunks.ATAC.sh <lib>` | seqkit split into 60 chunk pairs | SM2 only | 2 h, 10 CPU, 20 GB |
| 4 | `2_1_1_fastq_sample_assign_and_fixTn5.ATAC.sh <lib>` (+ `.py`) | per chunk: Tn5 barcode correction (one mismatch per half, N = mismatch, closest well), assignment to `SM2_B73` / `SM2_MUDR` / `SM2_At` by `config/PlateDesign_SM2_ATAC.legacy_demux.txt`, corrected barcode written into the read name | SM2 only | array 1 to 120, 3 h, 1 CPU, 20 GB |
| 5 | `3_merge_fastq.ATAC.sh <lib>` | concatenate the chunks per sample: `SM2_{B73,At,MUDR}_R{1,3}.bc1.bc2.fastq.gz` | SM2 only | 40 min, 4 CPU, 5 GB |

Root1 is a 10x scATAC library without Tn5 barcodes: its `Root1_rep1_R{1,3}.bc1.fastq.gz` from step 1 go
straight to mapping (`workflows/02_mapping/legacy_sm2_root1/`). MuDR wells are demultiplexed into
`SM2_MUDR_*` and never mapped.

Notes on the port: the archived scripts carried the sample arrays of the runs they were last used for
(other libraries); they are set to the paper's libraries here (`SM2_ATAC`, `Root1_rep1`), the raw SM2
FASTQs being `1_RawData/SMs/SM2_ATAC_R{1,2,3}.fastq.gz`. The scifi-demux merge script also listed an
`At_B73` run next to the Zhang libraries; that entry was dropped (see the unresolved note above).
