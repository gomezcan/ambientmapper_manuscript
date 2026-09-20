# Supplementary tables

Generators for Supplementary Tables S1 to S4. Each script reads the pipeline's own configuration
and output files, self-checks the result against the verified values, and writes a `.txt` (TSV
plus caption) and a formatted `.xlsx` (needs `openpyxl`; skipped with a note otherwise) into
`figures/tables/`. All four are read-only with respect to pipeline data and take `--outdir` to
write elsewhere. Run them from the repo root (paths are resolved relative to the script, so any
working directory works).

| Table | Script | What it is | Inputs (repo paths) | Cited by |
|---|---|---|---|---|
| S1 | `make_TableS1_library_design.py` | Library and plate design of the scifi-ATAC-seq libraries: plate index, pools, wells, whether the plate index identifies the sample, whether the design was supplied to AmbientMapper | `config/PlateDesign_SM2_ATAC.txt`, `config/PlateDesign_scifi_B73Mo17_rep{1,2}.txt`, `config/Well_to_Genotype_multiGenotypes_rep1.txt`, and the `*_decontam_params.json` of the five decontamination runs under `data/processed/` (listed in `DECONTAM_PARAMS`) | Methods, "scifi-ATAC-seq library construction" |
| S2 | `make_TableS2_dataset_inventory.py` | Dataset inventory: assay, material, reference panel, wells, cleaning mode, role, source, PubMed ID, accession, figures | `config/{SM2v2,Root1_rep1,B73Mo17_rep1,B73Mo17_rep2,multiGenotypes_rep1}.ambientmapper.json`, `config/scifi_Metadata_sra.clean.txt`, `config/srr_root.txt`, the four plate designs above plus `config/PlateDesign_scifi_multi_genotypes.txt`, the decontamination run directories under `data/processed/` (existence only), `data/processed/synthetic/{synthetic,synthetic_disc}/` (the `alpha_*` datasets and `barcoded/templates.tsv`) | Methods, "scifi-ATAC-seq datasets" (twice) |
| S3 | `make_TableS3_cross_plate_hopping.py` | Cross-plate barcode hopping by plate of origin and read depth | `data/processed/scifiATAC_B73_Arabidopsis/SM2v2/decontam_with_design_alpha05_v2/{SM2v2_cells_calls.decontam.tsv.gz, SM2v2_pre_barcode_composition.tsv.gz}` | Results section 3 |
| S4 | `make_TableS4_weak_doublet_pooling.py` | Evidence for pooling weak-doublet barcodes with singlets in the allele-purity analysis: per-dataset counts, share of the singlet class affected, median dominant-genome share of the assigned pair | `data/processed/zhang2024/<s>/diagnostics/06_48_barcode_purity/<s>_weak_doublet_diag.tsv` for `B73Mo17_rep1`, `B73Mo17_rep2`, `multiGenotypes_rep1` (from the 06_48 step of `workflows/03b_variant_based_comparison/`) | Methods, "Allele-level validation of decontamination"; Fig. 4L |

## Run

```bash
python3 analysis/tables/make_TableS1_library_design.py
python3 analysis/tables/make_TableS2_dataset_inventory.py
python3 analysis/tables/make_TableS3_cross_plate_hopping.py
python3 analysis/tables/make_TableS4_weak_doublet_pooling.py

# optional: also check the cited maize-root SRA run against the read headers of the
# demultiplexed FASTQ that was mapped (about 100 GB, not part of this repository)
python3 analysis/tables/make_TableS2_dataset_inventory.py --root-fastq /path/to/Root1_rep1_R1.bc1.fastq.gz
```

Outputs: `figures/tables/TableS1_library_and_plate_design.{txt,xlsx}`,
`TableS2_dataset_inventory.{txt,xlsx}`, `TableS3_cross_plate_barcode_hopping.{txt,xlsx}`,
`TableS4_weak_doublet_pooling.{txt,xlsx}`.

## Presentation rules

- **No internal identifiers in a published table.** `SM2`, `SM2v2`, `B73Mo17_rep1`,
  `multiGenotypes_rep1` and `Root1` must not appear in a table body or caption. Use "Maize and
  Arabidopsis", "B73/Mo17, replicate 1", "B73/Mo17, replicate 2", "Multi-genotype", "Maize root".
- **Plate index positions are shown, well by well is not.** "Columns 1 to 4, A1 to H4" rather than a
  32-item list. Nothing outside the study is named or referred to.
- **Per-well demultiplexing is not a design.** The B73/Mo17 libraries were split per well only to
  parallelise mapping; the analysis treats the plate as one pooled sample, and the tables say so.

## Rules

1. **Never hand-edit the `.txt` or `.xlsx`.** Edit the script and re-run. The output files are
   build products.
2. **The legend in the manuscript and the `CAPTION` string in the script must stay identical.** The
   caption is written into both output files, so an edit in only one place puts three copies of
   the same legend out of sync. Change both in the same pass.
3. **The self-checks are load-bearing.** Each script hard-fails if the design files, reference
   panels, accessions, run directories or measured totals drift from the verified structure. A
   failure means something upstream changed, not that the check is wrong. Investigate before
   updating the expected values.

## Traps these scripts encode

- **`config/srr_root.txt` spans two GEO samples and cannot be cited wholesale.** SRR12331466 and
  SRR12331467 are GSM4696884 ("Root1"), SRR12331468 and SRR12331469 are GSM4696885 ("Root2"), all
  under PRJNA648930. **Only SRR12331466 was mapped and analysed** (the BAM header shows bwa was
  given the `Root1_rep1` reads only), so that is the sole run Table S2 cites.
- **The root accession can be re-verified from the reads.** The demultiplexed reads keep their SRA
  accession in the read name (the `Root1_rep1` R1 FASTQ begins `@SRR12331466.1_...`). Table S2
  performs that check only when `--root-fastq` points at an existing file, prints a note when it
  skips, and hard-fails on a mismatch. An unreachable file is not evidence of a wrong accession.
- **PubMed IDs are resolved from NCBI, not from a reference manager.** Zhang et al. 2024 is
  38589969 (Genome Biol 25(1):90) and Marand et al. 2021 is 33964211 (Cell 184(11):3041-3055.e21).
- **Two plate designs exist for the interspecies library.** Table S1 and S2 use
  `config/PlateDesign_SM2_ATAC.txt`, the two-block design (maize, Arabidopsis) that was supplied to
  AmbientMapper and matches the released data. The three-block scifi-demux design
  (`PlateDesign_scifi_At_B73_rep1.txt`) carries a third block that is not part of this study, and
  Table S1 hard-fails if handed it.
- **`p_top2` from `cells_calls` is not the expected-genome fraction.** Table S3 uses
  `expected_frac` from `SM2v2_pre_barcode_composition.tsv.gz`. The wrong source gives 113
  pure-wrong barcodes instead of 5.
- **Cleaning mode is read from the run tree, not assumed from a directory name.** Table S1 reads
  `design_file` out of each `*_decontam_params.json`, and Table S2 checks which decontamination
  directories exist. Only the interspecies library was run both WD and ND.

## Open

Table S2 carries one accession placeholder, printed at the end of every run: the maize and
Arabidopsis library, pending the NCBI BioProject, SRA and GEO submission. Fill it in the `PROSE`
dict of `make_TableS2_dataset_inventory.py` and re-run.
