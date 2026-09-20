# Data

Nothing large is tracked here. `metadata/` holds the small tables the pipeline and the table
scripts read. `processed/` is the layout the `analysis/` scripts expect; its files are outputs of
`workflows/` and are not in git. Each dataset folder has a README listing the exact files.

## Sources

| Library | Raw reads | Processed tables |
|---|---|---|
| scifi-ATAC maize B73 and Arabidopsis (one plate, separate wells) | GEO, accession to be added | produced by `workflows/03` to `05`; the GEO record also carries the Socrates objects and Tn5 insertion BEDs |
| scifi-ATAC B73 and Mo17, replicates 1 and 2 (Zhang et al. 2024) | SRA SRR25320545 to SRR25320547, SRR25320539 to SRR25320541 | produced by `workflows/01` to `04` and `03b` |
| scifi-ATAC seven-genotype pool (Zhang et al. 2024) | SRA SRR25320542 to SRR25320544 | same |
| 10x scATAC-seq maize B73 root (Marand et al. 2021) | SRA SRR12331466, one run | produced by `workflows/02` to `03` |
| Synthetic benchmark | none, simulated | produced by `workflows/03_genotyping/synthetic/` |

## Layout of `processed/`

```
processed/
  scifiATAC_B73_Arabidopsis/
    SM2/        AmbientMapper run of 2026-01 (Fig. 1B to E)
    SM2v2/      AmbientMapper run of 2026-05, with_design (WD) and without_design (ND) (Fig. 2, 3, 5, Tables S1, S3)
    socrates/   Socrates objects and tables: SM2/, SM2v2_clean/, SM2v2_indep/, SM2v2_plate/, compare/, _data/
    bed/        concatenated-reference Tn5 BEDs (Fig. S1)
  marand2021_B73_root/Root1_rep1/   genotyping runs, evaluation tables, sub1k_B subsample
  synthetic/{synthetic,synthetic_disc}/   the two benchmark tracks
  zhang2024/{B73Mo17_rep1,B73Mo17_rep2,multiGenotypes_rep1}/   genotyping runs, souporcell, diagnostics
```

Subdirectory and file names mirror the pipeline outputs, so files can be dropped in under their
original names.
