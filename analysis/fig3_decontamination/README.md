# Figure 3: decontamination against the plate-design ground truth

Both scripts run from the repository root and source `analysis/_helpers/plotting.R`.
The manuscript figure is assembled by hand from the individual panel exports.

| Script | Panels | Inputs (repo paths, under `data/processed/scifiATAC_B73_Arabidopsis/SM2v2/`) | Run |
|---|---|---|---|
| `fig3.R` (+ `fig3.sh`) | A, B, C, E | `decontam_with_design_alpha05_v2/SM2v2_{barcode_policy, pre_barcode_genome_counts, post_barcode_genome_counts, pre_barcode_composition, post_barcode_composition, cells_calls.decontam, barcode_postclean}.tsv.gz` | `Rscript analysis/fig3_decontamination/fig3.R`, or `mkdir -p _logs && sbatch analysis/fig3_decontamination/fig3.sh` |
| `fig3D_rescue_anatomy.R` (+ `.sh`) | D | WD: `decontam_with_design_alpha05_v2/SM2v2_{pre_barcode_composition, post_barcode_composition, cells_calls.decontam}.tsv.gz`. ND, read only with `FIG3D_EXPLORATORY=1`: `decontam_without_design_alpha05_v2/SM2v2_post_barcode_genome_counts.tsv.gz` | `Rscript analysis/fig3_decontamination/fig3D_rescue_anatomy.R`, or `mkdir -p _logs && sbatch analysis/fig3_decontamination/fig3D_rescue_anatomy.sh` |

## Outputs (`figures/main/fig3/`)

| File | Manuscript panel |
|---|---|
| `Fig3a_Banyard_restoration.pdf` | A |
| `Fig3b_ContaminationShift_ecdf.pdf` | B (the density version `Fig3b_ContaminationShift.pdf` is exported for comparison only) |
| `Fig3c_ReadsRemoved.pdf` | C |
| `Fig3d_RescueReadLoss.{pdf,png}` | D |
| `Fig3E_Specificity.pdf` | E |
| `Fig3_preview_ABCE.{png,pdf}` | composite preview of A, B, C, E, not the shipped figure |

With `FIG3D_EXPLORATORY=1`, `fig3D_rescue_anatomy.R` also writes `FigS_rescue_*` (panels R1 to R4
and the ND read-fate panel). None of those is in the manuscript.

## Where the inputs come from

WD = the design-aware AmbientMapper decontam run on SM2v2, ND = the design-free run, both in
`workflows/04_decontamination/sm2v2/`, preceded by `workflows/03_genotyping/sm2v2/`.

## Notes

- Cleaning outcome per barcode (panel C, colours in `status_cols`): Rescue = pre-clean expected-genome fraction below 0.50, at least 100 post-clean reads and post-clean expected-genome fraction of at least 0.90; Heavily Contaminated = pre-clean contamination of at least 0.50 and at least 50% of reads removed; Clean (Preserved) = fewer than 10% of reads removed; Mixed = the rest. Only barcodes with more than 10 pre-clean reads are plotted.
- Panels B, C and E exclude the `Low reads` and `Indistinguishable` call groups.
- `fig3D_rescue_anatomy.R` rebuilds the rescue class exactly as `fig3.R` does and hard-stops (`stopifnot`) if the class is not 400 barcodes or the medians drift from the values quoted in the manuscript legend (1,607.5 reads before, 153.5 after, 88.1% removed).
- Panel E: in the design-aware run every barcode carries a single-genome `allowed_set`, so no expected-genome read is removed by construction. The panel reports the off-target reads removed per well type.
- `fig3.R` also constructs a depth-bin panel for inspection; it is not exported and is not a manuscript panel.
- R packages beyond the tidyverse: `patchwork`, `scales`, `viridis`, `ggpubr`, `scattermore` (panel C only).
