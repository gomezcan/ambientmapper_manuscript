# Figure 2: noise anatomy of the AmbientMapper genotyping model

Panel A is a schematic drawn in Illustrator and has no script. The script runs from the
repository root and sources `analysis/_helpers/plotting.R`.

| Script | Panels | Inputs (repo paths, under `data/processed/scifiATAC_B73_Arabidopsis/`) | Run |
|---|---|---|---|
| `fig2.R` | B, C, D, E | `SM2v2/decontam_with_design_alpha05_v2/SM2v2_cells_calls.decontam.tsv.gz` (B to E), `SM2v2/decontam_with_design_alpha05_v2/SM2v2_pre_barcode_genome_counts.tsv.gz` (E). Optional: `SM2v2/final/SM2v2_{eta_final,empty_jsd_tau}.json`, absent for this run, so the JSD tau is derived from the data | `Rscript analysis/fig2_ambientmapper_framework/fig2.R` |

## Outputs

`figures/main/fig2/`: `Fig2B_EmptyGate_NoiseRegimes.{png,pdf}`, `Fig2C_DeltaBIC_NoiseRegimes.{png,pdf}`,
`Fig2D_EmptyGate_vs_reads.pdf`, `Fig2E_Barnyard_per_call.pdf`.

## Where the inputs come from

The design-aware AmbientMapper decontam run on SM2v2 (`workflows/04_decontamination/sm2v2/`),
preceded by the SM2v2 genotyping run (`workflows/03_genotyping/sm2v2/`).

## Notes

- `FIG2_SAMPLE=SM2 Rscript analysis/fig2_ambientmapper_framework/fig2.R` rebuilds the panels from the run in `SM2/` into `figures/main/fig2_SM2/`. The tables of that arm use a different call vocabulary; the schema shim in the script maps both vocabularies and panel E's call levels are data-driven.
- The regime thresholds in the CONFIG block (`EMPTY_BIC_MARGIN = 10`, `BIC_MARGIN_SD = 6`, `DOUBLET_MINOR_MIN = 0.20`, `LOW_EVIDENCE_READS_MAX = 10`) are the figure's own display and classification thresholds for the dashed gate lines and the regime colouring. They are not read from the run configuration.
- The JSD tau drawn in panel B is derived from the data (90th percentile of `jsd_to_eta` among barcodes with `delta_empty >= EMPTY_BIC_MARGIN`); the panel subtitle says so.
- Panels B to D use barcodes with more than 50 reads (D: more than 5); panel E joins calls of barcodes with more than 50 reads to the pre-clean winner counts.
