# synthetic

Two tracks, `synthetic/` (all orthologous peaks) and `synthetic_disc/` (discriminative peaks), each
with 15 datasets `alpha_000` and `alpha_{002..050}_{Il14H,Ki11}` [03_genotyping/synthetic].

Per dataset: `cell_map_ref_chunks/*_filtered.tsv.gz` (Fig. 4A, S5D),
`genotyping_runs/factorial_phase2_2026-04-09/<config>/<ds>_cells_calls.tsv.gz` (Fig. 4B, S5),
`decontam_without_design_alpha05_C0/<ds>_{pre,post}_barcode_genome_counts.tsv.gz` (Fig. 4C),
`barcoded/<ds>/truth_table.tsv` (Fig. S5).
`synthetic/eval/phase2_2026-04-09/phase2_summary_metrics.tsv` [eval_phase_factorial.R phase2]: Fig. S5.
