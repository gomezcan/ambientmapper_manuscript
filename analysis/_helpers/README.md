# _helpers

Shared code sourced by the figure scripts, with repository-root-relative paths
(`source("analysis/_helpers/<file>")`), so every script is run from the repository root.

| File | Contents | Sourced by |
|---|---|---|
| `plotting.R` | species, stage and cleaning-outcome palettes, call-class colours, the shared ggplot theme, label helpers | Fig. 1, 2, 3 scripts |
| `data_utils.R` | small data helpers (`chop()` for barcode strings, loaders) | Fig. 1 and S1 scripts where needed |
| `fig4_helpers.R` | configuration order, labels, colours and shapes for the genotyping configurations (C0, S04, S10, ...), heatmap palette and theme | Fig. 4 parts 1 to 4, Fig. S5, Fig. S6 |
| `fig5_part2_helpers.R` | metadata readers for the per-genome plate objects, stage tokens (`SM2` = PreClean, `Clean.SM2v2wd` = WD, `Clean.SM2v2` = ND), UMAP and flow helpers | Fig. 5G to L, Fig. S8 |
