# Software versions used for the manuscript

Recovered from run logs, BAM `@PG` headers, JSON sidecars, installed package metadata and the
pinned environment export, then spot-checked (2026-08-23 and 2026-08-24). "Route" distinguishes the
legacy processing of the maize and Arabidopsis library and the root library from the scifi-demux
processing of the Zhang et al. 2024 libraries.

| Tool | Version | Where used |
|---|---|---|
| UMI-tools | 1.1.6 | preprocessing, both routes |
| Cutadapt | 5.1 | preprocessing, both routes |
| seqkit | 2.9.0 | FASTQ chunking, scifi-demux route |
| pigz | 2.8 | compression in clean-bams and BED steps |
| scifi-demux | 0.1.1 (step 1, pre-`cli/` build), 0.1.2 at commit e34e47b (step 2); repository release 0.1.3 | Zhang libraries |
| BWA-MEM | 0.7.17-r1188 | mapping, both routes, and WASP remapping |
| samtools | 1.18 (filter and sort) and 1.22 (view and merge), legacy route; 1.22, scifi-demux route | mapping |
| Picard MarkDuplicates | 2.18.29, legacy route; 3.4.0, scifi-demux route | mapping |
| AmbientMapper | 0.1.0; assignment and genotyping at commit 362a8d8, decontamination at d7cb3b2, clean-bams at 77acee9 | genotyping and decontamination |
| Python | 3.12.7 | base environment |
| pysam | 0.23.3 | AmbientMapper, WASP |
| pandas | 2.1.4 | AmbientMapper |
| numpy | 1.26.4 | base environment and SEACells environment |
| DuckDB | 1.5.1 | AmbientMapper |
| ART | 2.5.8 (2016-06-05), 75 bp paired-end reads, HS25 profile | synthetic benchmark |
| minimap2 | 2.14-r883 | synthetic benchmark, ortholog alignment |
| MACS3 | 3.0.2 | peak calling for the synthetic benchmark and the concatenated-reference pre-pass |
| MACS2 | 2.2.4 | Socrates `callACRs`, per-genome QC |
| BEDTools | 2.31.1 | benchmark, browser tracks |
| R | 4.3.3 (cluster) | Socrates, purity model, Fig. 3 and Fig. 4 scripts |
| Socrates | 0.0.1, github.com/plantformatics/Socrates at commit 65077b2 | QC, SVD, UMAP, clustering |
| Seurat | 5.2.1 | Leiden clustering via `FindClusters(algorithm = 4)` |
| uwot | 0.2.2 | UMAP |
| igraph (R) | 2.1.0 or later | Leiden |
| leidenalg | 0.10.2 | Leiden |
| SEACells | 0.3.3 with Python 3.10, scanpy 1.11.5, anndata 0.11.4, scipy 1.15.3 | meta-cells |
| Souporcell | 2.1 (container built 2024-04-16), vartrix 1.1.22 inside | variant-based genotyping comparison |
| bcftools | 1.12 | variant panel preparation |
| WASP | github.com/bmvdgeijn/WASP master at commit d3b8447 (after release 0.3.4) | allele-specific remapping |
| PyTables | 3.10.1 | WASP |
| wigToBigWig | v4 | browser tracks |
| pyGenomeTracks | 3.9 | browser tracks |
| R packages, cluster | data.table 1.15.4, ggplot2 3.5.1, patchwork 1.3.0, ggpubr 0.6.0, glmmTMB 1.1.9 | Fig. 3, Fig. 4, purity model |
| R packages, workstation | ggplot2 4.0.3, patchwork 1.3.2, data.table 1.18.2.1, ggpubr 1.0.0, rstatix 1.1.0 (R 4.6.1 at the time of the final render) | Fig. 1, 2, 5, S1, S2, S7, S8 |
| openpyxl | 3.1.5 | supplementary tables |
