# Workflows

The pipeline that produced every table read by `analysis/`, in execution order. Each stage README
lists its scripts per dataset, the order to run them, the manuscript inputs they produce and the
resources the SLURM headers request. Scripts are HPC job scripts. Paths are relative to
`PROJECT_ROOT`, the root of the working tree that held the raw data.

```
01_preprocessing            demultiplexing, barcode correction, trimming
02_mapping                  BWA-MEM to every candidate reference, BAM cleaning, Tn5 insertion BEDs
03_genotyping               AmbientMapper extract, filter, assign, genotyping; synthetic benchmark; evaluation tables
03b_variant_based_comparison Souporcell v2.1 supervised genotyping, WASP remapping, allele purity
04_decontamination          AmbientMapper decontam (WD and ND) and clean-bams; combined-genome cleaning
05_qc_and_embedding         Socrates QC, SVD, UMAP, Leiden, co-embedding, SEACells meta-cells, marker annotation
```

Two routes exist upstream of AmbientMapper. The Zhang et al. 2024 libraries were processed with
scifi-demux. The maize and Arabidopsis library and the root library were processed with the legacy
UMI-tools, cutadapt and Perl chain kept under `legacy_sm2_root1/` in stages 01 and 02.
