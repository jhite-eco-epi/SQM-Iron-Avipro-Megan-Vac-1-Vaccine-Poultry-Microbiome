# Bioinformatics Pipeline

This document details the bioinformatics pipeline used to process the raw 16S rRNA amplicon sequencing data into a feature table for downstream analysis. The analysis was conducted using QIIME 2 (v2024.2.0).

## 1. Denoising and ASV Table Generation

The raw, demultiplexed, paired-end FASTQ files were processed using the DADA2 plugin within QIIME 2 via the `denoise-paired` command. This single step performs quality filtering, denoising, merging of paired-end reads, and chimera removal to generate a feature table of Amplicon Sequence Variants (ASVs).

The specific parameters used for this step were:
- `--p-trim-left-f 8`
- `--p-trim-left-r 29`
- `--p-trunc-len-f 232`
- `--p-trunc-len-r 156`

**Results:** Quality filtering and denoising yielded a total of **617,853 paired-end reads** across **33 samples**, with a median of **13,846 reads per sample** (IQR: 9,137–26,942 reads). The initial processing generated **1,049 unique ASVs**.

The output of this process was the primary, un-rarefied feature table.

## 2. Core Diversity Analysis and Rarefaction

All core alpha and beta diversity metrics were calculated using the `qiime diversity core-metrics-phylogenetic` pipeline. This pipeline automatically rarefies the feature table before calculating diversity metrics.

The key parameter used for this step was:
- `--p-sampling-depth 2906`

This sampling depth was chosen to retain the maximum number of samples (31) while excluding two low-count technical replicates.

**Results:** After rarefaction, the final dataset contained **90,086 high-quality reads** classified into **893 unique ASVs** across **31 samples**. Each sample was standardized to exactly **2,906 sequences**.

One of the outputs of this pipeline was the file `data/rarefied_table.qza`, which is the feature table containing **31 samples** rarefied to an even depth of **2,906 sequences**. This table was used as the input for all subsequent statistical analyses and plots. On average, samples in this rarefied table contained approximately **195 ASVs**.

## 3. Post-QIIME 2 Analysis in R

Following the processing in QIIME 2, further data handling was performed in R for statistical analysis and visualization.

A key step was the cleaning and filtering of the sample metadata, performed by the `utils/metadata_utils.R` script. This script loads the full metadata file and removes specific samples before they are used in downstream analyses. The following four samples were excluded:
- `CN223` (removed as a suspected outlier)
- `CN414`, `CN415`, `CN416` (removed as technical replicates) 