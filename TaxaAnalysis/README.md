# Taxonomic Composition Analysis

This directory contains comprehensive taxonomic analysis including differential abundance testing, taxonomic validation, and composition visualization for the poultry vaccination × iron supplementation experiment.

## Overview

The analysis examines treatment effects on microbial taxa at both phylum (level-2) and genus (level-6) taxonomic levels using compositionally-aware statistical methods and rigorous taxonomic validation procedures.

## Analysis Pipeline

### 1. Data Export and Initial Processing
- Raw taxonomic data exported from QIIME2 at phylum and genus levels
- Initial quality filtering and sample alignment with metadata

### 2. Taxonomic Name Validation (`taxonomic_name_mapping_and_validation.R`)
- **Interactive validation scripts**: `interactive_genus_validation.R`, `interactive_phylum_validation.R`
- Cross-reference against NCBI Taxonomy Database for current nomenclature
- Manual review and correction of ambiguous classifications
- **Output**: `genus_name_mapping_complete.csv`, `phylum_name_mapping_complete.csv`

### 3. Final Dataset Creation
- **Taxonomic validation output**: Interactive validation scripts create validated taxonomic datasets
  - `final_validated_genera_complete.csv` (from `interactive_genus_validation.R`)
  - `final_validated_phyla_complete.csv` (from `interactive_phylum_validation.R`)

### 4. Data Quality Control and Filtering
- **Sample filtering**: Minimum 1,000 sequence counts per sample (applied in analysis scripts)
- **Taxa filtering**: Present in ≥10% of samples (prevalence cutoff)
- **Implementation**: 
  - `ancombc2_shared_analysis.R` (ANCOMBC2 built-in filtering: `prv_cut = 0.10, lib_cut = 1000`)
  - `taxa_shared.R` (matching filters for visualization consistency)

### 5. Differential Abundance Testing (`ancombc2_shared_analysis.R`)
- **Method**: ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction 2)
- **Model**: `vaccine_status * iron_treatment` (factorial design with interactions)
- **Multiple testing correction**: Benjamini-Hochberg (FDR < 0.05)

### 6. Visualization and Summary (`taxa_shared.R`, `run_taxa_all.R`)
- Taxonomic composition plots by treatment group
- Significant results visualization
- Summary tables and statistical outputs

## Key Files and Outputs

### Analysis Scripts
- `ancombc2_shared_analysis.R` - Main differential abundance analysis
- `taxonomic_name_mapping_and_validation.R` - Taxonomic validation utilities
- `interactive_genus_validation.R` - Interactive genus name validation
- `interactive_phylum_validation.R` - Interactive phylum name validation
- `taxa_shared.R` - Taxonomic composition visualization
- `run_taxa_all.R` - Complete analysis pipeline runner
- `create_significant_results_table.R` - Summary of significant findings

### Validated Datasets
- `final_validated_genera_complete.csv` - Quality-filtered genus-level abundance table (84 genera, 29 samples)
- `final_validated_phyla_complete.csv` - Quality-filtered phylum-level abundance table
- `genus_name_mapping_complete.csv` - NCBI-validated genus nomenclature mappings
- `phylum_name_mapping_complete.csv` - NCBI-validated phylum nomenclature mappings

### Statistical Results

#### Genus-Level ANCOM-BC2 Results
- `ancombc2_genus_results.md` - Comprehensive statistical summary
- `ancombc2_genus_vaccine_statusVaccinated.csv` - Vaccination main effects
- `ancombc2_genus_iron_treatmentFeSO4.csv` - FeSO₄ treatment effects  
- `ancombc2_genus_iron_treatmentSQM_Iron.csv` - SQM Iron treatment effects
- `ancombc2_genus_interaction_vaccine_FeSO4.csv` - Vaccination × FeSO₄ interactions
- `ancombc2_genus_interaction_vaccine_SQM.csv` - Vaccination × SQM interactions

#### Phylum-Level ANCOM-BC2 Results
- `ancombc2_phylum_results.md` - Comprehensive statistical summary
- `ancombc2_phylum_vaccine_statusVaccinated.csv` - Vaccination main effects
- `ancombc2_phylum_iron_treatmentFeSO4.csv` - FeSO₄ treatment effects
- `ancombc2_phylum_iron_treatmentSQM_Iron.csv` - SQM Iron treatment effects
- `ancombc2_phylum_interaction_vaccine_FeSO4.csv` - Vaccination × FeSO₄ interactions
- `ancombc2_phylum_interaction_vaccine_SQM.csv` - Vaccination × SQM interactions

### Significant Results Summary
- `ancombc2_significant_feso4.csv` - Significant taxa affected by FeSO₄ treatment
- `ancombc2_significant_feso4.png` - Visualization of FeSO₄ effects
- `ancombc2_significant_sqm.csv` - Significant taxa affected by SQM Iron treatment  
- `ancombc2_significant_sqm.png` - Visualization of SQM Iron effects

### Taxonomic Composition Tables
- `taxa_table_level2_full_classified_complete.csv` - Complete phylum-level composition
- `taxa_table_level2_facet.csv` - Phylum composition summary by treatment
- `taxa_table_level6_full_classified_complete.csv` - Complete genus-level composition
- `taxa_table_level6_facet_complete.csv` - Genus composition summary by treatment

### Visualization Outputs
- `taxa_level2.png` - Phylum-level composition plot
- `taxa_level6.png` - Genus-level composition plot

## Statistical Parameters

- **Prevalence cutoff**: 10% (taxa must be present in ≥10% of samples)
- **Library size cutoff**: 1,000 total counts per sample
- **Taxonomic validation**: NCBI-validated taxa only
- **Multiple testing correction**: Benjamini-Hochberg (FDR)
- **Significance threshold**: q < 0.05
- **Effect size categories**: Small (|log₂FC| < 0.5), Medium (0.5 ≤ |log₂FC| < 1.0), Large (|log₂FC| ≥ 1.0)

## Final Dataset Specifications

- **Samples analyzed**: 29 (after quality filtering)
- **Genera analyzed**: 84 (after prevalence and validation filtering)
- **Phyla analyzed**: Variable (after prevalence and validation filtering)
- **Treatment groups**: 6 total (3 iron treatments × 2 vaccination statuses)

## Reproducing the Analysis

1. **Run taxonomic validation**: `source("interactive_genus_validation.R")` and `source("interactive_phylum_validation.R")`
2. **Execute main analysis**: `source("ancombc2_shared_analysis.R")`
3. **Generate visualizations**: `source("taxa_shared.R")`
4. **Create summary tables**: `source("create_significant_results_table.R")`
5. **Run complete pipeline**: `source("run_taxa_all.R")`

All analyses conducted in R version 4.4.2 using ANCOMBC package (version 2.8.1) and tidyverse ecosystem.
