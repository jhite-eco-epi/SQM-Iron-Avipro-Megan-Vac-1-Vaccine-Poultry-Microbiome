# SalmonellaIron: Gut Microbiome Analysis

## Project Overview

This study investigates how vaccination against *Salmonella* and different iron supplementation regimens affect the gut microbial communities in poultry. Using 16S rRNA amplicon sequencing, we analyzed gut microbiome samples from birds across different treatment combinations:

- **Vaccination Status**: Unvaccinated vs. Vaccinated
- **Iron Supplementation**: Control vs. FeSO₄ vs. SQM® Iron

**Dataset**: 31 high-quality samples, rarefied to 2,906 sequences per sample, containing 893 unique ASVs after quality filtering.

## Key Findings

- **Beta Diversity**: Both vaccination and iron supplementation significantly affect community structure, with significant interaction effects
- **Alpha Diversity**: No significant effects on Shannon diversity or observed richness
- **Taxonomic Composition**: Iron supplementation significantly affects specific bacterial genera, with 11 significant changes detected (no significant phylum-level changes)


## Repository Structure

### 📁 **AlphaAnalysis/** - Within-Sample Diversity
Alpha diversity (within-sample) analyses and visualizations:
- **Main script**: `composite_alpha_diversity_plot.R`
- **Composite figure**: Observed richness and Shannon diversity plots
- **Statistical analysis**: Scheirer-Ray-Hare tests (non-parametric two-way ANOVA)
- **Key finding**: No significant effects of vaccination or iron supplementation
- 📖 **[Complete Methods & Results](AlphaAnalysis/ALPHA_README.md)**

### 📁 **BetaAnalysis/** - Between-Sample Community Structure  
Beta diversity (between-sample) community structure analyses:
- **Bray-Curtis analysis**: Abundance-based community differences ([Details](BetaAnalysis/BRAY_CURTIS_README.md))
- **Jaccard analysis**: Presence/absence-based community differences ([Details](BetaAnalysis/JACCARD_README.md))
- **PERMANOVA**: Significant effects of vaccination, iron supplementation, and their interaction
- **Composite figures**: Publication-ready combined visualizations
- 📖 **Key finding**: Both vaccination and iron treatments significantly affect community structure

### 📁 **TaxaAnalysis/** - Taxonomic Composition & Differential Abundance
Comprehensive taxonomic analysis with rigorous validation:
- **ANCOM-BC2 analysis**: Differential abundance testing at genus and phylum levels
- **Taxonomic validation**: NCBI-validated nomenclature with interactive validation scripts
- **Quality control**: Prevalence and library size filtering
- **Visualization**: Treatment-specific composition plots and significant taxa summaries
- **84 validated genera** analyzed across **29 samples** after filtering
- 📖 **[Complete Pipeline & Methods](TaxaAnalysis/README.md)**

### 📁 **Appendix/** - Supplementary Materials
Supplementary figures and quality control:
- **Rarefaction curves**: `supplementary_figure_1_rarefaction.png`
- **Analysis script**: `supplementary_figure_1_simple.R`

### 📁 **utils/** - Core Analysis Functions
Utility functions and data processing scripts:
- `metadata_utils.R` - Sample filtering and metadata cleaning functions
- `taxa_core.R` - Core taxonomic analysis functions  
- `taxonomic_name_mapping_and_validation.R` - NCBI validation utilities
- Reusable functions ensuring consistent data processing across all analyses

### 📁 **scripts/** - Bioinformatics Pipeline
External bioinformatics processing scripts:
- `classify_taxonomy_2024_10.sh` - QIIME2 taxonomic classification script

### 📁 **renv/** - R Environment Management
R environment management for reproducibility:
- `renv.lock` - Exact package versions used in analysis
- `activate.R` - Environment activation script
- `settings.json` - Project-specific renv settings
- Ensures consistent analysis environment across different systems

## Analysis Pipeline

1. **Data Processing** (QIIME2): Raw 16S sequences → Quality filtered ASVs
2. **Diversity Metrics**: Alpha diversity (Shannon, Richness) and Beta diversity (Bray-Curtis, Jaccard)
3. **Statistical Analysis**: PERMANOVA, ANCOM-BC2, Non-parametric tests
4. **Taxonomic Validation**: NCBI-based nomenclature validation
5. **Visualization**: Treatment comparisons, NMDS ordinations, taxonomic composition
6. **Integration**: Composite figures combining multiple analysis types

## Reproducing the Complete Analysis

This repository is designed for full reproducibility. Follow these steps to recreate the entire analysis:

### 1. Environment Setup
```r
# Install renv if not already installed
install.packages("renv")

# Activate the project environment
renv::activate()

# Restore exact package versions
renv::restore()
```

### 2. Run Analysis Modules (in order)

#### **Alpha Diversity Analysis**
```r
setwd("AlphaAnalysis/")
source("composite_alpha_diversity_plot.R")
```

#### **Beta Diversity Analysis**  
```r
setwd("BetaAnalysis/")
source("stats_bray_curtis.R")         # Bray-Curtis NMDS + PERMANOVA
source("stats_jaccard.R")            # Jaccard NMDS + PERMANOVA
source("composite_beta_diversity_plot.R")  # Combined visualizations
```

#### **Taxonomic Analysis**
```r
setwd("TaxaAnalysis/")
# Complete pipeline (includes validation, analysis, and visualization)
source("run_taxa_all.R")

# Or run components individually:
source("interactive_genus_validation.R")    # Taxonomic validation
source("interactive_phylum_validation.R")   
source("ancombc2_shared_analysis.R")        # Differential abundance
source("taxa_shared.R")                     # Composition visualization
source("create_significant_results_table.R") # Results summary
```

#### **Supplementary Figures**
```r
setwd("Appendix/")
source("supplementary_figure_1_simple.R")   # Rarefaction curves
```

### 3. Expected Outputs

After running the complete pipeline, you should have:
- **Alpha diversity plots** with statistical tables (AlphaAnalysis/)
- **Beta diversity NMDS plots** with PERMANOVA results (BetaAnalysis/)
- **Taxonomic composition plots** and differential abundance results (TaxaAnalysis/)
- **Supplementary rarefaction curves** (Appendix/)

## Key Methods

- **Sequencing**: 16S rRNA V4 region, paired-end Illumina
- **Processing**: QIIME2 DADA2 pipeline with phylogenetic diversity (see [BIOINFORMATICS_README.md](BIOINFORMATICS_README.md))
- **Statistics**: R-based analysis with vegan, ANCOMBC, and tidyverse
- **Validation**: NCBI Taxonomy Database cross-referencing
- **Visualization**: Publication-quality figures with consistent color schemes

## Dependencies

This project uses `renv` for package management. Key packages include:
- **qiime2R** (v1.0.13) - QIIME2 data import
- **vegan** (v2.7.1) - Ecological statistics and ordination
- **ANCOMBC** (v2.8.1) - Differential abundance testing
- **ggplot2, patchwork** - Visualization and composite figures
- **dplyr, tidyverse** - Data manipulation and analysis

Complete dependency list with exact versions: [`renv.lock`](renv.lock)

## Data Requirements

The analysis expects the following data files (typically in a `data/` directory):
- **Rarefied ASV table** from QIIME2 (`.qza` format)
- **Sample metadata** with treatment group information
- **Taxonomic classifications** at phylum and genus levels
- **Distance matrices** (Bray-Curtis and Jaccard) from QIIME2

For data processing details, see: 📖 [**BIOINFORMATICS_README.md**](BIOINFORMATICS_README.md)

---s
