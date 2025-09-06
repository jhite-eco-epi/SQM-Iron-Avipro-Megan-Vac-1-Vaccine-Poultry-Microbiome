# ANCOM-BC2 Phylum-Level Differential Abundance Analysis

**Analysis Date:** 2025-08-22  
**Validated phylum:** 12  
**Samples:** 29  

## Taxonomic Name Mappings Applied

The following taxonomic names were updated to current NCBI scientific names:

| Original Name | Current Scientific Name |
|---------------|-------------------------|
| Cyanobacteria | Cyanobacteriota |
| Actinobacteriota | Actinomycetota |
| Desulfobacterota | Thermodesulfobacteriota |
| Proteobacteria | Pseudomonadota |
| Euryarchaeota | Methanobacteriota |

## Statistical Results

**Vaccination Effects:** 8 phylum tested, 0 significant (q < 0.05), 1 large effects (|LFC| > 1.0)

### Vaccination Effects (Vaccinated vs Unvaccinated)

**Top Effects (including large non-significant):**

| Taxon | Log2 Fold Change | Effect Size | q-value | Significant | Direction |
|-------|------------------|-------------|---------|-------------|----------|
| Verrucomicrobiota | -1.613 | Large | 0.248 |  | Higher in Unvaccinated |

**Iron Treatment Effects (FeSO4 vs Control):** 8 phylum tested, 0 significant (q < 0.05), 1 large effects (|LFC| > 1.0)

### Iron Treatment Effects: FeSO4 vs Control

**Top Effects (including large non-significant):**

| Taxon | Log2 Fold Change | Effect Size | q-value | Significant | Direction |
|-------|------------------|-------------|---------|-------------|----------|
| Verrucomicrobiota | 2.351 | Large | 0.062 |  | Higher in FeSO4 |

**Iron Treatment Effects (SQM Iron vs Control):** 8 phylum tested, 0 significant (q < 0.05), 2 large effects (|LFC| > 1.0)

### Iron Treatment Effects: SQM Iron vs Control

**Top Effects (including large non-significant):**

| Taxon | Log2 Fold Change | Effect Size | q-value | Significant | Direction |
|-------|------------------|-------------|---------|-------------|----------|
| Verrucomicrobiota | 1.725 | Large | 0.432 |  | Higher in SQM_Iron |
| Thermodesulfobacteriota | 1.054 | Large | 0.619 |  | Higher in SQM_Iron |

**Interaction Effects (Vaccination × FeSO4):** 8 phylum tested, 0 significant (q < 0.05), 0 large effects (|LFC| > 1.0)

### Interaction Effects: Vaccination × FeSO4

No significant or large effects detected.

**Interaction Effects (Vaccination × SQM):** 8 phylum tested, 0 significant (q < 0.05), 3 large effects (|LFC| > 1.0)

### Interaction Effects: Vaccination × SQM

**Top Effects (including large non-significant):**

| Taxon | Log2 Fold Change | Effect Size | q-value | Significant | Direction |
|-------|------------------|-------------|---------|-------------|----------|
| Verrucomicrobiota | 1.438 | Large | 0.809 |  | Vaccination effect stronger in SQM |
| Cyanobacteriota | -1.391 | Large | 0.809 |  | Vaccination effect weaker in SQM |
| Campylobacterota | -1.004 | Large | 0.809 |  | Vaccination effect weaker in SQM |

## Methods

- **Statistical Method:** ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction)
- **Model:** Interaction model (`vaccine_status * iron_treatment`)
- **Multiple Testing Correction:** Benjamini-Hochberg (BH) procedure
- **Significance Threshold:** q < 0.05 (FDR-adjusted p-values)
- **Effect Size Categories:** Small (|LFC| < 0.5), Medium (0.5 ≤ |LFC| < 1.0), Large (|LFC| ≥ 1.0)
- **Taxonomic Validation:** All taxa validated against NCBI Taxonomy Database

## Data Files

Detailed results are available in CSV format:

**Main Effects:**
- [`ancombc2_phylum_vaccine_statusVaccinated.csv`](ancombc2_phylum_vaccine_statusVaccinated.csv) - Vaccination effects
- [`ancombc2_phylum_iron_treatmentFeSO4.csv`](ancombc2_phylum_iron_treatmentFeSO4.csv) - FeSO4 iron treatment effects
- [`ancombc2_phylum_iron_treatmentSQM_Iron.csv`](ancombc2_phylum_iron_treatmentSQM_Iron.csv) - SQM Iron treatment effects

**Interaction Effects:**
- [`ancombc2_phylum_interaction_vaccine_FeSO4.csv`](ancombc2_phylum_interaction_vaccine_FeSO4.csv) - Vaccination × FeSO4 interactions
- [`ancombc2_phylum_interaction_vaccine_SQM.csv`](ancombc2_phylum_interaction_vaccine_SQM.csv) - Vaccination × SQM interactions

