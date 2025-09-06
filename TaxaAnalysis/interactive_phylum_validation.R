#!/usr/bin/env Rscript

# Interactive Phylum Validation
# Run this script directly in your terminal to respond to prompts

source("utils/taxonomic_name_mapping_and_validation.R")

phylum_taxa <- extract_taxa_from_level(2)
result <- map_and_validate_taxa(phylum_taxa, "phylum")
write_csv(result$mapping, "TaxaAnalysis/phylum_name_mapping_complete.csv")
write_csv(result$final_valid, "TaxaAnalysis/final_validated_phyla_complete.csv")

