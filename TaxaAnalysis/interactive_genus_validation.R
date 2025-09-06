#!/usr/bin/env Rscript

# Interactive Genus Validation
# Run this script directly in your terminal to respond to prompts

source("utils/taxonomic_name_mapping_and_validation.R")

genus_taxa <- extract_taxa_from_level(6)
result <- map_and_validate_taxa(genus_taxa, "genus")
write_csv(result$mapping, "TaxaAnalysis/genus_name_mapping_complete.csv")
write_csv(result$final_valid, "TaxaAnalysis/final_validated_genera_complete.csv")
write_csv(result$manual_review, "TaxaAnalysis/genus_manual_review_required.csv")
