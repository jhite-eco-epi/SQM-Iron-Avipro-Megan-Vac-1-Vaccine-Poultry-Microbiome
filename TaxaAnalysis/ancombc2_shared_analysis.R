#!/usr/bin/env Rscript

# ANCOM-BC2 Differential Abundance Analysis - Shared Function for Phylum and Genus Levels
# Uses NCBI-validated and renamed taxa from validation pipeline

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ANCOMBC)
library(phyloseq)
library(tibble)

# Load environment variables and utility functions
readRenviron(".Renviron")
source("utils/metadata_utils.R")

# Shared ANCOM-BC2 analysis function
run_ancombc2_analysis <- function(taxonomic_level = c("phylum", "genus"), verbose = TRUE) {
  taxonomic_level <- match.arg(taxonomic_level)
  
  if (verbose) cat(sprintf("Loading microbiome data for %s-level analysis...\n", taxonomic_level))
  metadata <- load_clean_metadata()
  
  # Determine data source and validation files based on taxonomic level
  if (taxonomic_level == "phylum") {
    data_level <- "level-2"
    mapping_file <- "TaxaAnalysis/phylum_name_mapping_complete.csv"
    validated_file <- "TaxaAnalysis/final_validated_phyla_complete.csv"
    col_name <- "phylum"
  } else {
    data_level <- "level-6" 
    mapping_file <- "TaxaAnalysis/genus_name_mapping_complete.csv"
    validated_file <- "TaxaAnalysis/final_validated_genera_complete.csv"
    col_name <- "genus"
  }
  
  if (verbose) cat(sprintf("Loading original %s data for statistical analysis...\n", data_level))
  base_dir <- Sys.getenv("BASE_DATA_PATH", unset = "")
  if (base_dir == "") stop("BASE_DATA_PATH is not set. Please export BASE_DATA_PATH to the directory containing exported_taxa_barplot_data.")
  
  csv_path <- file.path(base_dir, "exported_taxa_barplot_data", paste0(data_level, ".csv"))
  if (!file.exists(csv_path)) stop(sprintf("%s data file not found: %s", data_level, csv_path))
  data <- read_csv(csv_path, show_col_types = FALSE)
  
  # Sync with metadata
  keep_samples <- intersect(metadata$sampleid, data$index)
  if (length(keep_samples) == 0) stop("No overlapping samples between counts and cleaned metadata.")
  data <- data[data$index %in% keep_samples, ]
  data_with_treatments <- data %>%
    left_join(metadata[, c("sampleid", "vaccine", "FeSO4", "SQMFe")],
              by = c("index" = "sampleid")) %>%
    mutate(
      vaccine_status = ifelse(vaccine == 1, "Vaccinated", "Unvaccinated"),
      iron_treatment = case_when(
        FeSO4 == 1 & SQMFe == 0 ~ "FeSO4",
        FeSO4 == 0 & SQMFe == 1 ~ "SQM_Iron", 
        FeSO4 == 0 & SQMFe == 0 ~ "Control",
        TRUE ~ "Other"
      )
    ) %>%
    filter(iron_treatment != "Other")
  
  if (verbose) cat(sprintf("Data loaded: %d samples\n", nrow(data_with_treatments)))
  
  # Load validated taxa list with name mappings
  if (verbose) cat(sprintf("Loading validated %s with name mappings...\n", taxonomic_level))
  taxa_mapping <- read_csv(mapping_file, show_col_types = FALSE)
  validated_taxa <- read_csv(validated_file, show_col_types = FALSE)[[col_name]]
  if (verbose) cat(sprintf("Loaded %d validated %s (NCBI-validated with name mappings)\n", length(validated_taxa), taxonomic_level))
  
  # Show key name mappings
  name_changes <- taxa_mapping %>% 
    filter(original_name != current_name & validation_status == "valid") %>%
    select(original_name, current_name)
  if (nrow(name_changes) > 0 && verbose) {
    cat(sprintf("Key %s name mappings applied:\n", taxonomic_level))
    for (i in 1:nrow(name_changes)) {
      cat(sprintf("  %s -> %s\n", name_changes$original_name[i], name_changes$current_name[i]))
    }
  }
  
  # Create name mapping
  name_mapping <- taxa_mapping %>%
    filter(validation_status == "valid") %>%
    select(original_name, current_name) %>%
    deframe()
  
  # Process count data with taxonomic filtering and name mapping
  if (verbose) cat("Processing count data with taxonomic filtering and name mapping...\n")
  
  if (taxonomic_level == "phylum") {
    # Phylum-level processing
    count_table <- data_with_treatments %>%
      pivot_longer(cols = -c(index, vaccine_status, iron_treatment),
                   names_to = "taxon_full", values_to = "count") %>%
      mutate(
        taxon_original = str_extract(taxon_full, "p__[^;]+") %>% str_remove("p__")
      ) %>%
      filter(
        str_detect(taxon_full, "p__"),
        !is.na(taxon_original)
      )
  } else {
    # Genus-level processing  
    count_table <- data_with_treatments %>%
      pivot_longer(cols = -c(index, vaccine_status, iron_treatment),
                   names_to = "taxon_full", values_to = "count") %>%
      mutate(
        taxon_original = str_extract(taxon_full, "g__[^;]+$") %>% str_remove("g__")
      ) %>%
      filter(
        str_detect(taxon_full, "g__"),
        !is.na(taxon_original)
      )
  }
  
  # Apply name mappings and filter to validated taxa
  count_table <- count_table %>%
    # Apply classifier spelling corrections
    mutate(
      taxon_corrected = case_when(
        taxonomic_level == "phylum" & taxon_original == "Halobacterota" ~ "Halobacteriota",
        TRUE ~ taxon_original
      )
    ) %>%
    # Apply NCBI name mappings
    left_join(taxa_mapping %>% select(original_name, current_name, validation_status), 
              by = c("taxon_corrected" = "original_name")) %>%
    mutate(
      taxon_final = ifelse(!is.na(current_name), current_name, taxon_corrected)
    ) %>%
    # Filter to only validated taxa
    filter(
      validation_status == "valid",
      taxon_final %in% validated_taxa
    ) %>%
    # Ensure count is numeric and filter out zero counts
    mutate(count = as.numeric(count)) %>%
    filter(count > 0) %>%
    # Create unique sample identifiers
    mutate(sample_id = paste0(vaccine_status, "_", iron_treatment, "_", index)) %>%
    select(taxon = taxon_final, sample_id, count) %>%
    group_by(taxon, sample_id) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = sample_id, values_from = count, values_fill = 0) %>%
    column_to_rownames("taxon")
  
  if (verbose) cat(sprintf("Final count table: %d %s x %d samples\n", nrow(count_table), taxonomic_level, ncol(count_table)))
  
  # Create sample metadata for phyloseq
  sample_names <- colnames(count_table)
  sample_metadata <- data.frame(
    sample_id = sample_names,
    stringsAsFactors = FALSE
  ) %>%
    separate(sample_id, into = c("vaccine_status", "iron_treatment", "original_id"), sep = "_", remove = FALSE) %>%
    column_to_rownames("sample_id")
  
  # Create phyloseq object
  otu_mat <- as.matrix(count_table)
  sample_df <- sample_data(sample_metadata)
  phyloseq_obj <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE), sample_df)
  
  if (verbose) cat("Running ANCOM-BC2 analysis...\n")
  
  # Convert to count-like data for ANCOM-BC2
  phyloseq_transformed <- transform_sample_counts(phyloseq_obj, function(x) round(x * 10000))
  
  # Run ANCOM-BC2 with interaction model
  ancombc2_result <- ancombc2(
    data = phyloseq_transformed,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = "vaccine_status * iron_treatment",
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut = 0.10,
    lib_cut = 1000,
    s0_perc = 0.05,
    group = "vaccine_status",
    struc_zero = TRUE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE,
    global = TRUE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = NULL,
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
    trend_control = NULL
  )
  
  # Extract and process results
  results_df <- ancombc2_result$res
  if (verbose) cat("Processing results...\n")
  
  # Debug: check what columns are available
  if (verbose) {
    cat("Available result columns:\n")
    result_cols <- names(results_df)
    cat(paste(result_cols, collapse = ", "), "\n\n")
    
    # Show vaccine columns
    vaccine_cols <- result_cols[grepl("vaccine", result_cols)]
    cat("Vaccine-related columns:\n")
    cat(paste(vaccine_cols, collapse = ", "), "\n\n")
    
    # Show iron treatment columns
    iron_cols <- result_cols[grepl("iron", result_cols)]
    cat("Iron treatment-related columns:\n")
    cat(paste(iron_cols, collapse = ", "), "\n\n")
  }
  
  # Process vaccination effect results
  vaccine_results <- results_df %>%
    select(
      taxon,
      lfc_vaccine = lfc_vaccine_statusVaccinated,
      se_vaccine = se_vaccine_statusVaccinated,
      p_vaccine = p_vaccine_statusVaccinated,
      q_vaccine = q_vaccine_statusVaccinated
    ) %>%
    filter(!is.na(lfc_vaccine)) %>%
    mutate(
      significant_vaccine = q_vaccine < 0.05,
      effect_size_vaccine = case_when(
        abs(lfc_vaccine) < 0.5 ~ "Small",
        abs(lfc_vaccine) < 1.0 ~ "Medium",
        TRUE ~ "Large"
      ),
      direction_vaccine = ifelse(lfc_vaccine > 0, "Higher in Vaccinated", "Higher in Unvaccinated")
    ) %>%
    arrange(desc(abs(lfc_vaccine)))
  
  # Find iron treatment columns dynamically
  iron_cols <- names(results_df)[grepl("iron_treatment", names(results_df))]
  feso4_cols <- iron_cols[grepl("FeSO4", iron_cols)]
  sqm_cols <- iron_cols[grepl("SQM", iron_cols)]
  
  # Process iron treatment effects - FeSO4
  if (length(feso4_cols) > 0) {
    lfc_col <- feso4_cols[grepl("^lfc_", feso4_cols)][1]
    se_col <- feso4_cols[grepl("^se_", feso4_cols)][1]
    p_col <- feso4_cols[grepl("^p_", feso4_cols)][1]
    q_col <- feso4_cols[grepl("^q_", feso4_cols)][1]
    
    iron_feso4_results <- results_df %>%
      select(
        taxon,
        lfc_iron = !!lfc_col,
        se_iron = !!se_col,
        p_iron = !!p_col,
        q_iron = !!q_col
      ) %>%
      filter(!is.na(lfc_iron)) %>%
      mutate(
        significant_iron = q_iron < 0.05,
        effect_size_iron = case_when(
          abs(lfc_iron) < 0.5 ~ "Small",
          abs(lfc_iron) < 1.0 ~ "Medium",
          TRUE ~ "Large"
        ),
        direction_iron = ifelse(lfc_iron > 0, "Higher in FeSO4", "Higher in Control")
      ) %>%
      arrange(desc(abs(lfc_iron)))
  } else {
    iron_feso4_results <- data.frame()
    if (verbose) cat("No FeSO4 iron treatment results found\n")
  }
  
  # Process iron treatment effects - SQM
  if (length(sqm_cols) > 0) {
    lfc_col <- sqm_cols[grepl("^lfc_", sqm_cols)][1]
    se_col <- sqm_cols[grepl("^se_", sqm_cols)][1]
    p_col <- sqm_cols[grepl("^p_", sqm_cols)][1]
    q_col <- sqm_cols[grepl("^q_", sqm_cols)][1]
    
    iron_sqm_results <- results_df %>%
      select(
        taxon,
        lfc_iron = !!lfc_col,
        se_iron = !!se_col,
        p_iron = !!p_col,
        q_iron = !!q_col
      ) %>%
      filter(!is.na(lfc_iron)) %>%
      mutate(
        significant_iron = q_iron < 0.05,
        effect_size_iron = case_when(
          abs(lfc_iron) < 0.5 ~ "Small",
          abs(lfc_iron) < 1.0 ~ "Medium",
          TRUE ~ "Large"
        ),
        direction_iron = ifelse(lfc_iron > 0, "Higher in SQM_Iron", "Higher in Control")
      ) %>%
      arrange(desc(abs(lfc_iron)))
  } else {
    iron_sqm_results <- data.frame()
    if (verbose) cat("No SQM iron treatment results found\n")
  }
  
  # Return results list
  # Process interaction effects
  # Interaction: Vaccination effect in FeSO4 group
  interaction_feso4_results <- results_df %>%
    select(
      taxon,
      lfc_interaction = `lfc_vaccine_statusVaccinated:iron_treatmentFeSO4`,
      se_interaction = `se_vaccine_statusVaccinated:iron_treatmentFeSO4`,
      p_interaction = `p_vaccine_statusVaccinated:iron_treatmentFeSO4`,
      q_interaction = `q_vaccine_statusVaccinated:iron_treatmentFeSO4`
    ) %>%
    filter(!is.na(lfc_interaction)) %>%
    mutate(
      significant_interaction = q_interaction < 0.05,
      effect_size_interaction = case_when(
        abs(lfc_interaction) < 0.5 ~ "Small",
        abs(lfc_interaction) < 1.0 ~ "Medium",
        TRUE ~ "Large"
      ),
      direction_interaction = ifelse(lfc_interaction > 0, 
                                   "Vaccination effect stronger in FeSO4", 
                                   "Vaccination effect weaker in FeSO4")
    ) %>%
    arrange(desc(abs(lfc_interaction)))
  
  # Interaction: Vaccination effect in SQM group
  interaction_sqm_results <- results_df %>%
    select(
      taxon,
      lfc_interaction = `lfc_vaccine_statusVaccinated:iron_treatmentSQM`,
      se_interaction = `se_vaccine_statusVaccinated:iron_treatmentSQM`,
      p_interaction = `p_vaccine_statusVaccinated:iron_treatmentSQM`,
      q_interaction = `q_vaccine_statusVaccinated:iron_treatmentSQM`
    ) %>%
    filter(!is.na(lfc_interaction)) %>%
    mutate(
      significant_interaction = q_interaction < 0.05,
      effect_size_interaction = case_when(
        abs(lfc_interaction) < 0.5 ~ "Small",
        abs(lfc_interaction) < 1.0 ~ "Medium",
        TRUE ~ "Large"
      ),
      direction_interaction = ifelse(lfc_interaction > 0, 
                                   "Vaccination effect stronger in SQM", 
                                   "Vaccination effect weaker in SQM")
    ) %>%
    arrange(desc(abs(lfc_interaction)))

  list(
    taxonomic_level = taxonomic_level,
    sample_metadata = sample_metadata,
    validated_taxa = validated_taxa,
    name_changes = name_changes,
    vaccine_results = vaccine_results,
    iron_feso4_results = iron_feso4_results,
    iron_sqm_results = iron_sqm_results,
    interaction_feso4_results = interaction_feso4_results,
    interaction_sqm_results = interaction_sqm_results
  )
}

# Function to create markdown summary with statistical tables
create_ancombc2_markdown_summary <- function(results) {
  level <- results$taxonomic_level
  level_title <- str_to_title(level)
  
  # Helper function to create markdown table
  create_results_table <- function(data, title, treatment_type) {
    if (nrow(data) == 0) {
      return(paste0("### ", title, "\n\nNo results available.\n\n"))
    }
    
    # Determine which columns to use based on treatment type
    if (treatment_type == "vaccine") {
      sig_col <- "significant_vaccine"
      lfc_col <- "lfc_vaccine"
      se_col <- "se_vaccine"
      p_col <- "p_vaccine"
      q_col <- "q_vaccine"
      effect_col <- "effect_size_vaccine"
      direction_col <- "direction_vaccine"
    } else if (treatment_type == "interaction") {
      sig_col <- "significant_interaction"
      lfc_col <- "lfc_interaction"
      se_col <- "se_interaction"
      p_col <- "p_interaction"
      q_col <- "q_interaction"
      effect_col <- "effect_size_interaction"
      direction_col <- "direction_interaction"
    } else {
      sig_col <- "significant_iron"
      lfc_col <- "lfc_iron"
      se_col <- "se_iron"
      p_col <- "p_iron"
      q_col <- "q_iron"
      effect_col <- "effect_size_iron"
      direction_col <- "direction_iron"
    }
    
    # Get significant results
    sig_results <- data %>% 
      filter(!!sym(sig_col) == TRUE) %>%
      arrange(desc(abs(!!sym(lfc_col))))
    
    # Get top effects (significant + large effects)
    top_results <- data %>%
      filter(!!sym(sig_col) == TRUE | !!sym(effect_col) == "Large") %>%
      arrange(desc(abs(!!sym(lfc_col)))) %>%
      slice_head(n = 10)
    
    if (nrow(top_results) == 0) {
      return(paste0("### ", title, "\n\nNo significant or large effects detected.\n\n"))
    }
    
    # Create table
    table_content <- paste0("### ", title, "\n\n")
    
    if (nrow(sig_results) > 0) {
      table_content <- paste0(table_content, "**Significant Results (q < 0.05):**\n\n")
      table_content <- paste0(table_content, "| Taxon | Log2 Fold Change | Standard Error | p-value | q-value | Effect Size | Direction |\n")
      table_content <- paste0(table_content, "|-------|------------------|----------------|---------|---------|-------------|----------|\n")
      
      for (i in 1:nrow(sig_results)) {
        row <- sig_results[i, ]
        table_content <- paste0(table_content, sprintf("| %s | %.3f | %.3f | %.2e | %.3f | %s | %s |\n",
                                                       row$taxon, 
                                                       row[[lfc_col]], 
                                                       row[[se_col]], 
                                                       row[[p_col]], 
                                                       row[[q_col]], 
                                                       row[[effect_col]], 
                                                       row[[direction_col]]))
      }
      table_content <- paste0(table_content, "\n")
    }
    
    if (nrow(top_results) > nrow(sig_results)) {
      table_content <- paste0(table_content, "**Top Effects (including large non-significant):**\n\n")
      table_content <- paste0(table_content, "| Taxon | Log2 Fold Change | Effect Size | q-value | Significant | Direction |\n")
      table_content <- paste0(table_content, "|-------|------------------|-------------|---------|-------------|----------|\n")
      
      for (i in 1:nrow(top_results)) {
        row <- top_results[i, ]
        sig_symbol <- ifelse(row[[sig_col]], "✓", "")
        
        table_content <- paste0(table_content, sprintf("| %s | %.3f | %s | %.3f | %s | %s |\n",
                                                       row$taxon, 
                                                       row[[lfc_col]], 
                                                       row[[effect_col]], 
                                                       row[[q_col]], 
                                                       sig_symbol, 
                                                       row[[direction_col]]))
      }
      table_content <- paste0(table_content, "\n")
    }
    
    return(table_content)
  }
  
  # Build markdown content
  markdown <- paste0(
    "# ANCOM-BC2 ", level_title, "-Level Differential Abundance Analysis\n\n",
    "**Analysis Date:** ", Sys.Date(), "  \n",
    "**Validated ", level, ":** ", length(results$validated_taxa), "  \n",
    "**Samples:** ", nrow(results$sample_metadata), "  \n\n"
  )
  
  # Add name mappings section
  if (nrow(results$name_changes) > 0) {
    markdown <- paste0(markdown, "## Taxonomic Name Mappings Applied\n\n")
    markdown <- paste0(markdown, "The following taxonomic names were updated to current NCBI scientific names:\n\n")
    markdown <- paste0(markdown, "| Original Name | Current Scientific Name |\n")
    markdown <- paste0(markdown, "|---------------|-------------------------|\n")
    for (i in 1:nrow(results$name_changes)) {
      markdown <- paste0(markdown, sprintf("| %s | %s |\n", 
                                          results$name_changes$original_name[i], 
                                          results$name_changes$current_name[i]))
    }
    markdown <- paste0(markdown, "\n")
  }
  
  # Add results sections
  markdown <- paste0(markdown, "## Statistical Results\n\n")
  
  # Vaccination effects
  vaccine_summary <- sprintf("**Vaccination Effects:** %d %s tested, %d significant (q < 0.05), %d large effects (|LFC| > 1.0)\n\n",
                             nrow(results$vaccine_results),
                             level,
                             sum(results$vaccine_results$significant_vaccine),
                             sum(results$vaccine_results$effect_size_vaccine == "Large"))
  markdown <- paste0(markdown, vaccine_summary)
  
  if (nrow(results$vaccine_results) > 0) {
    markdown <- paste0(markdown, create_results_table(results$vaccine_results, "Vaccination Effects (Vaccinated vs Unvaccinated)", "vaccine"))
  }
  
  # Iron treatment effects - FeSO4
  feso4_summary <- sprintf("**Iron Treatment Effects (FeSO4 vs Control):** %d %s tested, %d significant (q < 0.05), %d large effects (|LFC| > 1.0)\n\n",
                           nrow(results$iron_feso4_results),
                           level,
                           sum(results$iron_feso4_results$significant_iron),
                           sum(results$iron_feso4_results$effect_size_iron == "Large"))
  markdown <- paste0(markdown, feso4_summary)
  
  if (nrow(results$iron_feso4_results) > 0) {
    markdown <- paste0(markdown, create_results_table(results$iron_feso4_results, "Iron Treatment Effects: FeSO4 vs Control", "iron"))
  }
  
  # Iron treatment effects - SQM
  sqm_summary <- sprintf("**Iron Treatment Effects (SQM Iron vs Control):** %d %s tested, %d significant (q < 0.05), %d large effects (|LFC| > 1.0)\n\n",
                         nrow(results$iron_sqm_results),
                         level,
                         sum(results$iron_sqm_results$significant_iron),
                         sum(results$iron_sqm_results$effect_size_iron == "Large"))
  markdown <- paste0(markdown, sqm_summary)
  
  if (nrow(results$iron_sqm_results) > 0) {
    markdown <- paste0(markdown, create_results_table(results$iron_sqm_results, "Iron Treatment Effects: SQM Iron vs Control", "iron"))
  }
  
  # Interaction effects - FeSO4
  feso4_interaction_summary <- sprintf("**Interaction Effects (Vaccination × FeSO4):** %d %s tested, %d significant (q < 0.05), %d large effects (|LFC| > 1.0)\n\n",
                                      nrow(results$interaction_feso4_results),
                                      level,
                                      sum(results$interaction_feso4_results$significant_interaction),
                                      sum(results$interaction_feso4_results$effect_size_interaction == "Large"))
  markdown <- paste0(markdown, feso4_interaction_summary)
  
  if (nrow(results$interaction_feso4_results) > 0) {
    markdown <- paste0(markdown, create_results_table(results$interaction_feso4_results, "Interaction Effects: Vaccination × FeSO4", "interaction"))
  }
  
  # Interaction effects - SQM
  sqm_interaction_summary <- sprintf("**Interaction Effects (Vaccination × SQM):** %d %s tested, %d significant (q < 0.05), %d large effects (|LFC| > 1.0)\n\n",
                                    nrow(results$interaction_sqm_results),
                                    level,
                                    sum(results$interaction_sqm_results$significant_interaction),
                                    sum(results$interaction_sqm_results$effect_size_interaction == "Large"))
  markdown <- paste0(markdown, sqm_interaction_summary)
  
  if (nrow(results$interaction_sqm_results) > 0) {
    markdown <- paste0(markdown, create_results_table(results$interaction_sqm_results, "Interaction Effects: Vaccination × SQM", "interaction"))
  }
  
  # Add methods section
  markdown <- paste0(markdown, "## Methods\n\n")
  markdown <- paste0(markdown, "- **Statistical Method:** ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction)\n")
  markdown <- paste0(markdown, "- **Model:** Interaction model (`vaccine_status * iron_treatment`)\n")
  markdown <- paste0(markdown, "- **Multiple Testing Correction:** Benjamini-Hochberg (BH) procedure\n")
  markdown <- paste0(markdown, "- **Significance Threshold:** q < 0.05 (FDR-adjusted p-values)\n")
  markdown <- paste0(markdown, "- **Effect Size Categories:** Small (|LFC| < 0.5), Medium (0.5 ≤ |LFC| < 1.0), Large (|LFC| ≥ 1.0)\n")
  markdown <- paste0(markdown, "- **Taxonomic Validation:** All taxa validated against NCBI Taxonomy Database\n\n")
  
  # Add data files section
  markdown <- paste0(markdown, "## Data Files\n\n")
  markdown <- paste0(markdown, "Detailed results are available in CSV format:\n\n")
  markdown <- paste0(markdown, "**Main Effects:**\n")
  markdown <- paste0(markdown, sprintf("- [`ancombc2_%s_vaccine_statusVaccinated.csv`](ancombc2_%s_vaccine_statusVaccinated.csv) - Vaccination effects\n", level, level))
  markdown <- paste0(markdown, sprintf("- [`ancombc2_%s_iron_treatmentFeSO4.csv`](ancombc2_%s_iron_treatmentFeSO4.csv) - FeSO4 iron treatment effects\n", level, level))
  markdown <- paste0(markdown, sprintf("- [`ancombc2_%s_iron_treatmentSQM_Iron.csv`](ancombc2_%s_iron_treatmentSQM_Iron.csv) - SQM Iron treatment effects\n\n", level, level))
  markdown <- paste0(markdown, "**Interaction Effects:**\n")
  markdown <- paste0(markdown, sprintf("- [`ancombc2_%s_interaction_vaccine_FeSO4.csv`](ancombc2_%s_interaction_vaccine_FeSO4.csv) - Vaccination × FeSO4 interactions\n", level, level))
  markdown <- paste0(markdown, sprintf("- [`ancombc2_%s_interaction_vaccine_SQM.csv`](ancombc2_%s_interaction_vaccine_SQM.csv) - Vaccination × SQM interactions\n", level, level))
  
  return(markdown)
}

# Function to save results and create summary
save_ancombc2_results <- function(results, verbose = TRUE) {
  level <- results$taxonomic_level
  
  if (verbose) cat(sprintf("Saving %s results...\n", level))
  
  # Save results
  write_csv(results$vaccine_results, sprintf("TaxaAnalysis/ancombc2_%s_vaccine_statusVaccinated.csv", level))
  write_csv(results$iron_feso4_results, sprintf("TaxaAnalysis/ancombc2_%s_iron_treatmentFeSO4.csv", level))
  write_csv(results$iron_sqm_results, sprintf("TaxaAnalysis/ancombc2_%s_iron_treatmentSQM_Iron.csv", level))
  write_csv(results$interaction_feso4_results, sprintf("TaxaAnalysis/ancombc2_%s_interaction_vaccine_FeSO4.csv", level))
  write_csv(results$interaction_sqm_results, sprintf("TaxaAnalysis/ancombc2_%s_interaction_vaccine_SQM.csv", level))
  
  # Create markdown summary with statistical results tables
  markdown_content <- create_ancombc2_markdown_summary(results)
  
  writeLines(markdown_content, sprintf("TaxaAnalysis/ancombc2_%s_results.md", level))
  if (verbose) cat(sprintf("Markdown results saved to: ancombc2_%s_results.md\n", level))
}

# Main execution function
run_both_ancombc2_analyses <- function(verbose = TRUE) {
  if (verbose) {
    cat("ANCOM-BC2 DIFFERENTIAL ABUNDANCE ANALYSIS\n")
    cat("========================================\n")
    cat("Running analysis for both phylum and genus levels\n\n")
  }
  
  # Run phylum analysis
  if (verbose) cat("=== PHYLUM LEVEL ANALYSIS ===\n")
  phylum_results <- run_ancombc2_analysis("phylum", verbose)
  save_ancombc2_results(phylum_results, verbose)
  
  if (verbose) cat("\n=== GENUS LEVEL ANALYSIS ===\n")
  genus_results <- run_ancombc2_analysis("genus", verbose)
  save_ancombc2_results(genus_results, verbose)
  
  if (verbose) {
    cat("\n=== ANALYSIS COMPLETE ===\n")
    cat("Both phylum and genus analyses completed successfully!\n")
  }
  
  list(phylum = phylum_results, genus = genus_results)
}


run_both_ancombc2_analyses()
