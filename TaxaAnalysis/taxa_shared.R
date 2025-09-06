#!/usr/bin/env Rscript
# taxa_shared.R
# Shared, rank-agnostic pipeline to create donuts and tables for Level 2 (phylum) and Level 6 (genus)

readRenviron(".Renviron")

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

source("utils/metadata_utils.R")
# Legacy taxonomy utils replaced by NCBI validation pipeline
source("utils/taxa_core.R")

# Filter data to only validated taxa before computation with ANCOMBC2-style filtering
filter_to_validated_taxa <- function(data, validated_taxa, rank, name_mapping = NULL, apply_ancombc2_filters = TRUE) {
  if (rank == "phylum") {
    # Get phylum columns
    taxa_cols <- grep("^d__", names(data), value = TRUE)
    
    # Create reverse mapping to find original names that map to validated taxa
    original_names_to_keep <- names(name_mapping)[name_mapping %in% validated_taxa]
    
    # Add classifier spelling corrections (these were corrected at original_name level)
    classifier_corrections <- c("Halobacterota" = "Halobacteriota")  # Add other spelling errors here if needed
    corrected_names_to_keep <- names(classifier_corrections)[classifier_corrections %in% validated_taxa]
    
    all_valid_names <- c(validated_taxa, original_names_to_keep, corrected_names_to_keep)
    
    # Extract phylum names and check against validated list (including original names)
    valid_columns <- c()
    for (col in taxa_cols) {
      # Extract phylum name from full taxonomy string
      phylum_name <- str_extract(col, "p__[^;]+") %>% str_remove("p__")
      
      # Keep if it's a validated phylum (current or original name) or unclassified
      if (!is.na(phylum_name) && (phylum_name %in% all_valid_names || phylum_name == "")) {
        valid_columns <- c(valid_columns, col)
      } else if (str_detect(col, "d__Bacteria;__$|d__Archaea;__$")) {
        # Keep unclassified bacteria/archaea
        valid_columns <- c(valid_columns, col)
      }
    }
    
  } else if (rank == "genus") {
    # Get genus columns  
    taxa_cols <- grep("^d__", names(data), value = TRUE)
    
    # Create reverse mapping to find original names that map to validated taxa
    original_names_to_keep <- if (!is.null(name_mapping)) names(name_mapping)[name_mapping %in% validated_taxa] else c()
    all_valid_names <- c(validated_taxa, original_names_to_keep)
    
    # Extract genus names and check against validated list (including original names)
    valid_columns <- c()
    for (col in taxa_cols) {
      # Extract genus name from full taxonomy string
      genus_name <- str_extract(col, "g__[^;]+$") %>% str_remove("g__")
      
      # Keep if it's a validated genus (current or original name) or unclassified
      if (!is.na(genus_name) && (genus_name %in% all_valid_names || genus_name == "")) {
        valid_columns <- c(valid_columns, col)
      } else if (str_detect(col, ";__$")) {
        # Keep unclassified at any level
        valid_columns <- c(valid_columns, col)
      }
    }
  }
  
  # Return data with only validated taxa columns (plus metadata columns)
  metadata_cols <- setdiff(names(data), taxa_cols)
  filtered_data <- data[, c(metadata_cols, valid_columns)]
  
  # Apply ANCOMBC2-style filtering if requested
  if (apply_ancombc2_filters) {
    # 1. Library size cutoff: 1000 total counts per sample (like ANCOMBC2)
    sample_totals <- rowSums(filtered_data[, valid_columns], na.rm = TRUE)
    filtered_data <- filtered_data[sample_totals >= 1000, ]
    
    # 2. Prevalence cutoff: taxa must appear in at least 10% of samples (like ANCOMBC2)
    if (nrow(filtered_data) > 0) {
      prevalence_threshold <- 0.10
      n_samples <- nrow(filtered_data)
      
      # Calculate prevalence for each taxon (presence/absence)
      prevalence <- colSums(filtered_data[, valid_columns] > 0, na.rm = TRUE) / n_samples
      
      # Keep taxa that meet prevalence threshold
      keep_taxa <- names(prevalence)[prevalence >= prevalence_threshold]
      
      # Filter to keep only prevalent taxa
      filtered_data <- filtered_data[, c(metadata_cols, keep_taxa)]
      
      cat(sprintf("ANCOMBC2-style filtering applied (%s):\n", rank))
      cat(sprintf("  - Samples after library size filter (≥1000): %d\n", nrow(filtered_data)))
      cat(sprintf("  - Taxa after prevalence filter (≥10%%): %d\n", length(keep_taxa)))
    }
  }
  
  filtered_data
}

analyze_level <- function(
  data,
  rank = c("phylum", "genus"),
  output_img_prefix,           # e.g., "TaxaAnalysis/taxa_level2" or "TaxaAnalysis/taxa_level6"
  table_prefix = NULL          # e.g., "TaxaAnalysis/taxa_table_level2"; if NULL, no base table is written
) {
  rank <- match.arg(rank)

  # 1) Compute group means using NCBI-validated taxonomic data
  if (rank == "phylum") {
    # Load validated phyla with name mappings
    phylum_mapping_df <- read_csv("TaxaAnalysis/phylum_name_mapping_complete.csv", show_col_types = FALSE)
    validated_phyla <- read_csv("TaxaAnalysis/final_validated_phyla_complete.csv", show_col_types = FALSE)$phylum
    
    # Create mapping from original names to current names
    name_mapping <- phylum_mapping_df %>%
      filter(validation_status == "valid") %>%
      select(original_name, current_name) %>%
      deframe()
    
    # Filter input data to only validated phyla before computation
    filtered_data <- filter_to_validated_taxa(data, validated_phyla, rank = "phylum", name_mapping = name_mapping)
    
    # Use NCBI-validated mappings only on pre-filtered data
    means_collapsed <- compute_group_means_phylum(filtered_data, mapping = name_mapping)
    
    classified_col <- "is_classified_phylum"
    unclassified_labels <- c("Bacteria - unclassified phylum", "Archaea - unclassified phylum")
  } else {
    # Load validated genera
    genus_mapping_df <- read_csv("TaxaAnalysis/genus_name_mapping_complete.csv", show_col_types = FALSE)
    validated_genera <- read_csv("TaxaAnalysis/final_validated_genera_complete.csv", show_col_types = FALSE)$genus
    
    # Create mapping from original names to current names
    name_mapping <- genus_mapping_df %>%
      filter(validation_status == "valid") %>%
      select(original_name, current_name) %>%
      deframe()
    
    # Filter input data to only validated genera before computation
    filtered_data <- filter_to_validated_taxa(data, validated_genera, rank = "genus", name_mapping = name_mapping)
    
    # Compute on pre-filtered data
    means_collapsed <- compute_group_means_genus(filtered_data)
    
    classified_col <- "is_classified_genus"
    unclassified_labels <- c("Bacteria - unclassified genus", "Archaea - unclassified genus")
  }

  # For genus level, keep original data for tables and create filtered version for donut
  means_collapsed_original <- means_collapsed  # Keep unfiltered data for full tables
  
  if (rank == "genus") {
    # Calculate overall prevalence (mean abundance across all samples)
    taxa_prevalence <- means_collapsed %>%
      dplyr::filter(.data[[classified_col]]) %>%
      dplyr::group_by(Taxon_clean) %>%
      dplyr::summarise(mean_abundance = mean(raw_pct, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(desc(mean_abundance))
    
    # Get top 15 most prevalent taxa (excluding unclassified)
    top_taxa <- taxa_prevalence %>%
      dplyr::filter(!Taxon_clean %in% unclassified_labels) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::pull(Taxon_clean)
    
    # Keep only top 15 taxa plus unclassified for donut plot
    keep_taxa <- c(unclassified_labels, top_taxa)
    means_collapsed <- means_collapsed %>%
      dplyr::filter(Taxon_clean %in% keep_taxa)
    
    cat("Genus-level analysis: Filtering to top 15 most prevalent genera for donut plot.\n")
    cat("Top genera by mean abundance:\n")
    print(taxa_prevalence %>% dplyr::slice_head(n = 15))
  }
  
  # Order factor levels with unclassified first
  all_taxa <- sort(setdiff(unique(means_collapsed$Taxon_clean), unclassified_labels))
  ordered_levels <- c(unclassified_labels, all_taxa)
  means_collapsed <- means_collapsed %>%
    mutate(Taxon_clean = factor(Taxon_clean, levels = ordered_levels))

  # 2) Build donut data using the exact same base
  donut_data <- donut_from_means(means_collapsed, classified_only = TRUE, threshold = NULL, classified_col = classified_col)

  # 3) Plot helper
  plot_donut <- function(df) {
    ggplot(df %>% dplyr::filter(!is.na(Supplement), !is.na(vaccine_label)),
           aes(x = 2, y = Pct, fill = Taxon_clean)) +
      geom_col(width = 1, colour = "white") +
      coord_polar(theta = "y") +
      facet_grid(vaccine_label ~ Supplement, drop = TRUE, switch = "y") +
      xlim(0.5, 2.5) +
      scale_fill_viridis_d(name = "Taxa", option = "turbo", drop = TRUE) +
      theme_void() +
      theme(
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11),
        strip.text   = element_text(size = 14),
        strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside"
      )
  }

  # 4) Save donut
  dir.create("TaxaAnalysis", showWarnings = FALSE, recursive = TRUE)

  if (nrow(donut_data) > 0) {
    cat("Saving ", rank, " donut...\n", sep = "")
    ggsave(paste0(output_img_prefix, ".tiff"), plot = plot_donut(donut_data),
           width = 12, height = 10, dpi = 600, bg = "white")
    
    # Also save as PNG for preview
    ggsave(paste0(output_img_prefix, ".png"), plot = plot_donut(donut_data),
           width = 12, height = 10, dpi = 150, bg = "white")
    cat("Preview PNG saved as: ", paste0(output_img_prefix, ".png"), "\n", sep = "")
  } else {
    cat("No data for ", rank, " donut; skipping.\n", sep = "")
  }

  # 5) Save table data
  if (!is.null(table_prefix)) {
    main_levels <- c("Control", "FeSO4", "SQM® Iron")
    
    # Generate tables from original unfiltered data (all taxa)
    # Apply name mapping for genus data
    means_for_complete_tables <- means_collapsed_original
    if (rank == "genus") {
      means_for_complete_tables <- means_collapsed_original %>%
        mutate(Taxon_clean = update_taxon_names(
          as.character(Taxon_clean),
          rank = "genus",
          mapping = name_mapping
        ))
    }
    
    raw_wide_all_complete <- means_for_complete_tables %>%
      dplyr::filter(.data[[classified_col]]) %>%
      dplyr::select(vaccine_label, Supplement, Taxon_clean, raw_pct) %>%
      dplyr::filter(Supplement %in% main_levels) %>%
      dplyr::mutate(Supplement = forcats::fct_drop(Supplement)) %>%
      tidyr::pivot_wider(names_from = Supplement, values_from = raw_pct, values_fill = 0, names_expand = TRUE)

    for (nm in main_levels) {
      if (!nm %in% names(raw_wide_all_complete)) raw_wide_all_complete[[nm]] <- 0
    }

    full_classified_complete <- raw_wide_all_complete %>%
      dplyr::group_by(vaccine_label) %>%
      dplyr::mutate(
        FeSO4_abs_change = FeSO4 - Control,
        SQM_abs_change   = `SQM® Iron` - Control
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(vaccine_label, Taxon_clean, Control, FeSO4, `SQM® Iron`, FeSO4_abs_change, SQM_abs_change) %>%
      dplyr::arrange(vaccine_label, Taxon_clean)

    # Save complete (unfiltered) tables
    readr::write_csv(full_classified_complete, paste0(table_prefix, "_full_classified_complete.csv"))
    
    # Generate tables from filtered data (matches donut plot)
    raw_wide_all <- means_collapsed %>%
      dplyr::filter(.data[[classified_col]]) %>%
      dplyr::select(vaccine_label, Supplement, Taxon_clean, raw_pct) %>%
      dplyr::filter(Supplement %in% main_levels) %>%
      dplyr::mutate(Supplement = forcats::fct_drop(Supplement)) %>%
      tidyr::pivot_wider(names_from = Supplement, values_from = raw_pct, values_fill = 0, names_expand = TRUE)

    for (nm in main_levels) {
      if (!nm %in% names(raw_wide_all)) raw_wide_all[[nm]] <- 0
    }

    full_classified <- raw_wide_all %>%
      dplyr::group_by(vaccine_label) %>%
      dplyr::mutate(
        FeSO4_abs_change = FeSO4 - Control,
        SQM_abs_change   = `SQM® Iron` - Control
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(vaccine_label, Taxon_clean, Control, FeSO4, `SQM® Iron`, FeSO4_abs_change, SQM_abs_change) %>%
      dplyr::arrange(vaccine_label, Taxon_clean)

    # Save filtered tables (for genus: top 15 + unclassified; for phylum: all)
    readr::write_csv(full_classified, paste0(table_prefix, "_full_classified.csv"))
    
    # Facet table that matches donut exactly
    facet_long <- donut_data %>%
      dplyr::filter(Supplement %in% main_levels) %>%
      dplyr::select(vaccine_label, Supplement, Taxon_clean, raw_pct, Pct) %>%
      dplyr::arrange(vaccine_label, Supplement, Taxon_clean)
    readr::write_csv(facet_long, paste0(table_prefix, "_facet.csv"))
    
    # Also create complete facet table from original data
    if (rank == "genus") {
      # Apply name mapping to the complete data
      means_collapsed_complete_mapped <- means_collapsed_original %>%
        mutate(Taxon_clean = update_taxon_names(
          as.character(Taxon_clean),
          rank = "genus",
          mapping = name_mapping
        ))
      
      # Order factor levels for complete data with unclassified first
      all_taxa_complete <- sort(setdiff(unique(means_collapsed_complete_mapped$Taxon_clean), unclassified_labels))
      ordered_levels_complete <- c(unclassified_labels, all_taxa_complete)
      means_collapsed_complete_mapped <- means_collapsed_complete_mapped %>%
        mutate(Taxon_clean = factor(Taxon_clean, levels = ordered_levels_complete))
      
      # Generate complete donut data
      donut_data_complete <- donut_from_means(means_collapsed_complete_mapped, classified_only = TRUE, threshold = NULL, classified_col = classified_col)
      
      facet_long_complete <- donut_data_complete %>%
        dplyr::filter(Supplement %in% main_levels) %>%
        dplyr::select(vaccine_label, Supplement, Taxon_clean, raw_pct, Pct) %>%
        dplyr::arrange(vaccine_label, Supplement, Taxon_clean)
      readr::write_csv(facet_long_complete, paste0(table_prefix, "_facet_complete.csv"))
    }
  }

  invisible(list(
    means = means_collapsed,
    donut_data = donut_data
  ))
}


