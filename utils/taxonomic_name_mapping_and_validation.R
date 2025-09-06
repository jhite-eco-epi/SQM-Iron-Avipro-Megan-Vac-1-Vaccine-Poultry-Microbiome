# Taxonomic Name Mapping and Validation Pipeline
# 1. Name mapping with synonym checking (old name -> current name)
# 2. Validation of mapped names against NCBI
# 3. Final validated lists for analysis

library(readr)
library(dplyr)
library(stringr)

#' Map taxonomic names to current synonyms, then validate against NCBI
#' 
#' @param taxa_vector Character vector of taxonomic names
#' @param taxonomic_level Character, "genus" or "phylum"
#' @return List with mapping and validation results
#' @export
map_and_validate_taxa <- function(taxa_vector, taxonomic_level = "genus") {
  
  # STEP 1: Name mapping with synonym checking
  # Initialize mapping results
  mapping_df <- data.frame(
    original_name = taxa_vector,
    current_name = taxa_vector,  # Default to same
    mapping_status = "no_change",
    mapping_notes = "",
    stringsAsFactors = FALSE
  )
  
  # Pattern-based pre-filtering (skip API calls for obvious artifacts)
  technical_patterns <- c(
    "^UCG-[0-9]+$", "^GCA-[0-9-]+$", "^RF[0-9]+$", "^CHKCI[0-9]+$", 
    "^CAG-[0-9]+$", "^DTU[0-9]+$", "^HT[0-9]+$", "^ASF[0-9]+$"
  )
  
  placeholder_patterns <- c(
    "uncultured", "unclassified", "Unassigned", "Incertae_Sedis"
  )
  
  environmental_patterns <- c(
    "_group$", "_cluster$", "_clade$", "vadinBE", "vadinCA", "vadinHA"
  )
  
  # Apply known classifier spelling corrections first
  classifier_corrections <- if (taxonomic_level == "phylum") {
    list(
      "Halobacterota" = "Halobacteriota"
    )
  } else {
    list(
      # Add genus-specific corrections here if needed
      # "Misspelled_genus" = "Correct_genus"
    )
  }
  
  for (old_name in names(classifier_corrections)) {
    if (old_name %in% mapping_df$original_name) {
      new_name <- classifier_corrections[[old_name]]
      idx <- which(mapping_df$original_name == old_name)
      # Replace the original name itself since it's a spelling error
      mapping_df$original_name[idx] <- new_name
      mapping_df$current_name[idx] <- new_name
      mapping_df$mapping_status[idx] <- "classifier_correction"
      mapping_df$mapping_notes[idx] <- paste("Classifier spelling correction:", old_name, "->", new_name)
    }
  }
  
  for (i in 1:nrow(mapping_df)) {
    taxon <- mapping_df$original_name[i]
    
    # Skip API calls for technical artifacts
    if (any(str_detect(taxon, technical_patterns))) {
      mapping_df$mapping_status[i] <- "technical_artifact"
      mapping_df$mapping_notes[i] <- "Technical code - no mapping needed"
      
    } else if (any(str_detect(taxon, placeholder_patterns))) {
      mapping_df$mapping_status[i] <- "placeholder"
      mapping_df$mapping_notes[i] <- "Taxonomic placeholder - no mapping needed"
      
    } else if (any(str_detect(taxon, environmental_patterns))) {
      mapping_df$mapping_status[i] <- "environmental"
      mapping_df$mapping_notes[i] <- "Environmental group - no mapping needed"
      
    } else if (nchar(taxon) < 3 || !str_detect(taxon, "[A-Za-z]")) {
      mapping_df$mapping_status[i] <- "invalid_format"
      mapping_df$mapping_notes[i] <- "Invalid format - no mapping needed"
      
    } else {
      # Potentially mappable - try synonym check
      mapping_df$mapping_status[i] <- "needs_mapping"
    }
  }
  
  # Count what needs mapping (exclude classifier corrections - they're already mapped)
  needs_mapping <- sum(mapping_df$mapping_status == "needs_mapping")
  pre_filtered <- nrow(mapping_df) - needs_mapping
  
  
  cat(sprintf("  Pre-filtered (no mapping needed): %d\n", pre_filtered))
  cat(sprintf("  Needs synonym checking: %d\n", needs_mapping))
  
  # NCBI synonym checking for mappable taxa
  if (needs_mapping > 0) {
    if (!requireNamespace("taxize", quietly = TRUE)) {
      stop("taxize package required. Install with: renv::install('taxize')")
    }
    library(taxize)
    
    # Only process names that actually need API calls (not classifier corrections)
    mapping_candidates <- mapping_df$original_name[mapping_df$mapping_status == "needs_mapping"]
    
    # Process in small batches with delays
    batch_size <- 3
    for (i in seq(1, length(mapping_candidates), by = batch_size)) {
      end_idx <- min(i + batch_size - 1, length(mapping_candidates))
      batch <- mapping_candidates[i:end_idx]
      
      # Small delay between batches with API key (10 req/s limit)
      if (i > 1) Sys.sleep(1)
      
      for (taxon in batch) {
        row_idx <- which(mapping_df$original_name == taxon)
        
        # Try synonym checking with error handling
        tryCatch({
          # Try to get UID - allow interactive selection for multiple matches
          uid_result <- taxize::get_uid(taxon, messages = TRUE, ask = TRUE)
          
          # Process the UID we got (user will have selected if multiple)
          if (!is.na(uid_result)) {
            target_rank <- if(taxonomic_level == "genus") "genus" else "phylum"
            
            # Get classification for the selected UID
            classification_result <- taxize::classification(uid_result, db = "ncbi", messages = FALSE)
            
            cat(sprintf("    Classification result for UID %s: %s\n", uid_result, 
                                   ifelse(is.null(classification_result[[1]]), "NULL", 
                                          paste("rows:", nrow(classification_result[[1]])))))
            
            if (!is.null(classification_result[[1]]) && nrow(classification_result[[1]]) > 0) {
              # Get the scientific name at the target taxonomic level
              rank_row <- classification_result[[1]][classification_result[[1]]$rank == target_rank, ]
              
              if (nrow(rank_row) > 0) {
                current_scientific_name <- rank_row$name[1]
                
                if (current_scientific_name == taxon) {
                  # Input name matches current scientific name
                  mapping_df$mapping_status[row_idx] <- "current_valid"
                  mapping_df$mapping_notes[row_idx] <- paste("Name is currently valid (", target_rank, "level)")
                } else {
                  # Input name is outdated - map to current name
                  mapping_df$current_name[row_idx] <- current_scientific_name
                  mapping_df$mapping_status[row_idx] <- "mapped_to_current"
                  mapping_df$mapping_notes[row_idx] <- paste("Updated to current scientific name:", taxon, "->", current_scientific_name, "(", target_rank, "level)")
                }
              } else {
                # Couldn't find the target taxonomic level - this means it's not at the right rank
                mapping_df$mapping_status[row_idx] <- "wrong_rank"
                mapping_df$mapping_notes[row_idx] <- paste("Name found in NCBI but not at", target_rank, "level")
              }
            } else {
              # Couldn't get classification - mark for manual review instead of assuming valid
              mapping_df$mapping_status[row_idx] <- "classification_failed"
              mapping_df$mapping_notes[row_idx] <- "Found in NCBI but classification failed - needs manual review"
            }
            
          } else {
            # Name not found - mark as no mapping found
            mapping_df$mapping_status[row_idx] <- "no_mapping"
            mapping_df$mapping_notes[row_idx] <- "Not found in NCBI (no UID returned)"
          }
          
          # Short delay between calls (10 req/s with API key)
          Sys.sleep(0.1)
          
        }, error = function(e) {
          mapping_df$mapping_status[row_idx] <- "mapping_error"
          mapping_df$mapping_notes[row_idx] <- paste("API error:", e$message)
          cat(sprintf("    ERROR processing %s: %s\n", taxon, e$message))
        })
        
        
        cat(sprintf("    %s -> %s (status: %s)\n", 
                    taxon, 
                    mapping_df$current_name[row_idx], 
                    mapping_df$mapping_status[row_idx]))
      }
    }
  }
  
  # STEP 2: Validation of mapped names
  
  # Initialize validation results
  validation_df <- mapping_df %>%
    mutate(
      validation_status = "pending",
      validation_notes = ""
    )
  
  # Only validate names that could potentially be valid (exclude wrong_rank)
  validatable_statuses <- c("no_change", "current_valid", "mapped_synonym", "mapped_to_current", "classifier_correction")
  needs_validation <- validation_df$mapping_status %in% validatable_statuses
  
  # Debug: show what needs validation
  cat(sprintf("  Mapping statuses: %s\n", paste(table(validation_df$mapping_status), collapse = ", ")))
  cat(sprintf("  Names that need validation: %s\n", paste(validation_df$current_name[needs_validation], collapse = ", ")))
  
  if (sum(needs_validation) > 0) {
    validation_candidates <- validation_df$current_name[needs_validation]
    
    # Validate in batches
    batch_size <- 5
    for (i in seq(1, length(validation_candidates), by = batch_size)) {
      end_idx <- min(i + batch_size - 1, length(validation_candidates))
      batch <- validation_candidates[i:end_idx]
      
      if (i > 1) Sys.sleep(1)
      
      for (current_name in batch) {
        row_indices <- which(validation_df$current_name == current_name & needs_validation)
        
        tryCatch({
          # Allow interactive selection for validation too
          uid <- taxize::get_uid(current_name, messages = TRUE, ask = TRUE)
          
          if (!is.na(uid)) {
            validation_df$validation_status[row_indices] <- "valid"
            validation_df$validation_notes[row_indices] <- "NCBI validated"
          } else {
            validation_df$validation_status[row_indices] <- "invalid"
            validation_df$validation_notes[row_indices] <- "Not found in NCBI"
          }
        }, error = function(e) {
          validation_df$validation_status[row_indices] <- "validation_error"
          validation_df$validation_notes[row_indices] <- paste("Validation error:", e$message)
        })
      }
    }
  }
  
  # Mark pre-filtered items as invalid
  pre_filtered_statuses <- c("technical_artifact", "placeholder", "environmental", "invalid_format", "no_mapping", "mapping_error")
  validation_df$validation_status[validation_df$mapping_status %in% pre_filtered_statuses] <- "invalid"
  validation_df$validation_notes[validation_df$mapping_status %in% pre_filtered_statuses] <- "Pre-filtered as invalid"
  
  # Mark classification failures as needing manual review
  validation_df$validation_status[validation_df$mapping_status == "classification_failed"] <- "manual_review"
  validation_df$validation_notes[validation_df$mapping_status == "classification_failed"] <- "Classification failed - needs manual review"
  
  # STEP 3: Create final results
  
  # Final validated list (only valid names)
  final_valid <- validation_df %>%
    filter(validation_status == "valid" & mapping_status != "wrong_rank") %>%
    select(!!paste0(taxonomic_level) := current_name) %>%
    distinct() %>%
    arrange(!!sym(paste0(taxonomic_level)))
  
  # Manual review list (classification failures and other uncertain cases)
  manual_review <- validation_df %>%
    filter(validation_status == "manual_review" | mapping_status == "wrong_rank") %>%
    select(original_name, current_name, mapping_status, mapping_notes, validation_status, validation_notes) %>%
    arrange(original_name)
  
  # Summary
  
  cat(sprintf("\nFINAL PIPELINE SUMMARY:\n"))
  cat(sprintf("=======================\n"))
  cat(sprintf("Original taxa: %d\n", length(taxa_vector)))
  cat(sprintf("Valid after mapping & validation: %d\n", nrow(final_valid)))
  cat(sprintf("Requiring manual review: %d\n", nrow(manual_review)))
  cat(sprintf("Mapping status counts:\n"))
  print(table(validation_df$mapping_status))
  cat(sprintf("Validation status counts:\n"))
  print(table(validation_df$validation_status))

  return(list(
    mapping = validation_df,
    final_valid = final_valid,
    manual_review = manual_review
  ))
}

#' Extract taxa from QIIME2 level data
#' 
#' @param level Integer, taxonomic level (2 = phylum, 6 = genus)
#' @return Character vector of taxa names
extract_taxa_from_level <- function(level = 6) {
  
  readRenviron(".Renviron")
  base_dir <- Sys.getenv("BASE_DATA_PATH", unset = "")
  csv_file <- file.path(base_dir, "exported_taxa_barplot_data", paste0("level-", level, ".csv"))
  
  if (!file.exists(csv_file)) {
    stop("Level ", level, " data file not found: ", csv_file)
  }
  
  data <- read_csv(csv_file, show_col_types = FALSE)
  taxa_cols <- setdiff(colnames(data), "index")
  
  if (level == 6) {
    # Extract genus names
    genus_resolved <- taxa_cols[str_detect(taxa_cols, "g__")]
    taxa_names <- str_extract(genus_resolved, "g__[^;]+$") %>% 
                  str_remove("g__") %>%
                  unique()
  } else if (level == 2) {
    # Extract phylum names  
    phylum_resolved <- taxa_cols[str_detect(taxa_cols, "p__")]
    taxa_names <- str_extract(phylum_resolved, "p__[^;]+") %>%
                  str_remove("p__") %>%
                  unique()
  } else {
    stop("Level ", level, " not supported. Use 2 (phylum) or 6 (genus)")
  }
  
  return(taxa_names)
}
