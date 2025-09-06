#!/usr/bin/env Rscript
# taxa_core.R
# Shared core logic for taxa aggregation used by donuts and tables

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# Build treatment factors consistently
build_treatment_factors <- function(df) {
  df %>%
    mutate(
      vaccine_label = case_when(
        vaccination == "unvaccinated" ~ "Unvaccinated",
        vaccination == "vaccinated" ~ "Vaccinated",
        TRUE ~ "Unknown"
      ),
      vaccine_label = factor(vaccine_label, levels = c("Unvaccinated", "Vaccinated", "Unknown")),
      Supplement = case_when(
        FeSO4 == 1 & SQMFe == 0 ~ "FeSO4",
        FeSO4 == 0 & SQMFe == 1 ~ "SQM® Iron",
        FeSO4 == 0 & SQMFe == 0 ~ "Control",
        TRUE ~ "Unknown"
      ),
      Supplement = factor(Supplement, levels = c("Control", "FeSO4", "SQM® Iron", "Unknown"))
    )
}

# Extract a clean genus-level label and whether it is truly classified to genus
extract_genus_label <- function(taxon) {
  genus <- if_else(
    str_detect(taxon, ";g__[^;]*$"),
    str_remove(str_extract(taxon, ";g__[^;]*$"), "^;g__"),
    ""
  )
  # Treat placeholder genus labels (e.g., "uncultured", "unclassified", "unknown", "sp.") as unclassified
  genus_lower <- tolower(genus)
  is_placeholder <- genus != "" & stringr::str_detect(genus_lower, "uncultured|unclassified|unknown|unassigned|metagenome|^sp\\.?$")
  if (is_placeholder) {
    genus <- ""
  }
  is_classified <- genus != ""
  label <- case_when(
    is_classified ~ genus,
    str_detect(taxon, ";f__[^;]*$") & str_remove(str_extract(taxon, ";f__[^;]*$"), "^;f__") != "" ~
      paste0(str_remove(str_extract(taxon, ";f__[^;]*$"), "^;f__"), " (family)"),
    str_detect(taxon, ";o__[^;]*$") & str_remove(str_extract(taxon, ";o__[^;]*$"), "^;o__") != "" ~
      paste0(str_remove(str_extract(taxon, ";o__[^;]*$"), "^;o__"), " (order)"),
    str_detect(taxon, ";c__[^;]*$") & str_remove(str_extract(taxon, ";c__[^;]*$"), "^;c__") != "" ~
      paste0(str_remove(str_extract(taxon, ";c__[^;]*$"), "^;c__"), " (class)"),
    str_detect(taxon, ";p__[^;]*$") & str_remove(str_extract(taxon, ";p__[^;]*$"), "^;p__") != "" ~
      paste0(str_remove(str_extract(taxon, ";p__[^;]*$"), "^;p__"), " (phylum)"),
    str_detect(taxon, "^d__Bacteria") ~ "Bacteria - unclassified genus",
    str_detect(taxon, "^d__Archaea") ~ "Archaea - unclassified genus",
    TRUE ~ "Unclassified genus"
  )
  tibble(Taxon_clean = label, is_classified_genus = is_classified)
}

# Compute per-sample relative abundances and group means collapsed to clean genus labels
# Returns a tibble: vaccine_label, Supplement, Taxon_clean, is_classified_genus, raw_pct
compute_group_means_genus <- function(data, taxa_cols_pattern = "^d__") {
  taxa_cols <- grep(taxa_cols_pattern, names(data), value = TRUE)
  if (length(taxa_cols) == 0) stop("No taxon columns detected with pattern ", taxa_cols_pattern)

  long <- data %>%
    pivot_longer(cols = all_of(taxa_cols), names_to = "Taxon", values_to = "Count") %>%
    group_by(index) %>%
    mutate(RelAbundance = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup() %>%
    build_treatment_factors()

  group_means <- long %>%
    group_by(vaccine_label, Supplement, Taxon) %>%
    summarize(mean_pct = mean(RelAbundance * 100, na.rm = TRUE), .groups = "drop")

  label_df <- map_dfr(group_means$Taxon, extract_genus_label)
  means_labeled <- bind_cols(group_means, label_df)

  collapsed <- means_labeled %>%
    group_by(vaccine_label, Supplement, Taxon_clean) %>%
    summarize(
      raw_pct = sum(mean_pct, na.rm = TRUE),
      is_classified_genus = any(is_classified_genus),
      .groups = "drop"
    )

  # Keep only main supplements if present
  main_levels <- c("Control", "FeSO4", "SQM® Iron")
  if (any(collapsed$Supplement %in% main_levels)) {
    collapsed <- collapsed %>% filter(as.character(Supplement) %in% main_levels)
  }

  collapsed
}

# Compute per-sample relative abundances and group means collapsed to clean phylum labels
# Clean a QIIME phylum-level taxon string to a display label
clean_phylum_label <- function(taxon_string) {
  dplyr::case_when(
    taxon_string == "d__Bacteria;__" ~ "Bacteria - unclassified phylum",
    taxon_string == "d__Archaea;__" ~ "Archaea - unclassified phylum",
    stringr::str_detect(taxon_string, "^d__Bacteria;p__") ~ stringr::str_replace(taxon_string, "^d__Bacteria;p__", ""),
    stringr::str_detect(taxon_string, "^d__Archaea;p__") ~ stringr::str_replace(taxon_string, "^d__Archaea;p__", ""),
    TRUE ~ taxon_string
  )
}

# Apply taxonomic name mappings (from NCBI validation)
update_taxon_names <- function(names, rank = "phylum", mapping = NULL) {
  cleaned <- as.character(names)
  
  # Trim and collapse spaces
  cleaned <- stringr::str_trim(cleaned)
  cleaned <- stringr::str_replace_all(cleaned, "\\s+", " ")
  
  # Apply classifier spelling corrections first (same as in validation pipeline)
  if (rank == "phylum") {
    classifier_corrections <- c("Halobacterota" = "Halobacteriota")
    for (old_name in names(classifier_corrections)) {
      cleaned <- stringr::str_replace_all(cleaned, paste0("^", old_name, "$"), classifier_corrections[old_name])
    }
  }
  
  # Apply NCBI name mapping if supplied
  if (!is.null(mapping) && length(mapping) > 0) {
    idx <- match(cleaned, names(mapping))
    replace_mask <- !is.na(idx)
    cleaned[replace_mask] <- unname(mapping[idx[replace_mask]])
  }
  
  cleaned
}
# Returns a tibble: vaccine_label, Supplement, Taxon_clean, is_classified_phylum, raw_pct
compute_group_means_phylum <- function(data, taxa_cols_pattern = "^d__", mapping = NULL) {
  taxa_cols <- grep(taxa_cols_pattern, names(data), value = TRUE)
  if (length(taxa_cols) == 0) stop("No taxon columns detected with pattern ", taxa_cols_pattern)

  long <- data %>%
    pivot_longer(cols = all_of(taxa_cols), names_to = "Taxon", values_to = "Count") %>%
    group_by(index) %>%
    mutate(RelAbundance = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup() %>%
    build_treatment_factors()

  group_means <- long %>%
    group_by(vaccine_label, Supplement, Taxon) %>%
    summarize(mean_pct = mean(RelAbundance * 100, na.rm = TRUE), .groups = "drop")

  # Clean and modernize labels
  means_labeled <- group_means %>%
    mutate(Taxon_clean = clean_phylum_label(Taxon)) %>%
    mutate(Taxon_clean = stringr::str_replace(Taxon_clean, "Campilobacterota", "Campylobacterota")) %>%
    mutate(Taxon_clean = update_taxon_names(
      as.character(Taxon_clean),
      rank = "phylum",
      mapping = mapping
    )) %>%
    mutate(is_classified_phylum = !str_detect(Taxon_clean, "unclassified phylum"))

  collapsed <- means_labeled %>%
    group_by(vaccine_label, Supplement, Taxon_clean, is_classified_phylum) %>%
    summarize(raw_pct = sum(mean_pct, na.rm = TRUE), .groups = "drop")

  main_levels <- c("Control", "FeSO4", "SQM® Iron")
  if (any(collapsed$Supplement %in% main_levels)) {
    collapsed <- collapsed %>% filter(as.character(Supplement) %in% main_levels)
  }

  collapsed
}

# Build donut-ready normalized percentages (optionally classified-only and/or thresholded)
donut_from_means <- function(means_collapsed, classified_only = FALSE, threshold = NULL, classified_col = "is_classified_genus") {
  df <- means_collapsed
  if (classified_only) df <- df %>% filter(.data[[classified_col]])
  if (!is.null(threshold)) df <- df %>% group_by(vaccine_label, Supplement) %>% filter(raw_pct >= threshold) %>% ungroup()
  df %>%
    group_by(vaccine_label, Supplement) %>%
    filter(sum(raw_pct) > 0) %>%
    mutate(Pct = raw_pct / sum(raw_pct) * 100) %>%
    ungroup()
}

# Build a wide table of per-group raw_pct values from collapsed means (classified-only)
wide_table_from_means <- function(means_collapsed_classified, classified_col = "is_classified_genus") {
  means_collapsed_classified %>%
    filter(.data[[classified_col]]) %>%
    select(vaccine_label, Supplement, Taxon_clean, raw_pct) %>%
    tidyr::pivot_wider(
      names_from = Supplement, values_from = raw_pct,
      values_fill = 0, names_expand = TRUE
    ) %>%
    arrange(vaccine_label, Taxon_clean)
}


