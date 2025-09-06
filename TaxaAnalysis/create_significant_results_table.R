#!/usr/bin/env Rscript
# create_significant_results_table.R
# Creates separate tables for each significant ANCOMBC2 effect

library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)

# Read all ANCOMBC2 result files
vaccine_results <- read_csv("TaxaAnalysis/ancombc2_genus_vaccine_statusVaccinated.csv", show_col_types = FALSE)
feso4_results <- read_csv("TaxaAnalysis/ancombc2_genus_iron_treatmentFeSO4.csv", show_col_types = FALSE)
sqm_results <- read_csv("TaxaAnalysis/ancombc2_genus_iron_treatmentSQM_Iron.csv", show_col_types = FALSE)
interaction_feso4_results <- read_csv("TaxaAnalysis/ancombc2_genus_interaction_vaccine_FeSO4.csv", show_col_types = FALSE)
interaction_sqm_results <- read_csv("TaxaAnalysis/ancombc2_genus_interaction_vaccine_SQM.csv", show_col_types = FALSE)

# Function to extract significant results
extract_significant <- function(data, effect_name) {
  if ("significant_vaccine" %in% names(data)) {
    sig_col <- "significant_vaccine"
    lfc_col <- "lfc_vaccine"
    se_col <- "se_vaccine"
    p_col <- "p_vaccine"
    q_col <- "q_vaccine"
    effect_col <- "effect_size_vaccine"
    direction_col <- "direction_vaccine"
  } else if ("significant_iron" %in% names(data)) {
    sig_col <- "significant_iron"
    lfc_col <- "lfc_iron"
    se_col <- "se_iron"
    p_col <- "p_iron"
    q_col <- "q_iron"
    effect_col <- "effect_size_iron"
    direction_col <- "direction_iron"
  } else {
    sig_col <- "significant_interaction"
    lfc_col <- "lfc_interaction"
    se_col <- "se_interaction"
    p_col <- "p_interaction"
    q_col <- "q_interaction"
    effect_col <- "effect_size_interaction"
    direction_col <- "direction_interaction"
  }
  
  data %>%
    filter(!!sym(sig_col) == TRUE) %>%
    select(
      taxon,
      lfc = !!sym(lfc_col),
      se = !!sym(se_col),
      p_value = !!sym(p_col),
      q_value = !!sym(q_col),
      effect_size = !!sym(effect_col),
      direction = !!sym(direction_col)
    ) %>%
    mutate(effect = effect_name) %>%
    select(effect, taxon, lfc, se, p_value, q_value, effect_size, direction)
}

# Function to create table image using ggplot
create_table_image <- function(data, title, filename_base) {
  if (nrow(data) == 0) {
    cat("No significant results for", title, "\n")
    return(NULL)
  }
  
  # Prepare data for display
  display_data <- data %>%
    mutate(
      lfc = round(lfc, 3),
      se = round(se, 3),
      p_value = format(p_value, scientific = TRUE, digits = 2),
      q_value = round(q_value, 3)
    ) %>%
    select(taxon, lfc, se, p_value, q_value, effect_size, direction)
  
  # Save CSV
  write_csv(display_data, paste0("TaxaAnalysis/", filename_base, ".csv"))
  cat("Saved CSV:", paste0("TaxaAnalysis/", filename_base, ".csv"), "\n")
  
  # Create publication-style table grob with proper header line
  table_grob <- tableGrob(
    display_data,
    rows = NULL,
    cols = c("Genus", "Log2 FC", "SE", "p-value", "q-value", "Effect Size", "Direction"),
    theme = ttheme_minimal(
      core = list(
        fg_params = list(cex = 0.9, fontface = "plain"),
        bg_params = list(fill = "white", col = "white", lwd = 0)
      ),
      colhead = list(
        fg_params = list(cex = 0.9, fontface = "bold"),
        bg_params = list(fill = "white", col = "white", lwd = 0)
      ),
      rowhead = list(
        fg_params = list(cex = 0.9),
        bg_params = list(fill = "white", col = "white", lwd = 0)
      )
    )
  )
  
  # Add horizontal line UNDER the header
  table_grob <- gtable_add_grob(
    table_grob,
    grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), 
                        y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                        gp = gpar(lwd = 1.5, col = "black")),
    t = 1, b = 1, l = 1, r = ncol(table_grob)
  )
  
  # Save TIFF
  tiff_height <- max(3, nrow(display_data) * 0.25 + 1)
  tiff(paste0("TaxaAnalysis/", filename_base, ".tiff"), 
       width = 10, height = tiff_height, units = "in", res = 300)
  grid.draw(table_grob)
  dev.off()
  
  cat("Saved TIFF:", paste0("TaxaAnalysis/", filename_base, ".tiff"), "\n")
  
  # Also save as PNG for preview
  png(paste0("TaxaAnalysis/", filename_base, ".png"), 
      width = 10, height = tiff_height, units = "in", res = 150)
  grid.draw(table_grob)
  dev.off()
  
  cat("Preview PNG saved as:", paste0("TaxaAnalysis/", filename_base, ".png"), "\n")
  
  return(nrow(display_data))
}

# Extract and process each effect separately
effects_list <- list(
  list(data = vaccine_results, name = "Vaccination Effects", file = "ancombc2_significant_vaccination"),
  list(data = feso4_results, name = "FeSO₄ vs Control Effects", file = "ancombc2_significant_feso4"),
  list(data = sqm_results, name = "SQM Iron vs Control Effects", file = "ancombc2_significant_sqm"),
  list(data = interaction_feso4_results, name = "Vaccination × FeSO₄ Interactions", file = "ancombc2_significant_interactions_feso4"),
  list(data = interaction_sqm_results, name = "Vaccination × SQM Interactions", file = "ancombc2_significant_interactions_sqm")
)

# Process each effect
total_significant <- 0
for (effect in effects_list) {
  sig_data <- extract_significant(effect$data, effect$name)
  n_sig <- create_table_image(sig_data, effect$name, effect$file)
  if (!is.null(n_sig)) total_significant <- total_significant + n_sig
}

# Print summary
cat("\n=== Summary ===\n")
cat("Total significant results across all effects:", total_significant, "\n")

# Count by effect type
effect_counts <- c()
for (effect in effects_list) {
  sig_data <- extract_significant(effect$data, effect$name)
  if (nrow(sig_data) > 0) {
    effect_counts <- c(effect_counts, paste(effect$name, ":", nrow(sig_data)))
  }
}

if (length(effect_counts) > 0) {
  cat("Breakdown by effect:\n")
  for (count in effect_counts) {
    cat(" -", count, "\n")
  }
} else {
  cat("No significant results found across any effects (q < 0.05)\n")
}

cat("\nSignificant results table creation complete!\n")
