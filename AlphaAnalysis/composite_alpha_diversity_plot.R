#!/usr/bin/env Rscript
# composite_alpha_diversity_plot.R
#
# Creates a composite figure with two subplots and shared legend:
# - Left panel: Observed Richness (ASV/OTU counts)  
# - Right panel: Shannon Diversity Index
# Both panels show vaccination status (x-axis) vs iron supplementation (fill color)
# -------------------------------------------------------------------------

library(qiime2R)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)  # for get_legend function

# Load env vars and utility
readRenviron(".Renviron")

# Source metadata utility
source("utils/metadata_utils.R")

# -------------------------------
# Load Shannon Diversity Data
# -------------------------------

# Read the exported Shannon diversity TSV file
shannon_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_shannon.tsv")
shannon_data <- read.table(shannon_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(shannon_data)[1] <- "sampleid"
names(shannon_data)[2] <- "shannon"

# -------------------------------
# Load Observed Features (Richness) Data
# -------------------------------

# Read the exported observed features TSV file
observed_features_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_observed_features.tsv")
observed_features_data <- read.table(observed_features_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(observed_features_data)[1] <- "sampleid"
names(observed_features_data)[2] <- "observed_features"

cat("Loaded observed features data for", nrow(observed_features_data), "samples\n")

# -------------------------------
# Load and Merge Metadata
# -------------------------------

# Load and clean metadata
treatment_data <- load_clean_metadata()

# Merge data
merged_shannon <- merge(treatment_data, shannon_data, by = "sampleid", all.x = TRUE)
merged_richness <- merge(treatment_data, observed_features_data, by = "sampleid", all.x = TRUE)

# Check that we have both datasets
cat("Shannon data:", nrow(merged_shannon), "samples\n")
cat("Richness data:", nrow(merged_richness), "samples\n")

# Remove any samples with missing data
merged_shannon <- merged_shannon[!is.na(merged_shannon$shannon), ]
merged_richness <- merged_richness[!is.na(merged_richness$observed_features), ]

# Find common samples between both datasets
common_samples <- intersect(merged_shannon$sampleid, merged_richness$sampleid)
merged_shannon <- merged_shannon[merged_shannon$sampleid %in% common_samples, ]
merged_richness <- merged_richness[merged_richness$sampleid %in% common_samples, ]

cat("Final dataset:", length(common_samples), "samples with both metrics\n")

# -------------------------------
# Prepare Treatment Variables (consistent for both datasets)
# -------------------------------

prepare_plot_data <- function(merged_data) {
  # Convert treatment variables to factors
  merged_data$vaccine <- factor(merged_data$vaccine, levels = c("0", "1"))
  merged_data$FeSO4 <- factor(merged_data$FeSO4, levels = c("0", "1"))
  merged_data$SQMFe <- factor(merged_data$SQMFe, levels = c("0", "1"))
  
  # Label vaccination status for readability
  merged_data$vaccine_label <- factor(merged_data$vaccine,
                                      levels = c("0", "1"),
                                      labels = c("Unvaccinated", "Vaccinated"))
  
  # Derive a single Supplement variable with Control category
  plot_data <- merged_data %>%
    mutate(Supplement = dplyr::case_when(
      FeSO4 == "1" & SQMFe == "0" ~ "FeSO4",
      FeSO4 == "0" & SQMFe == "1" ~ "SQM® Iron",
      FeSO4 == "0" & SQMFe == "0" ~ "Control",
      TRUE ~ "Other"  # fallback in case both supplements applied
    )) %>%
    # keep only the categories of interest (Control, FeSO4, SQMFe)
    filter(Supplement %in% c("Control", "FeSO4", "SQM® Iron"))
  
  plot_data$Supplement <- factor(plot_data$Supplement,
                                 levels = c("Control", "FeSO4", "SQM® Iron"))
  
  return(plot_data)
}

# Prepare both datasets
shannon_plot_data <- prepare_plot_data(merged_shannon)
richness_plot_data <- prepare_plot_data(merged_richness)

# -------------------------------
# Create Base Plot Theme (matching existing Shannon plot)
# -------------------------------

base_theme <- theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

# -------------------------------
# Create Observed Richness Plot (Left Panel)
# -------------------------------

richness_plot <- ggplot(richness_plot_data,
                        aes(x = vaccine_label, y = observed_features, fill = Supplement)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.8),
              alpha = 0.7, size = 2, shape = 21, colour = "black") +
  labs(x = "Vaccination Status",
       y = "Observed Richness (ASVs)",
       fill = "Supplementation") +
  base_theme +
  scale_fill_viridis_d(option = "plasma") +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# -------------------------------
# Create Shannon Diversity Plot (Right Panel)
# -------------------------------

shannon_plot <- ggplot(shannon_plot_data,
                       aes(x = vaccine_label, y = shannon, fill = Supplement)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.8),
              alpha = 0.7, size = 2, shape = 21, colour = "black") +
  labs(x = "Vaccination Status",
       y = "Shannon Diversity Index",
       fill = "Supplementation") +
  base_theme +
  scale_fill_viridis_d(option = "plasma") +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# -------------------------------
# Combine Plots with Shared Legend
# -------------------------------

# Extract legend from one of the plots
legend <- get_legend(shannon_plot)

# Remove legends from individual plots
richness_plot_no_legend <- richness_plot + theme(legend.position = "none")
shannon_plot_no_legend <- shannon_plot + theme(legend.position = "none")

# Combine plots using patchwork
composite_plot <- richness_plot_no_legend + shannon_plot_no_legend + 
  plot_layout(ncol = 2)

# Add shared legend on the right
final_plot <- composite_plot + legend + 
  plot_layout(ncol = 3, widths = c(4, 4, 1))

# -------------------------------
# Display and Save Plot
# -------------------------------

print(final_plot)

# Save high-quality PNG
cat("Saving composite alpha diversity plot...\n")
ggsave(
  filename = "AlphaAnalysis/composite_alpha_diversity_plot.tiff",
  plot = final_plot,
  width = 16,
  height = 10,
  dpi = 600,
  bg = "white"
)

cat("Composite plot saved as: AlphaAnalysis/composite_alpha_diversity_plot.tiff\n")

# Also save as PNG for preview
ggsave(
  filename = "AlphaAnalysis/composite_alpha_diversity_plot.png",
  plot = final_plot,
  width = 16,
  height = 10,
  dpi = 150,
  bg = "white"
)

cat("Preview PNG saved as: AlphaAnalysis/composite_alpha_diversity_plot.png\n")

# -------------------------------
# Print Summary Statistics
# -------------------------------

cat("\n=== Summary Statistics ===\n")
cat("Observed Richness by Treatment:\n")
richness_summary <- richness_plot_data %>%
  group_by(vaccine_label, Supplement) %>%
  summarise(
    n = n(),
    mean_richness = round(mean(observed_features, na.rm = TRUE), 1),
    sd_richness = round(sd(observed_features, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(richness_summary)

cat("\nShannon Diversity by Treatment:\n")
shannon_summary <- shannon_plot_data %>%
  group_by(vaccine_label, Supplement) %>%
  summarise(
    n = n(),
    mean_shannon = round(mean(shannon, na.rm = TRUE), 2),
    sd_shannon = round(sd(shannon, na.rm = TRUE), 2),
    .groups = "drop"
  )
print(shannon_summary)

cat("\nComposite figure created successfully!\n")
cat("Figure shows observed richness (ASVs/OTUs) on the left and Shannon diversity on the right.\n")
cat("Both panels use the same treatment structure with shared legend.\n")

# -------------------------------
# Statistical Analysis for Both Metrics
# -------------------------------

library(rcompanion)

cat("\n=== Statistical Analysis: Observed Richness ===\n")
cat("Scheirer-Ray-Hare test (non-parametric 2-way ANOVA):\n")
richness_stats <- scheirerRayHare(observed_features ~ vaccine_label * Supplement, data = richness_plot_data)
print(richness_stats)

cat("\n=== Statistical Analysis: Shannon Diversity ===\n")
cat("Scheirer-Ray-Hare test (non-parametric 2-way ANOVA):\n")
shannon_stats <- scheirerRayHare(shannon ~ vaccine_label * Supplement, data = shannon_plot_data)
print(shannon_stats)

# -------------------------------
# Save Statistical Results as CSV and Tables
# -------------------------------

library(gridExtra)
library(grid)
library(gtable)

# Format richness stats for output (include all rows including Residuals)
richness_results <- data.frame(
  Effect = rownames(richness_stats),
  Df = richness_stats$Df,
  Sum_Sq = round(richness_stats$`Sum Sq`, 2),
  H_statistic = round(richness_stats$H, 3),
  p_value = round(richness_stats$p.value, 5),
  stringsAsFactors = FALSE
)

# Add Total row for richness
total_df_richness <- sum(richness_stats$Df, na.rm = TRUE)
total_sum_sq_richness <- sum(richness_stats$`Sum Sq`, na.rm = TRUE)
richness_results <- rbind(richness_results, 
                         data.frame(Effect = "Total", 
                                   Df = total_df_richness,
                                   Sum_Sq = round(total_sum_sq_richness, 2),
                                   H_statistic = NA,
                                   p_value = NA,
                                   stringsAsFactors = FALSE))

# Format shannon stats for output (include all rows including Residuals)
shannon_results <- data.frame(
  Effect = rownames(shannon_stats),
  Df = shannon_stats$Df,
  Sum_Sq = round(shannon_stats$`Sum Sq`, 2),
  H_statistic = round(shannon_stats$H, 3),
  p_value = round(shannon_stats$p.value, 5),
  stringsAsFactors = FALSE
)

# Add Total row for shannon
total_df_shannon <- sum(shannon_stats$Df, na.rm = TRUE)
total_sum_sq_shannon <- sum(shannon_stats$`Sum Sq`, na.rm = TRUE)
shannon_results <- rbind(shannon_results,
                        data.frame(Effect = "Total", 
                                  Df = total_df_shannon,
                                  Sum_Sq = round(total_sum_sq_shannon, 2),
                                  H_statistic = NA,
                                  p_value = NA,
                                  stringsAsFactors = FALSE))

# Save CSV files
write.csv(richness_results, "AlphaAnalysis/observed_richness_stats.csv", row.names = FALSE)
write.csv(shannon_results, "AlphaAnalysis/shannon_diversity_stats.csv", row.names = FALSE)

cat("Statistical results saved as CSV files:\n")
cat("- AlphaAnalysis/observed_richness_stats.csv\n")
cat("- AlphaAnalysis/shannon_diversity_stats.csv\n")

# Create formatted tables matching other manuscript tables
create_stats_table <- function(stats_df, filename_base) {
  # Keep all rows including Residuals and Total
  display_stats <- stats_df
  
  # Clean up effect names for better display
  display_stats$Effect <- gsub("vaccine_label", "Vaccine", display_stats$Effect)
  display_stats$Effect <- gsub("Supplement", "Iron Supplement", display_stats$Effect)
  display_stats$Effect <- gsub("Vaccine:Iron Supplement", "Vaccine × Iron Supplement", display_stats$Effect)
  
  # Format column names for display (removed Significance column)
  colnames(display_stats) <- c("Effect", "df", "Sum Sq", "H", "p-value")
  
  # Create publication-style table grob with proper styling
  table_grob <- tableGrob(
    display_stats,
    rows = NULL,
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
  
  # Add horizontal line UNDER the header (matching other tables)
  table_grob <- gtable_add_grob(
    table_grob,
    grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), 
                        y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                        gp = gpar(lwd = 1.5, col = "black")),
    t = 1, b = 1, l = 1, r = ncol(table_grob)
  )
  
  return(list(table_grob = table_grob, display_data = display_stats))
}

# Generate table images
richness_table_result <- create_stats_table(richness_results, "observed_richness_stats_table")
shannon_table_result <- create_stats_table(shannon_results, "shannon_diversity_stats_table")

# Save richness table
richness_table_grob <- richness_table_result$table_grob
richness_display_data <- richness_table_result$display_data

tiff_height_richness <- max(3, nrow(richness_display_data) * 0.25 + 1)
tiff("AlphaAnalysis/observed_richness_stats_table.tiff", 
     width = 10, height = tiff_height_richness, units = "in", res = 300)
grid.draw(richness_table_grob)
dev.off()

png("AlphaAnalysis/observed_richness_stats_table.png", 
    width = 10, height = tiff_height_richness, units = "in", res = 150)
grid.draw(richness_table_grob)
dev.off()

# Save shannon table  
shannon_table_grob <- shannon_table_result$table_grob
shannon_display_data <- shannon_table_result$display_data

tiff_height_shannon <- max(3, nrow(shannon_display_data) * 0.25 + 1)
tiff("AlphaAnalysis/shannon_diversity_stats_table.tiff", 
     width = 10, height = tiff_height_shannon, units = "in", res = 300)
grid.draw(shannon_table_grob)
dev.off()

png("AlphaAnalysis/shannon_diversity_stats_table.png", 
    width = 10, height = tiff_height_shannon, units = "in", res = 150)
grid.draw(shannon_table_grob)
dev.off()

cat("Statistical tables saved as images:\n")
cat("- AlphaAnalysis/observed_richness_stats_table.tiff/.png\n")
cat("- AlphaAnalysis/shannon_diversity_stats_table.tiff/.png\n") 