#!/usr/bin/env Rscript
# supplementary_figure_1_simple_v2.R
#
# Create Supplementary Figure 1: Rarefaction curves using a robust approach
# -------------------------------------------------------------------------

# Load required libraries
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(tidyverse)
  library(RColorBrewer)
})

# ---- 1. Load and prepare data ---------------------------------------------

cat("Loading data...\n")

# Create temp directory for exports
dir.create("temp_export", showWarnings = FALSE)

# Export feature table from QIIME2
# We are temporarily expoerting this file within this script so we don't commit this file, since it is large.
cat("Exporting feature table from QIIME2...\n")
system("source ~/.zshrc && conda activate qiime2-amplicon-2024.10 &&qiime tools export --input-path data/table.qza --output-path temp_export/")
system("biom convert -i temp_export/feature-table.biom -o temp_export/feature-table.tsv --to-tsv")

# Load the exported feature table
otu_table <- read.table("temp_export/feature-table.tsv", 
                        sep = "\t", header = TRUE, row.names = 1, 
                        skip = 1, comment.char = "", check.names = FALSE)

# Transpose so samples are rows, ASVs are columns (for vegan)
otu_matrix <- t(as.matrix(otu_table))

# Load metadata
metadata <- read.table("data/PoulVaccMetadata.tsv", 
                       sep = "\t", header = TRUE, row.names = 1)

# Create treatment labels
metadata$Treatment <- case_when(
  metadata$Group == "GroupA" ~ "Unvaccinated Control",
  metadata$Group == "GroupE" ~ "Unvaccinated + FeSO4", 
  metadata$Group == "GroupF" ~ "Unvaccinated + SQM® Iron",
  metadata$Group == "GroupB" ~ "Vaccinated Control",
  metadata$Group == "GroupC" ~ "Vaccinated + FeSO4",
  metadata$Group == "GroupD" ~ "Vaccinated + SQM® Iron",
  TRUE ~ "Other"
)

# Match samples between OTU table and metadata
common_samples <- intersect(rownames(otu_matrix), rownames(metadata))
otu_matrix <- otu_matrix[common_samples, ]
metadata <- metadata[common_samples, ]

cat(paste("Analyzing", nrow(otu_matrix), "samples with", ncol(otu_matrix), "ASVs\n"))

# Get sample statistics
sample_sums <- rowSums(otu_matrix)
cat(paste("Sample depth range:", min(sample_sums), "to", max(sample_sums), "\n"))
cat(paste("Median sample depth:", median(sample_sums), "\n"))

# ---- 2. Generate rarefaction curves manually ------------------------------

# Function to calculate observed richness at different depths
calculate_rarefaction <- function(sample_vector, max_depth = NULL) {
  if (is.null(max_depth)) max_depth <- sum(sample_vector)
  
  # Create sequence of depths to test
  depths <- seq(100, min(max_depth, 4000), by = 100)
  
  richness_values <- numeric(length(depths))
  
  for (i in seq_along(depths)) {
    depth <- depths[i]
    if (depth > sum(sample_vector)) {
      richness_values[i] <- NA
    } else {
      # Use vegan's rarefy function for single rarefaction
      richness_values[i] <- rarefy(sample_vector, depth)
    }
  }
  
  # Return valid points only
  valid <- !is.na(richness_values)
  return(data.frame(Depth = depths[valid], Richness = richness_values[valid]))
}

# Generate rarefaction data for all samples
cat("Calculating rarefaction curves...\n")
rarefaction_data <- data.frame()

for (i in 1:nrow(otu_matrix)) {
  sample_name <- rownames(otu_matrix)[i]
  sample_data <- otu_matrix[i, ]
  
  # Calculate rarefaction curve for this sample
  rare_curve <- calculate_rarefaction(sample_data)
  
  if (nrow(rare_curve) > 0) {
    rare_curve$SampleID <- sample_name
    rare_curve$Treatment <- metadata[sample_name, "Treatment"]
    rarefaction_data <- rbind(rarefaction_data, rare_curve)
  }
}

cat(paste("Generated", nrow(rarefaction_data), "rarefaction data points\n"))

# ---- 3. Create the plot ---------------------------------------------------

if (nrow(rarefaction_data) == 0) {
  stop("No rarefaction data generated. Check data loading.")
}

# Set up colors for treatment groups
treatments <- unique(rarefaction_data$Treatment)
n_treatments <- length(treatments)
colors <- brewer.pal(max(3, min(8, n_treatments)), "Set2")[1:n_treatments]
names(colors) <- treatments

# Create the rarefaction plot
p_rare <- ggplot(rarefaction_data, aes(x = Depth, y = Richness, color = Treatment)) +
  # Individual sample curves (thin lines)
  geom_line(aes(group = SampleID), alpha = 0.6, linewidth = 0.4) +
  
  # Treatment group means (thick lines with confidence intervals)
  stat_smooth(aes(fill = Treatment), method = "loess", se = TRUE, 
              linewidth = 1.5, alpha = 0.2) +
  
  # Add vertical line at chosen rarefaction depth
  geom_vline(xintercept = 2906, 
             linetype = "dashed", 
             color = "red", 
             linewidth = 1.2, 
             alpha = 0.8) +
  
  # Customize theme with clean white background
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.width = unit(1.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  # Labels and colors
  labs(
    x = "Number of Sequences per Sample",
    y = "Observed ASV Richness",
    color = "Treatment Group",
    fill = "Treatment Group"
  ) +
  
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  
  # Add annotation for the rarefaction depth
  annotate("text", 
           x = 3200, 
           y = max(rarefaction_data$Richness) * 0.9, 
           label = "Chosen rarefaction\ndepth: 2,906", 
           color = "red", 
           size = 3.5, 
           fontface = "bold",
           hjust = 0)

# ---- 4. Save the figure ---------------------------------------------------

cat("Saving figure...\n")

# Save as PNG
ggsave("Appendix/supplementary_figure_1_rarefaction.png", 
       plot = p_rare, 
       width = 12, height = 8, units = "in", dpi = 150)


# ---- 5. Generate statistics and caption -----------------------------------

cat("\n=== Rarefaction Summary Statistics ===\n")
cat(paste("Total samples analyzed:", nrow(otu_matrix), "\n"))
cat(paste("Total ASVs detected:", ncol(otu_matrix), "\n"))
cat(paste("Maximum sequencing depth:", max(sample_sums), "\n"))
cat(paste("Minimum sequencing depth:", min(sample_sums), "\n"))
cat(paste("Median sequencing depth:", median(sample_sums), "\n"))
cat(paste("Chosen rarefaction depth: 2,906 sequences per sample\n"))

# Check plateau status at 2906
if (nrow(rarefaction_data) > 0) {
  plateau_data <- rarefaction_data %>%
    filter(Depth <= 2906) %>%
    group_by(Treatment, SampleID) %>%
    slice_max(Depth, n = 1) %>%
    group_by(Treatment) %>%
    summarise(
      n_samples = n(),
      mean_richness_at_2906 = mean(Richness, na.rm = TRUE),
      sd_richness = sd(Richness, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n=== Richness at 2,906 sequences by treatment ===\n")
  print(plateau_data)
}

# Clean up
cat("\nCleaning up temporary files...\n")
unlink("temp_export", recursive = TRUE)

cat("\n=== Files Generated ===\n")
cat("✓ Appendix/supplementary_figure_1_rarefaction.png (150 dpi)\n")
