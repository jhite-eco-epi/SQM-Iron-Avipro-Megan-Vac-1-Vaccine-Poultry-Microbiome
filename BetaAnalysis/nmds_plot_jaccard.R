#!/usr/bin/env Rscript
# jaccard_nmds_clean_plot.R
#
# Generate clean Jaccard NMDS ordination plot without envfit overlays
# Shows treatment effects with 95% confidence ellipses
# Uses 4D NMDS for better stress values but plots first two dimensions
# -------------------------------------------------------------------------

readRenviron(".Renviron")

# ---- 0.  Packages ---------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)    # metaMDS
  library(ggplot2)
})

# ---- 1.  Load metadata & remove outliers ----------------------------------

source("utils/metadata_utils.R")
md <- load_clean_metadata()

# ---- 2.  Load pre-computed Jaccard distance matrix ------------------------

jaccard_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_jaccard_matrix.tsv")
if (!file.exists(jaccard_file)) {
  stop("Cannot find Jaccard distance matrix at ", jaccard_file)
}

jaccard_matrix <- read.table(jaccard_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
jaccard_dist <- as.dist(jaccard_matrix)

# Match samples between distance matrix and cleaned metadata
common_samples <- intersect(rownames(jaccard_matrix), md$sampleid)
md <- md %>% filter(sampleid %in% common_samples)

# Subset and reorder to match distance matrix order
jaccard_matrix_subset <- as.matrix(jaccard_matrix)[common_samples, common_samples]

# ---- 3.  Run 4D NMDS on pre-computed Jaccard ------------------------------

cat("Running 4D NMDS on", nrow(jaccard_matrix_subset), "samples using QIIME2 Jaccard distances …\n")

set.seed(123)
nmds_result <- metaMDS(jaccard_matrix_subset, k = 4, trymax = 100, trace = FALSE)

cat("NMDS stress:", round(nmds_result$stress, 3), "\n")

# ---- 4.  Extract NMDS scores and merge with metadata ----------------------

scores_df <- as.data.frame(scores(nmds_result))
scores_df$sampleid <- rownames(scores_df)

# Merge with metadata
plot_data <- scores_df %>%
  left_join(md, by = "sampleid") %>%
  mutate(
    vaccine_label = factor(vaccination,
                           levels = c("unvaccinated", "vaccinated"),
                           labels = c("Unvaccinated", "Vaccinated")),
    Supplement = case_when(
      FeSO4 == 1 & SQMFe == 0 ~ "FeSO4",
      FeSO4 == 0 & SQMFe == 1 ~ "SQM® Iron",
      TRUE ~ "Control"),
    group = interaction(vaccine_label, Supplement, sep = "_")
  )

# Ensure proper factor levels
plot_data$Supplement <- factor(plot_data$Supplement,
                               levels = c("Control", "FeSO4", "SQM® Iron"))

# ---- 5.  Create clean NMDS plot --------------------------------------------

p1 <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, 
                            color = Supplement, 
                            shape = vaccine_label)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = interaction(vaccine_label, Supplement),
                   linetype = vaccine_label),
               type = "t", level = 0.95, linewidth = 1) +
  scale_color_viridis_d(name = "Iron Supplement", option = "plasma") +
  scale_shape_manual(name = "Vaccination Status",
                     values = c("Unvaccinated" = 16, "Vaccinated" = 17)) +
  scale_linetype_manual(name = "95% Confidence Ellipse",
                        values = c("Unvaccinated" = "dotted", "Vaccinated" = "solid"),
                        labels = c("Unvaccinated (dotted)", "Vaccinated (solid)")) +
  labs(x = "NMDS1",
       y = "NMDS2") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right"
  ) +
  guides(
    color = guide_legend(name = "Iron Supplement",
                         override.aes = list(shape = 15, size = 4, linetype = "blank"),
                         order = 1),
    shape = guide_legend(name = "Vaccination Status",
                         override.aes = list(size = 4),
                         order = 2),
    linetype = guide_legend(name = "95% Confidence Ellipse",
                           override.aes = list(color = "black", 
                                             linewidth = 2),
                           keywidth = unit(2, "cm"),
                           keyheight = unit(0.8, "cm"),
                           order = 3)
  )

# Display plot
print(p1)

# ---- 6.  Save plot ---------------------------------------------------------

cat("Saving Jaccard NMDS plot...\n")

# Save PNG version
ggsave(
  filename = "BetaAnalysis/nmds_plot_jaccard.png",
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Save TIFF version for publication
ggsave(
  filename = "BetaAnalysis/nmds_plot_jaccard.tiff",
  plot = p1,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)

cat("Jaccard NMDS plots saved:\n")
cat("  - BetaAnalysis/nmds_plot_jaccard.png\n")
cat("  - BetaAnalysis/nmds_plot_jaccard.tiff\n")
cat("4D NMDS stress =", round(nmds_result$stress, 3), "\n")
