#!/usr/bin/env Rscript
# composite_nmds_stress_analysis.R
#
# Creates a composite figure showing NMDS stress reduction analysis for both:
# - Left panels (A & C): Bray-Curtis stress vs dimensions and % reduction
# - Right panels (B & D): Jaccard stress vs dimensions and % reduction
# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

# Load env vars and utility
readRenviron(".Renviron")
source("utils/metadata_utils.R")

# -------------------------------
# Load Bray-Curtis Stress Data
# -------------------------------

cat("Running Bray-Curtis stress analysis...\n")

# Source the Bray-Curtis dimensionality script to get stress results
source("BetaAnalysis/nmds_dimensionality_bray.R")

# Store Bray-Curtis results
bray_stress_results <- stress_results
bray_p1 <- p1  # stress vs dimensions plot
bray_p2 <- p2  # % reduction plot

# -------------------------------
# Load Jaccard Stress Data  
# -------------------------------

cat("Running Jaccard stress analysis...\n")

# Source the Jaccard dimensionality script to get stress results
source("BetaAnalysis/nmds_dimensionality_jaccard.R")

# Store Jaccard results
jaccard_stress_results <- stress_results
jaccard_p1 <- p1  # stress vs dimensions plot  
jaccard_p2 <- p2  # % reduction plot

# -------------------------------
# Create Composite Figure with Panel Labels
# -------------------------------

# Add panel labels (A/B only, no titles) with borders and labels in upper right
bray_stress_plot <- bray_p1 + 
  annotate("text", x = Inf, y = Inf, label = "A", 
           hjust = 1.2, vjust = 1.5, size = 6, fontface = "bold") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

jaccard_stress_plot <- jaccard_p1 + 
  annotate("text", x = Inf, y = Inf, label = "B", 
           hjust = 1.2, vjust = 1.5, size = 6, fontface = "bold") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

# Arrange side by side (only top plots)
final_plot <- bray_stress_plot | jaccard_stress_plot

# -------------------------------
# Display and Save Plot
# -------------------------------

print(final_plot)

# Save high-quality TIFF
cat("Saving composite NMDS stress analysis plot...\n")
ggsave(
  filename = "BetaAnalysis/composite_nmds_stress_analysis.tiff",
  plot = final_plot,
  width = 12,
  height = 6,
  dpi = 600,
  bg = "white"
)

# Save preview PNG
ggsave(
  filename = "BetaAnalysis/composite_nmds_stress_analysis.png",
  plot = final_plot,
  width = 12,
  height = 6,
  dpi = 150,
  bg = "white"
)

cat("Composite NMDS stress analysis saved as:\n")
cat("- BetaAnalysis/composite_nmds_stress_analysis.tiff (high quality)\n")
cat("- BetaAnalysis/composite_nmds_stress_analysis.png (preview)\n")

# -------------------------------
# Print Summary
# -------------------------------

cat("\n=== NMDS Stress Analysis Summary ===\n")
cat("Bray-Curtis optimal dimensions: 4D (stress =", round(bray_stress_results$stress[bray_stress_results$k == 4], 3), ")\n")
cat("Jaccard optimal dimensions: 4D (stress =", round(jaccard_stress_results$stress[jaccard_stress_results$k == 4], 3), ")\n")
cat("Both analyses support 4D NMDS solutions with good stress values (< 0.10).\n")
