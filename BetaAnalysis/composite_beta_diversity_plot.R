#!/usr/bin/env Rscript
# composite_beta_diversity_plot.R
#
# Creates a composite figure with two beta diversity NMDS plots and shared legend:
# - Left panel: Bray-Curtis NMDS with taxa vectors
# - Right panel: Jaccard NMDS with taxa vectors  
# Both panels show treatment effects with shared legend
# Uses patchwork for combining plots similar to alpha diversity composite
# -------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(cowplot)  # for get_legend function

# -------------------------------
# Load Beta Diversity Plots from Source Scripts
# -------------------------------

cat("Loading Bray-Curtis NMDS plot from source script...\n")
source("BetaAnalysis/nmds_plot_bray.R")
bray_plot <- p1  # The plot object from bray script

cat("Loading Jaccard NMDS plot from source script...\n") 
source("BetaAnalysis/nmds_plot_jaccard.R")
jaccard_plot <- p1  # The plot object from jaccard script

# -------------------------------
# Prepare Plots for Composite Figure
# -------------------------------

# Add subplot labels inside the plot areas in upper left corner
bray_plot_labeled <- bray_plot + 
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

jaccard_plot_labeled <- jaccard_plot + 
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# -------------------------------
# Combine Plots with Shared Legend
# -------------------------------

# Extract legend from one of the plots (they should have the same legend)
shared_legend <- get_legend(bray_plot_labeled)

# Remove legends from individual plots
bray_plot_no_legend <- bray_plot_labeled + theme(legend.position = "none")
jaccard_plot_no_legend <- jaccard_plot_labeled + theme(legend.position = "none")

# Combine plots using patchwork - side by side
composite_plot <- bray_plot_no_legend + jaccard_plot_no_legend + 
  plot_layout(ncol = 2)

# Add shared legend on the right with more space for multiple legends
final_plot <- composite_plot + shared_legend + 
  plot_layout(ncol = 3, widths = c(3.5, 3.5, 2))

# -------------------------------
# Display and Save Plot
# -------------------------------

print(final_plot)

# Save high-quality PNG
cat("Saving composite beta diversity plot...\n")
ggsave(
  filename = "BetaAnalysis/composite_beta_diversity_plot.tiff",
  plot = final_plot,
  width = 18,
  height = 10,
  dpi = 600,
  bg = "white"
)

cat("Composite beta diversity plot saved as: BetaAnalysis/composite_beta_diversity_plot.tiff\n")

# Also save as PNG for preview
ggsave(
  filename = "BetaAnalysis/composite_beta_diversity_plot.png",
  plot = final_plot,
  width = 18,
  height = 10,
  dpi = 150,
  bg = "white"
)

cat("Preview PNG saved as: BetaAnalysis/composite_beta_diversity_plot.png\n")
cat("Figure shows Bray-Curtis NMDS on the left (A) and Jaccard NMDS on the right (B).\n")
cat("Both panels show treatment effects with shared legend.\n")