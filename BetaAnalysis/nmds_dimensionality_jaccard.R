#!/usr/bin/env Rscript
# nmds_dimensionality_jaccard.R
#
# Test NMDS stress across different numbers of dimensions to find the optimal
# dimensionality for the poultry vaccination × iron-supplement trial.
# Uses the same QIIME2 Jaccard distances as other beta analyses.
# -------------------------------------------------------------------------

readRenviron(".Renviron")

# ---- 0.  Packages ---------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# ---- 1.  Load metadata & Jaccard distance matrix --------------------------

source("utils/metadata_utils.R")
md <- load_clean_metadata()

jaccard_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_jaccard_matrix.tsv")
if (!file.exists(jaccard_file)) {
  stop("Cannot find Jaccard distance matrix at ", jaccard_file)
}

jaccard_matrix <- read.table(jaccard_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Match samples between distance matrix and cleaned metadata
common_samples <- intersect(rownames(jaccard_matrix), md$sampleid)
jaccard_matrix_subset <- as.matrix(jaccard_matrix)[common_samples, common_samples]

cat("Testing NMDS dimensionality on", nrow(jaccard_matrix_subset), "samples\n")

# ---- 2.  Test NMDS across multiple dimensions -----------------------------

dimensions <- 1:6  # Test 1D through 6D
stress_results <- data.frame(
  k = integer(),
  stress = numeric(),
  converged = logical()
)

set.seed(123)

for (k in dimensions) {
  cat("Testing k =", k, "dimensions...")
  
  # Run NMDS with current dimensionality
  nmds_test <- metaMDS(jaccard_matrix_subset, k = k, trymax = 100, trace = FALSE)
  
  stress_results <- rbind(stress_results, data.frame(
    k = k,
    stress = nmds_test$stress,
    converged = nmds_test$converged
  ))
  
  cat(" stress =", round(nmds_test$stress, 4), 
      ifelse(nmds_test$converged, " (converged)", " (not converged)"), "\n")
}

# ---- 3.  Display results table --------------------------------------------

cat("\n=== NMDS Stress by Dimensionality (Jaccard) ===\n")
print(stress_results, row.names = FALSE)

# ---- 4.  Create stress plot (scree plot) ----------------------------------

# Calculate stress reduction between consecutive dimensions
stress_results <- stress_results %>%
  mutate(
    stress_reduction = lag(stress) - stress,
    pct_reduction = (stress_reduction / lag(stress)) * 100
  )

# Plot stress vs dimensions
p1 <- ggplot(stress_results, aes(x = k, y = stress)) +
  geom_line(size = 1, color = "blue") +
  geom_point(size = 3, color = "red") +
  geom_hline(yintercept = c(0.05, 0.1, 0.2), 
             linetype = "dashed", color = "gray50", alpha = 0.7) +
  annotate("text", x = max(dimensions), y = 0.05, label = "Excellent (< 0.05)", 
           hjust = 1, vjust = -0.5, size = 3, color = "gray50") +
  annotate("text", x = max(dimensions), y = 0.1, label = "Good (< 0.10)", 
           hjust = 1, vjust = -0.5, size = 3, color = "gray50") +
  annotate("text", x = max(dimensions), y = 0.2, label = "Fair (< 0.20)", 
           hjust = 1, vjust = -0.5, size = 3, color = "gray50") +
  scale_x_continuous(breaks = dimensions) +
  labs(x = "Number of Dimensions", y = "Stress") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

print(p1)

# ---- 5.  Stress reduction plot --------------------------------------------

p2 <- ggplot(stress_results %>% filter(k > 1), aes(x = k, y = pct_reduction)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = paste0(round(pct_reduction, 1), "%")), 
            vjust = -0.5, size = 3) +
  scale_x_continuous(breaks = dimensions[dimensions > 1]) +
  labs(x = "Dimensions Added", y = "% Stress Reduction") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

print(p2) 