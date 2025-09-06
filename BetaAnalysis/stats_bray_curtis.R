library(vegan)
library(dplyr)
library(ggplot2)
# Load env vars and utility
readRenviron(".Renviron")
source("utils/metadata_utils.R")


# Load cleaned metadata
treatment_data <- load_clean_metadata()

# -------------------------------------------
# 3. Load Bray-Curtis Distance Matrix
# -------------------------------------------
bray_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_bray_curtis_distance.tsv")
bray_matrix <- read.table(bray_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
bray_dist <- as.dist(bray_matrix)

# -------------------------------------------
# 4. Match Metadata with Distance Matrix Samples
# -------------------------------------------
# Determine common sample IDs between distance matrix and cleaned metadata
common_samples <- intersect(rownames(bray_matrix), treatment_data$sampleid)

# Subset and reorder treatment_data to include only those common samples, in the order they appear in the distance matrix
treatment_data <- treatment_data[match(common_samples, treatment_data$sampleid), ]

# Also, subset the distance matrix to include only these common samples:
bray_matrix_subset <- as.matrix(bray_matrix)[common_samples, common_samples]


# -------------------------------------------
# 5. PERMANOVA: Vaccination × Supplement full factorial
# -------------------------------------------

meta_for_adonis <- treatment_data %>%
  mutate(
    vaccine_label = factor(vaccine, levels = c("0", "1"), labels = c("Unvaccinated", "Vaccinated")),
    Supplement = dplyr::case_when(
      FeSO4 == "1" & SQMFe == "0" ~ "FeSO4",
      FeSO4 == "0" & SQMFe == "1" ~ "SQM® Iron",
      FeSO4 == "0" & SQMFe == "0" ~ "Control",
      TRUE ~ "Other"
    )
  ) %>%
  select(sampleid, vaccine_label, Supplement) %>%
  rename(SampleID = sampleid)

meta_for_adonis$vaccine_label <- factor(meta_for_adonis$vaccine_label,
                                        levels = c("Unvaccinated", "Vaccinated"))
meta_for_adonis$Supplement <- factor(meta_for_adonis$Supplement,
                                     levels = c("Control", "FeSO4", "SQM® Iron", "Other"))

set.seed(123)
permanova_results <- vegan::adonis2(bray_matrix_subset ~ vaccine_label * Supplement,
                                    data = meta_for_adonis,
                                    permutations = 999,
                                    by = "terms")

cat("\n=====================\nPERMANOVA (Bray-Curtis)\n=====================\n")
print(permanova_results)

# -------------------------------
# Save PERMANOVA Results as CSV and Tables
# -------------------------------

library(gridExtra)
library(grid)
library(gtable)

# Format PERMANOVA results for output
bray_permanova_results <- data.frame(
  Effect = rownames(permanova_results),
  Df = permanova_results$Df,
  Sum_Sq = round(permanova_results$SumOfSqs, 4),
  R2 = round(permanova_results$R2, 3),
  F_statistic = round(permanova_results$F, 3),
  p_value = round(permanova_results$`Pr(>F)`, 5),
  stringsAsFactors = FALSE
)

# Clean up effect names for better display
bray_permanova_results$Effect <- gsub("vaccine_label", "Vaccine", bray_permanova_results$Effect)
bray_permanova_results$Effect <- gsub("Supplement", "Iron Supplement", bray_permanova_results$Effect)
bray_permanova_results$Effect <- gsub("Vaccine:Iron Supplement", "Vaccine × Iron Supplement", bray_permanova_results$Effect)

# Save CSV file
write.csv(bray_permanova_results, "BetaAnalysis/bray_curtis_permanova_stats.csv", row.names = FALSE)
cat("PERMANOVA results saved as CSV: BetaAnalysis/bray_curtis_permanova_stats.csv\n")

# Create formatted table
create_permanova_table <- function(stats_df, title) {
  # Format column names for display
  colnames(stats_df) <- c("Effect", "df", "Sum Sq", "R²", "F", "p-value")
  
  # Create publication-style table grob
  table_grob <- tableGrob(
    stats_df,
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
  
  # Add horizontal line UNDER the header
  table_grob <- gtable_add_grob(
    table_grob,
    grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), 
                        y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                        gp = gpar(lwd = 1.5, col = "black")),
    t = 1, b = 1, l = 1, r = ncol(table_grob)
  )
  
  return(list(table_grob = table_grob, display_data = stats_df))
}

# Generate table
bray_table_result <- create_permanova_table(bray_permanova_results, "Bray-Curtis PERMANOVA")
bray_table_grob <- bray_table_result$table_grob
bray_display_data <- bray_table_result$display_data

# Save tables
tiff_height <- max(3, nrow(bray_display_data) * 0.25 + 1)
tiff("BetaAnalysis/bray_curtis_permanova_table.tiff", 
     width = 10, height = tiff_height, units = "in", res = 300)
grid.draw(bray_table_grob)
dev.off()

png("BetaAnalysis/bray_curtis_permanova_table.png", 
    width = 10, height = tiff_height, units = "in", res = 150)
grid.draw(bray_table_grob)
dev.off()

cat("PERMANOVA table saved as: BetaAnalysis/bray_curtis_permanova_table.tiff/.png\n")


# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# vegan::adonis2(formula = bray_matrix_subset ~ vaccine_label * Supplement, data = meta_for_adonis, permutations = 999, by = "terms")
#                          Df SumOfSqs      R2      F Pr(>F)    
# vaccine_label             1   0.2328 0.07113 2.6475  0.001 ***
# Supplement                2   0.6032 0.18433 3.4305  0.001 ***
# vaccine_label:Supplement  2   0.4143 0.12659 2.3559  0.001 ***
# Residual                 23   2.0221 0.61794                  
# Total                    28   3.2723 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



library(pairwiseAdonis)
meta_for_adonis$Group <- interaction(meta_for_adonis$vaccine_label, meta_for_adonis$Supplement)
pairwise_results <- pairwise.adonis(bray_matrix_subset, meta_for_adonis$Group,
                                   p.adjust.m = "bonferroni")
cat("\n=====================\nPairwise PERMANOVA\n=====================\n")
print(pairwise_results)



#                                         pairs Df   SumsOfSqs
# 1  Unvaccinated.Control vs Vaccinated.Control  1 0.008361285
# 2    Unvaccinated.Control vs Vaccinated.FeSO4  1 0.017314831
# 3    Unvaccinated.Control vs Vaccinated.SQMFe  1 0.020900406
# 4  Unvaccinated.Control vs Unvaccinated.FeSO4  1 0.030177215
# 5  Unvaccinated.Control vs Unvaccinated.SQMFe  1 0.015673401
# 6      Vaccinated.Control vs Vaccinated.FeSO4  1 0.014148058
# 7      Vaccinated.Control vs Vaccinated.SQMFe  1 0.017273242
# 8    Vaccinated.Control vs Unvaccinated.FeSO4  1 0.023221570
# 9    Vaccinated.Control vs Unvaccinated.SQMFe  1 0.009327633
# 10       Vaccinated.FeSO4 vs Vaccinated.SQMFe  1 0.011151811
# 11     Vaccinated.FeSO4 vs Unvaccinated.FeSO4  1 0.012044748
# 12     Vaccinated.FeSO4 vs Unvaccinated.SQMFe  1 0.009006730
# 13     Vaccinated.SQMFe vs Unvaccinated.FeSO4  1 0.010330115
# 14     Vaccinated.SQMFe vs Unvaccinated.SQMFe  1 0.009738457
# 15   Unvaccinated.FeSO4 vs Unvaccinated.SQMFe  1 0.008891723
#      F.Model        R2 p.value p.adjusted sig
# 1   2.689254 0.2775502   0.026      0.390    
# 2   3.963972 0.3313257   0.010      0.150    
# 3   6.719919 0.4565187   0.014      0.210    
# 4  10.870482 0.5760575   0.005      0.075    
# 5   4.368160 0.3531778   0.007      0.105    
# 6   3.155243 0.3107009   0.011      0.165    
# 7   5.669940 0.4475112   0.003      0.045   .
# 8   8.714919 0.5545634   0.012      0.180    
# 9   2.596337 0.2705550   0.020      0.300    
# 10  2.585507 0.2442497   0.049      0.735    
# 11  3.027039 0.2745106   0.057      0.855    
# 12  1.879894 0.1902747   0.098      1.000    
# 13  3.796136 0.3218118   0.008      0.120    
# 14  2.756232 0.2562451   0.009      0.135    
# 15  2.779442 0.2578466   0.030      0.450    
