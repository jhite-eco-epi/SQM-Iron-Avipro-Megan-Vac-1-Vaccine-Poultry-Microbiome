library(vegan)
library(dplyr)
library(ggplot2)
# Load env vars and utility
readRenviron(".Renviron")
source("utils/metadata_utils.R")

# -------------------------------------------
# 1. Load Vaccination/Supplement Metadata
# -------------------------------------------
treatment_data <- load_clean_metadata()

# -------------------------------------------
# 2. Load Jaccard Distance Matrix
# -------------------------------------------
jaccard_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_jaccard_matrix.tsv")
jaccard_matrix <- read.table(jaccard_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# -------------------------------------------
# 3. Match Metadata with Distance Matrix Samples
# -------------------------------------------
# Determine common sample IDs between distance matrix and cleaned metadata
common_samples <- intersect(rownames(jaccard_matrix), treatment_data$sampleid)

# Reorder metadata to match distance matrix order
treatment_data <- treatment_data[match(common_samples, treatment_data$sampleid), ]

# Subset the distance matrix to include only these samples
jaccard_matrix_subset <- as.matrix(jaccard_matrix)[common_samples, common_samples]


# -------------------------------------------
# 4. Perform NMDS
# -------------------------------------------
nmds_global <- metaMDS(jaccard_matrix_subset, k = 2, trymax = 100)
scores_global <- as.data.frame(scores(nmds_global))
scores_global$SampleID <- rownames(scores_global)
scores_global <- merge(scores_global, treatment_data, by.x = "SampleID", by.y = "sampleid")

# -------------------------------------------
# 5. Prepare factors and create box-plot
# -------------------------------------------

# Factor conversions
scores_global$vaccine <- factor(scores_global$vaccine, levels = c("0", "1"))
scores_global$FeSO4   <- factor(scores_global$FeSO4,   levels = c("0", "1"))
scores_global$SQMFe   <- factor(scores_global$SQMFe,   levels = c("0", "1"))

# Readable vaccination label
scores_global$vaccine_label <- factor(scores_global$vaccine,
                                      levels = c("0", "1"),
                                      labels = c("Unvaccinated", "Vaccinated"))

# Single Supplement variable
plot_data <- scores_global %>%
  mutate(Supplement = dplyr::case_when(
    FeSO4 == "1" & SQMFe == "0" ~ "FeSO4",
    FeSO4 == "0" & SQMFe == "1" ~ "SQM® Iron",
    FeSO4 == "0" & SQMFe == "0" ~ "Control",
    TRUE ~ "Other"
  )) %>%
  filter(Supplement %in% c("Control", "FeSO4", "SQM® Iron"))

plot_data$Supplement <- factor(plot_data$Supplement,
                               levels = c("Control", "FeSO4", "SQM® Iron"))

# Plot
beta_boxplot <- ggplot(plot_data,
                       aes(x = vaccine_label, y = NMDS1, fill = Supplement)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.8),
              alpha = 0.7, size = 2, shape = 21, colour = "black") +
  labs(x = "Vaccination Status",
       y = "Ordination Score (NMDS 1)",
       fill = "Supplementation") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  scale_fill_viridis_d(option = "plasma")

# beta_boxplot + scale_fill_lancet()

# Show 1-D box-plot
print(beta_boxplot)

# -------------------------------------------
# 5b. 2-D NMDS scatter plot
# -------------------------------------------

# library(ggsci)   # Lancet palette

nmds_scatter <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2,
                                      colour = Supplement,
                                      shape  = vaccine_label)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(linetype = vaccine_label,
                   group = interaction(vaccine_label, Supplement)),
               type = "t", size = 1) +
  scale_shape_manual(name = "Vaccination",
                     values = c("Unvaccinated" = 16, "Vaccinated" = 17)) +
  scale_linetype_manual(name = "Vaccination",
                        values = c("Unvaccinated" = "dashed", "Vaccinated" = "solid")) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  guides(colour = guide_legend(override.aes = list(shape = 15, linetype = "blank", size = 5)),
         linetype = guide_legend(override.aes = list(colour = "black")),
         shape = guide_legend(override.aes = list(linetype = "blank")))

nmds_scatter + scale_colour_viridis_d(option = "plasma")
# print(nmds_scatter)

# -------------------------------------------
# 6. PERMANOVA: Vaccination × Supplement full factorial
# -------------------------------------------

# Build metadata frame aligned with distance matrix labels
meta_for_adonis <- scores_global %>%
  mutate(Supplement = dplyr::case_when(
    FeSO4 == "1" & SQMFe == "0" ~ "FeSO4",
    FeSO4 == "0" & SQMFe == "1" ~ "SQM® Iron",
    FeSO4 == "0" & SQMFe == "0" ~ "Control",
    TRUE ~ "Other"
  )) %>%
  select(SampleID, vaccine_label, Supplement)

# Ensure factors
meta_for_adonis$vaccine_label <- factor(meta_for_adonis$vaccine_label,
                                        levels = c("Unvaccinated", "Vaccinated"))
meta_for_adonis$Supplement <- factor(meta_for_adonis$Supplement,
                                     levels = c("Control", "FeSO4", "SQM® Iron", "Other"))

set.seed(123)
permanova_results <- vegan::adonis2(jaccard_matrix_subset ~ vaccine_label * Supplement,
                                    data = meta_for_adonis,
                                    permutations = 999,
                                    by = "terms")

cat("\n=====================\nPERMANOVA (Jaccard)\n=====================\n")
print(permanova_results)

# -------------------------------
# Save PERMANOVA Results as CSV and Tables
# -------------------------------

library(gridExtra)
library(grid)
library(gtable)

# Format PERMANOVA results for output
jaccard_permanova_results <- data.frame(
  Effect = rownames(permanova_results),
  Df = permanova_results$Df,
  Sum_Sq = round(permanova_results$SumOfSqs, 4),
  R2 = round(permanova_results$R2, 3),
  F_statistic = round(permanova_results$F, 3),
  p_value = round(permanova_results$`Pr(>F)`, 5),
  stringsAsFactors = FALSE
)

# Clean up effect names for better display
jaccard_permanova_results$Effect <- gsub("vaccine_label", "Vaccine", jaccard_permanova_results$Effect)
jaccard_permanova_results$Effect <- gsub("Supplement", "Iron Supplement", jaccard_permanova_results$Effect)
jaccard_permanova_results$Effect <- gsub("Vaccine:Iron Supplement", "Vaccine × Iron Supplement", jaccard_permanova_results$Effect)

# Save CSV file
write.csv(jaccard_permanova_results, "BetaAnalysis/jaccard_permanova_stats.csv", row.names = FALSE)
cat("PERMANOVA results saved as CSV: BetaAnalysis/jaccard_permanova_stats.csv\n")

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
jaccard_table_result <- create_permanova_table(jaccard_permanova_results, "Jaccard PERMANOVA")
jaccard_table_grob <- jaccard_table_result$table_grob
jaccard_display_data <- jaccard_table_result$display_data

# Save tables
tiff_height <- max(3, nrow(jaccard_display_data) * 0.25 + 1)
tiff("BetaAnalysis/jaccard_permanova_table.tiff", 
     width = 10, height = tiff_height, units = "in", res = 300)
grid.draw(jaccard_table_grob)
dev.off()

png("BetaAnalysis/jaccard_permanova_table.png", 
    width = 10, height = tiff_height, units = "in", res = 150)
grid.draw(jaccard_table_grob)
dev.off()

cat("PERMANOVA table saved as: BetaAnalysis/jaccard_permanova_table.tiff/.png\n")

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# vegan::adonis2(formula = jaccard_matrix_subset ~ vaccine_label * Supplement, data = meta_for_adonis, permutations = 999, by = "terms")
#                          Df SumOfSqs      R2      F Pr(>F)    
# vaccine_label             1   0.3414 0.05378 1.7320  0.002 ** 
# Supplement                2   0.7985 0.12579 2.0256  0.001 ***
# vaccine_label:Supplement  2   0.6749 0.10630 1.7119  0.001 ***
# Residual                 23   4.5336 0.71413                  
# Total                    28   6.3485 1.00000                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# -------------------------------------------
# 7. Pairwise PERMANOVA
# -------------------------------------------

library(pairwiseAdonis)
# Combine factors to a single group variable
meta_for_adonis$Group <- interaction(meta_for_adonis$vaccine_label, meta_for_adonis$Supplement)
pairwise_results <- pairwise.adonis(jaccard_matrix_subset, meta_for_adonis$Group,
                                   p.adjust.m = "bonferroni")
cat("\n=====================\nPairwise PERMANOVA\n=====================\n")
print(pairwise_results)

#                                        pairs Df   SumsOfSqs  F.Model        R2 p.value p.adjusted
# 1  Unvaccinated.Control vs Vaccinated.Control  1 0.005905551 1.810619 0.2055042   0.182      1.000
# 2    Unvaccinated.Control vs Vaccinated.FeSO4  1 0.011122345 3.800366 0.3220549   0.021      0.315
# 3    Unvaccinated.Control vs Vaccinated.SQMFe  1 0.009136870 2.993146 0.2722738   0.052      0.780
# 4  Unvaccinated.Control vs Unvaccinated.FeSO4  1 0.017183855 6.002371 0.4286682   0.005      0.075
# 5  Unvaccinated.Control vs Unvaccinated.SQMFe  1 0.012616366 3.646113 0.3130755   0.033      0.495
# 6      Vaccinated.Control vs Vaccinated.FeSO4  1 0.006067849 3.050763 0.3035355   0.011      0.165
# 7      Vaccinated.Control vs Vaccinated.SQMFe  1 0.003967660 1.860218 0.2099517   0.041      0.615
# 8    Vaccinated.Control vs Unvaccinated.FeSO4  1 0.007475577 3.901578 0.3578911   0.032      0.480
# 9    Vaccinated.Control vs Unvaccinated.SQMFe  1 0.005215232 2.006817 0.2228109   0.099      1.000
# 10       Vaccinated.FeSO4 vs Vaccinated.SQMFe  1 0.002454062 1.265618 0.1365930   0.193      1.000
# 11     Vaccinated.FeSO4 vs Unvaccinated.FeSO4  1 0.005099824 2.915404 0.2670908   0.014      0.210
# 12     Vaccinated.FeSO4 vs Unvaccinated.SQMFe  1 0.004282835 1.825087 0.1857578   0.111      1.000
# 13     Vaccinated.SQMFe vs Unvaccinated.FeSO4  1 0.004025489 2.146680 0.2115648   0.045      0.675
# 14     Vaccinated.SQMFe vs Unvaccinated.SQMFe  1 0.002922443 1.181934 0.1287238   0.268      1.000
# 15   Unvaccinated.FeSO4 vs Unvaccinated.SQMFe  1 0.002476714 1.084926 0.1194205   0.366      1.000
#    sig
# 1     
# 2     
# 3     
# 4     
# 5     
# 6     
# 7     
# 8     
# 9     
# 10    
# 11    
# 12    
# 13    
# 14    
# 15    
