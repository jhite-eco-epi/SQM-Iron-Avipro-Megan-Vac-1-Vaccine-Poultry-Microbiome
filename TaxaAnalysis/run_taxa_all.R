#!/usr/bin/env Rscript
# run_taxa_all.R
# Orchestrates Level 2 (phylum) and Level 6 (genus) analyses through shared pipeline

readRenviron(".Renviron")

suppressPackageStartupMessages({
  library(readr)
})

source("TaxaAnalysis/taxa_shared.R")

base_dir <- Sys.getenv("BASE_DATA_PATH", unset = "")
if (base_dir == "") stop("BASE_DATA_PATH is not set. Please export BASE_DATA_PATH to the directory containing exported_taxa_barplot_data.")

csv_l2 <- file.path(base_dir, "exported_taxa_barplot_data", "level-2.csv")
csv_l6 <- file.path(base_dir, "exported_taxa_barplot_data", "level-6.csv")
stopifnot(file.exists(csv_l2), file.exists(csv_l6))

data_l2 <- read_csv(csv_l2, show_col_types = FALSE)
data_l6 <- read_csv(csv_l6, show_col_types = FALSE)

# Sync metadata for both
source("utils/metadata_utils.R")
md <- load_clean_metadata()
keep_l2 <- intersect(md$sampleid, data_l2$index)
keep_l6 <- intersect(md$sampleid, data_l6$index)
if (length(keep_l2) == 0 || length(keep_l6) == 0) stop("No overlapping samples between counts and cleaned metadata.")

data_l2 <- data_l2[data_l2$index %in% keep_l2, ]
data_l6 <- data_l6[data_l6$index %in% keep_l6, ]

data_l2 <- dplyr::left_join(data_l2, md[, c("sampleid", "vaccination", "FeSO4", "SQMFe")], by = c("index" = "sampleid"))
data_l6 <- dplyr::left_join(data_l6, md[, c("sampleid", "vaccination", "FeSO4", "SQMFe")], by = c("index" = "sampleid"))

dir.create("TaxaAnalysis", showWarnings = FALSE, recursive = TRUE)

cat("\n=== Running Level 2 (phylum) ===\n")
analyze_level(
  data = data_l2,
  rank = "phylum",
  output_img_prefix = "TaxaAnalysis/taxa_level2",
  table_prefix = "TaxaAnalysis/taxa_table_level2"
)

cat("\n=== Running Level 6 (genus) ===\n")
analyze_level(
  data = data_l6,
  rank = "genus",
  output_img_prefix = "TaxaAnalysis/taxa_level6",
  table_prefix = "TaxaAnalysis/taxa_table_level6"
)

cat("\nAll taxa outputs generated.\n")


