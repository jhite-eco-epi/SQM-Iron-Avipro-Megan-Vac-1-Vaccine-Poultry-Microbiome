# Utility: load and clean metadata for the poultry vaccination × iron-supplement experiment
#
# * Reads PoulVaccMetadata.tsv from the base data path
# * Removes user-specified sample IDs (e.g., technical replicates CN415/CN416)
# * Keeps only one row per Group × bird combination (first occurrence)
#
# Usage inside analysis scripts:
#   source("../utils/metadata_utils.R")         # path relative to script folder
#   md <- load_clean_metadata()                  # returns cleaned data.frame
#   223 is an outlier, 414/415/416 are technical replicates on bird 1 group D
# -----------------------------------------------------------------------------

load_clean_metadata <- function(base_path = Sys.getenv("BASE_DATA_PATH"),
                                remove_ids = c("CN223", "CN414", "CN415", "CN416")) {
  if (base_path == "") {
    stop("BASE_DATA_PATH environment variable is not set.")
  }

  metadata_file <- file.path(base_path, "PoulVaccMetadata.tsv")
  md <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  # Remove explicit sample IDs (duplicates/outliers)
  md <- md[!md$sampleid %in% remove_ids, ]

  md$sampleid <- as.character(md$sampleid)
  return(md)
} 