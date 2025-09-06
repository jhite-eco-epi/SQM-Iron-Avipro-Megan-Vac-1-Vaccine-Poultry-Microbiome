#!/usr/bin/env bash
set -euo pipefail

# QIIME 2 taxonomy classification + rarefied taxa-barplot export (2024.10)
#
# Usage:
#   bash scripts/classify_taxonomy_2024_10.sh \
#     -r data/rep-seqs.qza \
#     -R data/rarefied_table.qza \
#     [-m data/sample-metadata.tsv] \
#     [-c path/to/classifier.qza | -u https://classifier-url.qza] \
#     [-e qiime2-amplicon-2024.10] \
#     [-o data]
#
# Notes:
# - Provide either a local classifier (-c) or a download URL (-u). If both are given, -c wins.
# - The script activates the specified conda environment and runs QIIME 2 commands.
# - If -R is provided, uses an existing rarefied table and skips rarefaction.
# - Outputs taxonomy-updated.qza/.qzv and a taxa-barplot export under OUTPUT_DIR.

ENV_NAME="qiime2-amplicon-2024.10"
REP_SEQS=""
METADATA_TSV=""
CLASSIFIER_QZA=""
CLASSIFIER_URL=""
OUTPUT_DIR="data"

log() { printf "[%%s] %%s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die() { echo "Error: $*" >&2; exit 1; }

RA_TABLE_EXISTING=""

while getopts ":e:r:m:c:u:o:R:h" opt; do
  case $opt in
    e) ENV_NAME="$OPTARG" ;;
    r) REP_SEQS="$OPTARG" ;;
    m) METADATA_TSV="$OPTARG" ;;
    c) CLASSIFIER_QZA="$OPTARG" ;;
    u) CLASSIFIER_URL="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    R) RA_TABLE_EXISTING="$OPTARG" ;;
    h) sed -n '1,80p' "$0"; exit 0 ;;
    :) die "Option -$OPTARG requires an argument" ;;
    \?) die "Invalid option: -$OPTARG" ;;
  esac
done

[[ -n "$REP_SEQS" ]] || die "-r rep-seqs.qza is required"
[[ -n "$RA_TABLE_EXISTING" ]] || die "-R rarefied_table.qza is required"

mkdir -p "$OUTPUT_DIR"

# Prepare classifier
if [[ -z "$CLASSIFIER_QZA" && -z "$CLASSIFIER_URL" ]]; then
  die "Provide either -c classifier.qza or -u classifier URL"
fi
if [[ -z "$CLASSIFIER_QZA" ]]; then
  CLASSIFIER_QZA="$OUTPUT_DIR/classifier.qza"
  log "Downloading classifier from URL: $CLASSIFIER_URL"
  curl -L --fail --retry 3 -o "$CLASSIFIER_QZA" "$CLASSIFIER_URL" || die "Failed to download classifier"
fi
[[ -f "$CLASSIFIER_QZA" ]] || die "Classifier not found: $CLASSIFIER_QZA"

# Activate conda env with QIIME 2 (temporarily disable nounset for conda hooks)
if command -v conda >/dev/null 2>&1; then
  set +u
  # shellcheck disable=SC1090
  eval "$(conda shell.bash hook)"
  conda activate "$ENV_NAME" || die "Failed to activate conda env: $ENV_NAME"
  set -u
else
  die "conda not found in PATH"
fi

log "QIIME 2 info:"
qiime info | sed -n '1,40p' || die "qiime not available"

# Classify sequences
log "Classifying taxonomy with classifier: $CLASSIFIER_QZA"
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_QZA" \
  --i-reads "$REP_SEQS" \
  --o-classification "$OUTPUT_DIR/taxonomy-updated.qza"

qiime metadata tabulate \
  --m-input-file "$OUTPUT_DIR/taxonomy-updated.qza" \
  --o-visualization "$OUTPUT_DIR/taxonomy-updated.qzv"

[[ -f "$RA_TABLE_EXISTING" ]] || die "Provided -R rarefied table not found: $RA_TABLE_EXISTING"
log "Using existing rarefied table: $RA_TABLE_EXISTING"
RA_TABLE="$RA_TABLE_EXISTING"

# Create taxa barplot on the rarefied table
log "Building taxa barplot visualization"
BARPLOT_ARGS=(
  --i-table "$RA_TABLE"
  --i-taxonomy "$OUTPUT_DIR/taxonomy-updated.qza"
  --o-visualization "$OUTPUT_DIR/taxa-barplots-rarefied-updated.qzv"
)
if [[ -n "$METADATA_TSV" ]]; then
  BARPLOT_ARGS+=( --m-metadata-file "$METADATA_TSV" )
fi
qiime taxa barplot "${BARPLOT_ARGS[@]}"

# Export barplot data (level-*.csv) for downstream donut scripts
log "Exporting taxa barplot data"
mkdir -p "$OUTPUT_DIR/exported_taxa_barplot_data"
qiime tools export \
  --input-path "$OUTPUT_DIR/taxa-barplots-rarefied-updated.qzv" \
  --output-path "$OUTPUT_DIR/exported_taxa_barplot_data"

log "Done. Key outputs:"
echo "  Taxonomy:       $OUTPUT_DIR/taxonomy-updated.qza (viz: $OUTPUT_DIR/taxonomy-updated.qzv)"
echo "  Rarefied table: $RA_TABLE"
echo "  Barplots:       $OUTPUT_DIR/taxa-barplots-rarefied-updated.qzv"
echo "  Exported data:  $OUTPUT_DIR/exported_taxa_barplot_data (level-2.csv, level-7.csv, etc.)"


