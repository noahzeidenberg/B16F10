#!/bin/bash

# Set strict error handling
set -euo pipefail

# Define the base directory (current directory)
BASE_DIR="$(pwd)"

# Create logs directory if it doesn't exist
mkdir -p "${BASE_DIR}/logs"

# Function to log messages
log_message() {
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${message}" | tee -a "${BASE_DIR}/logs/process_all.log"
}

# Extract RNA-seq datasets from GDS table
log_message "Extracting RNA-seq datasets from GDS table..."
Rscript -e '
gds_df <- read.csv("gds_table.csv")
rna_seq_datasets <- gds_df[grep("sequencing", gds_df$Type, ignore.case=TRUE), "Accession"]
write.csv(rna_seq_datasets, "rna_seq_datasets.csv", row.names=FALSE)
'

# Process each RNA-seq dataset
while IFS= read -r accession; do
    if [[ -n "$accession" ]]; then
        log_message "Processing dataset: $accession"
        Rscript R/RNAseq_analysis.R "$accession" 2>&1 | tee -a "${BASE_DIR}/logs/RNAseq_${accession}.log"
    fi
done < rna_seq_datasets.csv

log_message "All datasets processed" 