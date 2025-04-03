#!/bin/bash
# Script to run a single GSE ID download without using SLURM array job

# Check if GSE ID is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <GSE_ID>"
    echo "Example: $0 GSE275540"
    exit 1
fi

GSE_ID=$1

# Create logs directory if it doesn't exist
mkdir -p logs/download_rnaseq

# Run the R script directly
Rscript download_data.R $GSE_ID

echo "Download and conversion complete for $GSE_ID" 