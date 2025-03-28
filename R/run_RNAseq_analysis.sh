#!/bin/bash
#SBATCH --account=def-username
#SBATCH --job-name=RNAseq_analysis
#SBATCH --output=RNAseq_output_%A_%a.log
#SBATCH --error=RNAseq_error_%A_%a.log
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:15:00
#SBATCH --mem=16000M

# Exit on error
set -e

# Load the required modules
module load gcc
module load r
module load sratoolkit

# Path to your R script and data directory
R_SCRIPT=./R/RNAseq_analysis.R
DATA_DIR=./data
RESULTS_DIR=./results

# Create necessary directories
mkdir -p "$DATA_DIR" "$RESULTS_DIR"

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "Error: R script not found at $R_SCRIPT"
    exit 1
fi

# Run the R script with the current array index
Rscript "$R_SCRIPT" --accession "${SLURM_ARRAY_TASK_ID}" || {
    echo "Error: R script failed for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
}
