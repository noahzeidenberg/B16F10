#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=00:30:00     # Less time needed for normalization
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G           # Less memory needed for normalization
#SBATCH --job-name=normalization
#SBATCH --output=logs/normalization/normalization_%A_%a.out
#SBATCH --error=logs/normalization/normalization_%A_%a.err
#SBATCH --array=1-60%5     # Will process up to 60 GSEs, 5 at a time.

# Enable debugging
set -x

# Load required modules
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get absolute paths
SUBMIT_DIR=$(realpath $SLURM_SUBMIT_DIR)
WORK_DIR=$SUBMIT_DIR

# Debug information
echo "=== Debug Information ==="
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "Current working directory: $(pwd)"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "=== End Debug Information ==="

# Create output directories
echo "Creating output directories..."
mkdir -p $SUBMIT_DIR/results/normalization
mkdir -p $SUBMIT_DIR/logs/normalization

# Get the GSE ID for this array task
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)
#GSE_ID=${1:-GSE90885}
echo "Processing GSE: $GSE_ID"

# Run the R script for this specific GSE
echo "=== Starting R script execution for $GSE_ID ==="
cd $SUBMIT_DIR
Rscript scripts/normalize_counts.R $GSE_ID
R_EXIT_CODE=$?
echo "R script exit code: $R_EXIT_CODE"

echo "Normalization complete for GSE: $GSE_ID" 