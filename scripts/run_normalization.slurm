#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=00:30:00     # Less time needed for normalization
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G           # Less memory needed for normalization
#SBATCH --job-name=normalization
#SBATCH --output=logs/normalization/normalization_%j.out
#SBATCH --error=logs/normalization/normalization_%j.err

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
echo "=== End Debug Information ==="

# Create output directories
echo "Creating output directories..."
mkdir -p $SUBMIT_DIR/results/normalization
mkdir -p $SUBMIT_DIR/logs/normalization

# Run the R script
echo "=== Starting R script execution ==="
cd $SUBMIT_DIR
Rscript scripts/normalize_counts.R
R_EXIT_CODE=$?
echo "R script exit code: $R_EXIT_CODE"

echo "Normalization complete for all datasets" 