#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=00:20:00     # Longer time for processing all datasets
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G          # More memory for processing all datasets
#SBATCH --job-name=batch_correction
#SBATCH --output=logs/batch_correction/batch_correction_%j.out
#SBATCH --error=logs/batch_correction/batch_correction_%j.err

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
mkdir -p $SUBMIT_DIR/results/batch_correction
mkdir -p $SUBMIT_DIR/logs/batch_correction

# Run the R script
echo "=== Starting R script execution ==="
cd $SUBMIT_DIR
Rscript scripts/batch_correction.R
R_EXIT_CODE=$?
echo "R script exit code: $R_EXIT_CODE"

echo "Batch correction complete" 