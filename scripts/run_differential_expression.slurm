#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --job-name=diff_exp
#SBATCH --output=logs/differential_expression/diff_exp_%A_%a.out
#SBATCH --error=logs/differential_expression/diff_exp_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-100%5     # Will process up to 100 GSEs, 5 at a time

# Enable debugging
set -x

# Load required modules
module load r

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
mkdir -p $SUBMIT_DIR/logs/differential_expression

# Get the GSE ID for this array task
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)

echo "Starting differential expression analysis for $GSE_ID"
echo "Current working directory: $(pwd)"
echo "Checking if GSE directory exists: $GSE_ID"
if [ -d "$GSE_ID" ]; then
  echo "GSE directory exists"
  echo "Checking if normalization directory exists: $GSE_ID/results/normalization"
  if [ -d "$GSE_ID/results/normalization" ]; then
    echo "Normalization directory exists"
    echo "Files in normalization directory:"
    ls -la $GSE_ID/results/normalization/
  else
    echo "Normalization directory does not exist"
    exit 1
  fi
else
  echo "GSE directory does not exist"
  exit 1
fi

# Run the R script
echo "Running differential_expression.R for $GSE_ID"
cd $SUBMIT_DIR
Rscript scripts/differential_expression.R $GSE_ID
R_EXIT_CODE=$?
echo "R script exit code: $R_EXIT_CODE"

# Check if the script ran successfully
if [ $R_EXIT_CODE -eq 0 ]; then
  echo "Differential expression analysis completed successfully for $GSE_ID"
else
  echo "Differential expression analysis failed for $GSE_ID"
  exit 1
fi 