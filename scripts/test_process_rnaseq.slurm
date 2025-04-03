#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=04:30:00     # could maybe work with 4 hours, but using 6 just in case
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=test_process_rnaseq
#SBATCH --output=logs/test_process_rnaseq/test_process_rnaseq_%j.out
#SBATCH --error=logs/test_process_rnaseq/test_process_rnaseq_%j.err

# Enable debugging
set -x

# Load required modules
module load star
module load fastqc
module load fastp
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Set the GSE ID to test
GSE_ID="GSE287957"

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
GSE_DIR=$SUBMIT_DIR/${GSE_ID}
mkdir -p $GSE_DIR/results/{fastqc,trimmed,alignment,counts}

# Run the R script
echo "=== Starting R script execution ==="
cd $SUBMIT_DIR
Rscript scripts/process_rnaseq.R $GSE_ID
R_EXIT_CODE=$?
echo "R script exit code: $R_EXIT_CODE"

echo "Processing complete for $GSE_ID" 