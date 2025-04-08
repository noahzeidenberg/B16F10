#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --job-name=diff_exp
#SBATCH --output=logs/differential_expression/diff_exp_%A_%a.out
#SBATCH --error=logs/differential_expression/diff_exp_%A_%a.err
#SBATCH --time=00:15:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-1

# Load required modules
module load r

# Set working directory
cd ~/scratch/B16F10

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get GSE ID from command line or use default
GSE_ID=${1:-GSE90885}

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
  fi
else
  echo "GSE directory does not exist"
fi

# Create output directory for logs if it doesn't exist
mkdir -p logs/differential_expression

# Run the R script
echo "Running differential_expression.R for $GSE_ID"
Rscript scripts/differential_expression.R $GSE_ID

# Check if the script ran successfully
if [ $? -eq 0 ]; then
  echo "Differential expression analysis completed successfully for $GSE_ID"
  
  # Check if results were generated
  if [ -f "$GSE_ID/results/differential_expression/${GSE_ID}_de_results.rds" ]; then
    echo "Results file created successfully"
  else
    echo "Warning: Results file not found. Analysis may have completed but no results were generated."
  fi
else
  echo "Differential expression analysis failed for $GSE_ID"
  exit 1
fi 