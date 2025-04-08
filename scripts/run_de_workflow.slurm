#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --job-name=de_workflow
#SBATCH --output=logs/de_workflow/de_workflow_%A_%a.out
#SBATCH --error=logs/de_workflow/de_workflow_%A_%a.err
#SBATCH --time=06:00:00
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
GSE_ID=${1:-GSE287957}

echo "Starting differential expression workflow for $GSE_ID"
echo "Current working directory: $(pwd)"

# Step 1: Create SRA to GEO mapping
echo "Step 1: Creating SRA to GEO mapping for $GSE_ID"
Rscript scripts/create_sra_geo_mapping.R $GSE_ID

# Check if the mapping creation was successful
if [ $? -ne 0 ]; then
  echo "Failed to create SRA to GEO mapping for $GSE_ID"
  exit 1
fi

# Step 2: Run differential expression analysis
echo "Step 2: Running differential expression analysis for $GSE_ID"
Rscript scripts/differential_expression.R $GSE_ID

# Check if the differential expression analysis was successful
if [ $? -ne 0 ]; then
  echo "Failed to run differential expression analysis for $GSE_ID"
  exit 1
fi

echo "Differential expression workflow completed successfully for $GSE_ID" 