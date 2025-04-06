#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --job-name=design_matrices
#SBATCH --output=logs/design_matrices_%A_%a.out
#SBATCH --error=logs/design_matrices_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1

# Load required modules
module load r

# Set working directory
cd ~/scratch/B16F10

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get GSE ID from command line or use default
GSE_ID=${1:-GSE254073}

# Run the R script
Rscript scripts/generate_design_matrices.R $GSE_ID

# Check if the script ran successfully
if [ $? -eq 0 ]; then
  echo "Design matrix generation completed successfully for $GSE_ID"
else
  echo "Design matrix generation failed for $GSE_ID"
  exit 1
fi 