#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --job-name=create_mapping
#SBATCH --output=logs/create_mapping/create_mapping_%A_%a.out
#SBATCH --error=logs/create_mapping/create_mapping_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-1

# Load required modules
module load r

# Set working directory
cd ~/scratch/B16F10

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get GSE ID from command line or use default
GSE_ID=${1:-GSE287957}

echo "Creating SRA to GEO mapping for $GSE_ID"
echo "Current working directory: $(pwd)"

# Run the R script
echo "Running create_sra_geo_mapping.R for $GSE_ID"
Rscript scripts/create_sra_geo_mapping.R $GSE_ID

# Check if the script ran successfully
if [ $? -eq 0 ]; then
  echo "SRA to GEO mapping creation completed successfully for $GSE_ID"
else
  echo "SRA to GEO mapping creation failed for $GSE_ID"
  exit 1
fi 