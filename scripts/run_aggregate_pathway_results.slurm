#!/bin/bash
#SBATCH --job-name=aggregate_pathways
#SBATCH --output=logs/pathway/aggregate_pathways_%j.out
#SBATCH --error=logs/pathway/aggregate_pathways_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Load required modules
module load r

# load environment
source .venv/bin/activate

# Set working directory
WORK_DIR="$HOME/scratch/B16F10"
cd $WORK_DIR

# Create logs directory if it doesn't exist
mkdir -p logs/pathway

# Run pathway aggregation analysis
echo "Starting pathway aggregation analysis at $(date)"
Rscript scripts/aggregate_pathway_results.R

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Pathway aggregation analysis completed successfully at $(date)"
else
    echo "Pathway aggregation analysis failed at $(date)"
    exit 1
fi

exit 0 