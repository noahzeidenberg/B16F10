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
echo "This script will skip GSEs that fail and exclude them from the final results."

# Run the R script and capture its exit status
Rscript scripts/aggregate_pathway_results.R
EXIT_STATUS=$?

# Check if the script ran successfully
if [ $EXIT_STATUS -eq 0 ]; then
    echo "Pathway aggregation analysis completed successfully at $(date)"
    
    # Check if any GSEs were skipped
    if [ -f "results/aggregated_pathway_analysis/skipped_gses.txt" ]; then
        SKIPPED_COUNT=$(wc -l < results/aggregated_pathway_analysis/skipped_gses.txt)
        echo "Note: $SKIPPED_COUNT GSEs were skipped due to errors or invalid results."
        echo "See results/aggregated_pathway_analysis/skipped_gses.txt for details."
    fi
    
    # Report the number of included GSEs
    if [ -f "results/aggregated_pathway_analysis/included_gses.txt" ]; then
        INCLUDED_COUNT=$(wc -l < results/aggregated_pathway_analysis/included_gses.txt)
        echo "Successfully included $INCLUDED_COUNT GSEs in the analysis."
    fi
else
    echo "Pathway aggregation analysis failed with exit code $EXIT_STATUS at $(date)"
    echo "Check the error log for details."
    exit 1
fi

exit 0 