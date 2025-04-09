#!/bin/bash
#SBATCH --job-name=pathway_enrichment
#SBATCH --output=logs/pathway/pathway_enrichment_%A_%a.out
#SBATCH --error=logs/pathway/pathway_enrichment_%A_%a.err
#SBATCH --time=0:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-60%5 # 60 is the number of GSE IDs, 5 is the number of tasks to run in parallel

# Exit on error
set -e

# Print debugging information
echo "Job Array ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Working Directory: $PWD"

# Check if GSE IDs file exists
if [ ! -f "rna_seq_gse_ids.txt" ]; then
    echo "Error: rna_seq_gse_ids.txt not found. Running create_gse_list.R..."
    module load r
    Rscript scripts/create_gse_list.R
fi

# Get the GSE ID for this array task
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)
if [ -z "$GSE_ID" ]; then
    echo "Error: Could not get GSE ID for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing GSE ID: $GSE_ID"

# Load required modules
module load r

# Set working directory to scratch space
SCRATCH_DIR="$SCRATCH/B16F10"
if [ ! -d "$SCRATCH_DIR" ]; then
    echo "Error: Scratch directory $SCRATCH_DIR does not exist"
    exit 1
fi

cd $SCRATCH_DIR

# Create logs directory if it doesn't exist
mkdir -p logs/pathway

# Check if differential expression results exist
DE_DIR="$SCRATCH_DIR/$GSE_ID/results/differential_expression"
if [ ! -d "$DE_DIR" ]; then
    echo "Error: Differential expression directory not found at $DE_DIR"
    exit 1
fi

# Run pathway enrichment analysis
echo "Starting pathway enrichment analysis for $GSE_ID at $(date)"
echo "Looking for differential expression results in: $DE_DIR"
ls -l "$DE_DIR"  # List contents for debugging

Rscript scripts/pathway_enrichment.R $GSE_ID

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Pathway enrichment analysis completed successfully for $GSE_ID at $(date)"
else
    echo "Pathway enrichment analysis failed for $GSE_ID at $(date)"
    exit 1
fi

exit 0 