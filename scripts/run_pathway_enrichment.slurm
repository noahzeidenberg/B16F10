#!/bin/bash
#SBATCH --job-name=pathway_enrichment
#SBATCH --output=logs/pathway/pathway_enrichment_%j.out
#SBATCH --error=logs/pathway/pathway_enrichment_%j.err
#SBATCH --time=0:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-60%5 # 60 is the number of GSE IDs, 5 is the number of tasks to run in parallel

# Get the GSE ID for this array task
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)

# Load required modules
module load r

# load environment
source .venv/bin/activate

# Set working directory
WORK_DIR="$HOME/scratch/B16F10"
cd $WORK_DIR

# Create logs directory if it doesn't exist
mkdir -p logs

# Run pathway enrichment analysis
echo "Starting pathway enrichment analysis for $GSE_ID at $(date)"
Rscript scripts/pathway_enrichment.R $GSE_ID

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Pathway enrichment analysis completed successfully for $GSE_ID at $(date)"
else
    echo "Pathway enrichment analysis failed for $GSE_ID at $(date)"
    exit 1
fi

exit 0 