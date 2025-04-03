#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=download_single_gse
#SBATCH --output=logs/download_rnaseq/download_single_gse_%j.out
#SBATCH --error=logs/download_rnaseq/download_single_gse_%j.err
#SBATCH --tmp=100G

# Script to run a single GSE ID download using SLURM

# Check if GSE ID is provided
if [ $# -ne 1 ]; then
    echo "Usage: sbatch $0 <GSE_ID>"
    echo "Example: sbatch $0 GSE275540"
    exit 1
fi

GSE_ID=$1

# Load required modules
module load sra-toolkit
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Create logs directory if it doesn't exist
mkdir -p logs/download_rnaseq

# Create a temporary directory for this job
TMP_DIR=$SLURM_TMPDIR/download_${GSE_ID}
mkdir -p $TMP_DIR
cd $TMP_DIR

# Copy the R script and .env file to the temporary directory
cp $SLURM_SUBMIT_DIR/download_data.R .
cp $SLURM_SUBMIT_DIR/.env .

# Run the R script
Rscript download_data.R $GSE_ID

# Create the proper directory structure in the permanent location
GSE_DIR=$SLURM_SUBMIT_DIR/${GSE_ID}
mkdir -p $GSE_DIR/samples

# Copy files to their correct locations
if [ -d "samples" ]; then
    for gsm_dir in samples/*; do
        if [ -d "$gsm_dir" ]; then
            gsm_id=$(basename "$gsm_dir")
            mkdir -p "$GSE_DIR/samples/$gsm_id"
            cp -r "$gsm_dir"/* "$GSE_DIR/samples/$gsm_id/"
        fi
    done
fi

# Copy any other files from the GSE directory
if [ -d "$GSE_ID" ]; then
    cp -r "$GSE_ID"/* "$GSE_DIR/"
fi

# Clean up only after successful copying
cd $SLURM_SUBMIT_DIR
rm -rf $TMP_DIR

echo "Download and conversion complete for $GSE_ID" 