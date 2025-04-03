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

# Load required modules again in case they were unloaded in temporary directory
module load sra-toolkit
module load r
module load gcc

# Run the R script
Rscript download_data.R $GSE_ID

# Create the proper directory structure in the permanent location
GSE_DIR=$SLURM_SUBMIT_DIR/${GSE_ID}
mkdir -p $GSE_DIR/samples

# Debug information
echo "Current directory: $(pwd)"
echo "Temporary directory: $TMP_DIR"
echo "GSE directory: $GSE_DIR"
echo "Contents of current directory:"
ls -la
echo "Contents of samples directory (if it exists):"
if [ -d "samples" ]; then
    ls -la samples
else
    echo "Samples directory not found in current directory"
fi

# Copy files to their correct locations
if [ -d "samples" ]; then
    echo "Copying sample directories..."
    for gsm_dir in samples/*; do
        if [ -d "$gsm_dir" ]; then
            gsm_id=$(basename "$gsm_dir")
            echo "Copying $gsm_id to $GSE_DIR/samples/$gsm_id"
            mkdir -p "$GSE_DIR/samples/$gsm_id"
            cp -rv "$gsm_dir"/* "$GSE_DIR/samples/$gsm_id/"
        fi
    done
else
    echo "Samples directory not found in $TMP_DIR"
fi

# Copy any other files from the GSE directory
if [ -d "$GSE_ID" ]; then
    echo "Copying GSE directory contents..."
    cp -rv "$GSE_ID"/* "$GSE_DIR/"
else
    echo "GSE directory not found in $TMP_DIR"
fi

# Verify the copy was successful
echo "Verifying files in permanent location:"
if [ -d "$GSE_DIR/samples" ]; then
    echo "Contents of $GSE_DIR/samples:"
    ls -la "$GSE_DIR/samples"
else
    echo "Samples directory not found in permanent location"
fi

# Clean up only after successful copying
cd $SLURM_SUBMIT_DIR
rm -rf $TMP_DIR

echo "Download and conversion complete for $GSE_ID" 