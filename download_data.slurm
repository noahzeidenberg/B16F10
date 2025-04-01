#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=download_rnaseq
#SBATCH --output=download_rnaseq_%A_%a.out
#SBATCH --error=download_rnaseq_%A_%a.err
#SBATCH --array=1-72 

# Load required modules
module load sra-toolkit
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get the GSE ID for this array job
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)

# Create a temporary directory for this job
TMP_DIR=$SLURM_TMPDIR/download_${GSE_ID}
mkdir -p $TMP_DIR
cd $TMP_DIR

# Copy the R script to the temporary directory
cp $SLURM_SUBMIT_DIR/download_data.R .

# Run the R script
Rscript download_data.R $GSE_ID

# Copy results back to the submission directory
RESULTS_DIR=$SLURM_SUBMIT_DIR/results/${GSE_ID}
mkdir -p $RESULTS_DIR
cp -r * $RESULTS_DIR/

# Clean up
cd $SLURM_SUBMIT_DIR
rm -rf $TMP_DIR 