#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=process_rnaseq
#SBATCH --output=process_rnaseq_%A_%a.out
#SBATCH --error=process_rnaseq_%A_%a.err
#SBATCH --array=1-72
#SBATCH --tmp=200G

# Load required modules
module load star/2.7.9a
module load fastqc/0.11.9
module load fastp/0.23.2
module load r/4.2.1
module load gcc/9.3.0

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get the GSE ID for this array job
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)

# Create a temporary directory for this job
TMP_DIR=$SLURM_TMPDIR/process_${GSE_ID}
mkdir -p $TMP_DIR
cd $TMP_DIR

# Copy the R script to the temporary directory
cp $SLURM_SUBMIT_DIR/scripts/process_rnaseq.R .

# Copy the FASTQ files from the previous step
cp -r $SLURM_SUBMIT_DIR/results/${GSE_ID}/* .

# Run the R script
Rscript process_rnaseq.R $GSE_ID

# Copy results back to the submission directory
RESULTS_DIR=$SLURM_SUBMIT_DIR/results/${GSE_ID}
mkdir -p $RESULTS_DIR/{fastqc,trimmed,alignment,counts}
cp -r fastqc/* $RESULTS_DIR/fastqc/
cp -r trimmed/* $RESULTS_DIR/trimmed/
cp -r alignment/* $RESULTS_DIR/alignment/
cp -r counts/* $RESULTS_DIR/counts/

# Clean up
cd $SLURM_SUBMIT_DIR
rm -rf $TMP_DIR 