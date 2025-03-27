#!/bin/bash
#SBATCH --account=def-username
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=RNAseq_analysis
#SBATCH --output=RNAseq_analysis_%j.out
#SBATCH --error=RNAseq_analysis_%j.err

# Load required modules
module load r/4.2.2
module load gcc/9.3.0

# Set working directory
cd $SLURM_SUBMIT_DIR

# Run R script
Rscript R/RNAseq_analysis.R 