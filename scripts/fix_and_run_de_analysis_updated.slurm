#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --job-name=fix_de_updated
#SBATCH --output=logs/fix_de_analysis/fix_de_updated_%j.out
#SBATCH --error=logs/fix_de_analysis/fix_de_updated_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Load required modules
module load r

# Set working directory
cd ~/scratch/B16F10

# Make the script executable
chmod +x scripts/fix_and_run_de_analysis.sh

# Run the script
./scripts/fix_and_run_de_analysis.sh 