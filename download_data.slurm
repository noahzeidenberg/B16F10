#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download_geo
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Load StdEnv/2023 and required modules
module load StdEnv/2023
module load python
module load r
module load edirect
module load sra-toolkit

# Activate virtual environment
source .venv/bin/activate

# Run only the download portion
python download_geo.py --all --parallel $SLURM_CPUS_PER_TASK

# If successful, submit the processing job
if [ $? -eq 0 ]; then
    echo "Download completed successfully. Submitting processing job..."
    sbatch process_data.slurm
else
    echo "Download failed. Check logs for details."
    exit 1
fi 