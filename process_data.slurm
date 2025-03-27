#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=16000M
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=process_geo
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Load StdEnv/2020 and required modules
module load StdEnv/2020
module load python
module load r
module load star
module load subread

# Activate virtual environment
source .venv/bin/activate

# Run only the processing portion
python process_geo.py --all --parallel $SLURM_CPUS_PER_TASK

# Check completion status
if [ $? -eq 0 ]; then
    echo "Processing completed successfully"
else
    echo "Processing failed. Check logs for details."
    exit 1
fi 