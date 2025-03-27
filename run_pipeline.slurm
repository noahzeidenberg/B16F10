#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=16000M
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=b16f10_pipeline
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Load modules in correct order
module purge  # Clear any previously loaded modules
module load StdEnv/2020
module load python
module load r
module load star
module load subread

# Load StdEnv/2023 modules separately
module load StdEnv/2023
module load sra-toolkit
module load edirect

# Activate virtual environment and verify R packages
source .venv/bin/activate

# Verify R packages are available
Rscript -e "if (!require('GEOquery')) install.packages('GEOquery', repos='https://cloud.r-project.org/')"
Rscript -e "if (!require('SRAdb')) install.packages('SRAdb', repos='https://cloud.r-project.org/')"

# Verify Python packages are available
python -c "
import sys
required_packages = ['numpy', 'pandas', 'scipy']
missing_packages = []
for package in required_packages:
    try:
        __import__(package)
    except ImportError:
        missing_packages.append(package)
if missing_packages:
    print(f'Missing packages: {missing_packages}')
    sys.exit(1)
print('All required Python packages are available')
"

# Run the pipeline with parallel processing
if ! python b16f10_pipeline.py --all --parallel $SLURM_CPUS_PER_TASK; then
    echo "Pipeline failed with an error. Check logs for details."
    exit 1
fi

# Check job state and completion
if [ "$SLURM_JOB_STATE" = "TIMEOUT" ]; then
    echo "Job timed out. Checking completion status..."
    if ! python check_completion.py; then
        echo "Job timed out and some accessions are incomplete. Resubmitting..."
        sbatch $0
    fi
else
    echo "Job completed or failed. Check logs for details."
fi