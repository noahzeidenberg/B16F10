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

# Function to check if all accessions are processed
check_completion() {
    python - <<EOF
import json
from pathlib import Path
import pandas as pd

metadata = pd.read_csv("${PATHS['metadata']}")
accessions = metadata['Accession'].unique()
all_completed = True

for acc in accessions:
    checkpoint_file = Path(f"checkpoints/{acc}_status.json")
    if not checkpoint_file.exists():
        print(f"Checkpoint file not found for {acc}")
        all_completed = False
        continue
    with open(checkpoint_file) as f:
        status = json.load(f)
        if not status.get("completed", False):
            print(f"Processing incomplete for {acc}")
            all_completed = False

exit(0 if all_completed else 1)
EOF
}

# Load required modules
module load python         # SLURM doesn't like specifying versions
module load r              # ""
module load sra-toolkit
module load star
module load subread
module load edirect        # for esearch, efetch, etc.

# Remove the temporary virtualenv creation and replace with existing .venv activation
source .venv/bin/activate

# Verify Python packages are available
python -c "import numpy; import pandas; print(f'NumPy version: {numpy.__version__}\nPandas version: {pandas.__version__}')"

# Run the pipeline with parallel processing
if ! python b16f10_pipeline.py --all --parallel $SLURM_CPUS_PER_TASK; then
    # If the pipeline failed with an error (non-zero exit code)
    echo "Pipeline failed with an error. Check logs for details."
    exit 1
fi

# Only check completion and resubmit if we hit the time limit
if [[ $SLURM_JOB_STATE == "TIMEOUT" ]]; then
    if ! check_completion; then
        echo "Job timed out and some accessions are incomplete. Resubmitting..."
        sbatch $0
        exit 0
    fi
fi

echo "Pipeline completed successfully for all accessions"