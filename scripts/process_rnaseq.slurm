#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=02:30:00     # 4.5 hours should be sufficient for most datasets
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  # 16 threads for fastp and STAR
#SBATCH --mem=64G           # 64GB should be sufficient for most datasets
#SBATCH --job-name=process_rnaseq
#SBATCH --output=logs/process_rnaseq/process_rnaseq_%A_%a.out
#SBATCH --error=logs/process_rnaseq/process_rnaseq_%A_%a.err
#SBATCH --array=1-72
#SBATCH --tmp=200

# Set to run in working directory instead of temporary directory
#export USE_TMPDIR=0 # for now just testing with this

# Enable debugging
set -x

# Load required modules
module load star
module load fastqc
module load fastp
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get the GSE ID for this array job
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)

# Get absolute paths
SUBMIT_DIR=$(realpath $SLURM_SUBMIT_DIR)
TMP_DIR=$(realpath $SLURM_TMPDIR)

# Debug information
echo "=== Debug Information ==="
echo "SLURM_TMPDIR: $SLURM_TMPDIR"
echo "TMP_DIR: $TMP_DIR"
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "Current working directory: $(pwd)"
echo "SLURM_TMPDIR permissions:"
ls -ld $SLURM_TMPDIR
echo "SLURM_TMPDIR contents:"
ls -la $SLURM_TMPDIR
echo "=== End Debug Information ==="

# Check if we should use temporary directory
if [ "$USE_TMPDIR" = "0" ]; then
    echo "Using permanent directory for processing"
    cd $SUBMIT_DIR
    WORK_DIR=$SUBMIT_DIR
else
    echo "Using temporary directory for processing"
    # Create a temporary directory for this job
    WORK_DIR=$TMP_DIR/process_${GSE_ID}
    echo "Creating WORK_DIR: $WORK_DIR"
    if ! mkdir -p $WORK_DIR; then
        echo "Error: Failed to create WORK_DIR"
        echo "TMP_DIR permissions:"
        ls -ld $TMP_DIR
        echo "TMP_DIR contents:"
        ls -la $TMP_DIR
        exit 1
    fi
    cd $WORK_DIR

    # Copy the R script to the temporary directory
    echo "Copying files to WORK_DIR"
    if ! cp $SUBMIT_DIR/scripts/process_rnaseq.R .; then
        echo "Error: Failed to copy R script to WORK_DIR"
        exit 1
    fi
fi

# Function to copy files from temporary to permanent storage
copy_files() {
    echo "=== Starting copy_files function ==="
    if [ "$USE_TMPDIR" = "0" ]; then
        echo "Using permanent directory, no need to copy files"
        return 0
    fi

    echo "=== Starting file copy operation ==="
    echo "Temporary directory: $WORK_DIR"
    echo "Permanent directory: $SUBMIT_DIR"
    echo "Current working directory: $(pwd)"
    
    # Verify WORK_DIR exists and is accessible
    if [ ! -d "$WORK_DIR" ]; then
        echo "Error: WORK_DIR does not exist: $WORK_DIR"
        return 1
    fi
    
    # Create the proper directory structure in the permanent location
    GSE_DIR=$SUBMIT_DIR/${GSE_ID}/results
    echo "Creating GSE_DIR: $GSE_DIR"
    if ! mkdir -p $GSE_DIR/{fastqc,trimmed,alignment,counts}; then
        echo "Error: Failed to create GSE_DIR"
        return 1
    fi
    
    # Copy results
    for dir in fastqc trimmed alignment counts; do
        if [ -d "$dir" ]; then
            echo "Copying $dir directory..."
            if ! cp -rv "$dir"/* "$GSE_DIR/$dir/" 2>&1; then
                echo "Warning: Error copying $dir directory"
                echo "Contents of destination after failed copy:"
                ls -la "$GSE_DIR/$dir"
            fi
        else
            echo "Warning: $dir directory not found in $WORK_DIR"
        fi
    done
    
    echo "=== File copy operation complete ==="
    return 0
}

# Function to check if output files exist
check_output_files() {
    echo "=== Checking for existing output files ==="
    GSE_DIR=$SUBMIT_DIR/${GSE_ID}/results
    
    # Check for feature counts file (final output for this script)
    if [ -f "$GSE_DIR/counts/feature_counts.rds" ]; then
        echo "Found feature counts file, skipping all processing"
        return 0
    fi
    
    echo "No feature counts file found, proceeding with processing"
    return 1
}

# Set up trap to handle signals
trap 'echo "Received signal, copying files..."; copy_files; exit 1' SIGTERM SIGINT SIGUSR1

# Check for existing output files
if check_output_files; then
    echo "Output files already exist for $GSE_ID. Skipping processing."
    exit 0
fi

# Copy the FASTQ files from the previous step
if [ "$USE_TMPDIR" = "1" ]; then
    echo "Copying input files to temporary directory"
    if ! cp -r $SUBMIT_DIR/${GSE_ID}/results/* .; then
        echo "Error: Failed to copy input files to temporary directory"
        exit 1
    fi
fi

# Run the R script
echo "=== Starting R script execution ==="
Rscript scripts/process_rnaseq.R $GSE_ID
R_EXIT_CODE=$?
echo "R script exit code: $R_EXIT_CODE"

# Copy files to permanent storage if using temporary directory
if [ "$USE_TMPDIR" = "1" ]; then
    echo "=== Starting final file copy ==="
    if ! copy_files; then
        echo "Error: Failed to copy files to permanent storage"
        exit 1
    fi

    # Clean up only after successful copying
    cd $SUBMIT_DIR
    echo "Cleaning up WORK_DIR: $WORK_DIR"
    if ! rm -rf $WORK_DIR; then
        echo "Warning: Failed to clean up WORK_DIR"
        ls -la $WORK_DIR
    fi
    echo "=== Cleanup complete ==="
fi

echo "Processing complete for $GSE_ID" 