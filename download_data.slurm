#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=download_rnaseq
#SBATCH --output=logs/download_rnaseq/download_rnaseq_%A_%a.out
#SBATCH --error=logs/download_rnaseq/download_rnaseq_%A_%a.err
#SBATCH --array=1-60%5
#SBATCH --tmp=200G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your.email@example.com

# Enable debugging
set -x

# Load required modules
module load sra-toolkit
module load r
module load gcc
module load parallel

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
    echo "Using permanent directory for downloads"
    cd $SUBMIT_DIR
    WORK_DIR=$SUBMIT_DIR
    
    # Ensure the GSE directory exists and has proper permissions
    GSE_DIR=$SUBMIT_DIR/${GSE_ID}
    echo "Ensuring GSE directory exists: $GSE_DIR"
    mkdir -p $GSE_DIR/samples
    chmod -R 755 $GSE_DIR
    echo "GSE directory permissions:"
    ls -ld $GSE_DIR
    echo "Samples directory permissions:"
    ls -ld $GSE_DIR/samples
else
    echo "Using temporary directory for downloads"
    # Create a temporary directory for this job
    WORK_DIR=$TMP_DIR/download_${GSE_ID}
    echo "Creating WORK_DIR: $WORK_DIR"
    if ! mkdir -p $WORK_DIR; then
        echo "Error: Failed to create WORK_DIR"
        echo "TMP_DIR permissions:"
        ls -ld $TMP_DIR
        echo "TMP_DIR contents:"
        ls -la $TMP_DIR
        exit 1
    fi
    
    # Set proper permissions
    chmod -R 755 $WORK_DIR
    echo "WORK_DIR permissions:"
    ls -ld $WORK_DIR
    
    cd $WORK_DIR

    # Copy the R script and .env file to the temporary directory
    echo "Copying files to WORK_DIR"
    if ! cp $SUBMIT_DIR/download_data.R . || ! cp $SUBMIT_DIR/.env .; then
        echo "Error: Failed to copy required files to WORK_DIR"
        exit 1
    fi
    
    # Create the GSE directory structure in the temporary location
    GSE_DIR=$WORK_DIR/${GSE_ID}
    echo "Creating GSE directory structure in temporary location: $GSE_DIR"
    mkdir -p $GSE_DIR/samples
    chmod -R 755 $GSE_DIR
    echo "Temporary GSE directory permissions:"
    ls -ld $GSE_DIR
    echo "Temporary samples directory permissions:"
    ls -ld $GSE_DIR/samples
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
    GSE_DIR=$SUBMIT_DIR/${GSE_ID}
    echo "Creating GSE_DIR: $GSE_DIR"
    if ! mkdir -p $GSE_DIR/samples; then
        echo "Error: Failed to create GSE_DIR"
        return 1
    fi
    
    # Set proper permissions
    chmod -R 755 $GSE_DIR
    echo "GSE_DIR permissions:"
    ls -ld $GSE_DIR
    echo "Samples directory permissions:"
    ls -ld $GSE_DIR/samples
    
    # Copy sample directories
    if [ -d "samples" ]; then
        echo "Copying sample directories from $WORK_DIR/samples..."
        echo "Contents of samples directory:"
        ls -la "samples"
        
        for gsm_dir in samples/*; do
            if [ -d "$gsm_dir" ]; then
                gsm_id=$(basename "$gsm_dir")
                echo "Copying $gsm_id to $GSE_DIR/samples/$gsm_id"
                echo "Contents of source directory:"
                ls -la "$gsm_dir"
                
                # Create destination directory
                if ! mkdir -p "$GSE_DIR/samples/$gsm_id"; then
                    echo "Error: Failed to create destination directory for $gsm_id"
                    continue
                fi
                
                # Set proper permissions
                chmod -R 755 "$GSE_DIR/samples/$gsm_id"
                
                # Copy with verbose output and error checking
                echo "Starting copy operation..."
                if ! cp -rv "$gsm_dir"/* "$GSE_DIR/samples/$gsm_id/" 2>&1; then
                    echo "Warning: Error copying $gsm_id directory, continuing with next sample"
                    echo "Contents of destination directory after failed copy:"
                    ls -la "$GSE_DIR/samples/$gsm_id"
                    continue
                fi
                
                # Verify the copy
                if [ ! -d "$GSE_DIR/samples/$gsm_id" ]; then
                    echo "Warning: Failed to create $gsm_id directory in permanent location, continuing with next sample"
                    continue
                fi
                
                # Verify contents were copied
                echo "Contents of destination directory after copy:"
                ls -la "$GSE_DIR/samples/$gsm_id"
                
                if [ ! "$(ls -A "$GSE_DIR/samples/$gsm_id")" ]; then
                    echo "Warning: $gsm_id directory is empty after copy, continuing with next sample"
                    continue
                fi
            fi
        done
    else
        echo "Warning: Samples directory not found in $WORK_DIR"
    fi
    
    # Copy any other files from the GSE directory
    if [ -d "$GSE_ID" ]; then
        echo "Copying GSE directory contents from $WORK_DIR/$GSE_ID..."
        echo "Contents to be copied:"
        ls -la "$GSE_ID"/*
        
        if ! cp -rv "$GSE_ID"/* "$GSE_DIR/" 2>&1; then
            echo "Warning: Error copying GSE directory contents"
            echo "Contents of destination after failed copy:"
            ls -la "$GSE_DIR"
        fi
    fi
    
    # Final verification of all copied files
    echo "=== Verifying files in permanent location ==="
    echo "Checking $GSE_DIR"
    if [ -d "$GSE_DIR/samples" ]; then
        echo "Contents of $GSE_DIR/samples:"
        ls -la "$GSE_DIR/samples"
        
        # Check each sample directory
        for gsm_dir in "$GSE_DIR/samples"/*; do
            if [ -d "$gsm_dir" ]; then
                gsm_id=$(basename "$gsm_dir")
                echo "Contents of $gsm_id:"
                ls -la "$gsm_dir"
                if [ -d "$gsm_dir/SRA" ]; then
                    echo "Contents of $gsm_id/SRA:"
                    ls -la "$gsm_dir/SRA"
                    if [ -d "$gsm_dir/SRA/FASTQ" ]; then
                        echo "Contents of $gsm_id/SRA/FASTQ:"
                        ls -la "$gsm_dir/SRA/FASTQ"
                    fi
                fi
            fi
        done
    else
        echo "Error: Samples directory not found in permanent location"
        return 1
    fi
    
    echo "=== File copy operation complete ==="
    return 0
}

# Function to check if FASTQ files exist
check_fastq_files() {
    echo "=== Checking for existing FASTQ files ==="
    GSE_DIR=$SUBMIT_DIR/${GSE_ID}
    
    if [ ! -d "$GSE_DIR/samples" ]; then
        echo "No samples directory found, proceeding with download"
        return 1
    fi
    
    for gsm_dir in "$GSE_DIR/samples"/*; do
        if [ -d "$gsm_dir" ]; then
            fastq_dir="$gsm_dir/SRA/FASTQ"
            if [ -d "$fastq_dir" ]; then
                fastq_files=$(find "$fastq_dir" -name "*.fastq.gz" | wc -l)
                if [ "$fastq_files" -gt 0 ]; then
                    echo "Found existing FASTQ files in $fastq_dir"
                    return 0
                fi
            fi
        fi
    done
    
    echo "No existing FASTQ files found, proceeding with download"
    return 1
}

# Set up trap to handle signals
trap 'echo "Received signal, copying files..."; copy_files; exit 1' SIGTERM SIGINT SIGUSR1

# Check for existing FASTQ files
if check_fastq_files; then
    echo "FASTQ files already exist for $GSE_ID. Skipping download and conversion."
    exit 0
fi

# Run the R script
echo "=== Starting R script execution ==="
Rscript download_data.R $GSE_ID
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

echo "Download and conversion complete for $GSE_ID" 