#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --job-name=download_rnaseq
#SBATCH --output=logs/download_rnaseq/download_rnaseq_%A_%a.out
#SBATCH --error=logs/download_rnaseq/download_rnaseq_%A_%a.err
#SBATCH --array=1-72
#SBATCH --tmp=100G

# Load required modules
module load sra-toolkit
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Get the GSE ID for this array job
GSE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rna_seq_gse_ids.txt)

# Get absolute paths
SUBMIT_DIR=$(realpath $SLURM_SUBMIT_DIR)
TMP_DIR=$(realpath $SLURM_TMPDIR)

# Check if we should use temporary directory
if [ "$USE_TMPDIR" = "0" ]; then
    echo "Using permanent directory for downloads"
    cd $SUBMIT_DIR
    WORK_DIR=$SUBMIT_DIR
else
    echo "Using temporary directory for downloads"
    # Create a temporary directory for this job
    WORK_DIR=$TMP_DIR/download_${GSE_ID}
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # Copy the R script and .env file to the temporary directory
    cp $SUBMIT_DIR/download_data.R .
    cp $SUBMIT_DIR/.env .
fi

# Run the R script
Rscript download_data.R $GSE_ID

# Copy files to permanent storage if using temporary directory
if [ "$USE_TMPDIR" = "1" ]; then
    echo "=== Starting file copy operation ==="
    echo "Temporary directory: $WORK_DIR"
    echo "Permanent directory: $SUBMIT_DIR"
    echo "Current working directory: $(pwd)"
    
    # Create the proper directory structure in the permanent location
    GSE_DIR=$SUBMIT_DIR/${GSE_ID}
    mkdir -p $GSE_DIR/samples
    
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
                mkdir -p "$GSE_DIR/samples/$gsm_id"
                
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
        exit 1
    fi
    
    echo "=== File copy operation complete ==="
    
    # Clean up only after successful copying
    cd $SUBMIT_DIR
    rm -rf $WORK_DIR
    echo "=== Cleanup complete ==="
fi

echo "Download and conversion complete for $GSE_ID" 