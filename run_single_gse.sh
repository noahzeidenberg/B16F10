#!/bin/bash
#SBATCH --account=def-lukens
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=download_single_gse
#SBATCH --output=logs/download_rnaseq/download_single_gse_%j.out
#SBATCH --error=logs/download_rnaseq/download_single_gse_%j.err
#SBATCH --tmp=100G

# Script to run a single GSE ID download using SLURM

# Check if GSE ID is provided
if [ $# -ne 1 ]; then
    echo "Usage: sbatch $0 <GSE_ID>"
    echo "Example: sbatch $0 GSE275540"
    exit 1
fi

GSE_ID=$1

# Load required modules
module load sra-toolkit
module load r
module load gcc

# Set up environment
source ~/scratch/B16F10/.venv/bin/activate

# Create logs directory if it doesn't exist
mkdir -p logs/download_rnaseq

# Check if we should use temporary directory
if [ "$USE_TMPDIR" = "0" ]; then
    echo "Using permanent directory for downloads"
    cd $SLURM_SUBMIT_DIR
else
    echo "Using temporary directory for downloads"
    # Create a temporary directory for this job
    TMP_DIR=$SLURM_TMPDIR/download_${GSE_ID}
    mkdir -p $TMP_DIR
    cd $TMP_DIR

    # Copy the R script and .env file to the temporary directory
    cp $SLURM_SUBMIT_DIR/download_data.R .
    cp $SLURM_SUBMIT_DIR/.env .
fi

# Load required modules again in case they were unloaded in temporary directory
module load sra-toolkit
module load r
module load gcc

# Function to copy files from temporary to permanent storage
copy_files() {
    if [ "$USE_TMPDIR" = "0" ]; then
        echo "Using permanent directory, no need to copy files"
        return 0
    fi

    echo "Copying files from temporary to permanent storage..."
    local TMP_GSE_DIR=$(find $TMP_DIR -type d -name "$GSE_ID" | head -n 1)
    local GSE_DIR=$SLURM_SUBMIT_DIR/${GSE_ID}
    
    if [ -z "$TMP_GSE_DIR" ]; then
        echo "Error: Could not find GSE directory in temporary location"
        return 1
    fi
    
    # Create the proper directory structure in the permanent location
    mkdir -p $GSE_DIR/samples
    
    # Copy sample directories
    if [ -d "$TMP_GSE_DIR/samples" ]; then
        echo "Copying sample directories..."
        for gsm_dir in "$TMP_GSE_DIR/samples"/*; do
            if [ -d "$gsm_dir" ]; then
                gsm_id=$(basename "$gsm_dir")
                echo "Copying $gsm_id to $GSE_DIR/samples/$gsm_id"
                mkdir -p "$GSE_DIR/samples/$gsm_id"
                
                # Copy with verbose output and error checking
                if ! cp -rv "$gsm_dir"/* "$GSE_DIR/samples/$gsm_id/" 2>&1; then
                    echo "Error copying $gsm_id directory"
                    return 1
                fi
                
                # Verify the copy
                if [ ! -d "$GSE_DIR/samples/$gsm_id" ]; then
                    echo "Error: Failed to create $gsm_id directory in permanent location"
                    return 1
                fi
            fi
        done
    else
        echo "Samples directory not found in $TMP_GSE_DIR"
        return 1
    fi
    
    # Copy any other files from the GSE directory
    echo "Copying GSE directory contents..."
    if ! cp -rv "$TMP_GSE_DIR"/* "$GSE_DIR/" 2>&1; then
        echo "Error copying GSE directory contents"
        return 1
    fi
    
    # Verify the copy was successful
    echo "Verifying files in permanent location:"
    if [ -d "$GSE_DIR/samples" ]; then
        echo "Contents of $GSE_DIR/samples:"
        ls -la "$GSE_DIR/samples"
        
        # Check each sample directory
        for gsm_dir in "$GSE_DIR/samples"/*; do
            if [ -d "$gsm_dir" ]; then
                echo "Contents of $(basename "$gsm_dir"):"
                ls -la "$gsm_dir"
                if [ -d "$gsm_dir/SRA" ]; then
                    echo "Contents of $(basename "$gsm_dir")/SRA:"
                    ls -la "$gsm_dir/SRA"
                fi
            fi
        done
    else
        echo "Samples directory not found in permanent location"
        return 1
    fi
    
    return 0
}

# Set up trap to handle signals
trap 'echo "Received signal, copying files..."; copy_files; exit 1' SIGTERM SIGINT SIGUSR1

# Run the R script
Rscript download_data.R $GSE_ID

# Copy files to permanent storage if using temporary directory
if [ "$USE_TMPDIR" = "1" ]; then
    if ! copy_files; then
        echo "Error: Failed to copy files to permanent storage"
        exit 1
    fi

    # Clean up only after successful copying
    cd $SLURM_SUBMIT_DIR
    rm -rf $TMP_DIR
fi

echo "Download and conversion complete for $GSE_ID"