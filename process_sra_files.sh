#!/bin/bash

# Set strict error handling
set -euo pipefail

# Define the base directory (current directory)
BASE_DIR="$(pwd)"

# Create logs directory if it doesn't exist
mkdir -p "${BASE_DIR}/logs"

# Function to log messages
log_message() {
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${message}" | tee -a "${BASE_DIR}/logs/process_sra.log"
}

# Function to convert .sra files to .fastq and run FastQC
process_srr_folder() {
    local srr_dir="$1"
    
    # Change to the SRR folder
    cd "$srr_dir" || {
        log_message "ERROR: Failed to change to directory $srr_dir"
        return 1
    }

    # Find .sra files in the current directory
    for sra_file in *.sra; do
        if [[ -f "$sra_file" ]]; then
            # Extract the base name (without extension)
            sra_base="${sra_file%.sra}"
            
            # Check if fastq files already exist
            if [[ -f "${sra_base}_1.fastq" && -f "${sra_base}_2.fastq" ]]; then
                log_message "FASTQ files already exist for $sra_file, skipping conversion"
                continue
            fi
            
            # Convert .sra to .fastq using fasterq-dump
            log_message "Converting $sra_file to FASTQ..."
            fasterq-dump --split-files "$sra_file" || {
                log_message "ERROR: Conversion of $sra_file failed"
                continue
            }
            
            # Run FastQC on each fastq file
            for fq_file in "${sra_base}_1.fastq" "${sra_base}_2.fastq"; do
                if [[ -f "$fq_file" ]]; then
                    log_message "Running FastQC on $fq_file..."
                    fastqc "$fq_file" || {
                        log_message "ERROR: FastQC failed on $fq_file"
                    }
                fi
            done
        fi
    done

    # Change back to the base directory
    cd "$BASE_DIR" || {
        log_message "ERROR: Failed to return to base directory"
        return 1
    }
}

# Main processing
log_message "Starting SRA file processing"

# Check if fasterq-dump is available
if ! command -v fasterq-dump &> /dev/null; then
    log_message "ERROR: fasterq-dump is not installed or not in PATH"
    exit 1
fi

# Check if fastqc is available
if ! command -v fastqc &> /dev/null; then
    log_message "ERROR: fastqc is not installed or not in PATH"
    exit 1
fi

# Iterate through all directories starting with "SRR"
for dir in SRR*; do
    if [[ -d "$dir" ]]; then
        log_message "Processing folder $dir..."
        process_srr_folder "$dir"
    fi
done

log_message "All folders processed successfully" 