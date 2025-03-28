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
    echo "[${timestamp}] ${message}" | tee -a "${BASE_DIR}/logs/differential_expression.log"
}

# Function to check if a command exists
check_command() {
    local cmd="$1"
    if ! command -v "$cmd" &> /dev/null; then
        log_message "ERROR: $cmd is not installed or not in PATH"
        exit 1
    fi
}

# Function to check if a file exists
check_file() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file $file not found"
        exit 1
    fi
}

# Check required commands
log_message "Checking required commands..."
for cmd in STAR featureCounts Rscript; do
    check_command "$cmd"
done

# Create necessary directories
mkdir -p "${BASE_DIR}/results/star_index"
mkdir -p "${BASE_DIR}/results/star_aligned"
mkdir -p "${BASE_DIR}/results/counts"
mkdir -p "${BASE_DIR}/results/deseq2"

# Download and prepare reference genome if not exists
if [[ ! -f "${BASE_DIR}/reference/mm10.fa" ]]; then
    log_message "Downloading mouse reference genome..."
    mkdir -p "${BASE_DIR}/reference"
    wget -O "${BASE_DIR}/reference/mm10.fa.gz" https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
    gunzip "${BASE_DIR}/reference/mm10.fa.gz"
fi

# Download and prepare gene annotation if not exists
if [[ ! -f "${BASE_DIR}/reference/mm10.gtf" ]]; then
    log_message "Downloading mouse gene annotation..."
    wget -O "${BASE_DIR}/reference/mm10.gtf.gz" https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm38.109.gtf.gz
    gunzip "${BASE_DIR}/reference/mm10.gtf.gz"
fi

# Generate STAR index if not exists
if [[ ! -f "${BASE_DIR}/results/star_index/Genome" ]]; then
    log_message "Generating STAR index..."
    STAR --runMode genomeGenerate \
         --genomeDir "${BASE_DIR}/results/star_index" \
         --genomeFastaFiles "${BASE_DIR}/reference/mm10.fa" \
         --sjdbGTFfile "${BASE_DIR}/reference/mm10.gtf" \
         --runThreadN 4
fi

# Function to process each sample
process_sample() {
    local sample_dir="$1"
    local sample_name=$(basename "$sample_dir")
    
    log_message "Processing sample: $sample_name"
    
    # Create output directory for this sample
    mkdir -p "${BASE_DIR}/results/star_aligned/${sample_name}"
    
    # Find fastq files
    local fastq1="${sample_dir}/${sample_name}_1.fastq"
    local fastq2="${sample_dir}/${sample_name}_2.fastq"
    
    if [[ ! -f "$fastq1" || ! -f "$fastq2" ]]; then
        log_message "ERROR: FASTQ files not found for $sample_name"
        return 1
    fi
    
    # Run STAR alignment
    log_message "Running STAR alignment for $sample_name..."
    STAR --runMode alignReads \
         --genomeDir "${BASE_DIR}/results/star_index" \
         --readFilesIn "$fastq1" "$fastq2" \
         --outFileNamePrefix "${BASE_DIR}/results/star_aligned/${sample_name}/" \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 4
    
    # Create symbolic link to BAM file
    ln -sf "${BASE_DIR}/results/star_aligned/${sample_name}/Aligned.sortedByCoord.out.bam" \
           "${BASE_DIR}/results/star_aligned/${sample_name}.bam"
}

# Process all samples
log_message "Processing all samples..."
for dir in SRR*; do
    if [[ -d "$dir" ]]; then
        process_sample "$dir"
    fi
done

# Run featureCounts
log_message "Running featureCounts..."
featureCounts -T 4 \
             -p -t exon -g gene_id \
             -a "${BASE_DIR}/reference/mm10.gtf" \
             -o "${BASE_DIR}/results/counts/all_samples_counts.txt" \
             "${BASE_DIR}/results/star_aligned/"*.bam

# Create sample information file for DESeq2
log_message "Creating sample information file..."
cat > "${BASE_DIR}/results/counts/sample_info.csv" << EOF
sample,condition
EOF

# Add sample information (you'll need to modify this based on your experimental design)
for dir in SRR*; do
    if [[ -d "$dir" ]]; then
        sample_name=$(basename "$dir")
        # Add your condition here (e.g., "control" or "treatment")
        echo "${sample_name},control" >> "${BASE_DIR}/results/counts/sample_info.csv"
    fi
done

# Run DESeq2 analysis
log_message "Running DESeq2 analysis..."
Rscript -e '
library(DESeq2)
library(tidyverse)

# Read count data
counts <- read.table("results/counts/all_samples_counts.txt", header=TRUE, row.names=1)
counts <- counts[,-c(1:5)]  # Remove annotation columns

# Read sample information
sample_info <- read.csv("results/counts/sample_info.csv", row.names=1)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                            colData = sample_info,
                            design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Save results
write.csv(res, "results/deseq2/differential_expression_results.csv")

# Create MA plot
pdf("results/deseq2/MA_plot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

# Create PCA plot
vsd <- vst(dds, blind=FALSE)
pdf("results/deseq2/PCA_plot.pdf")
plotPCA(vsd, intgroup="condition")
dev.off()
'

log_message "Differential expression analysis completed successfully" 