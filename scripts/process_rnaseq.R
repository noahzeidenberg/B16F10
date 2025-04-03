# RNA-seq Processing Pipeline
# This script performs QC, trimming, alignment, and counting of RNA-seq data

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "Rsubread",
  "sva",
  "edgeR",
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Function to run FastQC
run_fastqc <- function(fastq_dir, output_dir) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Run FastQC on all FASTQ files
  cmd <- sprintf("fastqc -o %s -t %d %s/*.fastq.gz",
                output_dir,
                parallel::detectCores() - 1,
                fastq_dir)
  system(cmd)
}

# Function to run fastp
run_fastp <- function(input_dir, output_dir) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all R1 files
  r1_files <- list.files(input_dir, pattern = "_1.fastq.gz$", full.names = TRUE)
  
  for (r1_file in r1_files) {
    r2_file <- sub("_1.fastq.gz$", "_2.fastq.gz", r1_file)
    if (!file.exists(r2_file)) next
    
    # Generate output filenames
    sample_name <- basename(sub("_1.fastq.gz$", "", r1_file))
    out_r1 <- file.path(output_dir, paste0(sample_name, "_1.trimmed.fastq.gz"))
    out_r2 <- file.path(output_dir, paste0(sample_name, "_2.trimmed.fastq.gz"))
    json_report <- file.path(output_dir, paste0(sample_name, ".fastp.json"))
    html_report <- file.path(output_dir, paste0(sample_name, ".fastp.html"))
    
    # Run fastp
    cmd <- sprintf("fastp --in1 %s --in2 %s --out1 %s --out2 %s --json %s --html %s --thread %d",
                  r1_file, r2_file, out_r1, out_r2, json_report, html_report,
                  parallel::detectCores() - 1)
    system(cmd)
  }
}

# Function to build STAR index
build_star_index <- function(fasta_file, gtf_file, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cmd <- sprintf("STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --runThreadN %d --sjdbGTFtagExonParentTranscript Parent",
                output_dir, fasta_file, gtf_file, parallel::detectCores() - 1)
  system(cmd)
}

# Function to run STAR alignment
run_star_alignment <- function(input_dir, genome_dir, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all trimmed R1 files
  r1_files <- list.files(input_dir, pattern = "_1.trimmed.fastq.gz$", full.names = TRUE)
  
  for (r1_file in r1_files) {
    r2_file <- sub("_1.trimmed.fastq.gz$", "_2.trimmed.fastq.gz", r1_file)
    if (!file.exists(r2_file)) next
    
    sample_name <- basename(sub("_1.trimmed.fastq.gz$", "", r1_file))
    out_prefix <- file.path(output_dir, sample_name)
    
    cmd <- sprintf("STAR --genomeDir %s --readFilesIn %s %s --readFilesCommand zcat --outFileNamePrefix %s_ --outSAMtype BAM SortedByCoordinate --runThreadN %d",
                  genome_dir, r1_file, r2_file, out_prefix, parallel::detectCores() - 1)
    system(cmd)
  }
}

# Function to run featureCounts
run_feature_counts <- function(bam_dir, gtf_file, output_file) {
  # Get all BAM files
  bam_files <- list.files(bam_dir, pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
  
  # Run featureCounts
  fc <- featureCounts(bam_files,
                     annot.ext = gtf_file,
                     isGTFAnnotationFile = TRUE,
                     GTF.featureType = "gene",
                     GTF.attrType = "gene_id",
                     nthreads = parallel::detectCores() - 1)
  
  # Save results
  saveRDS(fc, output_file)
  return(fc)
}

# Function to perform batch correction and normalization
normalize_counts <- function(counts, batch_info, output_file) {
  # Remove genes with zero counts across all samples
  keep <- rowSums(counts) > 0
  counts <- counts[keep, ]
  
  # ComBat-seq batch correction
  corrected_counts <- ComBat_seq(counts = as.matrix(counts),
                                batch = batch_info$batch)
  
  # Calculate TMM normalization factors
  y <- DGEList(counts = corrected_counts)
  y <- calcNormFactors(y, method = "TMM")
  
  # Get normalized counts (GeTMM)
  norm_counts <- cpm(y, log = TRUE)
  
  # Save results
  results <- list(
    raw_counts = counts,
    corrected_counts = corrected_counts,
    normalized_counts = norm_counts
  )
  saveRDS(results, output_file)
  return(results)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID")
    }
  }
  
  # Set up directories
  base_dir <- file.path(getwd(), "results", gse_id)
  fastqc_dir <- file.path(base_dir, "fastqc")
  trimmed_dir <- file.path(base_dir, "trimmed")
  star_index_dir <- file.path(base_dir, "star_index")
  alignment_dir <- file.path(base_dir, "alignment")
  counts_dir <- file.path(base_dir, "counts")
  
  # Reference files
  ref_fasta <- file.path(getwd(), "mm39", "GCF_000001635.27_GRCm39_genomic.fasta")
  ref_gff <- file.path(getwd(), "mm39", "GCF_000001635.27_GRCm39_genomic.gff")
  
  # Verify reference files exist
  if (!file.exists(ref_fasta) || !file.exists(ref_gff)) {
    stop("Reference files not found in mm39/ directory")
  }
  
  # Create directories
  dirs <- c(fastqc_dir, trimmed_dir, star_index_dir, alignment_dir, counts_dir)
  sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  
  # Process each sample
  sample_dirs <- list.dirs(file.path(base_dir, "samples"), recursive = FALSE)
  for (sample_dir in sample_dirs) {
    fastq_dir <- file.path(sample_dir, "SRA", "FASTQ")
    if (!dir.exists(fastq_dir)) next
    
    # Run FastQC on raw data
    run_fastqc(fastq_dir, fastqc_dir)
    
    # Run fastp
    run_fastp(fastq_dir, trimmed_dir)
  }
  
  # Build STAR index if not already built
  if (!file.exists(file.path(star_index_dir, "Genome"))) {
    build_star_index(ref_fasta, ref_gff, star_index_dir)
  }
  
  # Run STAR alignment
  run_star_alignment(trimmed_dir, star_index_dir, alignment_dir)
  
  # Run featureCounts
  counts <- run_feature_counts(
    alignment_dir,
    ref_gff,
    file.path(counts_dir, "feature_counts.rds")
  )
  
  # Get batch information (you'll need to modify this based on your metadata)
  batch_info <- data.frame(
    sample = colnames(counts$counts),
    batch = factor(rep(1, ncol(counts$counts)))  # Modify this based on your batches
  )
  
  # Perform normalization
  norm_results <- normalize_counts(
    counts$counts,
    batch_info,
    file.path(counts_dir, "normalized_counts.rds")
  )
  
  cat("Processing complete!\n")
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 