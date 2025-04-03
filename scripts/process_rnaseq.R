# RNA-seq Processing Pipeline
# This script performs QC, trimming, alignment, and counting of RNA-seq data

# Load required packages
cat("=== Starting RNA-seq Processing Pipeline ===\n")
cat("Loading required packages...\n")

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
  cat(sprintf("Loading package: %s\n", pkg))
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Function to run FastQC
run_fastqc <- function(fastq_dir, output_dir) {
  cat(sprintf("Running FastQC on directory: %s\n", fastq_dir))
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Run FastQC on all FASTQ files
  cmd <- sprintf("fastqc -o %s -t %d %s/*.fastq.gz",
                output_dir,
                parallel::detectCores() - 1,
                fastq_dir)
  cat(sprintf("Executing FastQC command: %s\n", cmd))
  system(cmd)
  cat("FastQC completed\n")
}

# Function to run fastp
run_fastp <- function(input_dir, output_dir) {
  cat(sprintf("Running fastp on directory: %s\n", input_dir))
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all R1 files
  r1_files <- list.files(input_dir, pattern = "_1.fastq.gz$", full.names = TRUE)
  cat(sprintf("Found %d R1 files to process\n", length(r1_files)))
  
  for (r1_file in r1_files) {
    r2_file <- sub("_1.fastq.gz$", "_2.fastq.gz", r1_file)
    if (!file.exists(r2_file)) {
      cat(sprintf("Skipping %s - R2 file not found\n", r1_file))
      next
    }
    
    # Generate output filenames
    sample_name <- basename(sub("_1.fastq.gz$", "", r1_file))
    out_r1 <- file.path(output_dir, paste0(sample_name, "_1.trimmed.fastq.gz"))
    out_r2 <- file.path(output_dir, paste0(sample_name, "_2.trimmed.fastq.gz"))
    json_report <- file.path(output_dir, paste0(sample_name, ".fastp.json"))
    html_report <- file.path(output_dir, paste0(sample_name, ".fastp.html"))
    
    cat(sprintf("Processing sample: %s\n", sample_name))
    # Run fastp. Note it allows up to 16 threads
    cmd <- sprintf("fastp --in1 %s --in2 %s --out1 %s --out2 %s --json %s --html %s --thread 16",
                  r1_file, r2_file, out_r1, out_r2, json_report, html_report,
                  parallel::detectCores() - 1)
    cat(sprintf("Executing fastp command: %s\n", cmd))
    system(cmd)
  }
  cat("fastp processing completed\n")
}

# Function to build STAR index
build_star_index <- function(fasta_file, gtf_file, output_dir) {
  cat(sprintf("Building STAR index in directory: %s\n", output_dir))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cmd <- sprintf("STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --runThreadN %d --sjdbGTFtagExonParentTranscript Parent",
                output_dir, fasta_file, gtf_file, parallel::detectCores() - 1)
  cat(sprintf("Executing STAR index command: %s\n", cmd))
  system(cmd)
  cat("STAR index building completed\n")
}

# Function to run STAR alignment
run_star_alignment <- function(input_dir, genome_dir, output_dir) {
  cat(sprintf("Running STAR alignment on directory: %s\n", input_dir))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all trimmed R1 files
  r1_files <- list.files(input_dir, pattern = "_1.trimmed.fastq.gz$", full.names = TRUE)
  cat(sprintf("Found %d trimmed R1 files to align\n", length(r1_files)))
  
  for (r1_file in r1_files) {
    r2_file <- sub("_1.trimmed.fastq.gz$", "_2.trimmed.fastq.gz", r1_file)
    if (!file.exists(r2_file)) {
      cat(sprintf("Skipping %s - R2 file not found\n", r1_file))
      next
    }
    
    sample_name <- basename(sub("_1.trimmed.fastq.gz$", "", r1_file))
    out_prefix <- file.path(output_dir, sample_name)
    
    cat(sprintf("Aligning sample: %s\n", sample_name))
    cmd <- sprintf("STAR --genomeDir %s --readFilesIn %s %s --readFilesCommand zcat --outFileNamePrefix %s_ --outSAMtype BAM SortedByCoordinate --runThreadN %d",
                  genome_dir, r1_file, r2_file, out_prefix, parallel::detectCores() - 1)
    cat(sprintf("Executing STAR alignment command: %s\n", cmd))
    system(cmd)
  }
  cat("STAR alignment completed\n")
}

# Function to run featureCounts
run_feature_counts <- function(bam_dir, gtf_file, output_file) {
  cat(sprintf("Running featureCounts on BAM files in: %s\n", bam_dir))
  # Get all BAM files
  bam_files <- list.files(bam_dir, pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
  cat(sprintf("Found %d BAM files to process\n", length(bam_files)))
  
  # Run featureCounts
  cat("Starting featureCounts analysis...\n")
  fc <- featureCounts(bam_files,
                     annot.ext = gtf_file,
                     isGTFAnnotationFile = TRUE,
                     GTF.featureType = "gene",
                     GTF.attrType = "gene_id",
                     nthreads = parallel::detectCores() - 1)
  
  # Save results
  cat(sprintf("Saving featureCounts results to: %s\n", output_file))
  saveRDS(fc, output_file)
  cat("featureCounts completed\n")
  return(fc)
}

# Function to perform batch correction and normalization
normalize_counts <- function(counts, batch_info, output_file) {
  cat("Starting count normalization...\n")
  # Remove genes with zero counts across all samples
  keep <- rowSums(counts) > 0
  counts <- counts[keep, ]
  cat(sprintf("Removed %d genes with zero counts\n", sum(!keep)))
  
  # ComBat-seq batch correction
  cat("Performing ComBat-seq batch correction...\n")
  corrected_counts <- ComBat_seq(counts = as.matrix(counts),
                                batch = batch_info$batch)
  
  # Calculate TMM normalization factors
  cat("Calculating TMM normalization factors...\n")
  y <- DGEList(counts = corrected_counts)
  y <- calcNormFactors(y, method = "TMM")
  
  # Get normalized counts (GeTMM)
  cat("Computing normalized counts...\n")
  norm_counts <- cpm(y, log = TRUE)
  
  # Save results
  cat(sprintf("Saving normalized counts to: %s\n", output_file))
  results <- list(
    raw_counts = counts,
    corrected_counts = corrected_counts,
    normalized_counts = norm_counts
  )
  saveRDS(results, output_file)
  cat("Normalization completed\n")
  return(results)
}

# Function to check if output files exist
check_output_files <- function(gse_id) {
  base_dir <- file.path(getwd(), gse_id, "results")
  
  # Check for normalized counts file (final output)
  if (file.exists(file.path(base_dir, "counts", "normalized_counts.rds"))) {
    cat("Found normalized counts file, skipping processing\n")
    return(TRUE)
  }
  
  # Check for feature counts file
  if (file.exists(file.path(base_dir, "counts", "feature_counts.rds"))) {
    cat("Found feature counts file, skipping processing\n")
    return(TRUE)
  }
  
  # Check for alignment files
  alignment_dir <- file.path(base_dir, "alignment")
  if (dir.exists(alignment_dir)) {
    bam_files <- list.files(alignment_dir, pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
    if (length(bam_files) > 0) {
      cat("Found alignment files, skipping processing\n")
      return(TRUE)
    }
  }
  
  # Check for trimmed files
  trimmed_dir <- file.path(base_dir, "trimmed")
  if (dir.exists(trimmed_dir)) {
    trimmed_files <- list.files(trimmed_dir, pattern = "\\.trimmed\\.fastq\\.gz$", full.names = TRUE)
    if (length(trimmed_files) > 0) {
      cat("Found trimmed files, skipping processing\n")
      return(TRUE)
    }
  }
  
  # Check for FastQC results
  fastqc_dir <- file.path(base_dir, "fastqc")
  if (dir.exists(fastqc_dir)) {
    fastqc_files <- list.files(fastqc_dir, pattern = "\\.html$", full.names = TRUE)
    if (length(fastqc_files) > 0) {
      cat("Found FastQC results, skipping processing\n")
      return(TRUE)
    }
  }
  
  cat("No existing output files found, proceeding with processing\n")
  return(FALSE)
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
  
  cat(sprintf("Processing GSE ID: %s\n", gse_id))
  
  # Check for existing output files
  if (check_output_files(gse_id)) {
    return(0)
  }
  
  # Set up directories
  base_dir <- file.path(getwd(), gse_id, "results")
  fastqc_dir <- file.path(base_dir, "fastqc")
  trimmed_dir <- file.path(base_dir, "trimmed")
  star_index_dir <- file.path(base_dir, "star_index")
  alignment_dir <- file.path(base_dir, "alignment")
  counts_dir <- file.path(base_dir, "counts")
  
  cat("Setting up directories...\n")
  cat(sprintf("Base directory: %s\n", base_dir))
  
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
  sample_dirs <- list.dirs(file.path(getwd(), gse_id, "samples"), recursive = FALSE)
  cat(sprintf("Found %d sample directories to process\n", length(sample_dirs)))
  
  for (sample_dir in sample_dirs) {
    fastq_dir <- file.path(sample_dir, "SRA", "FASTQ")
    if (!dir.exists(fastq_dir)) {
      cat(sprintf("Skipping %s - FASTQ directory not found\n", sample_dir))
      next
    }
    
    cat(sprintf("Processing sample directory: %s\n", sample_dir))
    
    # Run FastQC on raw data
    run_fastqc(fastq_dir, fastqc_dir)
    
    # Run fastp
    run_fastp(fastq_dir, trimmed_dir)
  }
  
  # Build STAR index if not already built
  if (!file.exists(file.path(star_index_dir, "Genome"))) {
    cat("Building STAR index...\n")
    build_star_index(ref_fasta, ref_gff, star_index_dir)
  } else {
    cat("Using existing STAR index\n")
  }
  
  # Run STAR alignment
  run_star_alignment(trimmed_dir, star_index_dir, alignment_dir)
  
  # Run featureCounts
  cat("Starting featureCounts analysis...\n")
  counts <- run_feature_counts(
    alignment_dir,
    ref_gff,
    file.path(counts_dir, "feature_counts.rds")
  )
  
  # Get batch information
  cat("Preparing batch information...\n")
  batch_info <- data.frame(
    sample = colnames(counts$counts),
    batch = factor(rep(1, ncol(counts$counts)))  # Modify this based on your batches
  )
  
  # Perform normalization
  cat("Starting count normalization...\n")
  norm_results <- normalize_counts(
    counts$counts,
    batch_info,
    file.path(counts_dir, "normalized_counts.rds")
  )
  
  cat("=== RNA-seq Processing Pipeline Completed Successfully ===\n")
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 