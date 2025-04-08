# RNA-seq Normalization Pipeline
# This script performs gene-length corrected TMM normalization

# Load required packages
cat("=== Starting RNA-seq Normalization Pipeline ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
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

# Function to load gene length information
load_gene_lengths <- function(gtf_file) {
  cat("Loading gene length information from GTF file...\n")
  
  # Read GTF file
  gtf <- read.table(gtf_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  
  # Extract gene IDs and lengths
  gene_lengths <- data.frame(
    gene_id = character(),
    length = integer(),
    stringsAsFactors = FALSE
  )
  
  # Process each gene
  for (i in 1:nrow(gtf)) {
    if (gtf$feature[i] == "gene") {
      # Extract gene ID from attributes
      attrs <- strsplit(gtf$attributes[i], ";")[[1]]
      gene_id <- gsub(".*gene_id \"(.*)\".*", "\\1", attrs[grep("gene_id", attrs)])
      
      # Calculate gene length
      gene_length <- gtf$end[i] - gtf$start[i] + 1
      
      # Add to data frame
      gene_lengths <- rbind(gene_lengths, data.frame(
        gene_id = gene_id,
        length = gene_length,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Remove duplicates
  gene_lengths <- gene_lengths[!duplicated(gene_lengths$gene_id), ]
  
  # Remove any "gene_id " prefix from gene IDs if present
  gene_lengths$gene_id <- gsub("^gene_id ", "", gene_lengths$gene_id)
  
  cat(sprintf("Loaded %d gene lengths\n", nrow(gene_lengths)))
  return(gene_lengths)
}

# Function to perform gene-length corrected TMM normalization
normalize_counts <- function(counts, gene_lengths, output_dir) {
  cat("Starting gene-length corrected TMM normalization...\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if counts is NULL or empty
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
    cat("Error: Counts data is NULL or empty\n")
    return(NULL)
  }
  
  # Print counts dimensions and first few rows for debugging
  cat(sprintf("Counts dimensions: %d rows x %d columns\n", nrow(counts), ncol(counts)))
  cat("First few row names of counts:\n")
  print(head(rownames(counts)))
  
  # Print gene lengths dimensions and first few rows for debugging
  cat(sprintf("Gene lengths dimensions: %d rows\n", nrow(gene_lengths)))
  cat("First few gene IDs from gene lengths:\n")
  print(head(gene_lengths$gene_id))
  
  # Match gene IDs
  common_genes <- intersect(rownames(counts), gene_lengths$gene_id)
  cat(sprintf("Found %d common genes\n", length(common_genes)))
  
  if (length(common_genes) == 0) {
    cat("Error: No common genes found between counts and gene lengths\n")
    cat("This might be due to different gene ID formats\n")
    
    # Try to clean up gene IDs in both datasets
    cat("Attempting to clean up gene IDs...\n")
    
    # Clean up counts row names
    clean_counts_rownames <- gsub("^gene_id ", "", rownames(counts))
    rownames(counts) <- clean_counts_rownames
    
    # Clean up gene lengths gene IDs
    clean_gene_lengths_ids <- gsub("^gene_id ", "", gene_lengths$gene_id)
    gene_lengths$gene_id <- clean_gene_lengths_ids
    
    # Try matching again
    common_genes <- intersect(rownames(counts), gene_lengths$gene_id)
    cat(sprintf("After cleanup: Found %d common genes\n", length(common_genes)))
    
    if (length(common_genes) == 0) {
      return(NULL)
    }
  }
  
  # Filter counts and gene lengths
  counts <- counts[common_genes, ]
  gene_lengths <- gene_lengths[gene_lengths$gene_id %in% common_genes, ]
  
  # Sort gene lengths to match counts
  gene_lengths <- gene_lengths[match(rownames(counts), gene_lengths$gene_id), ]
  
  # Create DGEList object
  y <- DGEList(counts = counts)
  
  # Calculate TMM normalization factors
  cat("Calculating TMM normalization factors...\n")
  y <- calcNormFactors(y, method = "TMM")
  
  # Get normalized counts (GeTMM)
  cat("Computing gene-length corrected TMM normalized counts...\n")
  
  # Calculate counts per million (CPM)
  cpm_counts <- cpm(y, log = TRUE)
  
  # Apply gene length correction
  # Divide by gene length (in kb) to get counts per million per kilobase
  gene_lengths_kb <- gene_lengths$length / 1000
  norm_counts <- sweep(cpm_counts, 1, gene_lengths_kb, "/")
  
  # Save results
  cat("Saving normalized counts...\n")
  results <- list(
    raw_counts = counts,
    normalized_counts = norm_counts,
    gene_lengths = gene_lengths
  )
  
  output_file <- file.path(output_dir, "normalized_counts.rds")
  saveRDS(results, output_file)
  cat(sprintf("Saved normalized counts to: %s\n", output_file))
  
  return(results)
}

# Main workflow
main <- function(gse_id = NULL) {
  # Set up directories
  base_dir <- getwd()
  
  # Check if GSE ID is provided
  if (is.null(gse_id)) {
    cat("Error: GSE ID is required\n")
    cat("Usage: Rscript normalize_counts.R <GSE_ID>\n")
    return(1)
  }
  
  # Create GSE-specific directories
  gse_dir <- file.path(base_dir, gse_id)
  batch_correction_dir <- file.path(base_dir, "results", "batch_correction")
  output_dir <- file.path(gse_dir, "results", "normalization")
  
  cat("Setting up directories...\n")
  cat(sprintf("Base directory: %s\n", base_dir))
  cat(sprintf("GSE directory: %s\n", gse_dir))
  cat(sprintf("Output directory: %s\n", output_dir))
  
  # Check if normalization has already been performed
  if (file.exists(file.path(output_dir, "normalized_counts.rds"))) {
    cat("Normalization has already been performed, skipping\n")
    return(0)
  }
  
  # Load batch-corrected counts 
  batch_corrected_file <- file.path(batch_correction_dir, "batch_corrected_counts.rds")
  if (!file.exists(batch_corrected_file)) {
    cat("Batch-corrected counts not found, cannot proceed\n")
    return(1)
  }
  
  cat("Loading batch-corrected counts...\n")
  feature_counts <- readRDS(batch_corrected_file)
  
  # Check if the feature counts object is valid
  if (is.null(feature_counts)) {
    cat("Error: Feature counts object is NULL\n")
    return(1)
  }
  
  # Extract the counts matrix from the feature counts object
  # The structure depends on how batch_correction.R saved it
  if (is.list(feature_counts) && !is.null(feature_counts$corrected_counts)) {
    # If it's a list with a 'corrected_counts' field (batch-corrected counts)
    counts_matrix <- feature_counts$corrected_counts
    cat("Using counts from feature_counts$corrected_counts\n")
  } else if (is.list(feature_counts) && !is.null(feature_counts$counts)) {
    # If it's a list with a 'counts' field
    counts_matrix <- feature_counts$counts
    cat("Using counts from feature_counts$counts\n")
  } else if (is.matrix(feature_counts)) {
    # If it's already a matrix
    counts_matrix <- feature_counts
    cat("Using counts directly from feature_counts (matrix)\n")
  } else {
    # Try to extract counts from the object
    cat("Attempting to extract counts from feature_counts object...\n")
    cat("Structure of feature_counts:\n")
    str(feature_counts)
    
    # Try to find a matrix in the object
    if (is.list(feature_counts)) {
      for (name in names(feature_counts)) {
        if (is.matrix(feature_counts[[name]])) {
          counts_matrix <- feature_counts[[name]]
          cat(sprintf("Using counts from feature_counts$%s\n", name))
          break
        }
      }
    }
    
    if (!exists("counts_matrix") || is.null(counts_matrix)) {
      cat("Error: Could not extract counts matrix from feature_counts object\n")
      return(1)
    }
  }
  
  # Filter counts matrix to only include samples from the current GSE
  cat(sprintf("Filtering counts matrix to only include samples from %s...\n", gse_id))
  
  # Get all sample names
  all_samples <- colnames(counts_matrix)
  cat(sprintf("Total samples in batch-corrected counts: %d\n", length(all_samples)))
  
  # Filter samples for the current GSE
  # Check for both old format (GSEID_SRRID) and new format (GSEID_GSMID_SRRID)
  gse_samples <- all_samples[grepl(paste0("^", gse_id, "_"), all_samples)]
  cat(sprintf("Samples for %s: %d\n", gse_id, length(gse_samples)))
  
  if (length(gse_samples) == 0) {
    cat(sprintf("Error: No samples found for %s in batch-corrected counts\n", gse_id))
    return(1)
  }
  
  # Subset the counts matrix to only include samples from the current GSE
  counts_matrix <- counts_matrix[, gse_samples]
  cat(sprintf("Filtered counts matrix dimensions: %d rows x %d columns\n", nrow(counts_matrix), ncol(counts_matrix)))
  
  # Load gene length information
  ref_gtf <- file.path("~/scratch/B16F10/mm39/GCF_000001635.27_GRCm39_genomic.gtf")
  if (!file.exists(ref_gtf)) {
    cat("Reference GTF file not found, cannot proceed\n")
    return(1)
  }
  
  gene_lengths <- load_gene_lengths(ref_gtf)
  
  # Perform normalization
  results <- normalize_counts(
    counts_matrix,
    gene_lengths,
    output_dir
  )
  
  if (is.null(results)) {
    cat("Normalization failed\n")
    return(1)
  }
  
  cat("=== RNA-seq Normalization Pipeline Completed Successfully ===\n")
  return(0)
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check if GSE ID is provided
  if (length(args) < 1) {
    cat("Error: GSE ID is required\n")
    cat("Usage: Rscript normalize_counts.R <GSE_ID>\n")
    quit(status = 1)
  }
  
  # Extract GSE ID from arguments
  gse_id <- args[1]
  
  # Validate GSE ID format (should start with "GSE")
  if (!grepl("^GSE[0-9]+$", gse_id)) {
    cat(sprintf("Error: Invalid GSE ID format: %s\n", gse_id))
    cat("GSE ID should start with 'GSE' followed by numbers\n")
    quit(status = 1)
  }
  
  # Run main function with GSE ID
  quit(status = main(gse_id))
} 