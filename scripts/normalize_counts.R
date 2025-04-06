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
  
  cat(sprintf("Loaded %d gene lengths\n", nrow(gene_lengths)))
  return(gene_lengths)
}

# Function to perform gene-length corrected TMM normalization
normalize_counts <- function(counts, gene_lengths, output_dir) {
  cat("Starting gene-length corrected TMM normalization...\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Match gene IDs
  common_genes <- intersect(rownames(counts), gene_lengths$gene_id)
  cat(sprintf("Found %d common genes\n", length(common_genes)))
  
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
main <- function() {
  # Set up directories
  base_dir <- getwd()
  batch_correction_dir <- file.path(base_dir, "results", "batch_correction")
  output_dir <- file.path(base_dir, "results", "normalization")
  
  cat("Setting up directories...\n")
  cat(sprintf("Base directory: %s\n", base_dir))
  cat(sprintf("Output directory: %s\n", output_dir))
  
  # Check if normalization has already been performed
  if (file.exists(file.path(output_dir, "normalized_counts.rds"))) {
    cat("Normalization has already been performed, skipping\n")
    return(0)
  }
  
  # Load batch-corrected counts # temporarily using non-batch corrected counts
  # batch_corrected_file <- file.path(batch_correction_dir, "batch_corrected_counts.rds")
  # if (!file.exists(batch_corrected_file)) {
  #   cat("Batch-corrected counts not found, cannot proceed\n")
  #   return(1)
  # }

  # Load non-batch corrected counts from feature counts
  batch_corrected_file <- file.path(base_dir, "results", "counts", "feature_counts.rds")
  if (!file.exists(batch_corrected_file)) {
    cat("Feature counts not found, cannot proceed\n")
    return(1)
  }
  
  cat("Loading batch-corrected counts...\n")
  batch_corrected <- readRDS(batch_corrected_file)
  
  # Load gene length information
  ref_gtf <- file.path(base_dir, "mm39", "GCF_000001635.27_GRCm39_genomic.gtf")
  if (!file.exists(ref_gtf)) {
    cat("Reference GTF file not found, cannot proceed\n")
    return(1)
  }
  
  gene_lengths <- load_gene_lengths(ref_gtf)
  
  # Perform normalization
  results <- normalize_counts(
    batch_corrected$corrected_counts,
    gene_lengths,
    output_dir
  )
  
  cat("=== RNA-seq Normalization Pipeline Completed Successfully ===\n")
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 