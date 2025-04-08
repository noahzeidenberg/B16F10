# Script to update the design matrix for GSE287957
# This script will create a design matrix that includes all samples from the normalized counts

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Set up directories
base_dir <- path.expand("~/scratch/B16F10")
gse_id <- "GSE287957"
gse_dir <- file.path(base_dir, gse_id)
normalization_dir <- file.path(gse_dir, "results", "normalization")
design_matrices_dir <- file.path(base_dir, "sample_design", "sample_design", "design_matrices")

# Create design matrices directory if it doesn't exist
dir.create(design_matrices_dir, showWarnings = FALSE, recursive = TRUE)

# Load normalized counts
normalized_file <- file.path(normalization_dir, "normalized_counts.rds")
cat(sprintf("Loading normalized counts from: %s\n", normalized_file))

if (!file.exists(normalized_file)) {
  cat(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
  quit(status = 1)
}

normalized_data <- readRDS(normalized_file)

# Extract counts
if (!"normalized_counts" %in% names(normalized_data)) {
  cat("Error: 'normalized_counts' not found in normalized data. Available keys:\n")
  cat(paste(names(normalized_data), collapse = ", "), "\n")
  quit(status = 1)
}

counts <- normalized_data$normalized_counts

# Get sample names from the counts matrix
samples <- colnames(counts)
cat(sprintf("Found %d samples in normalized counts\n", length(samples)))

# Load existing design matrix if it exists
design_file <- file.path(design_matrices_dir, paste0(gse_id, "_design_matrix.txt"))
existing_design <- data.frame()

if (file.exists(design_file)) {
  cat(sprintf("Loading existing design matrix: %s\n", design_file))
  existing_design <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("Loaded existing design matrix with %d samples\n", nrow(existing_design)))
}

# Load mapping file
mapping_file <- file.path(base_dir, "sample_design", "sample_design", "sra_to_geo_mapping.txt")
mapping <- data.frame()

if (file.exists(mapping_file)) {
  cat(sprintf("Loading mapping file: %s\n", mapping_file))
  mapping <- read.table(mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("Loaded mapping file with %d entries\n", nrow(mapping)))
} else {
  cat("Mapping file not found. Cannot proceed.\n")
  quit(status = 1)
}

# Extract SRA IDs from the sample names
sra_ids <- gsub("_Aligned.sortedByCoord.out.bam", "", samples)
cat(sprintf("Extracted %d SRA IDs\n", length(sra_ids)))

# Create a new design matrix
# First, get the GEO IDs for all SRA IDs
geo_ids <- mapping$GEO_ID[mapping$SRA_ID %in% sra_ids]
cat(sprintf("Found %d GEO IDs for the SRA IDs\n", length(geo_ids)))

# Create a new design matrix with all samples
# For simplicity, we'll assign all samples to a single group
new_design <- data.frame(
  Sample_geo_accession = geo_ids,
  Group = rep("all_samples", length(geo_ids)),
  stringsAsFactors = FALSE
)

# If we have an existing design matrix, try to preserve the group assignments
if (nrow(existing_design) > 0) {
  # For each GEO ID in the existing design matrix, update the group in the new design matrix
  for (i in 1:nrow(existing_design)) {
    geo_id <- existing_design$Sample_geo_accession[i]
    group <- existing_design$Group[i]
    
    # Find the index of this GEO ID in the new design matrix
    idx <- which(new_design$Sample_geo_accession == geo_id)
    
    # If found, update the group
    if (length(idx) > 0) {
      new_design$Group[idx] <- group
    }
  }
}

# Write the design matrix
cat(sprintf("Writing design matrix to: %s\n", design_file))
write.table(new_design, design_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Design matrix updated successfully.\n") 