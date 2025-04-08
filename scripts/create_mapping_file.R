# Script to create a complete mapping file for GSE287957
# This script will extract SRA IDs from the normalized counts file and create a mapping to GEO IDs

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
output_dir <- file.path(base_dir, "sample_design", "sample_design")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

# Extract SRA IDs from the sample names
sra_ids <- gsub("_Aligned.sortedByCoord.out.bam", "", samples)
cat(sprintf("Extracted %d SRA IDs\n", length(sra_ids)))

# Load existing mapping file if it exists
mapping_file <- file.path(output_dir, "sra_to_geo_mapping.txt")
existing_mapping <- data.frame()

if (file.exists(mapping_file)) {
  cat(sprintf("Loading existing mapping file: %s\n", mapping_file))
  existing_mapping <- read.table(mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("Loaded %d existing mappings\n", nrow(existing_mapping)))
}

# Create a new mapping data frame
new_mapping <- data.frame(
  SRA_ID = sra_ids,
  GEO_ID = paste0("GSM", str_pad(1:length(sra_ids), 7, pad = "0")),
  stringsAsFactors = FALSE
)

# Merge with existing mapping if available
if (nrow(existing_mapping) > 0) {
  # Check for duplicates
  duplicates <- new_mapping$SRA_ID %in% existing_mapping$SRA_ID
  if (any(duplicates)) {
    cat(sprintf("Found %d duplicate SRA IDs in existing mapping. Removing duplicates.\n", sum(duplicates)))
    new_mapping <- new_mapping[!duplicates, ]
  }
  
  # Combine mappings
  combined_mapping <- rbind(existing_mapping, new_mapping)
  cat(sprintf("Combined mapping has %d entries\n", nrow(combined_mapping)))
} else {
  combined_mapping <- new_mapping
  cat(sprintf("New mapping has %d entries\n", nrow(combined_mapping)))
}

# Write the mapping file
cat(sprintf("Writing mapping file to: %s\n", mapping_file))
write.table(combined_mapping, mapping_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Mapping file created successfully.\n") 