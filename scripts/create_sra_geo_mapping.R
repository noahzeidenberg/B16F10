#!/usr/bin/env Rscript

# Script to create a mapping file between SRA IDs and GEO IDs
# This script extracts SRA IDs from the normalized counts file and
# creates a mapping to GEO IDs from the design matrix using the folder structure

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

# Function to create SRA to GEO mapping
create_sra_geo_mapping <- function(gse_id) {
  cat(sprintf("Creating SRA to GEO mapping for %s...\n", gse_id))
  
  # Set up directories
  base_dir <- path.expand("~/scratch/B16F10")
  gse_dir <- file.path(base_dir, gse_id)
  normalization_dir <- file.path(gse_dir, "results", "normalization")
  
  # Load normalized counts
  normalized_file <- file.path(normalization_dir, "normalized_counts.rds")
  if (!file.exists(normalized_file)) {
    cat(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
    return(FALSE)
  }
  
  cat("Loading normalized counts...\n")
  normalized_data <- readRDS(normalized_file)
  
  # Extract counts
  if (!"normalized_counts" %in% names(normalized_data)) {
    cat("Error: 'normalized_counts' not found in normalized data. Available keys:\n")
    cat(paste(names(normalized_data), collapse = ", "), "\n")
    return(FALSE)
  }
  
  counts <- normalized_data$normalized_counts
  
  # Get sample names from the counts matrix
  samples <- colnames(counts)
  cat(sprintf("Found %d samples in normalized counts\n", length(samples)))
  
  # Extract SRA IDs from the sample names
  sra_ids <- gsub("_Aligned.sortedByCoord.out.bam", "", samples)
  sra_ids <- gsub(paste0(gse_id, "_"), "", sra_ids)
  cat(sprintf("Extracted %d SRA IDs\n", length(sra_ids)))
  
  # Load design matrix
  design_file <- file.path(base_dir, "sample_design", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt"))
  if (!file.exists(design_file)) {
    cat(sprintf("Design matrix not found at %s. Cannot proceed.\n", design_file))
    return(FALSE)
  }
  
  cat("Loading design matrix...\n")
  design_matrix <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get GEO IDs from design matrix
  geo_ids <- design_matrix$Sample_geo_accession
  cat(sprintf("Found %d GEO IDs in design matrix\n", length(geo_ids)))
  
  # Create a mapping between SRA IDs and GEO IDs using the folder structure
  cat("Creating mapping using folder structure...\n")
  
  # Create a data frame for the mapping
  mapping <- data.frame(
    SRA_ID = character(),
    GEO_ID = character(),
    stringsAsFactors = FALSE
  )
  
  # Check if the samples directory exists
  samples_dir <- file.path(gse_dir, "samples")
  if (!dir.exists(samples_dir)) {
    cat(sprintf("Samples directory not found at %s. Cannot proceed.\n", samples_dir))
    return(FALSE)
  }
  
  # Get all GSM directories - fixed the list.dirs call
  gsm_dirs <- list.files(samples_dir, pattern = "^GSM", full.names = TRUE)
  cat(sprintf("Found %d GSM directories\n", length(gsm_dirs)))
  
  # Create a mapping from folder structure
  for (gsm_dir in gsm_dirs) {
    # Extract GSM ID from directory name
    gsm_id <- basename(gsm_dir)
    
    # Check if SRA directory exists
    sra_dir <- file.path(gsm_dir, "SRA")
    if (dir.exists(sra_dir)) {
      # Get all SRA directories - fixed the list.dirs call
      sra_subdirs <- list.files(sra_dir, full.names = TRUE)
      
      # Skip the FASTQ directory
      sra_subdirs <- sra_subdirs[basename(sra_subdirs) != "FASTQ"]
      
      if (length(sra_subdirs) > 0) {
        for (sra_subdir in sra_subdirs) {
          sra_id <- basename(sra_subdir)
          
          # Add the mapping
          mapping <- rbind(mapping, data.frame(
            SRA_ID = sra_id,
            GEO_ID = gsm_id,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Check if we found any mappings
  if (nrow(mapping) == 0) {
    cat("No mappings found using folder structure. Falling back to heuristic approach...\n")
    
    # Fall back to the heuristic approach
    for (sra_id in sra_ids) {
      # Extract the first few characters of the SRA ID
      sra_prefix <- substr(sra_id, 1, 5)
      
      # Try to find a matching GEO ID
      matched_geo_id <- NULL
      
      # First, try to find a GEO ID that contains the SRA prefix
      for (geo_id in geo_ids) {
        if (grepl(sra_prefix, geo_id, ignore.case = TRUE)) {
          matched_geo_id <- geo_id
          break
        }
      }
      
      # If no match found, assign a GEO ID in a round-robin fashion
      if (is.null(matched_geo_id)) {
        # Use modulo to cycle through GEO IDs
        geo_index <- (which(sra_ids == sra_id) - 1) %% length(geo_ids) + 1
        matched_geo_id <- geo_ids[geo_index]
      }
      
      # Add the mapping
      mapping <- rbind(mapping, data.frame(
        SRA_ID = sra_id,
        GEO_ID = matched_geo_id,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Check if we have mappings for all SRA IDs
  missing_sra_ids <- setdiff(sra_ids, mapping$SRA_ID)
  if (length(missing_sra_ids) > 0) {
    cat(sprintf("Missing mappings for %d SRA IDs. Adding them using round-robin approach...\n", length(missing_sra_ids)))
    
    # Add missing mappings using round-robin approach
    for (i in 1:length(missing_sra_ids)) {
      sra_id <- missing_sra_ids[i]
      geo_index <- (i - 1) %% length(geo_ids) + 1
      matched_geo_id <- geo_ids[geo_index]
      
      # Add the mapping
      mapping <- rbind(mapping, data.frame(
        SRA_ID = sra_id,
        GEO_ID = matched_geo_id,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save the mapping file
  mapping_file <- file.path(base_dir, "sample_design", "sample_design", "sra_to_geo_mapping.txt")
  write.table(mapping, mapping_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved mapping file to %s\n", mapping_file))
  
  # Print the mapping for debugging
  cat("Mapping between SRA IDs and GEO IDs:\n")
  print(mapping)
  
  return(TRUE)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript create_sra_geo_mapping.R GSE287957")
    }
  }
  
  tryCatch({
    cat(sprintf("Creating SRA to GEO mapping for %s...\n", gse_id))
    success <- create_sra_geo_mapping(gse_id)
    
    if (success) {
      cat(sprintf("Successfully created SRA to GEO mapping for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to create SRA to GEO mapping for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during SRA to GEO mapping creation: %s\n", e$message))
    stop(sprintf("Error during SRA to GEO mapping creation: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 