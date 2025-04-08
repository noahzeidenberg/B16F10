# RNA-seq Batch Correction Pipeline
# This script performs batch correction across all datasets

# Load required packages
cat("=== Starting RNA-seq Batch Correction Pipeline ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
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

# Function to create SRR to GSM mapping
create_srr_gsm_mapping <- function(base_dir) {
  cat("Creating SRR to GSM mapping...\n")
  
  # Read GSE IDs from the file
  gse_ids_file <- file.path(base_dir, "rna_seq_gse_ids.txt")
  if (!file.exists(gse_ids_file)) {
    cat(sprintf("Error: GSE IDs file not found at %s\n", gse_ids_file))
    return(NULL)
  }
  
  # Read the GSE IDs
  gse_ids <- readLines(gse_ids_file)
  cat(sprintf("Read %d GSE IDs from %s\n", length(gse_ids), gse_ids_file))
  
  # Create a mapping from SRR to GSM
  srr_to_gsm <- list()
  
  for (gse_id in gse_ids) {
    gse_dir <- file.path(base_dir, gse_id)
    
    # Skip if directory doesn't exist
    if (!dir.exists(gse_dir)) {
      cat(sprintf("Warning: Directory for %s does not exist, skipping\n", gse_id))
      next
    }
    
    # Find all sample directories
    sample_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    
    for (sample_dir in sample_dirs) {
      gsm_id <- basename(sample_dir)
      
      # Check if this is a GSM ID
      if (grepl("^GSM", gsm_id)) {
        # Look for SRA directory
        sra_dir <- file.path(sample_dir, "SRA")
        if (dir.exists(sra_dir)) {
          # Look for SRR directories (folders starting with SRR)
          srr_dirs <- list.dirs(sra_dir, recursive = FALSE)
          srr_dirs <- srr_dirs[grep("^SRR", basename(srr_dirs))]
          
          for (srr_dir in srr_dirs) {
            srr_id <- basename(srr_dir)
            
            # Add to mapping
            srr_to_gsm[[srr_id]] <- gsm_id
          }
        }
      }
    }
  }
  
  cat(sprintf("Created mapping for %d SRR IDs\n", length(srr_to_gsm)))
  return(srr_to_gsm)
}

# Function to load design matrices
load_design_matrices <- function(design_dir) {
  cat("Loading design matrices...\n")
  
  # Find all design matrix files
  design_files <- list.files(design_dir, pattern = "GSE.*_design_matrix.txt", full.names = TRUE)
  
  if (length(design_files) == 0) {
    cat("No design matrix files found\n")
    return(NULL)
  }
  
  cat(sprintf("Found %d design matrix files\n", length(design_files)))
  
  # Load all design matrices
  design_matrices <- list()
  
  for (file in design_files) {
    gse_id <- gsub("_design_matrix.txt", "", basename(file))
    cat(sprintf("Loading design matrix for %s\n", gse_id))
    
    # First, examine the file to understand its structure
    tryCatch({
      # Read the first few lines to understand the structure
      first_lines <- readLines(file, n = 5)
      cat("First few lines of the file:\n")
      print(first_lines)
      
      # Check for BOM character in the first line
      if (length(first_lines) > 0 && grepl("^﻿", first_lines[1])) {
        cat("Detected BOM character in the file\n")
        # Remove BOM from the first line
        first_lines[1] <- sub("^﻿", "", first_lines[1])
      }
      
      # Count the number of columns in the first data line
      if (length(first_lines) > 1) {
        data_line <- first_lines[2]
        num_cols <- length(strsplit(data_line, "\t")[[1]])
        cat(sprintf("Number of columns in data: %d\n", num_cols))
      }
      
      # Try to read the file with different approaches
      if (grepl("\t", first_lines[1])) {
        # Tab-separated file
        cat("Detected tab-separated file\n")
        # Use fileEncoding to handle BOM
        design_matrix <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                                   fileEncoding = "UTF-8-BOM")
      } else if (grepl(",", first_lines[1])) {
        # Comma-separated file
        cat("Detected comma-separated file\n")
        # Use fileEncoding to handle BOM
        design_matrix <- read.table(file, header = TRUE, sep = ",", stringsAsFactors = FALSE, 
                                   fileEncoding = "UTF-8-BOM")
      } else {
        # Space-separated file
        cat("Detected space-separated file\n")
        # Use fileEncoding to handle BOM
        design_matrix <- read.table(file, header = TRUE, stringsAsFactors = FALSE, 
                                   fileEncoding = "UTF-8-BOM")
      }
      
      # Print the structure of the loaded data
      cat("Structure of loaded design matrix:\n")
      str(design_matrix)
      
      # Store in list
      design_matrices[[gse_id]] <- design_matrix
    }, error = function(e) {
      cat(sprintf("Error loading design matrix for %s: %s\n", gse_id, e$message))
    })
  }
  
  return(design_matrices)
}

# Function to collect all feature counts files
collect_feature_counts <- function(base_dir, design_matrices, srr_to_gsm) {
  cat("Collecting all feature counts files...\n")
  
  # Read GSE IDs from the file
  gse_ids_file <- file.path(base_dir, "rna_seq_gse_ids.txt")
  if (!file.exists(gse_ids_file)) {
    cat(sprintf("Error: GSE IDs file not found at %s\n", gse_ids_file))
    return(NULL)
  }
  
  # Read the GSE IDs
  gse_ids <- readLines(gse_ids_file)
  cat(sprintf("Read %d GSE IDs from %s\n", length(gse_ids), gse_ids_file))
  
  # Collect all feature counts files
  all_counts <- list()
  all_samples <- c()
  all_gse_ids <- c()
  all_groups <- c()
  
  for (gse_id in gse_ids) {
    gse_dir <- file.path(base_dir, gse_id)
    
    # Skip if directory doesn't exist
    if (!dir.exists(gse_dir)) {
      cat(sprintf("Warning: Directory for %s does not exist, skipping\n", gse_id))
      next
    }
    
    counts_file <- file.path(gse_dir, "results", "counts", "feature_counts.rds")
    
    if (file.exists(counts_file)) {
      cat(sprintf("Loading counts from %s\n", gse_id))
      counts_obj <- readRDS(counts_file)
      
      # Extract counts and sample names
      # Make sure we're getting the counts matrix correctly
      if (is.null(counts_obj$counts)) {
        cat(sprintf("Warning: No counts matrix found in %s\n", counts_file))
        next
      }
      
      counts <- counts_obj$counts
      samples <- colnames(counts)
      
      # Skip GSEs with only one sample
      if (length(samples) <= 1) {
        cat(sprintf("Skipping %s: Only %d sample(s) found (ComBat-seq requires multiple samples per batch)\n", 
                   gse_id, length(samples)))
        next
      }
      
      # Add GSE ID to sample names to make them unique
      # Extract SRR IDs from the sample names
      srr_ids <- samples
      for (i in 1:length(samples)) {
        sample_id <- samples[i]
        
        # Extract SRR ID from the sample ID if it's a BAM file
        if (grepl("_Aligned\\.sortedByCoord\\.out\\.bam$", sample_id)) {
          srr_ids[i] <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", sample_id)
        } else if (grepl("^SRR", sample_id)) {
          srr_ids[i] <- sample_id
        }
      }
      
      # Get GSM IDs for each sample
      gsm_ids <- rep(NA, length(samples))
      for (i in 1:length(samples)) {
        srr_id <- srr_ids[i]
        if (srr_id %in% names(srr_to_gsm)) {
          gsm_ids[i] <- srr_to_gsm[[srr_id]]
        }
      }
      
      # Create new sample names with GSE ID, GSM ID, and SRR ID
      new_samples <- samples
      for (i in 1:length(samples)) {
        if (!is.na(gsm_ids[i])) {
          # If we have a GSM ID, use the new format
          new_samples[i] <- paste0(gse_id, "_", gsm_ids[i], "_", srr_ids[i])
        } else {
          # If we don't have a GSM ID, use the old format
          new_samples[i] <- paste0(gse_id, "_", samples[i])
        }
      }
      
      colnames(counts) <- new_samples
      
      # Get group information from design matrix if available
      groups <- rep(NA, length(samples))
      if (!is.null(design_matrices) && gse_id %in% names(design_matrices)) {
        design_matrix <- design_matrices[[gse_id]]
        
        # Match samples to their groups
        for (i in 1:length(samples)) {
          sample_id <- samples[i]
          
          # Extract SRR ID from the sample ID if it's a BAM file
          srr_id <- NULL
          if (grepl("_Aligned\\.sortedByCoord\\.out\\.bam$", sample_id)) {
            srr_id <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", sample_id)
          } else if (grepl("^SRR", sample_id)) {
            srr_id <- sample_id
          }
          
          # If we have an SRR ID, try to map it to a GSM ID
          if (!is.null(srr_id) && srr_id %in% names(srr_to_gsm)) {
            gsm_id <- srr_to_gsm[[srr_id]]
            
            # Find the row in the design matrix that matches this GSM ID
            match_idx <- which(design_matrix$Sample_geo_accession == gsm_id)
            if (length(match_idx) > 0) {
              groups[i] <- design_matrix$Group[match_idx]
            }
          } else {
            # Try direct matching with the sample ID
            match_idx <- which(design_matrix$Sample_geo_accession == sample_id)
            if (length(match_idx) > 0) {
              groups[i] <- design_matrix$Group[match_idx]
            }
          }
        }
        
        # Print matching information for debugging
        cat(sprintf("  Matched %d/%d samples to groups for %s\n", 
                   sum(!is.na(groups)), length(samples), gse_id))
        if (sum(!is.na(groups)) < length(samples)) {
          cat("  Unmatched samples: ")
          cat(paste(samples[is.na(groups)], collapse=", "))
          cat("\n")
        }
      } else {
        cat(sprintf("  No design matrix found for %s\n", gse_id))
      }
      
      # Add to lists
      all_counts[[gse_id]] <- counts
      all_samples <- c(all_samples, new_samples)
      all_gse_ids <- c(all_gse_ids, rep(gse_id, length(samples)))
      all_groups <- c(all_groups, groups)
    } else {
      cat(sprintf("Warning: No feature counts file found for %s\n", gse_id))
    }
  }
  
  if (length(all_counts) == 0) {
    cat("No feature counts files found, cannot proceed\n")
    return(NULL)
  }
  
  # Combine all counts
  cat("Combining all counts...\n")
  combined_counts <- do.call(cbind, all_counts)
  
  # Create batch information
  batch_info <- data.frame(
    sample = all_samples,
    gse_id = all_gse_ids,
    group = all_groups,
    batch = factor(all_gse_ids)  # Each GSE is a batch
  )
  
  # Print summary of group information
  cat("\nGroup information summary:\n")
  group_summary <- table(batch_info$group, useNA="ifany")
  print(group_summary)
  
  return(list(
    counts = combined_counts,
    batch_info = batch_info
  ))
}

# Function to perform batch correction
perform_batch_correction <- function(counts, batch_info, output_dir) {
  cat("Starting batch correction...\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Remove genes with zero counts across all samples
  keep <- rowSums(counts) > 0
  counts <- counts[keep, ]
  cat(sprintf("Removed %d genes with zero counts\n", sum(!keep)))
  
  # Check if we have multiple batches
  num_batches <- length(unique(batch_info$batch))
  cat(sprintf("Number of batches: %d\n", num_batches))
  
  # Check if we have at least 2 samples per batch
  samples_per_batch <- table(batch_info$batch)
  min_samples_per_batch <- min(samples_per_batch)
  cat(sprintf("Minimum samples per batch: %d\n", min_samples_per_batch))
  
  if (min_samples_per_batch < 2) {
    cat("Error: Some batches have fewer than 2 samples. ComBat-seq requires at least 2 samples per batch.\n")
    cat("Please ensure all batches have at least 2 samples before running batch correction.\n")
    return(NULL)
  }
  
  if (num_batches > 1) {
    # Check if we have group information
    has_groups <- !all(is.na(batch_info$group))
    
    if (has_groups) {
      cat("Performing ComBat-seq batch correction with model matrix...\n")
      
      # Check for NA groups
      na_groups <- is.na(batch_info$group)
      if (any(na_groups)) {
        cat(sprintf("Found %d samples with NA groups. Removing these samples before batch correction.\n", sum(na_groups)))
        
        # Remove samples with NA groups
        counts <- counts[, !na_groups]
        batch_info <- batch_info[!na_groups, ]
        
        # Recalculate batch information
        batch_info$batch <- factor(batch_info$batch)
        
        # Check if we still have enough samples per batch
        samples_per_batch <- table(batch_info$batch)
        min_samples_per_batch <- min(samples_per_batch)
        cat(sprintf("After removing NA groups, minimum samples per batch: %d\n", min_samples_per_batch))
        
        if (min_samples_per_batch < 2) {
          cat("Error: After removing NA groups, some batches have fewer than 2 samples.\n")
          cat("Cannot proceed with batch correction.\n")
          return(NULL)
        }
      }
      
      # Create a treatment factor that identifies control vs. treatment
      # First, identify all unique groups
      all_groups <- unique(batch_info$group)
      cat("All unique groups:\n")
      print(all_groups)
      
      # Create a treatment factor (control vs. treatment)
      # We'll identify control groups based on common naming patterns
      control_patterns <- c("control", "control[0-9]+", "ctrl", "untreated", "vehicle", "sham", "wt", "wild type")
      
      # Initialize treatment factor
      treatment_factor <- rep("treatment", nrow(batch_info))
      
      # Identify control groups - improved approach
      control_groups <- c()
      for (group in all_groups) {
        group_lower <- tolower(group)
        
        # Check for exact matches first
        if (group_lower %in% c("control", "ctrl", "untreated", "vehicle", "sham", "wt", "wild type")) {
          cat(sprintf("Exact match for control group: %s\n", group))
          control_groups <- c(control_groups, group)
          treatment_factor[batch_info$group == group] <- "control"
          next
        }
        
        # Check for pattern matches
        is_control <- FALSE
        for (pattern in control_patterns) {
          if (grepl(pattern, group_lower, perl = TRUE)) {
            cat(sprintf("Pattern match for control group: %s (matched pattern: %s)\n", group, pattern))
            control_groups <- c(control_groups, group)
            treatment_factor[batch_info$group == group] <- "control"
            is_control <- TRUE
            break
          }
        }
        
        if (!is_control) {
          cat(sprintf("Treatment group: %s\n", group))
        }
      }
      
      # If no control groups were identified, try a more aggressive approach
      if (length(control_groups) == 0) {
        cat("No control groups identified with standard patterns. Trying alternative approach...\n")
        
        # Try to identify control groups based on common experimental design patterns
        # Often the first group in each batch is a control
        for (batch in unique(batch_info$batch)) {
          batch_samples <- batch_info$batch == batch
          batch_groups <- unique(batch_info$group[batch_samples])
          
          if (length(batch_groups) > 0) {
            # Assume the first group is a control
            control_group <- batch_groups[1]
            cat(sprintf("Assuming first group in batch %s is control: %s\n", batch, control_group))
            control_groups <- c(control_groups, control_group)
            treatment_factor[batch_info$group == control_group] <- "control"
          }
        }
      }
      
      # If still no control groups, use a manual approach
      if (length(control_groups) == 0) {
        cat("Still no control groups identified. Using manual approach...\n")
        
        # For each batch, try to identify a control group
        for (batch in unique(batch_info$batch)) {
          batch_samples <- batch_info$batch == batch
          batch_groups <- unique(batch_info$group[batch_samples])
          
          # Ask the user to identify control groups
          cat(sprintf("Batch %s has the following groups:\n", batch))
          for (i in 1:length(batch_groups)) {
            cat(sprintf("  %d: %s\n", i, batch_groups[i]))
          }
          
          # For now, we'll use a heuristic: assume groups with "control" in the name are controls
          # In a real implementation, you might want to prompt the user
          for (group in batch_groups) {
            if (grepl("control", tolower(group))) {
              cat(sprintf("Manually identified control group: %s\n", group))
              control_groups <- c(control_groups, group)
              treatment_factor[batch_info$group == group] <- "control"
            }
          }
        }
      }
      
      # Print identified control groups
      cat("Identified control groups:\n")
      print(control_groups)
      
      # Standardize control group names
      cat("Standardizing control group names...\n")
      # Group "control", "control ", and "control1" together as just "control"
      for (i in 1:nrow(batch_info)) {
        if (batch_info$group[i] %in% c("control", "control ", "control1")) {
          batch_info$group[i] <- "control"
          cat(sprintf("Standardized group name: %s -> control\n", batch_info$group[i]))
        }
      }
      
      # Convert to factor
      treatment_factor <- factor(treatment_factor, levels = c("control", "treatment"))
      
      # Print treatment factor summary
      cat("Treatment factor summary:\n")
      print(table(treatment_factor))
      
      # Check if we have both control and treatment samples
      if (length(unique(treatment_factor)) < 2) {
        cat("Error: Could not identify both control and treatment groups. Cannot proceed with batch correction.\n")
        return(NULL)
      }
      
      # Create a model matrix for the treatment factor
      # This will preserve the contrast between control and treatment
      mod <- model.matrix(~treatment_factor)
      
      # Print model matrix
      cat("Model matrix:\n")
      print(mod)
      
      # Ensure counts is a matrix
      counts_matrix <- as.matrix(counts)
      
      # Ensure batch is a factor
      batch_factor <- as.factor(batch_info$batch)
      
      # Print dimensions for debugging
      cat(sprintf("Counts matrix dimensions: %d x %d\n", nrow(counts_matrix), ncol(counts_matrix)))
      cat(sprintf("Number of batches: %d\n", length(unique(batch_factor))))
      
      # ComBat-seq batch correction with model matrix
      # Set group=NULL as we're using mod to preserve treatment differences
      cat("Running ComBat-seq with model matrix...\n")
      corrected_counts <- ComBat_seq(counts = counts_matrix,
                                    batch = batch_factor,
                                    group = NULL,
                                    covar_mod = mod)
      
      # Print a sample of the raw and corrected counts for comparison
      cat("Sample comparison of raw vs. corrected counts (first 5 genes, first 5 samples):\n")
      sample_comparison <- data.frame(
        raw = counts_matrix[1:min(5, nrow(counts_matrix)), 1:min(5, ncol(counts_matrix))],
        corrected = corrected_counts[1:min(5, nrow(corrected_counts)), 1:min(5, ncol(corrected_counts))]
      )
      print(sample_comparison)
    } else {
      cat("Performing ComBat-seq batch correction without group information...\n")
      # Ensure counts is a matrix
      counts_matrix <- as.matrix(counts)
      
      # Ensure batch is a factor
      batch_factor <- as.factor(batch_info$batch)
      
      # Print dimensions for debugging
      cat(sprintf("Counts matrix dimensions: %d x %d\n", nrow(counts_matrix), ncol(counts_matrix)))
      cat(sprintf("Number of batches: %d\n", length(unique(batch_factor))))
      
      # ComBat-seq batch correction without group information
      corrected_counts <- ComBat_seq(counts = counts_matrix,
                                    batch = batch_factor)
    }
    
    # Check if correction made any changes
    if (identical(as.matrix(counts), corrected_counts)) {
      cat("Warning: Batch correction did not change the counts. This might indicate an issue with the correction process.\n")
      
      # Try an alternative approach if the first one didn't work
      cat("Trying alternative batch correction approach...\n")
      
      # Create a simple model matrix with just an intercept
      simple_mod <- matrix(1, nrow = ncol(counts_matrix), ncol = 1)
      
      # Run ComBat-seq with the simple model matrix
      corrected_counts <- ComBat_seq(counts = counts_matrix,
                                    batch = batch_factor,
                                    group = NULL,
                                    covar_mod = simple_mod)
      
      # Check again if correction made any changes
      if (identical(as.matrix(counts), corrected_counts)) {
        cat("Warning: Alternative batch correction also did not change the counts.\n")
      } else {
        # Calculate the difference between raw and corrected counts
        diff_counts <- sum(abs(as.matrix(counts) - corrected_counts))
        cat(sprintf("Total absolute difference between raw and corrected counts: %f\n", diff_counts))
      }
    } else {
      # Calculate the difference between raw and corrected counts
      diff_counts <- sum(abs(as.matrix(counts) - corrected_counts))
      cat(sprintf("Total absolute difference between raw and corrected counts: %f\n", diff_counts))
    }
  } else {
    # Skip batch correction if we only have one batch
    cat("Skipping batch correction (only one batch detected)...\n")
    corrected_counts <- as.matrix(counts)
  }
  
  # Save results
  cat("Saving batch-corrected counts...\n")
  results <- list(
    raw_counts = counts,
    corrected_counts = corrected_counts,
    batch_info = batch_info
  )
  
  output_file <- file.path(output_dir, "batch_corrected_counts.rds")
  saveRDS(results, output_file)
  cat(sprintf("Saved batch-corrected counts to: %s\n", output_file))
  
  return(results)
}

# Main workflow
main <- function() {
  # Set up directories
  base_dir <- getwd()
  output_dir <- file.path(base_dir, "results", "batch_correction")
  design_dir <- file.path(base_dir, "sample_design", "sample_design", "design_matrices")
  
  cat("Setting up directories...\n")
  cat(sprintf("Base directory: %s\n", base_dir))
  cat(sprintf("Output directory: %s\n", output_dir))
  cat(sprintf("Design matrices directory: %s\n", design_dir))
  
  # Check if batch correction has already been performed
  if (file.exists(file.path(output_dir, "batch_corrected_counts.rds"))) {
    cat("Batch correction has already been performed, skipping\n")
    return(0)
  }
  
  # Create SRR to GSM mapping
  srr_to_gsm <- create_srr_gsm_mapping(base_dir)
  
  # Load design matrices
  design_matrices <- load_design_matrices(design_dir)
  
  # Collect all feature counts
  all_data <- collect_feature_counts(base_dir, design_matrices, srr_to_gsm)
  
  if (is.null(all_data)) {
    cat("Failed to collect feature counts, cannot proceed\n")
    return(1)
  }
  
  # Perform batch correction
  results <- perform_batch_correction(
    all_data$counts,
    all_data$batch_info,
    output_dir
  )
  
  cat("=== RNA-seq Batch Correction Pipeline Completed Successfully ===\n")
  cat("Batch correction has been performed. Please run normalize_counts.R next.\n")
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 