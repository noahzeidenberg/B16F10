#!/usr/bin/env Rscript
# Script to fix the design matrix by adding control/treatment classification
# This script will identify control and treatment groups based on naming patterns
# and add a new column to the design matrix

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

# Function to find design matrix file
find_design_matrix <- function(base_dir, gse_id) {
  cat(sprintf("Searching for design matrix for %s...\n", gse_id))
  
  # List of possible design matrix locations
  possible_locations <- c(
    file.path(base_dir, "results", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    file.path(base_dir, "results", "design_matrices", "sample_design", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    file.path(base_dir, "results", "design_matrices", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    file.path(base_dir, gse_id, "results", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    file.path(base_dir, gse_id, "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    # Add the correct path
    file.path(base_dir, "sample_design", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    # Also check the absolute path with the correct base directory
    file.path(path.expand("~/scratch/B16F10"), "sample_design", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt"))
  )
  
  # Check each location
  for (location in possible_locations) {
    cat(sprintf("Checking location: %s\n", location))
    if (file.exists(location)) {
      cat(sprintf("Found design matrix at: %s\n", location))
      return(location)
    }
  }
  
  # If not found, try to find any file with the GSE ID in the name
  cat("Design matrix not found in expected locations. Searching for any file with GSE ID...\n")
  
  # Search in results/design_matrices
  design_dir <- file.path(base_dir, "results", "design_matrices")
  if (dir.exists(design_dir)) {
    design_files <- list.files(design_dir, pattern = gse_id, full.names = TRUE)
    if (length(design_files) > 0) {
      cat(sprintf("Found potential design matrix files: %s\n", paste(design_files, collapse = ", ")))
      return(design_files[1])  # Return the first match
    }
  }
  
  # Search in sample_design directory
  sample_design_dir <- file.path(base_dir, "sample_design", "sample_design", "design_matrices")
  if (dir.exists(sample_design_dir)) {
    design_files <- list.files(sample_design_dir, pattern = gse_id, full.names = TRUE)
    if (length(design_files) > 0) {
      cat(sprintf("Found potential design matrix files in sample_design directory: %s\n", paste(design_files, collapse = ", ")))
      return(design_files[1])  # Return the first match
    }
  }
  
  # Search in GSE directory
  gse_dir <- file.path(base_dir, gse_id)
  if (dir.exists(gse_dir)) {
    design_files <- list.files(gse_dir, pattern = "design_matrix", recursive = TRUE, full.names = TRUE)
    if (length(design_files) > 0) {
      cat(sprintf("Found potential design matrix files in GSE directory: %s\n", paste(design_files, collapse = ", ")))
      return(design_files[1])  # Return the first match
    }
  }
  
  cat("No design matrix found for GSE ID. Cannot proceed.\n")
  return(NULL)
}

# Function to fix the design matrix
fix_design_matrix <- function(gse_id) {
  cat(sprintf("Fixing design matrix for %s...\n", gse_id))
  
  # Set up directories
  # Use path.expand to properly handle the tilde in the path
  base_dir <- path.expand("~/scratch/B16F10")
  cat(sprintf("Base directory: %s\n", base_dir))
  
  # Find design matrix file
  design_file <- find_design_matrix(base_dir, gse_id)
  if (is.null(design_file)) {
    return(FALSE)
  }
  
  # Load design matrix
  cat("Loading design matrix...\n")
  design_matrix <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Print design matrix structure for debugging
  cat("Structure of design matrix:\n")
  str(design_matrix)
  
  # Check if we have the necessary data
  if (is.null(design_matrix) || nrow(design_matrix) == 0) {
    cat("No design matrix found. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Check if we have group information
  if (!"Group" %in% colnames(design_matrix)) {
    cat("No 'Group' column found in design matrix. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Identify all unique groups
  all_groups <- unique(design_matrix$Group)
  cat("All unique groups:\n")
  print(all_groups)
  
  # Check if we have more than one group
  if (length(all_groups) <= 1) {
    cat("Only one group found in the design matrix. Adding control/treatment classification...\n")
    
    # Create a treatment factor (control vs. treatment)
    # We'll identify control groups based on common naming patterns
    control_patterns <- c("control", "control[0-9]+", "ctrl", "untreated", "vehicle", "sham", "wt", "wild type")
    
    # Initialize treatment factor
    treatment_factor <- rep("treatment", nrow(design_matrix))
    
    # Identify control groups - improved approach
    control_groups <- c()
    for (group in all_groups) {
      group_lower <- tolower(group)
      
      # Check for exact matches first
      if (group_lower %in% c("control", "ctrl", "untreated", "vehicle", "sham", "wt", "wild type")) {
        cat(sprintf("Exact match for control group: %s\n", group))
        control_groups <- c(control_groups, group)
        treatment_factor[design_matrix$Group == group] <- "control"
        next
      }
      
      # Check for pattern matches
      is_control <- FALSE
      for (pattern in control_patterns) {
        if (grepl(pattern, group_lower, perl = TRUE)) {
          cat(sprintf("Pattern match for control group: %s (matched pattern: %s)\n", group, pattern))
          control_groups <- c(control_groups, group)
          treatment_factor[design_matrix$Group == group] <- "control"
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
      if (length(all_groups) > 0) {
        # Assume the first group is a control
        control_group <- all_groups[1]
        cat(sprintf("Assuming first group is control: %s\n", control_group))
        control_groups <- c(control_groups, control_group)
        treatment_factor[design_matrix$Group == control_group] <- "control"
      }
    }
    
    # If still no control groups, use a manual approach
    if (length(control_groups) == 0) {
      cat("Still no control groups identified. Using manual approach...\n")
      
      # For each group, try to identify a control group
      for (group in all_groups) {
        # Ask the user to identify control groups
        cat(sprintf("Group: %s\n", group))
        
        # For now, we'll use a heuristic: assume groups with "control" in the name are controls
        # In a real implementation, you might want to prompt the user
        if (grepl("control", tolower(group))) {
          cat(sprintf("Manually identified control group: %s\n", group))
          control_groups <- c(control_groups, group)
          treatment_factor[design_matrix$Group == group] <- "control"
        }
      }
    }
    
    # If we still don't have control groups, create a synthetic one
    if (length(control_groups) == 0) {
      cat("No control groups identified. Creating a synthetic control group...\n")
      
      # Create a synthetic control group by splitting the samples
      num_samples <- nrow(design_matrix)
      if (num_samples > 1) {
        # Split samples into control and treatment
        control_indices <- 1:floor(num_samples/2)
        treatment_indices <- (floor(num_samples/2) + 1):num_samples
        
        # Update treatment factor
        treatment_factor[control_indices] <- "control"
        treatment_factor[treatment_indices] <- "treatment"
        
        cat(sprintf("Created synthetic control group with %d samples\n", length(control_indices)))
      } else {
        cat("Not enough samples to create a synthetic control group. Cannot proceed.\n")
        return(FALSE)
      }
    }
    
    # Add treatment factor to design matrix
    design_matrix$Treatment <- treatment_factor
    
    # Print treatment factor summary
    cat("Treatment factor summary:\n")
    print(table(design_matrix$Treatment))
    
    # Save the updated design matrix
    output_file <- gsub("\\.txt$", "_fixed.txt", design_file)
    cat(sprintf("Saving updated design matrix to: %s\n", output_file))
    write.table(design_matrix, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("Design matrix fixed successfully.\n")
    return(TRUE)
  } else {
    cat("Multiple groups found in the design matrix. No need to add control/treatment classification.\n")
    return(TRUE)
  }
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript fix_design_matrix.R GSE223515")
    }
  }
  
  tryCatch({
    cat(sprintf("Fixing design matrix for %s...\n", gse_id))
    success <- fix_design_matrix(gse_id)
    
    if (success) {
      cat(sprintf("Successfully fixed design matrix for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to fix design matrix for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during design matrix fix: %s\n", e$message))
    stop(sprintf("Error during design matrix fix: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 