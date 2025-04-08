# Differential Expression Analysis Script
# This script performs differential expression analysis using edgeR

# Load required packages
cat("=== Starting Differential Expression Analysis ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "pheatmap",
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  cat(sprintf("Loading package: %s\n", pkg))
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

# Function to perform differential expression analysis
perform_de_analysis <- function(gse_id) {
  cat(sprintf("Performing differential expression analysis for %s...\n", gse_id))
  
  # Set up directories
  # Use path.expand to properly handle the tilde in the path
  base_dir <- path.expand("~/scratch/B16F10")
  cat(sprintf("Base directory: %s\n", base_dir))
  
  gse_dir <- file.path(base_dir, gse_id)
  cat(sprintf("GSE directory: %s\n", gse_dir))
  
  # Check if GSE directory exists
  if (!dir.exists(gse_dir)) {
    cat(sprintf("GSE directory %s does not exist. Cannot proceed.\n", gse_dir))
    return(FALSE)
  }
  
  normalization_dir <- file.path(gse_dir, "results", "normalization")
  cat(sprintf("Normalization directory: %s\n", normalization_dir))
  
  # Check if normalization directory exists
  if (!dir.exists(normalization_dir)) {
    cat(sprintf("Normalization directory %s does not exist. Cannot proceed.\n", normalization_dir))
    return(FALSE)
  }
  
  output_dir <- file.path(gse_dir, "results", "differential_expression")
  cat(sprintf("Output directory: %s\n", output_dir))
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if analysis has already been performed
  results_file <- file.path(output_dir, paste0(gse_id, "_de_results.rds"))
  if (file.exists(results_file)) {
    cat(sprintf("Differential expression analysis already performed for %s. Skipping...\n", gse_id))
    return(TRUE)
  }
  
  # List files in normalization directory for debugging
  cat("Files in normalization directory:\n")
  list.files(normalization_dir, full.names = TRUE)
  
  # Load normalized counts
  normalized_file <- file.path(normalization_dir, "normalized_counts.rds")
  cat(sprintf("Looking for normalized counts at: %s\n", normalized_file))
  
  if (!file.exists(normalized_file)) {
    cat(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
    return(FALSE)
  }
  
  cat("Loading normalized counts...\n")
  normalized_data <- readRDS(normalized_file)
  
  # Print structure of normalized data for debugging
  cat("Structure of normalized data:\n")
  str(normalized_data)
  
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
  cat("Sample names:\n")
  cat(paste(samples, collapse = ", "), "\n")
  
  # Find design matrix file
  design_file <- find_design_matrix(base_dir, gse_id)
  if (is.null(design_file)) {
    return(FALSE)
  }
  
  cat("Loading design matrix...\n")
  design_matrix <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Print design matrix structure for debugging
  cat("Structure of design matrix:\n")
  str(design_matrix)
  
  # Check if we have the necessary data
  if (is.null(counts) || nrow(counts) == 0) {
    cat("No counts found. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Check if we have group information
  if (!"Group" %in% colnames(design_matrix)) {
    cat("No 'Group' column found in design matrix. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Create a mapping between SRA IDs and GEO accession IDs
  cat("Creating mapping between SRA IDs and GEO accession IDs...\n")
  
  # Extract SRA IDs and GSM IDs from the sample names
  # Check for new format (GSEID_GSMID_SRRID) first
  sra_ids <- character(length(samples))
  gsm_ids <- character(length(samples))
  
  for (i in 1:length(samples)) {
    sample_name <- samples[i]
    
    # Check if the sample name follows the new format (GSEID_GSMID_SRRID)
    if (grepl(paste0("^", gse_id, "_GSM[0-9]+_SRR"), sample_name)) {
      # Extract GSM ID and SRR ID from the new format
      parts <- strsplit(sample_name, "_")[[1]]
      gsm_ids[i] <- parts[2]  # GSM ID is the second part
      sra_ids[i] <- parts[3]  # SRR ID is the third part
    } else {
      # Fall back to the old format (GSEID_SRRID)
      sra_ids[i] <- gsub("_Aligned.sortedByCoord.out.bam", "", sample_name)
      sra_ids[i] <- gsub(paste0(gse_id, "_"), "", sra_ids[i])
      gsm_ids[i] <- NA  # We'll need to look up the GSM ID in the mapping file
    }
  }
  
  cat("Extracted SRA IDs:\n")
  cat(paste(sra_ids, collapse = ", "), "\n")
  
  cat("Extracted GSM IDs:\n")
  cat(paste(gsm_ids, collapse = ", "), "\n")
  
  # Try to find a mapping file for samples with missing GSM IDs
  mapping_file <- file.path(base_dir, "sample_design", "sample_design", "sra_to_geo_mapping.txt")
  if (file.exists(mapping_file) && any(is.na(gsm_ids))) {
    cat(sprintf("Found mapping file at: %s\n", mapping_file))
    mapping <- read.table(mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Check if we have the necessary columns
    if ("SRA_ID" %in% colnames(mapping) && "GEO_ID" %in% colnames(mapping)) {
      cat("Mapping file has required columns.\n")
      
      # Filter mapping to only include SRA IDs in our samples that don't have GSM IDs
      missing_gsm_indices <- which(is.na(gsm_ids))
      missing_sra_ids <- sra_ids[missing_gsm_indices]
      mapping <- mapping[mapping$SRA_ID %in% missing_sra_ids, ]
      
      if (nrow(mapping) > 0) {
        cat(sprintf("Found %d mappings between SRA IDs and GEO IDs.\n", nrow(mapping)))
        
        # Update GSM IDs for samples with missing GSM IDs
        for (i in 1:nrow(mapping)) {
          sra_id <- mapping$SRA_ID[i]
          geo_id <- mapping$GEO_ID[i]
          
          # Find the index of this SRA ID in our samples
          sra_index <- which(sra_ids == sra_id)
          if (length(sra_index) > 0) {
            gsm_ids[sra_index] <- geo_id
          }
        }
        
        cat("Updated GSM IDs from mapping file:\n")
        cat(paste(gsm_ids, collapse = ", "), "\n")
      }
    }
  }
  
  # Create a mapping from sample indices to groups
  sample_to_group <- data.frame(
    Sample_Index = integer(),
    Group = character(),
    stringsAsFactors = FALSE
  )
  
  # Match samples to their groups using GSM IDs
  for (i in 1:length(samples)) {
    gsm_id <- gsm_ids[i]
    
    if (!is.na(gsm_id)) {
      # Find the row in the design matrix that matches this GSM ID
      match_idx <- which(design_matrix$Sample_geo_accession == gsm_id)
      if (length(match_idx) > 0) {
        group <- design_matrix$Group[match_idx]
        sample_to_group <- rbind(sample_to_group, data.frame(
          Sample_Index = i,
          Group = group,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Filter counts to only include samples with valid groups
  if (nrow(sample_to_group) > 0) {
    valid_indices <- sample_to_group$Sample_Index
    valid_samples <- samples[valid_indices]
    
    cat(sprintf("Filtering counts to %d samples with valid groups.\n", length(valid_samples)))
    counts <- counts[, valid_samples]
    
    # Simplify group names by removing numbers
    simplified_groups <- gsub("[0-9]", "", sample_to_group$Group)
    
    # Create a factor for the simplified groups
    group_factor <- factor(simplified_groups)
    
    # Print the group factor for debugging
    cat("Simplified group factor:\n")
    print(group_factor)
    cat("Simplified group factor levels:\n")
    print(levels(group_factor))
  } else {
    cat("No valid samples found with matching GSM IDs in design matrix. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Create DGEList object
  cat("Creating DGEList object...\n")
  
  # Check for negative values in counts
  if (any(counts < 0)) {
    cat("Warning: Negative values detected in counts. This is expected for log-transformed data.\n")
    cat("Converting to raw counts before creating DGEList...\n")
    
    # Print summary of negative values
    cat("Summary of negative values:\n")
    print(summary(as.vector(counts[counts < 0])))
    
    # Convert log-transformed counts back to raw counts
    # For log-transformed data, we need to reverse the transformation
    # If the data was log-transformed with log=TRUE in cpm(), we need to use exp()
    raw_counts <- exp(counts)
    
    # Print summary of raw counts
    cat("Summary of raw counts after conversion:\n")
    print(summary(as.vector(raw_counts)))
    
    # Check for extreme values
    if (any(raw_counts > 1e6)) {
      cat("Warning: Extreme values detected in raw counts. This might cause issues.\n")
      cat("Summary of extreme values:\n")
      print(summary(as.vector(raw_counts[raw_counts > 1e6])))
      
      # Cap extreme values
      raw_counts[raw_counts > 1e6] <- 1e6
      cat("Capped extreme values to 1e6.\n")
    }
    
    # Create DGEList with raw counts
    y <- DGEList(counts = raw_counts)
  } else {
    # Create DGEList with original counts
    y <- DGEList(counts = counts)
  }
  
  # Filter low count genes with very lenient criteria
  cat("Filtering low count genes with very lenient criteria...\n")
  # Use very lenient filtering criteria to keep more genes
  keep <- filterByExpr(y, min.count = 0, min.total.count = 0, large.n = 10, min.samples = 1)
  y <- y[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("Kept %d genes out of %d\n", sum(keep), length(keep)))
  
  # If we're still keeping too few genes, try an alternative approach
  if (sum(keep) < 1000) {
    cat("Keeping too few genes with filterByExpr. Using alternative approach...\n")
    # Keep genes with at least one count in at least one sample
    keep <- rowSums(y$counts > 0) >= 1
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat(sprintf("Kept %d genes out of %d using alternative approach\n", sum(keep), length(keep)))
  }
  
  # Create design matrix for edgeR
  cat("Creating design matrix for edgeR...\n")
  
  # Print group factor levels for debugging
  cat("Group factor levels:\n")
  print(levels(group_factor))
  
  # Create design matrix
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  
  # Print design matrix for debugging
  cat("Design matrix structure:\n")
  str(design)
  cat("Number of samples in DGEList:", ncol(y), "\n")
  cat("Number of rows in design matrix:", nrow(design), "\n")
  
  # Ensure the design matrix has the correct structure
  if (ncol(design) != length(levels(group_factor))) {
    cat("Warning: Design matrix columns do not match group factor levels. Adjusting...\n")
    # Recreate the design matrix with explicit levels
    design <- model.matrix(~0 + factor(group_factor, levels = levels(group_factor)))
    colnames(design) <- levels(group_factor)
  }
  
  # Estimate dispersion
  cat("Estimating dispersion...\n")
  y <- estimateDisp(y, design)
  
  # Plot dispersion
  cat("Plotting dispersion...\n")
  pdf(file.path(output_dir, paste0(gse_id, "_dispersion_plot.pdf")))
  plotBCV(y)
  dev.off()
  
  # Fit model with error handling
  cat("Fitting model...\n")
  tryCatch({
    # Check if we have enough samples per group
    group_counts <- table(group_factor)
    cat("Sample counts per group:\n")
    print(group_counts)
    
    # Check if we have at least one sample in each group
    if (min(group_counts) < 1) {
      cat("Error: Not enough samples in each group. Cannot fit model.\n")
      return(FALSE)
    }
    
    # Check if we have enough genes
    if (nrow(y) < 10) {
      cat("Error: Not enough genes after filtering. Cannot fit model.\n")
      return(FALSE)
    }
    
    # Fit the model
    fit <- glmQLFit(y, design)
    cat("Model fitting successful.\n")
  }, error = function(e) {
    cat(sprintf("Error during model fitting: %s\n", e$message))
    cat("Design matrix structure:\n")
    str(design)
    cat("DGEList structure:\n")
    str(y)
    return(FALSE)
  })
  
  # Create contrasts for all pairwise comparisons
  cat("Creating contrasts...\n")
  contrasts <- list()
  
  # Get the group levels from the design matrix column names
  groups <- colnames(design)
  cat(sprintf("Groups found: %s\n", paste(groups, collapse = ", ")))
  
  if (length(groups) > 1) {
    # Check if we have control and treatment groups
    if ("control" %in% groups && "treatment" %in% groups) {
      # Create contrast for control vs treatment
      contrast_name <- "treatment_vs_control"
      cat(sprintf("Creating contrast: %s\n", contrast_name))
      
      # Create the contrast using makeContrasts
      contrast_formula <- "treatment - control"
      cat(sprintf("Contrast formula: %s\n", contrast_formula))
      
      contrasts[[contrast_name]] <- makeContrasts(
        contrast_formula,
        levels = design
      )
      
      # Check if we have unknown groups
      if ("unknown" %in% groups) {
        # Create contrast for control vs treatment+unknown
        contrast_name <- "treatment_unknown_vs_control"
        cat(sprintf("Creating contrast: %s\n", contrast_name))
        
        # Create the contrast using makeContrasts
        contrast_formula <- "(treatment + unknown) - control"
        cat(sprintf("Contrast formula: %s\n", contrast_formula))
        
        contrasts[[contrast_name]] <- makeContrasts(
          contrast_formula,
          levels = design
        )
      }
    } else {
      # If we don't have control and treatment groups, create contrasts for each group compared to the first group
      control_group <- groups[1]  # Assuming the first group is the control
      for (j in 2:length(groups)) {
        contrast_name <- paste0(groups[j], "_vs_", control_group)
        cat(sprintf("Creating contrast: %s\n", contrast_name))
        
        # Create the contrast using makeContrasts
        contrast_formula <- paste0(groups[j], " - ", control_group)
        cat(sprintf("Contrast formula: %s\n", contrast_formula))
        
        contrasts[[contrast_name]] <- makeContrasts(
          contrast_formula,
          levels = design
        )
      }
    }
  } else {
    cat("Only one group found in the design matrix. Cannot create contrasts.\n")
  }
  
  # Perform differential expression analysis for each contrast
  results <- list()
  
  for (contrast_name in names(contrasts)) {
    cat(sprintf("Testing contrast: %s\n", contrast_name))
    
    # Test for differential expression
    qlf <- glmQLFTest(fit, contrast = contrasts[[contrast_name]])
    
    # Get results
    res <- topTags(qlf, n = Inf)
    
    # Add to results list
    results[[contrast_name]] <- res$table
    
    # Save results
    res_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_de_results.txt"))
    write.table(res$table, res_file, sep = "\t", quote = FALSE)
    cat(sprintf("Saved results to: %s\n", res_file))
    
    # Create MA plot
    cat(sprintf("Creating MA plot for %s...\n", contrast_name))
    pdf(file.path(output_dir, paste0(gse_id, "_", contrast_name, "_ma_plot.pdf")))
    plotMD(qlf)
    abline(h = c(-1, 1), col = "blue")
    dev.off()
    
    # Create volcano plot
    cat(sprintf("Creating volcano plot for %s...\n", contrast_name))
    pdf(file.path(output_dir, paste0(gse_id, "_", contrast_name, "_volcano_plot.pdf")))
    with(res$table, plot(logFC, -log10(FDR), pch = 20, main = paste("Volcano Plot:", contrast_name),
                        xlim = c(-5, 5), ylim = c(0, 10)))
    with(subset(res$table, FDR < 0.05 & abs(logFC) > 1), points(logFC, -log10(FDR), pch = 20, col = "red"))
    abline(h = -log10(0.05), col = "blue", lty = 2)
    abline(v = c(-1, 1), col = "blue", lty = 2)
    dev.off()
  }
  
  # Save all results
  saveRDS(results, results_file)
  cat(sprintf("Saved all results to: %s\n", results_file))
  
  # Create heatmap of differentially expressed genes
  if (length(results) > 0) {
    cat("Creating heatmap of differentially expressed genes...\n")
    
    # Get differentially expressed genes (FDR < 0.05 and |logFC| > 1)
    de_genes <- list()
    for (contrast_name in names(results)) {
      de_genes[[contrast_name]] <- rownames(results[[contrast_name]])[
        results[[contrast_name]]$FDR < 0.05 & abs(results[[contrast_name]]$logFC) > 1
      ]
    }
    
    # Get unique differentially expressed genes
    all_de_genes <- unique(unlist(de_genes))
    
    if (length(all_de_genes) > 0) {
      # Extract expression data for differentially expressed genes
      de_expr <- y$counts[all_de_genes, ]
      
      # Scale for better visualization
      # Note: Since we're working with log-transformed values, we can use scale() directly
      de_expr_scaled <- t(scale(t(de_expr)))
      
      # Create heatmap
      pdf(file.path(output_dir, paste0(gse_id, "_de_heatmap.pdf")))
      pheatmap(de_expr_scaled, 
               main = "Differentially Expressed Genes",
               scale = "none",
               clustering_method = "ward.D2",
               show_rownames = FALSE)
      dev.off()
      cat(sprintf("Created heatmap with %d differentially expressed genes\n", length(all_de_genes)))
    } else {
      cat("No differentially expressed genes found for any contrast\n")
    }
  }
  
  return(TRUE)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript differential_expression.R GSE223515")
    }
  }
  
  tryCatch({
    cat(sprintf("Performing differential expression analysis for %s...\n", gse_id))
    success <- perform_de_analysis(gse_id)
    
    if (success) {
      cat(sprintf("Successfully performed differential expression analysis for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to perform differential expression analysis for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during differential expression analysis: %s\n", e$message))
    stop(sprintf("Error during differential expression analysis: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 