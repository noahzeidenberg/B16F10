# Modified Differential Expression Analysis Script
# This script performs differential expression analysis using limma instead of edgeR
# It can handle negative counts by using log-transformed data

# Load required packages
cat("=== Starting Modified Differential Expression Analysis ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
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
  cat(sprintf("Performing modified differential expression analysis for %s...\n", gse_id))
  
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
  results_file <- file.path(output_dir, paste0(gse_id, "_de_results_modified.rds"))
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
  # This is a simplified approach - in a real scenario, you might need a more robust mapping
  cat("Creating mapping between SRA IDs and GEO accession IDs...\n")
  
  # Extract SRA IDs from the sample names (assuming format like SRR12345678_Aligned.sortedByCoord.out.bam)
  sra_ids <- gsub("_Aligned.sortedByCoord.out.bam", "", samples)
  cat("Extracted SRA IDs:\n")
  cat(paste(sra_ids, collapse = ", "), "\n")
  
  # Try to find a mapping file
  mapping_file <- file.path(base_dir, "sample_design", "sample_design", "sra_to_geo_mapping.txt")
  if (file.exists(mapping_file)) {
    cat(sprintf("Found mapping file at: %s\n", mapping_file))
    mapping <- read.table(mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Check if we have the necessary columns
    if ("SRA_ID" %in% colnames(mapping) && "GEO_ID" %in% colnames(mapping)) {
      cat("Mapping file has required columns.\n")
      
      # Filter mapping to only include SRA IDs in our samples
      mapping <- mapping[mapping$SRA_ID %in% sra_ids, ]
      
      if (nrow(mapping) > 0) {
        cat(sprintf("Found %d mappings between SRA IDs and GEO IDs.\n", nrow(mapping)))
        
        # Filter design matrix to only include GEO IDs in our mapping
        design_matrix <- design_matrix[design_matrix$Sample_geo_accession %in% mapping$GEO_ID, ]
        
        if (nrow(design_matrix) > 0) {
          cat(sprintf("Filtered design matrix to %d samples.\n", nrow(design_matrix)))
        } else {
          cat("No matching GEO IDs found in design matrix. Cannot proceed.\n")
          return(FALSE)
        }
      } else {
        cat("No matching SRA IDs found in mapping file. Cannot proceed.\n")
        return(FALSE)
      }
    } else {
      cat("Mapping file does not have required columns. Cannot proceed.\n")
      return(FALSE)
    }
  } else {
    cat("No mapping file found. Attempting to use all samples in design matrix...\n")
    
    # If we don't have a mapping file, we'll try to use all samples in the design matrix
    # This is a fallback approach and might not work in all cases
    if (ncol(counts) != nrow(design_matrix)) {
      cat(sprintf("Number of samples in counts (%d) does not match number of samples in design matrix (%d). Cannot proceed.\n", 
                  ncol(counts), nrow(design_matrix)))
      return(FALSE)
    }
  }
  
  # Check for negative values in counts
  negative_count <- sum(counts < 0)
  cat(sprintf("Found %d negative values in normalized counts\n", negative_count))
  
  if (negative_count > 0) {
    cat("Handling negative values by adding a constant to make all values non-negative...\n")
    # Add a constant to make all values non-negative
    min_value <- min(counts)
    if (min_value < 0) {
      counts <- counts - min_value + 1
      cat("Added constant to make all values non-negative.\n")
    }
  }
  
  # Log-transform the counts
  cat("Log-transforming the counts...\n")
  # Add a small constant to avoid log(0)
  log_counts <- log2(counts + 1)
  
  # Create design matrix for limma
  cat("Creating design matrix for limma...\n")
  
  # Create factor from group information
  group_factor <- factor(design_matrix$Group)
  
  # Create design matrix
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  
  # Fit linear model
  cat("Fitting linear model...\n")
  fit <- lmFit(log_counts, design)
  
  # Create contrasts for all pairwise comparisons
  cat("Creating contrasts...\n")
  groups <- levels(group_factor)
  contrasts <- list()
  
  if (length(groups) > 1) {
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        contrast_name <- paste0(groups[j], "_vs_", groups[i])
        contrasts[[contrast_name]] <- makeContrasts(
          paste0("group_factor", groups[j], " - group_factor", groups[i]),
          levels = design
        )
      }
    }
  }
  
  # Perform differential expression analysis for each contrast
  results <- list()
  
  for (contrast_name in names(contrasts)) {
    cat(sprintf("Testing contrast: %s\n", contrast_name))
    
    # Test for differential expression
    fit2 <- contrasts.fit(fit, contrasts[[contrast_name]])
    fit2 <- eBayes(fit2)
    
    # Get results
    res <- topTable(fit2, coef = 1, number = Inf)
    
    # Add to results list
    results[[contrast_name]] <- res
    
    # Save results
    res_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_de_results_modified.txt"))
    write.table(res, res_file, sep = "\t", quote = FALSE)
    cat(sprintf("Saved results to: %s\n", res_file))
    
    # Create MA plot
    cat(sprintf("Creating MA plot for %s...\n", contrast_name))
    pdf(file.path(output_dir, paste0(gse_id, "_", contrast_name, "_ma_plot_modified.pdf")))
    plotMA(fit2, coef = 1)
    dev.off()
    
    # Create volcano plot
    cat(sprintf("Creating volcano plot for %s...\n", contrast_name))
    pdf(file.path(output_dir, paste0(gse_id, "_", contrast_name, "_volcano_plot_modified.pdf")))
    with(res, plot(logFC, -log10(adj.P.Val), pch = 20, main = paste("Volcano Plot:", contrast_name),
                  xlim = c(-5, 5), ylim = c(0, 10)))
    with(subset(res, adj.P.Val < 0.05 & abs(logFC) > 1), points(logFC, -log10(adj.P.Val), pch = 20, col = "red"))
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
    
    # Get differentially expressed genes (adj.P.Val < 0.05 and |logFC| > 1)
    de_genes <- list()
    for (contrast_name in names(results)) {
      de_genes[[contrast_name]] <- rownames(results[[contrast_name]])[
        results[[contrast_name]]$adj.P.Val < 0.05 & abs(results[[contrast_name]]$logFC) > 1
      ]
    }
    
    # Get unique differentially expressed genes
    all_de_genes <- unique(unlist(de_genes))
    
    if (length(all_de_genes) > 0) {
      # Extract expression data for differentially expressed genes
      de_expr <- log_counts[all_de_genes, ]
      
      # Scale for better visualization
      de_expr_scaled <- t(scale(t(de_expr)))
      
      # Create heatmap
      pdf(file.path(output_dir, paste0(gse_id, "_de_heatmap_modified.pdf")))
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
      stop("Please provide a valid GSE ID. For example: Rscript differential_expression_modified.R GSE223515")
    }
  }
  
  tryCatch({
    cat(sprintf("Performing modified differential expression analysis for %s...\n", gse_id))
    success <- perform_de_analysis(gse_id)
    
    if (success) {
      cat(sprintf("Successfully performed modified differential expression analysis for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to perform modified differential expression analysis for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during modified differential expression analysis: %s\n", e$message))
    stop(sprintf("Error during modified differential expression analysis: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 