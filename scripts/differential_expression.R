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
  "data.table",
  "plotly",  # Add plotly for interactive plots
  "dplyr",   # Add dplyr for data manipulation
  "htmlwidgets"  # For saving interactive plots
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
    
    # Print the original group assignments for debugging
    cat("Original group assignments:\n")
    print(sample_to_group$Group)
    
    # Use the original group names without simplification
    group_factor <- factor(sample_to_group$Group)
    
    # Print the group factor for debugging
    cat("Group factor:\n")
    print(group_factor)
    cat("Group factor levels:\n")
    print(levels(group_factor))
    
    # Print a table of the group assignments
    cat("Group assignments table:\n")
    print(table(group_factor))
  } else {
    cat("No valid samples found with matching GSM IDs in design matrix. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Create DGEList object
  cat("Creating DGEList object...\n")
  
  # Based on the normalization approach in normalize_counts.R,
  # the counts are log2-transformed CPM values with gene length correction
  cat("Converting log2-transformed CPM values back to raw counts...\n")
  
  # Instead of trying to convert back to raw counts (which can lead to extremely large values),
  # we'll use a more direct approach for edgeR analysis
  
  # First, check if we have raw counts available in the normalized data
  if ("raw_counts" %in% names(normalized_data)) {
    cat("Using raw counts from normalized data...\n")
    raw_counts <- normalized_data$raw_counts
    
    # Print summary of raw counts
    cat("Raw counts summary:\n")
    print(summary(as.vector(raw_counts)))
    
    # Create DGEList with raw counts
    y <- DGEList(counts = raw_counts)
  } else {
    # If raw counts are not available, we'll use a more conservative approach
    cat("Raw counts not available. Using a more conservative approach...\n")
    
    # First, convert from log2(CPM/kb) to CPM/kb
    cpm_kb <- 2^counts
    
    # Then, multiply by gene length (in kb) to get CPM
    if ("gene_lengths" %in% names(normalized_data)) {
      gene_lengths_kb <- normalized_data$gene_lengths$length / 1000
      
      # Ensure gene lengths match the counts matrix
      if (length(gene_lengths_kb) == nrow(counts)) {
        # Convert from CPM/kb to CPM
        cpm <- sweep(cpm_kb, 1, gene_lengths_kb, "*")
        
        # Instead of converting to raw counts, we'll use the CPM values directly
        # This is a reasonable approach for edgeR when we don't have raw counts
        cat("Using CPM values directly for edgeR analysis...\n")
        
        # Print summary of CPM values
        cat("CPM values summary:\n")
        print(summary(as.vector(cpm)))
        
        # Create DGEList with CPM values
        y <- DGEList(counts = cpm)
        
        # Set the library sizes to 1e6 to indicate these are CPM values
        y$samples$lib.size <- rep(1e6, ncol(y))
      } else {
        cat("Error: Gene lengths do not match counts matrix dimensions\n")
        cat(sprintf("Gene lengths: %d, Counts rows: %d\n", length(gene_lengths_kb), nrow(counts)))
        return(FALSE)
      }
    } else {
      cat("Error: Gene lengths not found in normalized data\n")
      return(FALSE)
    }
  }
  
  # Filter low count genes with more stringent criteria
  cat("Filtering low count genes with more stringent criteria...\n")
  
  # First, check if we're using raw counts or CPM values
  if (all(y$samples$lib.size == 1e6)) {
    # We're using CPM values
    cat("Using CPM-specific filtering criteria...\n")
    
    # For CPM values, we'll use a minimum CPM threshold
    # Keep genes with at least 1 CPM in at least 2 samples
    min_cpm <- 1
    min_samples <- 2
    
    # Calculate CPM values (they should already be in CPM format)
    cpm_values <- y$counts
    
    # Filter based on CPM threshold
    keep <- rowSums(cpm_values >= min_cpm) >= min_samples
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat(sprintf("Kept %d genes out of %d (CPM >= %g in at least %d samples)\n", 
                sum(keep), length(keep), min_cpm, min_samples))
  } else {
    # We're using raw counts
    cat("Using raw count-specific filtering criteria...\n")
    
    # For raw counts, we'll use edgeR's filterByExpr with more stringent parameters
    # Keep genes with at least 10 counts in at least 2 samples
    min_count <- 10
    min_samples <- 2
    
    # Use filterByExpr with more stringent criteria
    keep <- filterByExpr(y, min.count = min_count, min.total.count = 0, 
                         large.n = 10, min.samples = min_samples)
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat(sprintf("Kept %d genes out of %d (count >= %d in at least %d samples)\n", 
                sum(keep), length(keep), min_count, min_samples))
  }
  
  # If we're still keeping too many genes, apply additional filtering
  if (sum(keep) > 20000) {
    cat("Still keeping too many genes. Applying additional filtering...\n")
    
    # Calculate coefficient of variation (CV) for each gene
    # Genes with high CV are more likely to be differentially expressed
    if (all(y$samples$lib.size == 1e6)) {
      # For CPM values
      cpm_values <- y$counts
      cv <- apply(cpm_values, 1, function(x) sd(x)/mean(x))
    } else {
      # For raw counts, convert to log-CPM first
      log_cpm <- cpm(y$counts, log = TRUE)
      cv <- apply(log_cpm, 1, function(x) sd(x)/mean(x))
    }
    
    # Keep genes with CV above the median
    cv_threshold <- median(cv, na.rm = TRUE)
    keep_cv <- cv >= cv_threshold
    y <- y[keep_cv, , keep.lib.sizes = FALSE]
    cat(sprintf("After CV filtering: Kept %d genes out of %d (CV >= %g)\n", 
                sum(keep_cv), sum(keep), cv_threshold))
  }
  
  # If we're still keeping too few genes, try a more lenient approach
  if (sum(keep) < 1000) {
    cat("Keeping too few genes. Using more lenient filtering criteria...\n")
    
    if (all(y$samples$lib.size == 1e6)) {
      # For CPM values, use a lower threshold
      min_cpm <- 0.5
      min_samples <- 1
      cpm_values <- y$counts
      keep <- rowSums(cpm_values >= min_cpm) >= min_samples
    } else {
      # For raw counts, use a lower count threshold
      min_count <- 5
      min_samples <- 1
      keep <- filterByExpr(y, min.count = min_count, min.total.count = 0, 
                           large.n = 10, min.samples = min_samples)
    }
    
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat(sprintf("With lenient criteria: Kept %d genes out of %d\n", sum(keep), length(keep)))
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
  
  # Find the control group (assuming it contains "control" in the name)
  control_groups <- groups[grep("control", groups, ignore.case = TRUE)]
  treatment_groups <- groups[grep("treatment", groups, ignore.case = TRUE)]
  
  cat("DEBUG: Control groups:", paste(control_groups, collapse=", "), "\n")
  cat("DEBUG: Treatment groups:", paste(treatment_groups, collapse=", "), "\n")
  
  if (length(control_groups) > 0 && length(treatment_groups) > 0) {
    # Use the first control group as the reference
    control_group <- control_groups[1]
    cat(sprintf("Using %s as the control group\n", control_group))
    
    # Create contrasts for each treatment group compared to the control group
    cat("DEBUG: Starting loop through treatment groups\n")
    
    # Create contrast for treatment1 vs control1
    if ("treatment1" %in% treatment_groups) {
      contrast_name <- "treatment1_vs_control1"
      cat(sprintf("Creating contrast: %s\n", contrast_name))
      
      tryCatch({
        contrasts[[contrast_name]] <- makeContrasts(
          treatment1 - control1,
          levels = design
        )
        cat(sprintf("DEBUG: Successfully created contrast for %s\n", contrast_name))
      }, error = function(e) {
        cat(sprintf("DEBUG: Error creating contrast: %s\n", e$message))
      })
    }
    
    # Create contrast for treatment2 vs control1
    if ("treatment2" %in% treatment_groups) {
      contrast_name <- "treatment2_vs_control1"
      cat(sprintf("Creating contrast: %s\n", contrast_name))
      
      tryCatch({
        contrasts[[contrast_name]] <- makeContrasts(
          treatment2 - control1,
          levels = design
        )
        cat(sprintf("DEBUG: Successfully created contrast for %s\n", contrast_name))
      }, error = function(e) {
        cat(sprintf("DEBUG: Error creating contrast: %s\n", e$message))
      })
    }
    
    # Create contrast for treatment3 vs control1
    if ("treatment3" %in% treatment_groups) {
      contrast_name <- "treatment3_vs_control1"
      cat(sprintf("Creating contrast: %s\n", contrast_name))
      
      tryCatch({
        contrasts[[contrast_name]] <- makeContrasts(
          treatment3 - control1,
          levels = design
        )
        cat(sprintf("DEBUG: Successfully created contrast for %s\n", contrast_name))
      }, error = function(e) {
        cat(sprintf("DEBUG: Error creating contrast: %s\n", e$message))
      })
    }
    
    # Create contrast for treatment4 vs control1
    if ("treatment4" %in% treatment_groups) {
      contrast_name <- "treatment4_vs_control1"
      cat(sprintf("Creating contrast: %s\n", contrast_name))
      
      tryCatch({
        contrasts[[contrast_name]] <- makeContrasts(
          treatment4 - control1,
          levels = design
        )
        cat(sprintf("DEBUG: Successfully created contrast for %s\n", contrast_name))
      }, error = function(e) {
        cat(sprintf("DEBUG: Error creating contrast: %s\n", e$message))
      })
    }
  } else {
    # If we don't have control and treatment groups, create contrasts for each group compared to the first group
    control_group <- groups[1]  # Assuming the first group is the control
    cat(sprintf("DEBUG: No control/treatment groups found. Using %s as control\n", control_group))
    
    # Create contrasts for each group compared to the control group
    for (j in 2:length(groups)) {
      contrast_name <- paste0(groups[j], "_vs_", control_group)
      cat(sprintf("Creating contrast: %s\n", contrast_name))
      
      # Create the contrast using makeContrasts
      contrast_formula <- paste0(groups[j], " - ", control_group)
      cat(sprintf("DEBUG: Contrast formula: %s\n", contrast_formula))
      
      tryCatch({
        contrasts[[contrast_name]] <- makeContrasts(
          contrast_formula,
          levels = design
        )
        cat(sprintf("DEBUG: Successfully created contrast for %s\n", contrast_name))
      }, error = function(e) {
        cat(sprintf("DEBUG: Error creating contrast: %s\n", e$message))
      })
    }
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
    
    # Create interactive MA plot with plotly
    ma_data <- data.frame(
      logCPM = qlf$table$logCPM,
      logFC = qlf$table$logFC,
      FDR = qlf$table$FDR,
      genes = rownames(qlf$table)
    )
    
    # Add significance information
    ma_data$Significance <- dplyr::case_when(
      ma_data$FDR < 0.05 & ma_data$logFC > 1 ~ "Upregulated",
      ma_data$FDR < 0.05 & ma_data$logFC < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
    
    # Define colors for the points
    color_map <- c("Upregulated" = "#c46666", "Downregulated" = "#1a80bb", "Not Significant" = "#b0b0b0")
    
    # Create interactive MA plot
    ma_plot <- plot_ly(
      data = ma_data,
      x = ~logCPM,
      y = ~logFC,
      text = ~paste("Gene:", genes,
                    "<br>logCPM:", round(logCPM, 2),
                    "<br>logFC:", round(logFC, 2),
                    "<br>FDR:", formatC(FDR, format = "e", digits = 2)),
      hoverinfo = "text",
      mode = "markers",
      marker = list(size = 6)
    ) %>%
      add_markers(
        color = ~Significance,
        colors = color_map
      ) %>%
      layout(
        title = paste("MA Plot:", contrast_name),
        xaxis = list(title = "logCPM"),
        yaxis = list(title = "logFC"),
        shapes = list(
          list(type = "line", x0 = min(ma_data$logCPM), x1 = max(ma_data$logCPM), 
               y0 = 1, y1 = 1, line = list(dash = "dash", color = "black"), opacity = 0.3),
          list(type = "line", x0 = min(ma_data$logCPM), x1 = max(ma_data$logCPM), 
               y0 = -1, y1 = -1, line = list(dash = "dash", color = "black"), opacity = 0.3)
        )
      )
    
    # Save the interactive plot
    htmlwidgets::saveWidget(ma_plot, file.path(output_dir, paste0(gse_id, "_", contrast_name, "_ma_plot.html")), selfcontained = TRUE)
    
    # Create volcano plot
    cat(sprintf("Creating volcano plot for %s...\n", contrast_name))
    
    # Add logP for y axis
    volcano_data <- data.frame(
      logFC = res$table$logFC,
      FDR = res$table$FDR,
      genes = rownames(res$table)
    )
    volcano_data$logP <- -log10(volcano_data$FDR)
    
    # Define significance thresholds
    volcano_data$Significance <- dplyr::case_when(
      volcano_data$FDR < 0.05 & volcano_data$logFC > 1 ~ "Upregulated",
      volcano_data$FDR < 0.05 & volcano_data$logFC < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
    
    # Create interactive volcano plot
    volcano_plot <- plot_ly(
      data = volcano_data,
      x = ~logFC,
      y = ~logP,
      text = ~paste("Gene:", genes,
                    "<br>logFC:", round(logFC, 2),
                    "<br>FDR:", formatC(FDR, format = "e", digits = 2)),
      hoverinfo = "text",
      mode = "markers",
      marker = list(size = 6)
    ) %>%
      add_markers(
        color = ~Significance,
        colors = color_map
      ) %>%
      layout(
        title = paste("Volcano Plot:", contrast_name),
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(FDR)"),
        shapes = list(
          list(type = "line", x0 = -1, x1 = -1, y0 = 0, y1 = max(volcano_data$logP), 
               line = list(dash = "dot"), opacity = 0.5),
          list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = max(volcano_data$logP), 
               line = list(dash = "dot"), opacity = 0.5),
          list(type = "line", x0 = min(volcano_data$logFC), x1 = max(volcano_data$logFC), 
               y0 = -log10(0.05), y1 = -log10(0.05), 
               line = list(dash = "dash", color = "black"), opacity = 0.3)
        )
      )
    
    # Save the interactive plot
    htmlwidgets::saveWidget(volcano_plot, file.path(output_dir, paste0(gse_id, "_", contrast_name, "_volcano_plot.html")), selfcontained = TRUE)
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
      
      # Handle extreme values
      cat("Processing expression data for heatmap...\n")
      cat("Original expression data summary:\n")
      print(summary(as.vector(de_expr)))
      
      # Log transform the data to handle extreme values
      de_expr_log <- log2(de_expr + 1)
      cat("Log-transformed expression data summary:\n")
      print(summary(as.vector(de_expr_log)))
      
      # Scale the log-transformed data
      tryCatch({
        # Replace any infinite values with NA before scaling
        de_expr_log[is.infinite(de_expr_log)] <- NA
        
        # Scale the data, handling NA values
        de_expr_scaled <- t(scale(t(de_expr_log), center = TRUE, scale = TRUE))
        cat("Scaled expression data summary:\n")
        print(summary(as.vector(de_expr_scaled)))
        
        # Calculate row means for sorting, ignoring NA values
        row_means <- rowMeans(de_expr_scaled, na.rm = TRUE)
        # Sort rows by mean expression
        de_expr_scaled <- de_expr_scaled[order(row_means, decreasing = TRUE), ]
        
        # Create annotation for samples
        sample_annotation <- data.frame(
          Group = group_factor,
          row.names = colnames(de_expr_scaled)
        )
        
        # Create heatmap with more robust parameters
        # Save as PDF for static view
        pdf(file.path(output_dir, paste0(gse_id, "_de_heatmap.pdf")))
        pheatmap(de_expr_scaled,
                main = "Differentially Expressed Genes",
                scale = "none",  # Already scaled
                clustering_method = "ward.D2",
                show_rownames = FALSE,
                annotation_col = sample_annotation,
                na_col = "grey")  # Color NA values in grey
        dev.off()
        
        # Create interactive heatmap with plotly
        # Convert to data frame for plotly
        heatmap_df <- as.data.frame(de_expr_scaled)
        heatmap_df$genes <- rownames(heatmap_df)
        
        # Create interactive heatmap
        heatmap_plot <- plot_ly(
          x = colnames(de_expr_scaled),
          y = rownames(de_expr_scaled),
          z = as.matrix(de_expr_scaled),
          type = "heatmap",
          colorscale = "RdBu",  # Red for upregulated, blue for downregulated
          reversescale = TRUE,
          colorbar = list(title = "Z-score")
        ) %>%
          layout(
            title = "Differentially Expressed Genes",
            xaxis = list(title = "Samples"),
            yaxis = list(title = "Genes")
          )
        
        # Save the interactive plot
        htmlwidgets::saveWidget(heatmap_plot, file.path(output_dir, paste0(gse_id, "_de_heatmap.html")), selfcontained = TRUE)
        
        cat(sprintf("Created heatmap with %d differentially expressed genes\n", length(all_de_genes)))
      }, error = function(e) {
        cat(sprintf("Error during heatmap creation: %s\n", e$message))
        cat("Expression data statistics:\n")
        cat("Number of genes:", nrow(de_expr_scaled), "\n")
        cat("Number of samples:", ncol(de_expr_scaled), "\n")
        cat("Number of NA values:", sum(is.na(de_expr_scaled)), "\n")
        cat("Number of infinite values:", sum(is.infinite(de_expr_scaled)), "\n")
        
        # Try an alternative approach for heatmap creation
        cat("Trying alternative approach for heatmap creation...\n")
        tryCatch({
          # Replace NA values with 0 for visualization purposes
          de_expr_scaled_alt <- de_expr_scaled
          de_expr_scaled_alt[is.na(de_expr_scaled_alt)] <- 0
          
          # Create heatmap with alternative approach
          pdf(file.path(output_dir, paste0(gse_id, "_de_heatmap_alt.pdf")))
          pheatmap(de_expr_scaled_alt,
                  main = "Differentially Expressed Genes (NA values set to 0)",
                  scale = "none",
                  clustering_method = "ward.D2",
                  show_rownames = FALSE,
                  annotation_col = sample_annotation)
          dev.off()
          cat("Created alternative heatmap with NA values set to 0\n")
        }, error = function(e2) {
          cat(sprintf("Error during alternative heatmap creation: %s\n", e2$message))
        })
        
        return(FALSE)
      })
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