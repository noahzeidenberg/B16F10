# Aggregate Pathway Enrichment Results
# This script aggregates pathway enrichment results from multiple GSE datasets
# to identify common pathways and patterns across different treatments in B16F10 cells

# Load required packages
cat("=== Starting Pathway Results Aggregation ===\n")
cat("Loading required packages...\n")

# Function to safely load packages with offline fallback
safe_load_package <- function(pkg, bioc_pkg = FALSE) {
  cat(sprintf("Loading package: %s\n", pkg))
  
  # First try to load the package directly
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    return(TRUE)
  }
  
  # If direct loading fails, try to install
  message(sprintf("Package %s not found. Attempting to install...\n", pkg))
  
  # Check if we're offline
  is_offline <- !curl::has_internet()
  if (is_offline) {
    message("No internet connection detected. Cannot install packages.")
    message("Please ensure all required packages are installed before running in offline mode.")
    return(FALSE)
  }
  
  # Try to install the package
  tryCatch({
    if (bioc_pkg) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error installing package %s: %s\n", pkg, e$message))
    return(FALSE)
  })
}

# List of required packages with their source
required_packages <- list(
  "dplyr" = FALSE,          # CRAN
  "tidyr" = FALSE,          # CRAN
  "ggplot2" = FALSE,        # CRAN
  "stringr" = FALSE,        # CRAN
  "curl" = FALSE,           # CRAN - for internet connectivity check
  "pheatmap" = FALSE,       # CRAN - for heatmap visualization
  "RColorBrewer" = FALSE,   # CRAN - for color palettes
  "org.Mm.eg.db" = TRUE,    # Bioconductor - for gene annotations
  "clusterProfiler" = TRUE  # Bioconductor - for pathway analysis
)

# Load curl first to check internet connectivity
if (!safe_load_package("curl", bioc_pkg = FALSE)) {
  cat("Warning: Could not load 'curl' package. Internet connectivity check will be skipped.\n")
}

# Load all packages
packages_loaded <- TRUE
for (pkg in names(required_packages)) {
  if (!safe_load_package(pkg, bioc_pkg = required_packages[[pkg]])) {
    cat(sprintf("Failed to load package: %s\n", pkg))
    packages_loaded <- FALSE
  }
}

# Check if all packages were loaded successfully
if (!packages_loaded) {
  cat("Some required packages could not be loaded. The script may not function correctly.\n")
  cat("Please ensure all required packages are installed before running the script.\n")
  cat("You can install them manually with:\n")
  cat("install.packages(c('dplyr', 'tidyr', 'ggplot2', 'stringr', 'curl', 'pheatmap', 'RColorBrewer'))\n")
  cat("BiocManager::install(c('org.Mm.eg.db', 'clusterProfiler'))\n")
  
  # Ask if the user wants to continue anyway
  cat("Do you want to continue anyway? (y/n): ")
  user_input <- readline()
  if (tolower(user_input) != "y") {
    stop("Script execution aborted due to missing packages.")
  }
}

# Function to find all GSE directories
find_gse_directories <- function(base_dir) {
  cat(sprintf("Searching for GSE directories in %s...\n", base_dir))
  
  # List all directories in the base directory
  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  
  # Filter for directories that match the GSE pattern
  gse_dirs <- all_dirs[grep("^GSE[0-9]+$", basename(all_dirs))]
  
  cat(sprintf("Found %d GSE directories\n", length(gse_dirs)))
  return(gse_dirs)
}

# Function to check if a GSE has pathway enrichment results
has_pathway_results <- function(gse_dir) {
  results_file <- file.path(gse_dir, "results", "pathway_enrichment", 
                           paste0(basename(gse_dir), "_pathway_enrichment_results.rds"))
  return(file.exists(results_file))
}

# Function to load pathway enrichment results for a GSE
load_pathway_results <- function(gse_dir) {
  gse_id <- basename(gse_dir)
  results_file <- file.path(gse_dir, "results", "pathway_enrichment", 
                           paste0(gse_id, "_pathway_enrichment_results.rds"))
  
  cat(sprintf("Loading pathway results for %s...\n", gse_id))
  
  tryCatch({
    results <- readRDS(results_file)
    
    # Validate the results structure
    if (is.null(results) || length(results) == 0) {
      message(sprintf("Warning: Empty results for %s. Skipping.\n", gse_id))
      return(NULL)
    }
    
    # Check if any contrasts have valid results
    valid_contrasts <- 0
    for (contrast_name in names(results)) {
      contrast_results <- results[[contrast_name]]
      if (!is.null(contrast_results) && 
          (("kegg" %in% names(contrast_results) && !is.null(contrast_results$kegg) && nrow(contrast_results$kegg@result) > 0) ||
           ("go_BP" %in% names(contrast_results) && !is.null(contrast_results$go_BP) && nrow(contrast_results$go_BP@result) > 0) ||
           ("go_MF" %in% names(contrast_results) && !is.null(contrast_results$go_MF) && nrow(contrast_results$go_MF@result) > 0) ||
           ("go_CC" %in% names(contrast_results) && !is.null(contrast_results$go_CC) && nrow(contrast_results$go_CC@result) > 0))) {
        valid_contrasts <- valid_contrasts + 1
      }
    }
    
    if (valid_contrasts == 0) {
      message(sprintf("Warning: No valid pathway results found for %s. Skipping.\n", gse_id))
      return(NULL)
    }
    
    cat(sprintf("Successfully loaded pathway results for %s with %d valid contrasts\n", gse_id, valid_contrasts))
    return(list(gse_id = gse_id, results = results))
  }, error = function(e) {
    message(sprintf("Error loading pathway results for %s: %s\n", gse_id, e$message))
    return(NULL)
  })
}

# Function to extract pathway data from enrichment results
extract_pathway_data <- function(enrichment_result, gse_id, contrast_name, analysis_type) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    return(NULL)
  }
  
  # Extract the result data frame
  result_df <- enrichment_result@result
  
  # Add metadata columns
  result_df$GSE_ID <- gse_id
  result_df$Contrast <- contrast_name
  result_df$Analysis_Type <- analysis_type
  
  # For KEGG pathways, extract the pathway name from the Description
  if (analysis_type == "KEGG") {
    result_df$Pathway_Name <- result_df$Description
  } else {
    # For GO terms, use the Description as is
    result_df$Pathway_Name <- result_df$Description
  }
  
  # Select and rename columns for consistency
  result_df <- result_df %>%
    select(
      GSE_ID, Contrast, Analysis_Type, Pathway_Name, 
      pvalue, p.adjust, qvalue, Count, GeneRatio, 
      Description, GeneID
    )
  
  return(result_df)
}

# Function to aggregate pathway data across all GSEs
aggregate_pathway_data <- function(pathway_results_list) {
  cat("Aggregating pathway data across all GSEs...\n")
  
  # Initialize empty data frames for each analysis type
  kegg_data <- data.frame()
  go_bp_data <- data.frame()
  go_mf_data <- data.frame()
  go_cc_data <- data.frame()
  
  # Process each GSE's results
  for (gse_result in pathway_results_list) {
    gse_id <- gse_result$gse_id
    results <- gse_result$results
    
    # Process each contrast
    for (contrast_name in names(results)) {
      contrast_results <- results[[contrast_name]]
      
      # Extract KEGG pathway data
      if ("kegg" %in% names(contrast_results)) {
        kegg_df <- extract_pathway_data(contrast_results$kegg, gse_id, contrast_name, "KEGG")
        if (!is.null(kegg_df)) {
          kegg_data <- rbind(kegg_data, kegg_df)
        }
      }
      
      # Extract GO BP data
      if ("go_BP" %in% names(contrast_results)) {
        go_bp_df <- extract_pathway_data(contrast_results$go_BP, gse_id, contrast_name, "GO_BP")
        if (!is.null(go_bp_df)) {
          go_bp_data <- rbind(go_bp_data, go_bp_df)
        }
      }
      
      # Extract GO MF data
      if ("go_MF" %in% names(contrast_results)) {
        go_mf_df <- extract_pathway_data(contrast_results$go_MF, gse_id, contrast_name, "GO_MF")
        if (!is.null(go_mf_df)) {
          go_mf_data <- rbind(go_mf_data, go_mf_df)
        }
      }
      
      # Extract GO CC data
      if ("go_CC" %in% names(contrast_results)) {
        go_cc_df <- extract_pathway_data(contrast_results$go_CC, gse_id, contrast_name, "GO_CC")
        if (!is.null(go_cc_df)) {
          go_cc_data <- rbind(go_cc_data, go_cc_df)
        }
      }
    }
  }
  
  # Combine all data
  all_pathway_data <- rbind(kegg_data, go_bp_data, go_mf_data, go_cc_data)
  
  # Print summary statistics
  cat("\nSummary of aggregated pathway data:\n")
  cat(sprintf("Total number of pathway entries: %d\n", nrow(all_pathway_data)))
  cat(sprintf("Number of unique GSE IDs: %d\n", length(unique(all_pathway_data$GSE_ID))))
  cat(sprintf("Number of unique contrasts: %d\n", length(unique(all_pathway_data$Contrast))))
  cat(sprintf("Number of unique pathways: %d\n", length(unique(all_pathway_data$Pathway_Name))))
  
  # Print breakdown by analysis type
  cat("\nBreakdown by analysis type:\n")
  analysis_counts <- table(all_pathway_data$Analysis_Type)
  for (analysis_type in names(analysis_counts)) {
    cat(sprintf("%s: %d entries\n", analysis_type, analysis_counts[analysis_type]))
  }
  
  return(all_pathway_data)
}

# Function to identify common pathways across contrasts
find_common_pathways <- function(pathway_data, min_occurrences = 2) {
  cat(sprintf("Finding pathways that appear in at least %d contrasts...\n", min_occurrences))
  
  # Count occurrences of each pathway
  pathway_counts <- pathway_data %>%
    group_by(Pathway_Name, Analysis_Type) %>%
    summarise(
      Occurrence_Count = n(),
      Contrast_Count = n_distinct(Contrast),
      GSE_Count = n_distinct(GSE_ID),
      .groups = "drop"
    ) %>%
    arrange(desc(Occurrence_Count))
  
  # Filter for pathways that appear in multiple contrasts
  common_pathways <- pathway_counts %>%
    filter(Occurrence_Count >= min_occurrences)
  
  cat(sprintf("Found %d common pathways\n", nrow(common_pathways)))
  
  return(common_pathways)
}

# Function to create a heatmap of pathway enrichment across contrasts
create_pathway_heatmap <- function(pathway_data, common_pathways, output_dir) {
  cat("Creating pathway enrichment heatmap...\n")
  
  # Filter pathway data to include only common pathways
  filtered_data <- pathway_data %>%
    filter(Pathway_Name %in% common_pathways$Pathway_Name)
  
  # Create a wide format data frame for the heatmap
  heatmap_data <- filtered_data %>%
    select(GSE_ID, Contrast, Pathway_Name, p.adjust) %>%
    mutate(
      Contrast_ID = paste(GSE_ID, Contrast, sep = "_"),
      log_pvalue = -log10(p.adjust)
    ) %>%
    select(Pathway_Name, Contrast_ID, log_pvalue) %>%
    pivot_wider(
      names_from = Contrast_ID,
      values_from = log_pvalue,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  # Set row names to pathway names
  rownames(heatmap_data) <- heatmap_data$Pathway_Name
  heatmap_data$Pathway_Name <- NULL
  
  # Limit to top 50 pathways if there are too many
  if (nrow(heatmap_data) > 50) {
    cat("Limiting heatmap to top 50 pathways by occurrence\n")
    top_pathways <- common_pathways %>%
      arrange(desc(Occurrence_Count)) %>%
      head(50) %>%
      pull(Pathway_Name)
    
    heatmap_data <- heatmap_data[rownames(heatmap_data) %in% top_pathways, ]
  }
  
  # Create the heatmap
  heatmap_file <- file.path(output_dir, "pathway_enrichment_heatmap.pdf")
  pdf(heatmap_file, width = 12, height = 10)
  
  # Use pheatmap for better visualization
  pheatmap::pheatmap(
    as.matrix(heatmap_data),
    scale = "row",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Pathway Enrichment Across B16F10 Experiments",
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  )
  
  dev.off()
  cat(sprintf("Saved pathway enrichment heatmap to: %s\n", heatmap_file))
}

# Function to identify treatment-specific pathways
find_treatment_specific_pathways <- function(pathway_data, output_dir) {
  cat("Identifying treatment-specific pathways...\n")
  
  # Extract treatment information from contrast names
  pathway_data <- pathway_data %>%
    mutate(
      Treatment = str_extract(Contrast, "treatment[0-9]+"),
      Treatment = ifelse(is.na(Treatment), "unknown", Treatment)
    )
  
  # Group pathways by treatment
  treatment_pathways <- pathway_data %>%
    group_by(Treatment, Pathway_Name, Analysis_Type) %>%
    summarise(
      Occurrence_Count = n(),
      Contrast_Count = n_distinct(Contrast),
      GSE_Count = n_distinct(GSE_ID),
      .groups = "drop"
    ) %>%
    arrange(Treatment, desc(Occurrence_Count))
  
  # Save treatment-specific pathways
  treatment_file <- file.path(output_dir, "treatment_specific_pathways.txt")
  write.table(treatment_pathways, treatment_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved treatment-specific pathways to: %s\n", treatment_file))
  
  return(treatment_pathways)
}

# Function to create a summary report
create_summary_report <- function(pathway_data, common_pathways, treatment_pathways, output_dir) {
  cat("Creating summary report...\n")
  
  # Create a summary data frame
  summary_data <- pathway_data %>%
    group_by(Analysis_Type) %>%
    summarise(
      Total_Pathways = n_distinct(Pathway_Name),
      Total_Entries = n(),
      .groups = "drop"
    )
  
  # Create a summary of common pathways by analysis type
  common_summary <- common_pathways %>%
    group_by(Analysis_Type) %>%
    summarise(
      Common_Pathways = n(),
      .groups = "drop"
    )
  
  # Create a summary of treatment-specific pathways
  treatment_summary <- treatment_pathways %>%
    group_by(Treatment) %>%
    summarise(
      Unique_Pathways = n_distinct(Pathway_Name),
      .groups = "drop"
    )
  
  # Create the report
  report_file <- file.path(output_dir, "pathway_analysis_summary.txt")
  sink(report_file)
  
  cat("=== B16F10 Pathway Enrichment Analysis Summary ===\n\n")
  cat(sprintf("Analysis Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  cat("=== Dataset Summary ===\n")
  cat(sprintf("Total number of GSE datasets: %d\n", length(unique(pathway_data$GSE_ID))))
  cat(sprintf("Total number of contrasts: %d\n", length(unique(pathway_data$Contrast))))
  cat(sprintf("Total number of unique pathways: %d\n\n", length(unique(pathway_data$Pathway_Name))))
  
  cat("=== Pathway Analysis Summary ===\n")
  print(summary_data)
  cat("\n")
  
  cat("=== Common Pathways Summary ===\n")
  print(common_summary)
  cat("\n")
  
  cat("=== Treatment-Specific Pathways Summary ===\n")
  print(treatment_summary)
  cat("\n")
  
  cat("=== Top Common Pathways ===\n")
  top_common <- common_pathways %>%
    arrange(desc(Occurrence_Count)) %>%
    head(20)
  print(top_common)
  cat("\n")
  
  cat("=== Top Treatment-Specific Pathways ===\n")
  for (treatment in unique(treatment_pathways$Treatment)) {
    cat(sprintf("\nTreatment: %s\n", treatment))
    top_treatment <- treatment_pathways %>%
      filter(Treatment == treatment) %>%
      arrange(desc(Occurrence_Count)) %>%
      head(10)
    print(top_treatment)
  }
  
  sink()
  cat(sprintf("Saved summary report to: %s\n", report_file))
}

# Main function to aggregate pathway results
aggregate_pathway_results <- function() {
  # Set up directories
  base_dir <- path.expand("~/scratch/B16F10")
  output_dir <- file.path(base_dir, "results", "aggregated_pathway_analysis")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find all GSE directories
  gse_dirs <- find_gse_directories(base_dir)
  
  # Filter for GSEs with pathway enrichment results
  gse_dirs_with_results <- gse_dirs[sapply(gse_dirs, has_pathway_results)]
  
  if (length(gse_dirs_with_results) == 0) {
    message("No GSE directories with pathway enrichment results found.")
    return(FALSE)
  }
  
  cat(sprintf("Found %d GSE directories with pathway enrichment results\n", length(gse_dirs_with_results)))
  
  # Load pathway results for each GSE
  pathway_results_list <- list()
  skipped_gses <- character(0)
  
  for (gse_dir in gse_dirs_with_results) {
    gse_id <- basename(gse_dir)
    result <- load_pathway_results(gse_dir)
    
    if (is.null(result)) {
      skipped_gses <- c(skipped_gses, gse_id)
    } else {
      pathway_results_list[[length(pathway_results_list) + 1]] <- result
    }
  }
  
  # Report skipped GSEs
  if (length(skipped_gses) > 0) {
    cat(sprintf("\nSkipped %d GSEs due to errors or invalid results:\n", length(skipped_gses)))
    cat(paste(skipped_gses, collapse = ", "), "\n")
  }
  
  # Filter out NULL results
  pathway_results_list <- pathway_results_list[!sapply(pathway_results_list, is.null)]
  
  if (length(pathway_results_list) == 0) {
    message("No pathway results could be loaded.")
    return(FALSE)
  }
  
  cat(sprintf("\nSuccessfully loaded pathway results for %d GSEs\n", length(pathway_results_list)))
  
  # Aggregate pathway data
  pathway_data <- aggregate_pathway_data(pathway_results_list)
  
  # Save aggregated pathway data
  data_file <- file.path(output_dir, "aggregated_pathway_data.rds")
  saveRDS(pathway_data, data_file)
  cat(sprintf("Saved aggregated pathway data to: %s\n", data_file))
  
  # Save list of included GSEs
  included_gses <- sapply(pathway_results_list, function(x) x$gse_id)
  included_file <- file.path(output_dir, "included_gses.txt")
  write.table(included_gses, included_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("Saved list of included GSEs to: %s\n", included_file))
  
  # Save list of skipped GSEs
  if (length(skipped_gses) > 0) {
    skipped_file <- file.path(output_dir, "skipped_gses.txt")
    write.table(skipped_gses, skipped_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(sprintf("Saved list of skipped GSEs to: %s\n", skipped_file))
  }
  
  # Find common pathways
  common_pathways <- find_common_pathways(pathway_data, min_occurrences = 2)
  
  # Save common pathways
  common_file <- file.path(output_dir, "common_pathways.txt")
  write.table(common_pathways, common_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved common pathways to: %s\n", common_file))
  
  # Create pathway enrichment heatmap
  create_pathway_heatmap(pathway_data, common_pathways, output_dir)
  
  # Find treatment-specific pathways
  treatment_pathways <- find_treatment_specific_pathways(pathway_data, output_dir)
  
  # Create summary report
  create_summary_report(pathway_data, common_pathways, treatment_pathways, output_dir)
  
  cat("\nPathway analysis aggregation completed successfully!\n")
  return(TRUE)
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  aggregate_pathway_results()
} 