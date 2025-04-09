# Pathway Enrichment Analysis Script
# This script performs KEGG pathway and Gene Ontology enrichment analysis using clusterProfiler

# Load required packages
cat("=== Starting Pathway Enrichment Analysis ===\n")
cat("Loading required packages...\n")

# Get the scratch directory from environment variable
get_base_dir <- function() {
  scratch_dir <- Sys.getenv("SCRATCH")
  if (scratch_dir == "") {
    # Fallback to home directory if not on Graham
    return(path.expand("~/scratch/B16F10"))
  }
  return(file.path(scratch_dir, "B16F10"))
}

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
  "clusterProfiler" = TRUE,  # Bioconductor
  "org.Mm.eg.db" = TRUE,     # Bioconductor
  "enrichplot" = TRUE,       # Bioconductor
  "ggplot2" = TRUE,          # CRAN
  "dplyr" = TRUE,            # CRAN
  "stringr" = TRUE,          # CRAN
  "curl" = TRUE,             # CRAN - for internet connectivity check
  "AnnotationDbi" = TRUE     # Bioconductor
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
  cat("install.packages(c('ggplot2', 'dplyr', 'stringr', 'curl'))\n")
  cat("BiocManager::install(c('clusterProfiler', 'org.Mm.eg.db', 'enrichplot'))\n")
  
  # Ask if the user wants to continue anyway
  cat("Do you want to continue anyway? (y/n): ")
  user_input <- readline()
  if (tolower(user_input) != "y") {
    stop("Script execution aborted due to missing packages.")
  }
}

# Function to perform pathway enrichment analysis
perform_pathway_enrichment <- function(gse_id) {
  cat(sprintf("Performing pathway enrichment analysis for %s...\n", gse_id))
  
  # Set up directories using the new function
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Print directory information for debugging
  cat(sprintf("Base directory: %s\n", base_dir))
  cat(sprintf("GSE directory: %s\n", gse_dir))
  
  # Check if GSE directory exists
  if (!dir.exists(gse_dir)) {
    message(sprintf("GSE directory %s does not exist. Cannot proceed.\n", gse_dir))
    return(FALSE)
  }
  
  # Check if differential expression analysis has been performed
  de_dir <- file.path(gse_dir, "results", "differential_expression")
  if (!dir.exists(de_dir)) {
    message(sprintf("Differential expression directory %s does not exist. Cannot proceed.\n", de_dir))
    return(FALSE)
  }
  
  # Create output directory for pathway enrichment results
  output_dir <- file.path(gse_dir, "results", "pathway_enrichment")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if analysis has already been performed
  results_file <- file.path(output_dir, paste0(gse_id, "_pathway_enrichment_results.rds"))
  if (file.exists(results_file)) {
    cat(sprintf("Pathway enrichment analysis already performed for %s. Skipping...\n", gse_id))
    return(TRUE)
  }
  
  # Load normalized counts to get gene IDs
  normalization_dir <- file.path(gse_dir, "results", "normalization")
  normalized_file <- file.path(normalization_dir, "normalized_counts.rds")
  
  if (!file.exists(normalized_file)) {
    message(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
    return(FALSE)
  }
  
  cat("Loading normalized counts...\n")
  normalized_data <- readRDS(normalized_file)
  
  # Print structure of normalized data for debugging
  cat("Structure of normalized data:\n")
  str(normalized_data)
  
  # Extract gene IDs from the normalized counts
  if (!"normalized_counts" %in% names(normalized_data)) {
    message("Error: 'normalized_counts' not found in normalized data. Available keys:\n")
    cat(paste(names(normalized_data), collapse = ", "), "\n")
    return(FALSE)
  }
  
  # Get all gene IDs from the normalized counts
  all_genes <- rownames(normalized_data$normalized_counts)
  cat(sprintf("Found %d genes in normalized counts\n", length(all_genes)))
  
  # Load differential expression results
  de_results_file <- file.path(de_dir, paste0(gse_id, "_de_results.rds"))
  if (!file.exists(de_results_file)) {
    message(sprintf("Differential expression results not found at %s. Cannot proceed.\n", de_results_file))
    return(FALSE)
  }
  
  cat("Loading differential expression results...\n")
  de_results <- readRDS(de_results_file)
  
  # Check if we have any differential expression results
  if (length(de_results) == 0) {
    message("No differential expression results found. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Initialize results list
  pathway_results <- list()
  
  # Process each contrast
  for (contrast_name in names(de_results)) {
    cat(sprintf("Processing contrast: %s\n", contrast_name))
    
    # Get differentially expressed genes (FDR < 0.05 and |logFC| > 1)
    de_genes <- rownames(de_results[[contrast_name]])[
      de_results[[contrast_name]]$FDR < 0.05 & abs(de_results[[contrast_name]]$logFC) > 1
    ]
    
    cat(sprintf("Found %d differentially expressed genes for %s\n", length(de_genes), contrast_name))
    
    # Print examples of differentially expressed genes
    cat("Examples of differentially expressed genes:\n")
    print(head(de_genes, n = 10))
    
    # Skip if no differentially expressed genes
    if (length(de_genes) == 0) {
      cat(sprintf("No differentially expressed genes found for %s. Skipping...\n", contrast_name))
      next
    }
    
    # Convert gene IDs from SYMBOL to ENTREZID
    cat("Converting gene IDs from SYMBOL to ENTREZID...\n")
    
    # Clean up gene symbols - remove version numbers and other common modifications
    clean_ids <- gsub("\\..*$", "", de_genes)  # Remove version numbers
    clean_ids <- gsub("_.*$", "", clean_ids)   # Remove suffixes
    
    # Try to convert gene symbols to ENTREZ IDs
    gene_list <- NULL
    tryCatch({
      gene_list <- bitr(clean_ids, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = FALSE)  # Don't drop genes that fail to map
      
      # Print conversion statistics
      cat("\nGene ID conversion results:\n")
      cat(sprintf("Original number of genes: %d\n", length(clean_ids)))
      cat(sprintf("Number of genes successfully converted: %d\n", nrow(gene_list)))
      cat(sprintf("Conversion rate: %.2f%%\n", nrow(gene_list) / length(clean_ids) * 100))
      
      # Print examples of successful conversions
      cat("\nExamples of successful conversions:\n")
      print(head(gene_list))
      
      # Find and print examples of failed conversions
      failed_ids <- setdiff(clean_ids, gene_list$SYMBOL)
      if (length(failed_ids) > 0) {
        cat(sprintf("\nExamples of failed conversions (%d total):\n", length(failed_ids)))
        print(head(failed_ids, n = 10))
      }
      
    }, error = function(e) {
      message(sprintf("Error during gene ID conversion: %s\n", e$message))
      return(NULL)
    })
    
    # Verify gene_list exists and has the correct structure
    if (is.null(gene_list) || nrow(gene_list) == 0) {
      message("Fatal error: No genes were successfully converted to ENTREZID")
      next
    }
    
    cat(sprintf("\nProceeding with pathway analysis using %d successfully converted genes\n", nrow(gene_list)))
    
    # Perform KEGG pathway enrichment
    cat("Performing KEGG pathway enrichment...\n")
    tryCatch({
      kegg_result <- clusterProfiler::enrichKEGG(
        gene = gene_list$ENTREZID,
        organism = "mmu",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      )
      
      # Check if we got any results
      if (is.null(kegg_result)) {
        cat("Warning: KEGG enrichment returned NULL result\n")
      } else if (nrow(kegg_result@result) == 0) {
        cat("Warning: No enriched KEGG pathways found\n")
      } else {
        cat(sprintf("Found %d enriched KEGG pathways\n", nrow(kegg_result@result)))
        
        # Print the first few results for debugging
        cat("\nTop KEGG pathways:\n")
        print(head(kegg_result@result))
        
        # Save KEGG results
        kegg_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_kegg_enrichment.txt"))
        write.table(kegg_result@result, kegg_file, sep = "\t", quote = FALSE)
        cat(sprintf("Saved KEGG results to: %s\n", kegg_file))
        
        # Create KEGG dotplot
        if (nrow(kegg_result@result) > 0) {
          kegg_plot_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_kegg_dotplot.pdf"))
          pdf(kegg_plot_file, width = 10, height = 8)
          print(enrichplot::dotplot(kegg_result, showCategory = 20, title = paste("KEGG Pathway Enrichment:", contrast_name)))
          dev.off()
          cat(sprintf("Saved KEGG dotplot to: %s\n", kegg_plot_file))
        }
        
        # Add to results list
        pathway_results[[contrast_name]]$kegg <- kegg_result
      }
    }, error = function(e) {
      message(sprintf("Error performing KEGG pathway enrichment: %s\n", e$message))
    })
    
    # Perform GO enrichment
    cat("Performing GO enrichment...\n")
    
    # Perform GO enrichment for each ontology (BP, MF, CC)
    for (ont in c("BP", "MF", "CC")) {
      cat(sprintf("Performing GO enrichment for %s ontology...\n", ont))
      
      tryCatch({
        go_result <- clusterProfiler::enrichGO(
          gene = gene_list$ENTREZID,
          OrgDb = "org.Mm.eg.db",
          ont = ont,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH"
        )
        
        # Check if we got any results
        if (is.null(go_result)) {
          cat(sprintf("Warning: GO %s enrichment returned NULL result\n", ont))
          next
        }
        
        if (nrow(go_result@result) == 0) {
          cat(sprintf("Warning: No enriched GO terms found for %s ontology\n", ont))
          next
        }
        
        cat(sprintf("Found %d enriched GO terms for %s ontology\n", nrow(go_result@result), ont))
        
        # Print the first few results for debugging
        cat(sprintf("\nTop GO %s terms:\n", ont))
        print(head(go_result@result))
        
        # Save GO results
        go_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_go_", ont, "_enrichment.txt"))
        write.table(go_result@result, go_file, sep = "\t", quote = FALSE)
        cat(sprintf("Saved GO %s results to: %s\n", ont, go_file))
        
        # Create GO dotplot
        go_plot_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_go_", ont, "_dotplot.pdf"))
        pdf(go_plot_file, width = 10, height = 8)
        print(enrichplot::dotplot(go_result, showCategory = 20, title = paste("GO", ont, "Enrichment:", contrast_name)))
        dev.off()
        cat(sprintf("Saved GO %s dotplot to: %s\n", ont, go_plot_file))
        
        # Add to results list
        pathway_results[[contrast_name]][[paste0("go_", ont)]] <- go_result
        
      }, error = function(e) {
        message(sprintf("Error performing GO %s enrichment: %s\n", ont, e$message))
      })
    }
  }
  
  # Check if we got any results before saving
  if (length(pathway_results) == 0) {
    message("Warning: No pathway enrichment results were generated")
    return(FALSE)
  }
  
  # Save all results
  saveRDS(pathway_results, results_file)
  cat(sprintf("\nSaved all pathway enrichment results to: %s\n", results_file))
  
  return(TRUE)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript pathway_enrichment.R GSE223515")
    }
  }
  
  # Set up error handling
  options(warn = 1)  # Print warnings as they occur
  
  tryCatch({
    cat(sprintf("Performing pathway enrichment analysis for %s...\n", gse_id))
    
    # Check if the GSE directory exists
    base_dir <- get_base_dir()
    gse_dir <- file.path(base_dir, gse_id)
    
    if (!dir.exists(gse_dir)) {
      message(sprintf("GSE directory %s does not exist. Cannot proceed.\n", gse_dir))
      return(FALSE)
    }
    
    # Perform the analysis
    success <- perform_pathway_enrichment(gse_id)
    
    if (success) {
      cat(sprintf("Successfully performed pathway enrichment analysis for %s\n", gse_id))
    } else {
      message(sprintf("Failed to perform pathway enrichment analysis for %s\n", gse_id))
    }
    
    return(success)
  }, error = function(e) {
    # Log the error to stderr
    message(sprintf("Error during pathway enrichment analysis: %s\n", e$message))
    message("Stack trace:")
    message(paste(capture.output(traceback()), collapse = "\n"))
    return(FALSE)
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 