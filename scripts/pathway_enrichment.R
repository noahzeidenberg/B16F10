# Pathway Enrichment Analysis Script
# This script performs KEGG pathway and Gene Ontology enrichment analysis using clusterProfiler

# Load required packages
cat("=== Starting Pathway Enrichment Analysis ===\n")
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
  cat(sprintf("Package %s not found. Attempting to install...\n", pkg))
  
  # Check if we're offline
  is_offline <- !curl::has_internet()
  if (is_offline) {
    cat("No internet connection detected. Cannot install packages.\n")
    cat("Please ensure all required packages are installed before running in offline mode.\n")
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
    cat(sprintf("Error installing package %s: %s\n", pkg, e$message))
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
  "curl" = TRUE              # CRAN - for internet connectivity check
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
  
  # Check if differential expression analysis has been performed
  de_dir <- file.path(gse_dir, "results", "differential_expression")
  if (!dir.exists(de_dir)) {
    cat(sprintf("Differential expression directory %s does not exist. Cannot proceed.\n", de_dir))
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
    cat(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
    return(FALSE)
  }
  
  cat("Loading normalized counts...\n")
  normalized_data <- readRDS(normalized_file)
  
  # Print structure of normalized data for debugging
  cat("Structure of normalized data:\n")
  str(normalized_data)
  
  # Extract gene IDs from the normalized counts
  if (!"normalized_counts" %in% names(normalized_data)) {
    cat("Error: 'normalized_counts' not found in normalized data. Available keys:\n")
    cat(paste(names(normalized_data), collapse = ", "), "\n")
    return(FALSE)
  }
  
  # Get all gene IDs from the normalized counts
  all_genes <- rownames(normalized_data$normalized_counts)
  cat(sprintf("Found %d genes in normalized counts\n", length(all_genes)))
  
  # Load differential expression results
  de_results_file <- file.path(de_dir, paste0(gse_id, "_de_results.rds"))
  if (!file.exists(de_results_file)) {
    cat(sprintf("Differential expression results not found at %s. Cannot proceed.\n", de_results_file))
    return(FALSE)
  }
  
  cat("Loading differential expression results...\n")
  de_results <- readRDS(de_results_file)
  
  # Check if we have any differential expression results
  if (length(de_results) == 0) {
    cat("No differential expression results found. Cannot proceed.\n")
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
    
    # Convert gene IDs from ENSEMBL to ENTREZID
    cat("Converting gene IDs from ENSEMBL to ENTREZID...\n")
    
    # Try to identify the gene ID format
    gene_id_format <- "ENSEMBL"  # Default to ENSEMBL
    
    # Print some statistics about the gene IDs
    cat("\nAnalyzing gene ID format:\n")
    ensembl_pattern <- sum(grepl("^ENS", de_genes))
    symbol_pattern <- sum(grepl("^[A-Za-z]", de_genes))
    cat(sprintf("Number of genes matching ENSEMBL pattern (^ENS): %d\n", ensembl_pattern))
    cat(sprintf("Number of genes matching SYMBOL pattern (^[A-Za-z]): %d\n", symbol_pattern))
    cat(sprintf("Total number of genes: %d\n", length(de_genes)))
    
    # Print some example gene IDs
    cat("\nFirst few gene IDs:\n")
    print(head(de_genes))
    cat("\nLast few gene IDs:\n")
    print(tail(de_genes))
    
    # Check if the gene IDs are in ENSEMBL format
    if (!all(grepl("^ENS", de_genes))) {
      # Check if the gene IDs are in SYMBOL format
      if (all(grepl("^[A-Za-z]", de_genes))) {
        gene_id_format <- "SYMBOL"
        cat("\nGene IDs appear to be in SYMBOL format.\n")
      } else {
        cat("\nWarning: Gene IDs do not appear to be in ENSEMBL or SYMBOL format.\n")
        cat("Attempting to identify format by checking a sample of IDs:\n")
        
        # Take a sample of gene IDs to analyze
        sample_size <- min(10, length(de_genes))
        sample_ids <- head(de_genes, n = sample_size)
        cat("Sample IDs:\n")
        print(sample_ids)
        
        # Try to detect any common patterns
        patterns <- list(
          ensembl = "^ENS[A-Z]*[0-9]+",
          symbol = "^[A-Za-z][A-Za-z0-9]*$",
          entrez = "^[0-9]+$"
        )
        
        pattern_matches <- sapply(patterns, function(p) sum(grepl(p, sample_ids)))
        cat("\nPattern matches in sample:\n")
        print(pattern_matches)
        
        # Use the most common pattern
        best_pattern <- names(which.max(pattern_matches))
        cat(sprintf("\nMost common pattern appears to be: %s\n", best_pattern))
        
        if (best_pattern == "symbol") {
          gene_id_format <- "SYMBOL"
        } else if (best_pattern == "entrez") {
          gene_id_format <- "ENTREZID"
        } else {
          gene_id_format <- "ENSEMBL"
        }
        
        cat(sprintf("Using %s as the gene ID format.\n", gene_id_format))
      }
    } else {
      cat("\nGene IDs appear to be in ENSEMBL format.\n")
    }
    
    # Function to print gene ID conversion results
    print_conversion_results <- function(result, original_ids) {
      cat("\nGene ID conversion results:\n")
      cat(sprintf("Original number of genes: %d\n", length(original_ids)))
      cat(sprintf("Number of genes successfully converted: %d\n", nrow(result)))
      cat(sprintf("Conversion rate: %.2f%%\n", nrow(result) / length(original_ids) * 100))
      
      # Print examples of successful conversions
      cat("\nExamples of successful conversions:\n")
      print(head(result))
      
      # Find and print examples of failed conversions
      failed_ids <- setdiff(original_ids, result[[1]])
      if (length(failed_ids) > 0) {
        cat(sprintf("\nExamples of failed conversions (%d total):\n", length(failed_ids)))
        print(head(failed_ids, n = 10))
      }
    }
    
    # Convert gene IDs with detailed error handling
    gene_list <- NULL
    conversion_error <- NULL
    
    cat(sprintf("\nAttempting to convert gene IDs from %s to ENTREZID...\n", gene_id_format))
    
    tryCatch({
      gene_list <- bitr(de_genes, 
                       fromType = gene_id_format, 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db")
      
      print_conversion_results(gene_list, de_genes)
      
    }, warning = function(w) {
      cat(sprintf("\nWarning during gene ID conversion: %s\n", w$message))
    }, error = function(e) {
      conversion_error <- e$message
      cat(sprintf("\nError during gene ID conversion: %s\n", e$message))
      cat("Trying alternative approach...\n")
      
      # Try with a different gene ID format
      alt_format <- if (gene_id_format == "ENSEMBL") "SYMBOL" else "ENSEMBL"
      cat(sprintf("Attempting conversion with %s format...\n", alt_format))
      
      tryCatch({
        gene_list <<- bitr(de_genes, 
                          fromType = alt_format, 
                          toType = "ENTREZID", 
                          OrgDb = "org.Mm.eg.db")
        
        print_conversion_results(gene_list, de_genes)
        
      }, warning = function(w2) {
        cat(sprintf("\nWarning during alternative gene ID conversion: %s\n", w2$message))
      }, error = function(e2) {
        cat(sprintf("\nError during alternative gene ID conversion: %s\n", e2$message))
        
        # Try one more time with SYMBOL format and more lenient matching
        cat("\nTrying one more time with more lenient matching...\n")
        
        # Clean up gene IDs - remove version numbers and other common modifications
        clean_ids <- gsub("\\..*$", "", de_genes)  # Remove version numbers
        clean_ids <- gsub("_.*$", "", clean_ids)   # Remove suffixes
        
        tryCatch({
          gene_list <<- bitr(clean_ids, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = "org.Mm.eg.db")
          
          print_conversion_results(gene_list, clean_ids)
          
        }, warning = function(w3) {
          cat(sprintf("\nWarning during final gene ID conversion attempt: %s\n", w3$message))
        }, error = function(e3) {
          cat(sprintf("\nError during final gene ID conversion attempt: %s\n", e3$message))
          return(NULL)
        })
      })
    })
    
    # Verify gene_list exists and has the correct structure
    if (is.null(gene_list)) {
      cat("\nFatal error: gene_list is NULL after all conversion attempts\n")
      next
    }
    
    if (!is.data.frame(gene_list)) {
      cat("\nFatal error: gene_list is not a data frame\n")
      cat("gene_list class:", class(gene_list), "\n")
      next
    }
    
    if (!"ENTREZID" %in% colnames(gene_list)) {
      cat("\nFatal error: ENTREZID column not found in gene_list\n")
      cat("Available columns:", paste(colnames(gene_list), collapse = ", "), "\n")
      next
    }
    
    if (nrow(gene_list) == 0) {
      cat("\nFatal error: No genes were successfully converted to ENTREZID\n")
      next
    }
    
    cat(sprintf("\nProceeding with pathway analysis using %d successfully converted genes\n", nrow(gene_list)))
    
    # Perform KEGG pathway enrichment
    cat("Performing KEGG pathway enrichment...\n")
    tryCatch({
      kegg_result <- enrichKEGG(
        gene = gene_list$ENTREZID,
        organism = "mmu",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      )
      
      # Save KEGG results
      kegg_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_kegg_enrichment.txt"))
      write.table(kegg_result@result, kegg_file, sep = "\t", quote = FALSE)
      cat(sprintf("Saved KEGG results to: %s\n", kegg_file))
      
      # Create KEGG dotplot
      kegg_plot_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_kegg_dotplot.pdf"))
      pdf(kegg_plot_file, width = 10, height = 8)
      print(dotplot(kegg_result, showCategory = 20, title = paste("KEGG Pathway Enrichment:", contrast_name)))
      dev.off()
      cat(sprintf("Saved KEGG dotplot to: %s\n", kegg_plot_file))
      
      # Create KEGG barplot
      kegg_bar_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_kegg_barplot.pdf"))
      pdf(kegg_bar_file, width = 10, height = 8)
      print(barplot(kegg_result, showCategory = 20, title = paste("KEGG Pathway Enrichment:", contrast_name)))
      dev.off()
      cat(sprintf("Saved KEGG barplot to: %s\n", kegg_bar_file))
      
      # Create KEGG emapplot (enrichment map)
      if (nrow(kegg_result@result) > 1) {
        kegg_emap_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_kegg_emapplot.pdf"))
        pdf(kegg_emap_file, width = 10, height = 8)
        print(emapplot(kegg_result, showCategory = 20))
        dev.off()
        cat(sprintf("Saved KEGG emapplot to: %s\n", kegg_emap_file))
      }
      
      # Add to results list
      pathway_results[[contrast_name]]$kegg <- kegg_result
    }, error = function(e) {
      cat(sprintf("Error performing KEGG pathway enrichment: %s\n", e$message))
    })
    
    # Perform GO enrichment
    cat("Performing GO enrichment...\n")
    
    # Perform GO enrichment for each ontology (BP, MF, CC)
    for (ont in c("BP", "MF", "CC")) {
      cat(sprintf("Performing GO enrichment for %s ontology...\n", ont))
      
      tryCatch({
        go_result <- enrichGO(
          gene = gene_list$ENTREZID,
          OrgDb = "org.Mm.eg.db",
          ont = ont,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH"
        )
        
        # Skip if no enriched terms
        if (nrow(go_result@result) == 0) {
          cat(sprintf("No enriched GO terms found for %s ontology. Skipping...\n", ont))
          next
        }
        
        # Save GO results
        go_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_go_", ont, "_enrichment.txt"))
        write.table(go_result@result, go_file, sep = "\t", quote = FALSE)
        cat(sprintf("Saved GO %s results to: %s\n", ont, go_file))
        
        # Create GO dotplot
        go_plot_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_go_", ont, "_dotplot.pdf"))
        pdf(go_plot_file, width = 10, height = 8)
        print(dotplot(go_result, showCategory = 20, title = paste("GO", ont, "Enrichment:", contrast_name)))
        dev.off()
        cat(sprintf("Saved GO %s dotplot to: %s\n", ont, go_plot_file))
        
        # Create GO barplot
        go_bar_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_go_", ont, "_barplot.pdf"))
        pdf(go_bar_file, width = 10, height = 8)
        print(barplot(go_result, showCategory = 20, title = paste("GO", ont, "Enrichment:", contrast_name)))
        dev.off()
        cat(sprintf("Saved GO %s barplot to: %s\n", ont, go_bar_file))
        
        # Create GO emapplot (enrichment map)
        if (nrow(go_result@result) > 1) {
          go_emap_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_go_", ont, "_emapplot.pdf"))
          pdf(go_emap_file, width = 10, height = 8)
          print(emapplot(go_result, showCategory = 20))
          dev.off()
          cat(sprintf("Saved GO %s emapplot to: %s\n", ont, go_emap_file))
        }
        
        # Add to results list
        pathway_results[[contrast_name]][[paste0("go_", ont)]] <- go_result
      }, error = function(e) {
        cat(sprintf("Error performing GO %s enrichment: %s\n", ont, e$message))
      })
    }
  }
  
  # Save all results
  saveRDS(pathway_results, results_file)
  cat(sprintf("Saved all pathway enrichment results to: %s\n", results_file))
  
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
  
  # Create a log file
  log_file <- file.path(getwd(), paste0(gse_id, "_pathway_enrichment.log"))
  cat(sprintf("Starting pathway enrichment analysis for %s at %s\n", gse_id, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_file)
  
  # Redirect output to both console and log file
  sink(log_file, append = TRUE)
  
  tryCatch({
    cat(sprintf("Performing pathway enrichment analysis for %s...\n", gse_id))
    
    # Check if the GSE directory exists
    base_dir <- path.expand("~/scratch/B16F10")
    gse_dir <- file.path(base_dir, gse_id)
    
    if (!dir.exists(gse_dir)) {
      cat(sprintf("GSE directory %s does not exist. Cannot proceed.\n", gse_dir))
      return(FALSE)
    }
    
    # Perform the analysis
    success <- perform_pathway_enrichment(gse_id)
    
    if (success) {
      cat(sprintf("Successfully performed pathway enrichment analysis for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to perform pathway enrichment analysis for %s\n", gse_id))
    }
    
    # Close the sink
    sink()
    
    # Print a summary to the console
    cat(sprintf("Analysis for %s completed. Check the log file at %s for details.\n", gse_id, log_file))
    
    return(success)
  }, error = function(e) {
    # Log the error
    cat(sprintf("Error during pathway enrichment analysis: %s\n", e$message))
    cat("Stack trace:\n")
    cat(paste(capture.output(traceback()), collapse = "\n"))
    
    # Close the sink
    sink()
    
    # Print a summary to the console
    cat(sprintf("Analysis for %s failed with error: %s\n", gse_id, e$message))
    cat(sprintf("Check the log file at %s for details.\n", log_file))
    
    return(FALSE)
  }, finally = {
    # Ensure the sink is closed even if there's an error
    if (sink.number() > 0) {
      sink()
    }
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 