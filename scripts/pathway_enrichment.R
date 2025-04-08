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
    
    # Skip if no differentially expressed genes
    if (length(de_genes) == 0) {
      cat(sprintf("No differentially expressed genes found for %s. Skipping...\n", contrast_name))
      next
    }
    
    # Convert gene IDs from ENSEMBL to ENTREZID
    cat("Converting gene IDs from ENSEMBL to ENTREZID...\n")
    
    # Try to identify the gene ID format
    gene_id_format <- "ENSEMBL"  # Default to ENSEMBL
    
    # Check if the gene IDs are in ENSEMBL format
    if (!all(grepl("^ENS", de_genes))) {
      # Check if the gene IDs are in SYMBOL format
      if (all(grepl("^[A-Za-z]", de_genes))) {
        gene_id_format <- "SYMBOL"
        cat("Gene IDs appear to be in SYMBOL format.\n")
      } else {
        cat("Gene IDs do not appear to be in ENSEMBL or SYMBOL format. Using ENSEMBL as default.\n")
      }
    } else {
      cat("Gene IDs appear to be in ENSEMBL format.\n")
    }
    
    # Convert gene IDs
    tryCatch({
      gene_list <- bitr(de_genes, fromType = gene_id_format, toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
      cat(sprintf("Successfully converted %d genes to ENTREZID\n", nrow(gene_list)))
    }, error = function(e) {
      cat(sprintf("Error converting gene IDs: %s\n", e$message))
      cat("Trying alternative approach...\n")
      
      # Try with a different gene ID format
      if (gene_id_format == "ENSEMBL") {
        gene_id_format <- "SYMBOL"
      } else {
        gene_id_format <- "ENSEMBL"
      }
      
      tryCatch({
        gene_list <- bitr(de_genes, fromType = gene_id_format, toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
        cat(sprintf("Successfully converted %d genes to ENTREZID using %s format\n", nrow(gene_list), gene_id_format))
      }, error = function(e2) {
        cat(sprintf("Error converting gene IDs with alternative approach: %s\n", e2$message))
        return(NULL)
      })
    })
    
    # Skip if gene ID conversion failed
    if (is.null(gene_list) || nrow(gene_list) == 0) {
      cat("Gene ID conversion failed. Skipping this contrast.\n")
      next
    }
    
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