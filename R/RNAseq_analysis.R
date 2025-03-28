#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("GEOquery", "DESeq2", "clusterProfiler", "org.Mm.eg.db", "tidyverse", "dotenv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Load environment variables from .env file
load_dot_env()

# Function to get next API key
get_next_api_key <- function() {
  # Get all API keys from environment
  api_keys <- c(
    Sys.getenv("NCBI_API_KEY_1"),
    Sys.getenv("NCBI_API_KEY_2")
  )
  
  # Remove any NULL or empty values
  api_keys <- api_keys[!is.null(api_keys) & api_keys != ""]
  
  if (length(api_keys) == 0) {
    stop("No NCBI API keys found in environment variables. Please check your .env file.")
  }
  
  # Get the current key from environment
  current_key <- Sys.getenv("ENTREZ_KEY")
  
  # If no current key, use the first one
  if (current_key == "") {
    return(api_keys[1])
  }
  
  # Find the current key's position
  current_pos <- which(api_keys == current_key)
  
  # If current key not found or it's the last one, return the first key
  if (length(current_pos) == 0 || current_pos == length(api_keys)) {
    return(api_keys[1])
  }
  
  # Otherwise, return the next key
  return(api_keys[current_pos + 1])
}

# Function to handle rate limiting with exponential backoff and API key rotation
handle_rate_limit <- function(fn, max_retries = 3, initial_delay = 1) {
  last_error <- NULL
  for (i in 1:max_retries) {
    tryCatch({
      return(fn())
    }, error = function(e) {
      last_error <- e
      if (grepl("429|403|404|Failed to perform HTTP request|cannot open the connection", e$message)) {
        delay <- initial_delay * (2^(i-1))  # Exponential backoff
        message(paste("Rate limit or connection error hit. Waiting", delay, "seconds before retry", i, "of", max_retries))
        
        # Rotate to next API key
        new_key <- get_next_api_key()
        Sys.setenv(ENTREZ_KEY = new_key)
        message(paste("Rotating to next API key:", substr(new_key, 1, 8), "..."))
        
        # Add a small random delay to prevent thundering herd
        Sys.sleep(delay + runif(1, 0, 1))
      } else {
        stop(e)
      }
    })
  }
  
  # If we get here, all retries failed
  if (!is.null(last_error)) {
    message(paste("All retries failed. Last error:", last_error$message))
    stop(last_error)
  }
  stop("Max retries reached")
}

# Initialize with first API key
Sys.setenv(ENTREZ_KEY = get_next_api_key())

# Set options for better error handling
options(warn = 1)  # Show warnings immediately
options(timeout = 300)  # Increase timeout to 5 minutes

# Read the GDS table subset
gds_df <- read.csv("gds_table_new.csv") # from create_gds_table.R

# Function to check if a GSE is RNA-seq
is_rnaseq <- function(gse_id) {
  tryCatch({
    gse <- handle_rate_limit(function() {
      getGEO(gse_id)
    })
    # Check if it's RNA-seq by looking at the technology type and library strategy
    tech <- unique(unlist(lapply(gse, function(x) x@experimentData@other$type)))
    lib_strategy <- unique(unlist(lapply(gse, function(x) x@phenoData@data$library_strategy)))
    return(any(grepl("RNA-seq|RNA sequencing", tech, ignore.case = TRUE)) ||
           any(grepl("RNA-Seq", lib_strategy, ignore.case = TRUE)))
  }, error = function(e) {
    message(paste("Error checking RNA-seq status for", gse_id, ":", e$message))
    return(FALSE)
  })
}

# Function to process a single RNA-seq dataset
process_rnaseq <- function(gse_id) {
  message(paste("Processing", gse_id))
  
  tryCatch({
    # Download the data with rate limiting
    gse <- handle_rate_limit(function() {
      getGEO(gse_id)
    })
    
    # Extract expression matrix and phenotype data
    expr_matrix <- exprs(gse[[1]])
    pheno_data <- pData(gse[[1]])
    
    # Check if expression matrix is empty or all zeros
    if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
      message("Expression matrix is empty")
      return(NULL)
    }
    
    if (all(expr_matrix == 0)) {
      message("All values in expression matrix are zero")
      # Try to get supplementary files with counts
      supp_files <- gse[[1]]@experimentData@other$supplementary_file
      if (!is.null(supp_files)) {
        message("Found supplementary files")
        message(paste("Files:", paste(supp_files, collapse="\n    ")))
        
        # Look for count files
        count_files <- supp_files[grep("counts|Counts|RAW", supp_files)]
        if (length(count_files) > 0) {
          message("Found count files")
          # Download and process count files
          for (file in count_files) {
            message(paste("Downloading:", file))
            tryCatch({
              download.file(file, destfile = file.path(results_dir, basename(file)), mode = "wb")
              # If it's a tar file, extract it
              if (grepl("\\.tar$", file)) {
                system(paste("tar -xf", file.path(results_dir, basename(file)), "-C", results_dir))
              }
            }, error = function(e) {
              message(paste("Failed to download:", file, "-", e$message))
            })
          }
          return(NULL)  # Skip DESeq2 analysis for now
        }
      }
      return(NULL)
    }
    
    # Get experiment metadata
    exp_metadata <- gse[[1]]@experimentData
    
    # Create results directory for this dataset
    results_dir <- file.path("results", gse_id)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save metadata
    metadata_list <- list(
      experiment = list(
        title = exp_metadata@title,
        abstract = exp_metadata@abstract,
        pubMedIds = exp_metadata@pubMedIds,
        type = exp_metadata@other$type,
        platform_id = exp_metadata@other$platform_id,
        overall_design = exp_metadata@other$overall_design
      ),
      samples = pheno_data
    )
    
    saveRDS(metadata_list, file = file.path(results_dir, paste0(gse_id, "_metadata.rds")))
    
    # Check for supplementary files with counts
    supp_files <- exp_metadata@other$supplementary_file
    if (!is.null(supp_files)) {
      message("Found supplementary files")
      message(paste("Files:", paste(supp_files, collapse="\n    ")))
      
      # Look for count files
      count_files <- supp_files[grep("counts|Counts|RAW", supp_files)]
      if (length(count_files) > 0) {
        message("Found count files")
        # Download and process count files
        for (file in count_files) {
          message(paste("Downloading:", file))
          tryCatch({
            download.file(file, destfile = file.path(results_dir, basename(file)), mode = "wb")
            # If it's a tar file, extract it
            if (grepl("\\.tar$", file)) {
              system(paste("tar -xf", file.path(results_dir, basename(file)), "-C", results_dir))
            }
          }, error = function(e) {
            message(paste("Failed to download:", file, "-", e$message))
          })
        }
      }
    }
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = round(expr_matrix),
      colData = pheno_data,
      design = ~ 1  # You'll need to adjust this based on your experimental design
    )
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds)
    
    # Convert to data frame
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    # Perform KEGG pathway analysis
    gene_list <- res_df$gene[!is.na(res_df$padj) & res_df$padj < 0.05]
    
    if (length(gene_list) > 0) {
      kegg_result <- enrichKEGG(
        gene = gene_list,
        organism = 'mmu',
        keyType = 'kegg',
        pvalueCutoff = 0.05
      )
      
      # Save results
      write.csv(res_df, file = file.path(results_dir, paste0(gse_id, "_DESeq2_results.csv")))
      write.csv(as.data.frame(kegg_result), file = file.path(results_dir, paste0(gse_id, "_KEGG_results.csv")))
      
      # Save expression matrix and phenotype data
      write.csv(expr_matrix, file = file.path(results_dir, paste0(gse_id, "_expression_matrix.csv")))
      write.csv(pheno_data, file = file.path(results_dir, paste0(gse_id, "_phenotype_data.csv")))
      
      # Save DESeq2 object
      saveRDS(dds, file = file.path(results_dir, paste0(gse_id, "_DESeq2_object.rds")))
      
      message(paste("Successfully processed and saved results for", gse_id))
    } else {
      message(paste("No significant genes found for", gse_id))
    }
    
    return(list(deseq_results = res_df, kegg_results = kegg_result))
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
    return(NULL)
  })
}

# Create main results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Process each dataset
results <- list()
for (gse_id in gds_df$Accession) {
  tryCatch({
    if (is_rnaseq(gse_id)) {
      message(paste("\nProcessing RNA-seq dataset:", gse_id))
      results[[gse_id]] <- process_rnaseq(gse_id)
    } else {
      message(paste("\nSkipping non-RNA-seq dataset:", gse_id))
    }
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
  })
}

# Save session info
sink("results/session_info.txt")
sessionInfo()
sink() 