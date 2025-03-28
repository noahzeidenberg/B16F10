#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("GEOquery", "DESeq2", "clusterProfiler", "org.Mm.eg.db", "tidyverse", "dotenv", "SRAdb")
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
    Sys.getenv("NCBI_API_KEY_2"),
    Sys.getenv("NCBI_API_KEY_3")  # Added third key
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
      if (grepl("429|403|404|Failed to perform HTTP request|cannot open the connection|Unrecognized option", e$message)) {
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

# Add a small delay between API calls
Sys.sleep(0.1)  # 100ms delay between calls

# Function to check if data is sequencing data
is_sequencing_data <- function(metadata) {
  # Check various metadata fields for sequencing indicators
  sequencing_indicators <- c(
    # Check type field
    any(grepl("sequencing|SRA|RNA-Seq|RNA_SEQ", toupper(metadata$type))),
    # Check library_strategy field
    any(grepl("RNA-Seq|RNA_SEQ|RNA|TRANSCRIPTOMIC", toupper(metadata$library_strategy))),
    # Check library_source field
    any(grepl("TRANSCRIPTOMIC|RNA", toupper(metadata$library_source))),
    # Check instrument_model field for sequencers
    any(grepl("ILLUMINA|NEXTSEQ|HISEQ|NOVASEQ|MISEQ", toupper(metadata$instrument_model)))
  )
  
  return(any(sequencing_indicators))
}

# Function to get SRR accessions from an SRX ID
get_srr_from_srx <- function(srx) {
  tryCatch({
    message(paste("  Getting SRR accessions for", srx))
    
    # Get current API key
    current_key <- Sys.getenv("ENTREZ_KEY")
    
    # Set API key as environment variable for the command
    # Search directly in sra database and get runinfo format
    cmd <- paste("export NCBI_API_KEY='", current_key, "' && ",
                "esearch -db sra -query ", srx, 
                " | efetch -format runinfo", 
                " | cut -d',' -f1", sep="")
    
    message(paste("  Running command:", cmd))
    
    # Execute the command with rate limit handling
    srr_ids <- handle_rate_limit(function() {
      system(cmd, intern = TRUE)
    })
    
    # Remove header row if present
    if (length(srr_ids) > 0 && srr_ids[1] == "Run") {
      srr_ids <- srr_ids[-1]
    }
    
    if (length(srr_ids) == 0) {
      message("  No SRR accessions found")
      return(NULL)
    }
    
    message(paste("  Found SRR accessions:", paste(srr_ids, collapse=", ")))
    return(srr_ids)
  }, error = function(e) {
    message(paste("Warning: Failed to get SRR info for", srx, ":", e$message))
    return(NULL)
  })
}

# Function to get SRX accessions from a GSM ID
get_srx_from_gsm <- function(gsm) {
  tryCatch({
    message(paste("  Fetching data for", gsm))
    
    # Get the GSM data with rate limit handling
    gsm_data <- handle_rate_limit(function() {
      getGEO(gsm)
    })
    
    # Extract the 'relation' field
    relations <- gsm_data@header["relation"][[1]]
    
    # Find the SRA link
    sra_link <- relations[grep("SRA:", relations)]
    
    if (length(sra_link) == 0) {
      message("  No SRA link found")
      return(NULL)
    }
    
    # Extract the SRA accession number (SRX ID)
    sra_acc <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
    message(paste("  Found SRA accession:", sra_acc))
    
    return(sra_acc)
  }, error = function(e) {
    message(paste("Warning: Failed to get SRA info for", gsm, ":", e$message))
    return(NULL)
  })
}

# Function to process a single RNA-seq dataset
process_rnaseq <- function(gse_id) {
  message(paste("Processing", gse_id))
  
  tryCatch({
    # Create results directory for this dataset
    results_dir <- file.path("results", gse_id)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Download the data with rate limiting
    gset <- handle_rate_limit(function() {
      getGEO(gse_id, GSEMatrix = TRUE)
    })
    
    # Save the GEO object
    saveRDS(gset, file = file.path(results_dir, paste0(gse_id, "_geo_object.rds")))
    
    # Extract metadata
    metadata <- pData(gset[[1]])
    write.csv(metadata, file = file.path(results_dir, paste0(gse_id, "_metadata.csv")))
    
    # Get GSM accessions
    gsm_accessions <- metadata$geo_accession
    message(paste("\nFound", length(gsm_accessions), "GSM accessions"))
    
    # Create samples directory
    samples_dir <- file.path(results_dir, "samples")
    dir.create(samples_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get SRX accessions for each GSM
    srx_info <- c()
    for (gsm in gsm_accessions) {
      message(paste("\nProcessing", gsm))
      
      # Create GSM-specific directory
      gsm_dir <- file.path(samples_dir, gsm)
      dir.create(gsm_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Get GSM data and save RDS
      gsm_data <- handle_rate_limit(function() {
        getGEO(gsm)
      })
      saveRDS(gsm_data, file = file.path(gsm_dir, paste0(gsm, ".rds")))
      
      # Get SRX accession
      srx_id <- get_srx_from_gsm(gsm)
      if (!is.null(srx_id)) {
        srx_info <- c(srx_info, srx_id)
      }
      
      # Add a small delay between GSM requests
      Sys.sleep(0.5)
    }
    
    srx_info <- unique(srx_info)
    
    if (length(srx_info) == 0) {
      message("Could not find any SRA run accessions for the GSM IDs")
      return(NULL)
    }
    
    message(paste("\nFound", length(srx_info), "unique SRA accessions"))
    
    # Process each SRX
    for (srx in srx_info) {
      message(paste("\nProcessing", srx))
      
      # Get SRR accessions for this SRX
      srr_ids <- get_srr_from_srx(srx)
      if (is.null(srr_ids)) {
        message(paste("  Skipping", srx, "as no SRR accessions were found"))
        next
      }
      
      # Create SRA directory
      sra_dir <- file.path(results_dir, "SRA")
      fastq_dir <- file.path(sra_dir, "FASTQ")
      dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Download and convert each SRR
      for (srr in srr_ids) {
        message(paste("  Processing SRR:", srr))
        
        # Use prefetch and fasterq-dump for better performance
        cmd <- paste("module load sra-toolkit &&",
                    "prefetch", srr, "&&",
                    "fasterq-dump", srr,
                    "--outdir", fastq_dir,
                    "--split-files",
                    "--threads 8")
        
        message(paste("  Running command:", cmd))
        
        # Execute the command and check for errors
        result <- system(cmd)
        if (result != 0) {
          message(paste("  Warning: Failed to download/convert", srr, "with exit code", result))
          message("  Trying alternative method...")
          
          # Try alternative method using fastq-dump
          alt_cmd <- paste("module load sra-toolkit &&",
                          "fastq-dump --split-files --gzip",
                          "--outdir", fastq_dir, srr)
          message(paste("  Running alternative command:", alt_cmd))
          alt_result <- system(alt_cmd)
          if (alt_result != 0) {
            message(paste("  Error: Both download methods failed for", srr))
          }
        } else {
          # Move the SRA file to the SRA directory
          sra_file <- file.path(".", paste0(srr, ".sra"))
          if (file.exists(sra_file)) {
            file.rename(sra_file, file.path(sra_dir, paste0(srr, ".sra")))
          }
          # Compress the FASTQ files
          message("  Compressing FASTQ files...")
          system(paste("gzip", file.path(fastq_dir, paste0(srr, "_*.fastq"))))
        }
      }
    }
    
    message(paste("\nSuccessfully processed", gse_id))
    return(TRUE)
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
    return(NULL)
  })
}

# Read the GDS table subset
gds_df <- read.csv("gds_table_new.csv")

# Create main results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Process each dataset
results <- list()
for (gse_id in gds_df$Accession) {
  tryCatch({
    message(paste("\nProcessing dataset:", gse_id))
    results[[gse_id]] <- process_rnaseq(gse_id)
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
  })
}

# Save session info
sink("results/session_info.txt")
sessionInfo()
sink()