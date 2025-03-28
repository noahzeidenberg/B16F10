#!/usr/bin/env Rscript

# Load required libraries
library(GEOquery)
library(SRAdb)
library(dotenv)

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
  for (i in 1:max_retries) {
    tryCatch({
      return(fn())
    }, error = function(e) {
      if (grepl("429", e$message)) {
        delay <- initial_delay * (2^(i-1))  # Exponential backoff
        message(paste("Rate limit hit. Waiting", delay, "seconds before retry", i, "of", max_retries))
        
        # Rotate to next API key
        new_key <- get_next_api_key()
        Sys.setenv(ENTREZ_KEY = new_key)
        message(paste("Rotating to next API key:", substr(new_key, 1, 8), "..."))
        
        Sys.sleep(delay)
      } else {
        stop(e)
      }
    })
  }
  stop("Max retries reached")
}

# Initialize with first API key
Sys.setenv(ENTREZ_KEY = get_next_api_key())

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript geo_download.R <accession> <output_dir>")
}

accession <- args[1]
output_dir <- args[2]

# Function to get SRR accessions from an SRX ID
get_srr_from_srx <- function(srx) {
  tryCatch({
    message(paste("  Getting SRR accessions for", srx))
    
    # Use esearch and elink with API key
    cmd <- paste("esearch -api_key $ENTREZ_KEY -db sra -query", srx, 
                "| elink -api_key $ENTREZ_KEY -target sra", 
                "| efetch -api_key $ENTREZ_KEY -format docsum", 
                "| xtract -pattern DocumentSummary -element Run@acc")
    
    message(paste("  Running command:", cmd))
    
    # Execute the command with rate limit handling
    srr_ids <- handle_rate_limit(function() {
      system(cmd, intern = TRUE)
    })
    
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

# Function to check if data is microarray data
is_microarray_data <- function(metadata) {
  # Check various metadata fields for microarray indicators
  microarray_indicators <- c(
    # Check type field
    any(grepl("array|microarray|AFFYMETRIX|ILLUMINA", toupper(metadata$type))),
    # Check platform field
    any(grepl("GPL|AFFYMETRIX|ILLUMINA", toupper(metadata$platform_id))),
    # Check instrument model field
    any(grepl("AFFYMETRIX|ILLUMINA|SCANNER", toupper(metadata$instrument_model))),
    # Check for specific Affymetrix platform identifiers
    any(grepl("HT_MG-430_PM|HT_MG-430A_PM|HT_MG-430B_PM|HT_MG-430_2_PM", toupper(metadata$platform_id))),
    # Check for specific Affymetrix instrument models
    any(grepl("AFFYMETRIX_GENE_ARRAY|AFFYMETRIX_GENE_CHIP|AFFYMETRIX_SCANNER", toupper(metadata$instrument_model)))
  )
  
  # Debug: Print the values being checked
  message("\nChecking microarray indicators:")
  message(paste("Type:", paste(unique(metadata$type), collapse=", ")))
  message(paste("Platform ID:", paste(unique(metadata$platform_id), collapse=", ")))
  message(paste("Instrument Model:", paste(unique(metadata$instrument_model), collapse=", ")))
  
  return(any(microarray_indicators))
}

# Function to download microarray data
download_microarray_data <- function(gsm_data, gsm_dir) {
  tryCatch({
    message(paste("  Downloading microarray data for", gsm_data@header$geo_accession))
    
    # Create necessary directories
    array_dir <- file.path(gsm_dir, "array_data")
    dir.create(array_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get supplementary files
    suppl_files <- gsm_data@header$supplementary_file
    if (!is.null(suppl_files)) {
      message("  Found supplementary files")
      message(paste("  Files:", paste(suppl_files, collapse="\n    ")))
      
      # Look for various array file types (case insensitive)
      cel_files <- suppl_files[grep("\\.(cel|CEL)(\\.gz)?$", suppl_files)]
      chp_files <- suppl_files[grep("\\.(chp|CHP)(\\.gz)?$", suppl_files)]
      idat_files <- suppl_files[grep("\\.(idat|IDAT)(\\.gz)?$", suppl_files)]
      txt_files <- suppl_files[grep("\\.(txt|TXT)(\\.gz)?$", suppl_files)]
      
      # Debug: Print found files
      message("\nFound files:")
      message(paste("CEL files:", paste(cel_files, collapse=", ")))
      message(paste("CHP files:", paste(chp_files, collapse=", ")))
      message(paste("IDAT files:", paste(idat_files, collapse=", ")))
      message(paste("TXT files:", paste(txt_files, collapse=", ")))
      
      files_downloaded <- FALSE
      
      # Try downloading files in order of preference
      if (length(cel_files) > 0) {
        message("  Found CEL files")
        for (file in cel_files) {
          message(paste("  Downloading:", file))
          tryCatch({
            download.file(file, destfile = file.path(array_dir, basename(file)), mode = "wb")
            files_downloaded <- TRUE
          }, error = function(e) {
            message(paste("    Failed to download:", file, "-", e$message))
          })
        }
      }
      
      if (!files_downloaded && length(chp_files) > 0) {
        message("  Found CHP files")
        for (file in chp_files) {
          message(paste("  Downloading:", file))
          tryCatch({
            download.file(file, destfile = file.path(array_dir, basename(file)), mode = "wb")
            files_downloaded <- TRUE
          }, error = function(e) {
            message(paste("    Failed to download:", file, "-", e$message))
          })
        }
      }
      
      if (!files_downloaded && length(idat_files) > 0) {
        message("  Found IDAT files")
        for (file in idat_files) {
          message(paste("  Downloading:", file))
          tryCatch({
            download.file(file, destfile = file.path(array_dir, basename(file)), mode = "wb")
            files_downloaded <- TRUE
          }, error = function(e) {
            message(paste("    Failed to download:", file, "-", e$message))
          })
        }
      }
      
      if (!files_downloaded && length(txt_files) > 0) {
        message("  Found TXT files")
        for (file in txt_files) {
          message(paste("  Downloading:", file))
          tryCatch({
            download.file(file, destfile = file.path(array_dir, basename(file)), mode = "wb")
            files_downloaded <- TRUE
          }, error = function(e) {
            message(paste("    Failed to download:", file, "-", e$message))
          })
        }
      }
      
      if (!files_downloaded) {
        message("  No array data files could be downloaded from supplementary files")
        message("  Checking for raw data in other fields...")
        
        # Try to find raw data in other fields
        raw_data <- gsm_data@header$raw_data
        if (!is.null(raw_data)) {
          message("  Found raw data field")
          message(paste("  Raw data:", paste(raw_data, collapse="\n    ")))
          # Try to download raw data files
          for (file in raw_data) {
            message(paste("  Downloading raw data:", file))
            tryCatch({
              download.file(file, destfile = file.path(array_dir, basename(file)), mode = "wb")
              files_downloaded <- TRUE
            }, error = function(e) {
              message(paste("    Failed to download:", file, "-", e$message))
            })
          }
        }
      }
      
      return(files_downloaded)
    } else {
      message("  No supplementary files found")
      return(FALSE)
    }
  }, error = function(e) {
    message(paste("Warning: Failed to download microarray data:", e$message))
    return(FALSE)
  })
}

# Download GEO data
download_geo_data <- function(accession, output_dir) {
  tryCatch({
    # Get GEO data with rate limit handling
    message("\nDownloading GEO metadata...")
    gset <- handle_rate_limit(function() {
      getGEO(accession, GSEMatrix = TRUE)
    })
    
    # Get platform information in multiple ways
    platform_info <- list()
    
    # Try getting platform from annotation
    if (!is.null(gset[[1]]@annotation) && gset[[1]]@annotation != "") {
      platform_info$annotation <- gset[[1]]@annotation
      message("Platform from annotation: ", platform_info$annotation)
    }
    
    # Try getting platform from platform_id
    if (!is.null(gset[[1]]$platform_id)) {
      platform_info$platform_id <- unique(gset[[1]]$platform_id)
      message("Platform from platform_id: ", paste(platform_info$platform_id, collapse=", "))
    }
    
    # Save platform information
    saveRDS(platform_info, file = file.path(output_dir, paste0(accession, "_platform_info.rds")))
    write.csv(data.frame(
      source = names(platform_info),
      platform = unlist(platform_info)
    ), file = file.path(output_dir, paste0(accession, "_platform_info.csv")))
    
    # Get detailed platform information using the first available platform ID
    platform_id <- NULL
    if (!is.null(platform_info$platform_id)) {
      platform_id <- platform_info$platform_id[1]
    } else if (!is.null(platform_info$annotation)) {
      platform_id <- platform_info$annotation
    }
    
    if (!is.null(platform_id)) {
      message("\nGetting detailed platform information...")
      tryCatch({
        platform_data <- getGEO(platform_id)
        saveRDS(platform_data, file = file.path(output_dir, paste0(platform_id, "_details.rds")))
        
        # Extract manufacturer and technology info
        manufacturer <- platform_data@header$manufacturer
        technology <- platform_data@header$technology
        message("Platform ", platform_id, ":")
        message("  Manufacturer: ", manufacturer)
        message("  Technology: ", technology)
      }, error = function(e) {
        message(paste("Warning: Failed to get detailed platform information:", e$message))
      })
    }
    
    # Continue with existing metadata processing...
    metadata <- pData(gset[[1]])
    write.csv(metadata, file = file.path(output_dir, paste0(accession, "_metadata.csv")))
    
    # Debug: Print metadata structure
    message("\nMetadata structure:")
    message(paste("Columns:", paste(colnames(metadata), collapse=", ")))
    message(paste("Number of samples:", nrow(metadata)))
    
    # Print relevant fields for debugging
    message("\nRelevant metadata fields:")
    relevant_fields <- c("type", "library_strategy", "library_source", "instrument_model", "platform_id")
    for (field in relevant_fields) {
      if (field %in% colnames(metadata)) {
        message(paste(field, ":", paste(unique(metadata[[field]]), collapse=", ")))
      }
    }
    
    # Create samples directory
    samples_dir <- file.path(output_dir, "samples")
    dir.create(samples_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get GSM accessions
    gsm_accessions <- metadata$geo_accession
    message(paste("\nFound", length(gsm_accessions), "GSM accessions:", paste(gsm_accessions, collapse=", ")))
    
    # Check data type and process accordingly
    if (is_sequencing_data(metadata)) {
      message("\nRNA-seq data detected. Fetching SRA information...")
      
      # Get SRX accessions for each GSM with rate limiting
      srx_info <- c()
      for (gsm in gsm_accessions) {
        message(paste("\nProcessing", gsm))
        
        # Create GSM-specific directory
        gsm_dir <- file.path(samples_dir, gsm)
        dir.create(gsm_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Get GSM data and save RDS with rate limiting
        gsm_data <- handle_rate_limit(function() {
          getGEO(gsm)
        })
        saveRDS(gsm_data, file = file.path(gsm_dir, paste0(gsm, ".rds")))
        
        # Get SRX accession with rate limiting
        srx_id <- get_srx_from_gsm(gsm)
        if (!is.null(srx_id)) {
          srx_info <- c(srx_info, srx_id)
        }
        
        # Add a small delay between GSM requests
        Sys.sleep(0.5)
      }
      
      srx_info <- unique(srx_info)
      
      if (length(srx_info) == 0) {
        stop("Could not find any SRA run accessions for the GSM IDs")
      }
      
      message(paste("\nFound", length(srx_info), "unique SRA accessions:", paste(srx_info, collapse=", ")))
      
      # Download and convert SRA files to FASTQ
      message("\nDownloading and converting SRA files to FASTQ...")
      for (srx in srx_info) {
        message(paste("\nProcessing", srx))
        
        # Get SRR accessions for this SRX
        srr_ids <- get_srr_from_srx(srx)
        if (is.null(srr_ids)) {
          message(paste("  Skipping", srx, "as no SRR accessions were found"))
          next
        }
        
        # Find the corresponding GSM for this SRX
        for (gsm in gsm_accessions) {
          gsm_dir <- file.path(samples_dir, gsm)
          sra_dir <- file.path(gsm_dir, "SRA")
          fastq_dir <- file.path(sra_dir, "FASTQ")
          
          # Create necessary directories
          dir.create(sra_dir, recursive = TRUE, showWarnings = FALSE)
          dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Download and convert each SRR
          for (srr in srr_ids) {
            message(paste("  Processing SRR:", srr))
            # Use prefetch and fasterq-dump for better performance
            cmd <- paste("module load sra-toolkit &&",
                        "prefetch", srr, "&&",           # don't use --type raw
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
      }
    } else if (is_microarray_data(metadata)) {
      message("\nMicroarray data detected. Processing array files...")
      
      # Save platform information to a separate file
      platform_info <- list(
        platform_id = metadata$platform_id[1],
        manufacturer = platform_data@header$manufacturer,
        technology = platform_data@header$technology
      )
      
      # Save platform information
      write.csv(
        data.frame(
          platform_id = platform_info$platform_id,
          manufacturer = platform_info$manufacturer,
          technology = platform_info$technology
        ),
        file = file.path(output_dir, "platform_info.csv")
      )
      
      # Process each GSM (just download files, no processing)
      for (gsm in gsm_accessions) {
        message(paste("\nProcessing", gsm))
        
        # Create GSM-specific directory
        gsm_dir <- file.path(samples_dir, gsm)
        dir.create(gsm_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Get GSM data and save RDS
        gsm_data <- getGEO(gsm)
        saveRDS(gsm_data, file = file.path(gsm_dir, paste0(gsm, ".rds")))
        
        # Download microarray data (without processing)
        download_microarray_data(gsm_data, gsm_dir)
      }
      
      message("\nDownload complete. Please install the following annotation package before processing:")
      message(paste("Platform ID:", platform_info$platform_id))
      message(paste("Manufacturer:", platform_info$manufacturer))
      message(paste("Technology:", platform_info$technology))
    } else {
      message("\nThis does not appear to be RNA-seq or microarray data. Available metadata:")
      print(metadata[1, ])
    }
    
    message("\nSuccessfully downloaded GEO data for ", accession)
  }, error = function(e) {
    stop(paste("Failed to download GEO data:", e$message))
  })
}

# Execute download
download_geo_data(accession, output_dir) 