#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("GEOquery", "SRAdb")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Function to get SRR accessions from an SRX ID
get_srr_from_srx <- function(srx) {
  tryCatch({
    message(paste("  Getting SRR accessions for", srx))
    # Use efetch with runinfo format to get SRR accessions from SRX
    cmd <- paste("efetch -db sra -id", srx, "-format runinfo | cut -d',' -f1")
    message(paste("  Running command:", cmd))
    
    # Execute the command and capture output
    srr_ids <- system(cmd, intern = TRUE)
    
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
    # Get the GSM data
    gsm_data <- getGEO(gsm)
    
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

# Function to safely download GEO data with retries
safe_getGEO <- function(accession, max_attempts = 3, delay = 2) {
  for (attempt in 1:max_attempts) {
    tryCatch({
      message(paste("Attempt", attempt, "of", max_attempts, "to download", accession))
      
      # Add delay between attempts
      if (attempt > 1) {
        message(paste("Waiting", delay, "seconds before retry..."))
        Sys.sleep(delay)
      }
      
      # Set options for GEOquery
      options('download.file.method.GEOquery'='curl')
      options('GEOquery.inmemory.gpl'=FALSE)
      
      # Try to download with different methods
      tryCatch({
        # First try with GSEMatrix=TRUE and destdir specified
        message("  Attempting download with GSEMatrix=TRUE...")
        gset <- getGEO(accession, GSEMatrix = TRUE, destdir = "temp")
      }, error = function(e) {
        message("  First attempt failed, trying without GSEMatrix...")
        # If that fails, try without GSEMatrix
        gset <- getGEO(accession, GSEMatrix = FALSE, destdir = "temp")
      })
      
      # Handle case where getGEO returns a list
      if (is.list(gset)) {
        message("  getGEO returned a list, extracting first element")
        gset <- gset[[1]]
      }
      
      # Verify we got valid data
      if (is.null(gset) || !inherits(gset, "GEOData")) {
        stop("Failed to get valid GEO data object")
      }
      
      return(gset)
      
    }, error = function(e) {
      if (attempt == max_attempts) {
        # On final attempt, try to get supplementary files directly
        message("  All attempts failed, trying to get supplementary files directly...")
        tryCatch({
          # Try to get supplementary files with curl method
          options('download.file.method.GEOquery'='curl')
          sfiles <- getGEOSuppFiles(accession, baseDir = "temp")
          if (!is.null(sfiles) && nrow(sfiles) > 0) {
            message("  Successfully retrieved supplementary files")
            # Create a minimal GEOData object with the supplementary files
            gset <- new("GEOData")
            gset@header <- list(
              geo_accession = accession,
              supplementary_file = rownames(sfiles)
            )
            return(gset)
          }
        }, error = function(e2) {
          message(paste("  Supplementary files retrieval failed:", e2$message))
        })
        
        # Try one last time with direct FTP access
        message("  Attempting direct FTP access...")
        tryCatch({
          ftp_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", 
                          substr(accession, 1, 6), "nnn/", accession, "/")
          message(paste("  Checking FTP URL:", ftp_url))
          
          # Try to get the series matrix file directly
          matrix_url <- paste0(ftp_url, accession, "_series_matrix.txt.gz")
          message(paste("  Attempting to download:", matrix_url))
          
          # Create temp directory if it doesn't exist
          dir.create("temp", showWarnings = FALSE, recursive = TRUE)
          
          # Download the file
          download.file(matrix_url, 
                       destfile = file.path("temp", paste0(accession, "_series_matrix.txt.gz")),
                       mode = "wb")
          
          # Try to read the downloaded file
          gset <- getGEO(accession, GSEMatrix = TRUE, destdir = "temp")
          if (!is.null(gset)) {
            message("  Successfully retrieved data through direct FTP")
            return(gset)
          }
        }, error = function(e3) {
          message(paste("  Direct FTP access failed:", e3$message))
        })
        
        stop(paste("Failed to download after", max_attempts, "attempts:", e$message))
      } else {
        message(paste("Attempt", attempt, "failed:", e$message))
        return(NULL)
      }
    })
  }
}

# Function to process a single dataset
process_dataset <- function(accession) {
  tryCatch({
    message(paste("Processing dataset:", accession))
    
    # Create temp directory if it doesn't exist
    dir.create("temp", showWarnings = FALSE, recursive = TRUE)
    
    # Path to your GDS table
    gds_df <- read.csv("gds_table.csv")
    
    # Get the current row based on the accession number
    row <- which(gds_df$Accession == accession)
    
    if (length(row) == 0) {
      stop(paste("Accession", accession, "not found in GDS table"))
    }
    
    # Get the dataset metadata
    gsm_data <- safe_getGEO(accession)
    if (is.null(gsm_data)) {
      stop(paste("Failed to retrieve data for accession:", accession))
    }
    
    # Create output directory
    output_dir <- file.path("results", accession)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get the header information
    header_info <- gsm_data@header
    
    # Check if we have supplementary files
    if (!is.null(header_info$supplementary_file)) {
      message("  Found supplementary files")
      message(paste("  Files:", paste(header_info$supplementary_file, collapse="\n    ")))
      
      # Create supplementary files directory
      suppl_dir <- file.path(output_dir, "supplementary_files")
      dir.create(suppl_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Download supplementary files
      for (file in header_info$supplementary_file) {
        message(paste("  Downloading:", file))
        tryCatch({
          # Try to download with curl method
          options('download.file.method.GEOquery'='curl')
          download.file(file, destfile = file.path(suppl_dir, basename(file)), mode = "wb")
        }, error = function(e) {
          message(paste("    Failed to download:", file, "-", e$message))
        })
      }
    }
    
    # Check if it's sequencing data
    if (is_sequencing_data(header_info)) {
      message(paste("Dataset", accession, "is sequencing data"))
      # Get SRX accessions
      srx_ids <- get_srx_from_gsm(accession)
      if (!is.null(srx_ids)) {
        # Get SRR accessions
        srr_ids <- lapply(srx_ids, get_srr_from_srx)
        srr_ids <- unlist(srr_ids[!sapply(srr_ids, is.null)])
        
        if (length(srr_ids) > 0) {
          # Save SRR IDs to file
          write.csv(data.frame(SRR = srr_ids), 
                   file.path(output_dir, "srr_ids.csv"), 
                   row.names = FALSE)
          message(paste("Saved", length(srr_ids), "SRR IDs"))
        }
      }
    } else if (is_microarray_data(header_info)) {
      message(paste("Dataset", accession, "is microarray data"))
      # Download microarray data
      success <- download_microarray_data(gsm_data, output_dir)
      if (success) {
        message("Successfully downloaded microarray data")
      } else {
        message("Failed to download microarray data")
      }
    } else {
      message(paste("Dataset", accession, "is neither sequencing nor microarray data"))
    }
    
    # Save metadata
    write.csv(as.data.frame(header_info), 
              file.path(output_dir, "metadata.csv"), 
              row.names = TRUE)
    
    message(paste("Completed processing dataset:", accession))
    return(TRUE)
    
  }, error = function(e) {
    message(paste("Error processing dataset", accession, ":", e$message))
    return(FALSE)
  })
}

# Main script logic
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check if we have an accession number
  if (length(args) < 1) {
    stop("No accession number provided")
  }
  
  # Get the accession number (first argument)
  accession <- args[1]
  
  # Call the function to process the dataset
  success <- process_dataset(accession)
  
  # Exit with appropriate status code
  if (!success) {
    quit(status = 1)
  }
}