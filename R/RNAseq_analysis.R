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

# Function to process a single dataset
process_dataset <- function(accession) {
  tryCatch({
    # Create results directory for this dataset
    results_dir <- file.path("results", accession)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get GEO data
    message("\nDownloading GEO metadata...")
    gset <- getGEO(accession, GSEMatrix = TRUE)
    saveRDS(gset, file = file.path(results_dir, paste0(accession, "_geo_object.rds")))
    metadata <- pData(gset[[1]])
    write.csv(metadata, file = file.path(results_dir, paste0(accession, "_metadata.csv")))
    
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
    samples_dir <- file.path(results_dir, "samples")
    dir.create(samples_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get GSM accessions
    gsm_accessions <- metadata$geo_accession
    message(paste("\nFound", length(gsm_accessions), "GSM accessions:", paste(gsm_accessions, collapse=", ")))
    
    # Check data type and process accordingly
    if (is_sequencing_data(metadata)) {
      message("\nRNA-seq data detected. Fetching SRA information...")
      
      # Get SRX accessions for each GSM
      srx_info <- c()
      for (gsm in gsm_accessions) {
        message(paste("\nProcessing", gsm))
        
        # Create GSM-specific directory
        gsm_dir <- file.path(samples_dir, gsm)
        dir.create(gsm_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Get GSM data and save RDS
        gsm_data <- getGEO(gsm)
        saveRDS(gsm_data, file = file.path(gsm_dir, paste0(gsm, ".rds")))
        
        # Get SRX accession
        srx_id <- get_srx_from_gsm(gsm)
        if (!is.null(srx_id)) {
          srx_info <- c(srx_info, srx_id)
        }
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
            
            # Check if SRA toolkit commands are available
            if (system("which prefetch", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
              message("  Error: prefetch command not found. Please ensure SRA toolkit is installed and in PATH")
              next
            }
            
            if (system("which fasterq-dump", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
              message("  Error: fasterq-dump command not found. Please ensure SRA toolkit is installed and in PATH")
              next
            }
            
            # Use prefetch and fasterq-dump for better performance
            cmd <- paste("prefetch", srr, "&&",           # don't use --type raw
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
              
              # Check if fastq-dump is available
              if (system("which fastq-dump", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
                message("  Error: fastq-dump command not found. Please ensure SRA toolkit is installed and in PATH")
                next
              }
              
              # Try alternative method using fastq-dump
              alt_cmd <- paste("fastq-dump --split-files --gzip",
                              "--outdir", fastq_dir, srr)
              message(paste("  Running alternative command:", alt_cmd))
              alt_result <- system(alt_cmd)
              if (alt_result != 0) {
                message(paste("  Error: Both download methods failed for", srr))
                message("  Please ensure SRA toolkit is properly installed and configured")
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
      message("\nMicroarray data detected. Downloading array files...")
      
      # Process each GSM
      for (gsm in gsm_accessions) {
        message(paste("\nProcessing", gsm))
        
        # Create GSM-specific directory
        gsm_dir <- file.path(samples_dir, gsm)
        dir.create(gsm_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Get GSM data and save RDS
        gsm_data <- getGEO(gsm)
        saveRDS(gsm_data, file = file.path(gsm_dir, paste0(gsm, ".rds")))
        
        # Download microarray data
        if (!download_microarray_data(gsm_data, gsm_dir)) {
          message(paste("  Warning: Failed to download microarray data for", gsm))
        }
      }
    } else {
      message("\nThis does not appear to be RNA-seq or microarray data. Available metadata:")
      print(metadata[1, ])
    }
    
    message("\nSuccessfully processed ", accession)
  }, error = function(e) {
    message(paste("Error processing", accession, ":", e$message))
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
for (accession in gds_df$Accession) {
  tryCatch({
    message(paste("\nProcessing dataset:", accession))
    results[[accession]] <- process_dataset(accession)
  }, error = function(e) {
    message(paste("Error processing", accession, ":", e$message))
  })
}

# Save session info
sink("results/session_info.txt")
sessionInfo()
sink() 