# SRA Download and Conversion Pipeline
# This script downloads SRA files and converts them to FASTQ format

# Verify required commands are available
tryCatch({
  # Check for required commands
  if (!check_command("prefetch") || !check_command("fasterq-dump")) {
    stop("Required commands not found. Please ensure sra-toolkit is loaded in the SLURM script.")
  }
  if (!check_command("pigz")) {
    stop("pigz command not found. Please ensure pigz is loaded in the SLURM script.")
  }
}, error = function(e) {
  cat("Error verifying commands:", conditionMessage(e), "\n")
  cat("Please ensure the required modules are loaded in the SLURM script.\n")
  quit(status = 1)
})

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "GEOquery",
  "rentrez",
  "xml2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Read and validate API keys from .env file
tryCatch({
  if (!file.exists(".env")) {
    stop("Cannot find .env file. Please ensure it exists in the current directory.")
  }
  
  env_lines <- readLines(".env")
  api_keys <- env_lines[grep("^NCBI_API_KEY_", env_lines)]
  
  if (length(api_keys) == 0) {
    stop("No NCBI API keys found in .env file. Keys should start with 'NCBI_API_KEY_'")
  }
  
  # Clean and extract keys
  api_keys <- gsub("'", "", api_keys)
  api_keys <- sapply(api_keys, function(x) {
    key <- strsplit(x, "=")[[1]]
    if (length(key) != 2) {
      stop("Invalid API key format in .env file. Expected format: NCBI_API_KEY_X='your-key-here'")
    }
    return(trimws(key[2]))
  })
  
  # Validate key format (basic check)
  invalid_keys <- api_keys == "" | is.na(api_keys) | nchar(api_keys) < 10
  if (any(invalid_keys)) {
    stop("Found invalid or empty API keys in positions: ", 
         paste(which(invalid_keys), collapse = ", "))
  }
  
  # Get array job ID and select key
  array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
  selected_key <- api_keys[((array_id - 1) %% length(api_keys)) + 1]
  
  # Set the selected API key for both rentrez and GEOquery
  rentrez::set_entrez_key(selected_key)
  options(GEOquery.api.key = selected_key)
  
  cat(sprintf("Successfully configured API key %d of %d for array job %d\n", 
              ((array_id - 1) %% length(api_keys)) + 1,
              length(api_keys),
              array_id))
  
}, error = function(e) {
  cat("Error setting up API keys:", conditionMessage(e), "\n")
  cat("Please check your .env file configuration.\n")
  quit(status = 1)
})

# Function to check available disk space (in GB)
check_disk_space <- function(path) {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      # Windows specific command
      df <- system(sprintf('wmic logicaldisk where "DeviceID=%s" get size', path))
      df <- strsplit(df, "\n")[[1]][2]
      df <- as.numeric(gsub("[^0-9.]", "", df))
      return(df)
    } else {
      # Linux specific command - use df -B1G to get size in GB
      df <- system(sprintf('df -B1G %s | awk \'NR==2 {print $4}\'', path), intern = TRUE)
      df <- as.numeric(df)
      return(df)
    }
  }, error = function(e) {
    cat("Error checking disk space:", conditionMessage(e), "\n")
    return(NA)
  })
}

# Function to determine the appropriate base directory
get_base_dir <- function() {
  # Get current directory
  current_dir <- getwd()
  cat(sprintf("Current working directory: %s\n", current_dir))
  
  # Check if we're running in SLURM
  if (Sys.getenv("SLURM_TMPDIR") != "") {
    cat("Running in SLURM environment\n")
    return(Sys.getenv("SLURM_TMPDIR"))
  }
  
  # When running manually, always use current directory
  cat("Running manually, using current directory\n")
  return(current_dir)
}

# Function to create directory structure for a GSE ID
create_gse_structure <- function(gse_id) {
  # Create main GSE directory in the appropriate directory
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  cat(sprintf("Checking GSE directory: %s\n", gse_dir))
  
  # If directory exists, check if it's writable
  if (dir.exists(gse_dir)) {
    cat("GSE directory already exists, checking permissions...\n")
    # Try to create a test file to check write permissions
    test_file <- file.path(gse_dir, ".test_write")
    if (file.create(test_file)) {
      file.remove(test_file)
      cat("Directory exists and is writable\n")
    } else {
      stop(sprintf("GSE directory exists but is not writable: %s", gse_dir))
    }
  } else {
    # Create directory with error checking
    if (!dir.create(gse_dir, showWarnings = TRUE, recursive = TRUE)) {
      cat(sprintf("Directory creation failed. Checking permissions...\n"))
      cat(sprintf("Parent directory exists: %s\n", dir.exists(dirname(gse_dir))))
      cat(sprintf("Parent directory permissions: %s\n", 
                  system(sprintf("ls -ld %s", dirname(gse_dir)), intern = TRUE)))
      stop(sprintf("Failed to create GSE directory: %s", gse_dir))
    }
    cat("Created new GSE directory\n")
  }
  
  # Create subdirectories
  dirs <- c(
    file.path(gse_dir, "samples"),
    file.path(gse_dir, "logs")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      if (!dir.create(dir, showWarnings = TRUE, recursive = TRUE)) {
        stop(sprintf("Failed to create directory: %s", dir))
      }
      cat(sprintf("Created directory: %s\n", dir))
    } else {
      cat(sprintf("Directory already exists: %s\n", dir))
    }
  }
  
  cat(sprintf("Directory structure ready at: %s\n", gse_dir))
  return(gse_dir)
}

# Function to create sample directory structure
create_sample_structure <- function(gse_dir, gsm_id) {
  # Create sample directory
  sample_dir <- file.path(gse_dir, "samples", gsm_id)
  
  # Create directory with error checking
  if (!dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)) {
    stop(sprintf("Failed to create sample directory: %s", sample_dir))
  }
  
  # Create subdirectories
  dirs <- c(
    file.path(sample_dir, "SRA"),
    file.path(sample_dir, "SRA", "FASTQ")
  )
  
  for (dir in dirs) {
    if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) {
      stop(sprintf("Failed to create directory: %s", dir))
    }
  }
  
  cat(sprintf("Created sample directory structure at: %s\n", sample_dir))
  return(sample_dir)
}

# Add rate limiting mechanism
last_request_time <- Sys.time()
request_count <- 0
MAX_REQUESTS_PER_SECOND <- 10

# Function to enforce rate limiting
enforce_rate_limit <- function() {
  current_time <- Sys.time()
  elapsed <- as.numeric(difftime(current_time, last_request_time, units = "secs"))
  
  if (elapsed < 1) {
    if (request_count >= MAX_REQUESTS_PER_SECOND) {
      # Wait until the next second
      Sys.sleep(1 - elapsed)
      request_count <<- 0
      last_request_time <<- Sys.time()
    }
  } else {
    request_count <<- 0
    last_request_time <<- current_time
  }
  
  request_count <<- request_count + 1
}

# Function to make API request with rate limiting
make_api_request <- function(fn, ...) {
  enforce_rate_limit()
  tryCatch({
    result <- fn(...)
    return(result)
  }, error = function(e) {
    if (grepl("HTTP 429|HTTP 403|Failed to perform HTTP request", e$message)) {
      # Rate limit hit, wait and retry
      Sys.sleep(1)
      enforce_rate_limit()
      return(fn(...))
    }
    stop(e)
  })
}

# Function to get SRX IDs from GSM IDs and save GSM objects
get_srx_ids <- function(gse_id) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Get the GEO object without matrix
  cat(sprintf("Fetching GSM information for %s...\n", gse_id))
  
  # Add retry logic for getGEO
  max_retries <- 3
  retry_count <- 0
  success <- FALSE
  
  while (!success && retry_count < max_retries) {
    tryCatch({
      # Add a delay before each attempt to avoid rate limiting
      if (retry_count > 0) {
        delay_time <- 5 * retry_count  # Increase delay with each retry
        cat(sprintf("Retry %d: Waiting %d seconds before trying again...\n", 
                   retry_count, delay_time))
        Sys.sleep(delay_time)
      }
      
      # Try to get the GEO object with rate limiting
      gse <- make_api_request(getGEO, gse_id, GSEMatrix = FALSE)
      success <- TRUE
      
      # Check if this is HTS data
      data_type <- Meta(gse)$type
      if (is.null(data_type) || !grepl("high throughput sequencing", data_type, ignore.case = TRUE)) {
        cat(sprintf("Dataset %s is not High-Throughput Sequencing data. Found type: %s. Skipping...\n", 
                   gse_id, ifelse(is.null(data_type), "unknown", data_type)))
        return(NULL)
      }
      
      # Get list of GSM objects
      gsm_list <- GSMList(gse)
      total_gsm <- length(gsm_list)
      cat(sprintf("Found %d GSM objects to process\n", total_gsm))
      
      # Save each GSM object and extract SRX IDs
      srx_ids <- sapply(seq_along(names(gsm_list)), function(i) {
        gsm_name <- names(gsm_list)[i]
        gsm <- gsm_list[[gsm_name]]
        
        cat(sprintf("Processing GSM %d of %d: %s\n", i, total_gsm, gsm_name))
        
        tryCatch({
          # Create sample directory structure
          sample_dir <- create_sample_structure(gse_dir, gsm_name)
          
          # Save GSM object
          rds_file <- file.path(sample_dir, paste0(gsm_name, ".rds"))
          saveRDS(gsm, rds_file)
          cat(sprintf("Saved GSM object to: %s\n", rds_file))
          
          # Extract SRX ID
          relations <- Meta(gsm)$relation
          sra_link <- relations[grep("SRA:", relations)]
          if (length(sra_link) > 0) {
            srx_id <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
            cat(sprintf("Found SRX ID: %s\n", srx_id))
            return(srx_id)
          }
          cat("No SRX ID found for this GSM\n")
          return(NA)
        }, error = function(e) {
          cat(sprintf("Error processing GSM %s: %s\n", gsm_name, e$message))
          return(NA)
        })
      })
      
      # Remove any NA values
      srx_ids <- srx_ids[!is.na(srx_ids)]
      
      if (length(srx_ids) == 0) {
        cat(sprintf("No SRX IDs found for GSE ID: %s. Skipping...\n", gse_id))
        return(NULL)
      }
      
      cat(sprintf("Successfully processed %d GSM objects and found %d SRX IDs\n", 
                  total_gsm, length(srx_ids)))
      return(srx_ids)
      
    }, error = function(e) {
      retry_count <- retry_count + 1
      
      if (grepl("HTTP 429|HTTP 403|Failed to perform HTTP request", e$message)) {
        cat(sprintf("Rate limit hit (attempt %d of %d). Waiting before retry...\n", 
                   retry_count, max_retries))
        Sys.sleep(10 * retry_count)  # Exponential backoff
      } else {
        cat(sprintf("Error fetching GEO data (attempt %d of %d): %s\n", 
                   retry_count, max_retries, e$message))
      }
      
      if (retry_count >= max_retries) {
        cat("Maximum retries reached. Skipping this GSE ID.\n")
        return(NULL)
      }
    })
  }
  
  return(NULL)
}

# Function to convert SRX IDs to SRA IDs using rentrez
convert_srx_to_sra <- function(srx_ids) {
  sra_ids <- list()
  
  for (i in seq_along(srx_ids)) {
    srx_id <- srx_ids[i]
    cat(sprintf("Processing SRX ID %d of %d: %s\n", i, length(srx_ids), srx_id))
    
    max_retries <- 3
    retry_count <- 0
    success <- FALSE
    
    while (!success && retry_count < max_retries) {
      tryCatch({
        if (retry_count > 0) {
          delay_time <- 5 * retry_count
          cat(sprintf("Retry %d: Waiting %d seconds before trying again...\n", 
                     retry_count, delay_time))
          Sys.sleep(delay_time)
        }
        
        # Use rate-limited API request
        xml_result <- make_api_request(rentrez::entrez_fetch, db = "sra", id = srx_id, rettype = "xml")
        xml_doc <- xml2::read_xml(xml_result)
        srr_ids <- xml2::xml_find_all(xml_doc, ".//RUN") |> xml2::xml_attr("accession")
        
        valid_srr_ids <- srr_ids[grepl("^SRR", srr_ids)]
        
        if (length(valid_srr_ids) > 0) {
          sra_ids[[srx_id]] <- valid_srr_ids
          cat(sprintf("Found %d SRA ID(s): %s\n", 
                     length(valid_srr_ids), 
                     paste(valid_srr_ids, collapse = ", ")))
          success <- TRUE
        } else {
          cat(sprintf("No valid SRA IDs found for SRX ID: %s\n", srx_id))
          success <- TRUE
        }
      }, error = function(e) {
        retry_count <- retry_count + 1
        
        if (grepl("HTTP 429|HTTP 403|Failed to perform HTTP request", e$message)) {
          cat(sprintf("Rate limit hit (attempt %d of %d). Waiting before retry...\n", 
                     retry_count, max_retries))
          Sys.sleep(10 * retry_count)
        } else {
          cat(sprintf("Error processing SRX ID %s (attempt %d of %d): %s\n", 
                     srx_id, retry_count, max_retries, e$message))
        }
        
        if (retry_count >= max_retries) {
          cat("Maximum retries reached. Skipping this SRX ID.\n")
          success <- TRUE
        }
      })
    }
    
    # Add a small delay between SRX IDs
    Sys.sleep(1)
  }
  
  all_sra_ids <- unlist(sra_ids)
  
  if (length(all_sra_ids) == 0) {
    cat("No SRA IDs found for any of the SRX IDs\n")
    return(NULL)
  }
  
  return(all_sra_ids)
}

# Function to check if a command is available
check_command <- function(cmd) {
  tryCatch({
    system(sprintf("which %s", cmd), intern = TRUE, ignore.stderr = TRUE)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# Function to download SRA files using prefetch
download_sra_files <- function(gse_id, sra_ids) {
  # Check if prefetch is available
  if (!check_command("prefetch")) {
    stop("prefetch command not found. Please ensure SRA toolkit is installed and in PATH")
  }
  
  # Use the existing GSE directory structure
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Check available space in the directory
  available_space <- check_disk_space(base_dir)
  if (!is.na(available_space) && available_space < 50) {  # Require at least 50GB free
    stop(sprintf("Insufficient disk space in directory. Only %.1fGB available. Need at least 50GB.", 
                 available_space))
  }
  
  # Create a mapping of SRA IDs to GSM directories
  gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
  gsm_ids <- basename(gsm_dirs)
  
  # Download SRA files for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory by checking which GSM has this SRA ID
    for (i in seq_along(gsm_dirs)) {
      gsm_dir <- gsm_dirs[i]
      gsm_id <- gsm_ids[i]
      sra_dir <- file.path(gsm_dir, "SRA")
      
      # Check space again before each download
      available_space <- check_disk_space(base_dir)
      if (!is.na(available_space) && available_space < 20) {  # Require at least 20GB free for each file
        stop(sprintf("Insufficient disk space for downloading %s. Only %.1fGB available.", 
                     sra_id, available_space))
      }
      
      # Check if this SRA ID belongs to this GSM
      gsm_rds <- file.path(gsm_dir, paste0(gsm_id, ".rds"))
      if (file.exists(gsm_rds)) {
        gsm_obj <- readRDS(gsm_rds)
        relations <- Meta(gsm_obj)$relation
        sra_link <- relations[grep("SRA:", relations)]
        if (length(sra_link) > 0) {
          srx_id <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
          # Check if this SRX ID maps to the current SRA ID
          xml_result <- make_api_request(rentrez::entrez_fetch, db = "sra", id = srx_id, rettype = "xml")
          xml_doc <- xml2::read_xml(xml_result)
          srr_ids <- xml2::xml_find_all(xml_doc, ".//RUN") |> xml2::xml_attr("accession")
          if (sra_id %in% srr_ids) {
            # This is the correct GSM for this SRA ID
            cat(sprintf("Downloading SRA %s for GSM %s\n", sra_id, gsm_id))
            
            # Create SRA directory if it doesn't exist
            if (!dir.exists(sra_dir)) {
              dir.create(sra_dir, recursive = TRUE)
            }
            
            # Download with prefetch
            cmd <- sprintf('prefetch --max-size 100G -O %s %s', sra_dir, sra_id)
            result <- tryCatch({
              system(cmd, intern = TRUE)
            }, error = function(e) {
              warning(sprintf("Error downloading %s: %s", sra_id, e$message))
              return(NULL)
            })
            
            # Check if download was successful
            if (!is.null(result)) {
              if (length(grep("error|failed", result, ignore.case = TRUE)) > 0) {
                warning(sprintf("Error downloading %s: %s", sra_id, paste(result, collapse="\n")))
              } else {
                # Verify the SRA file exists
                sra_file <- file.path(sra_dir, paste0(sra_id, ".sra"))
                if (!file.exists(sra_file)) {
                  warning(sprintf("SRA file not found after download: %s", sra_file))
                } else {
                  cat(sprintf("Successfully downloaded SRA file: %s\n", sra_file))
                }
              }
            }
            
            # Break the inner loop once we've found the correct GSM
            break
          }
        }
      }
    }
  }
}

# Function to convert SRA to FASTQ
convert_sra_to_fastq <- function(gse_id, sra_ids, threads = 8) {
  # Check if fasterq-dump is available
  if (!check_command("fasterq-dump")) {
    stop("fasterq-dump command not found. Please ensure SRA toolkit is installed and in PATH")
  }
  
  # Use the existing GSE directory structure
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Convert SRA to FASTQ for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_file <- file.path(gsm_dir, "SRA", paste0(sra_id, ".sra"))
      fastq_dir <- file.path(gsm_dir, "SRA", "FASTQ")
      
      if (file.exists(sra_file)) {
        cat(sprintf("Converting SRA %s to FASTQ in %s\n", sra_id, fastq_dir))
        
        # Create FASTQ directory if it doesn't exist
        if (!dir.exists(fastq_dir)) {
          dir.create(fastq_dir, recursive = TRUE)
        }
        
        # Convert to FASTQ
        cmd <- sprintf('fasterq-dump --split-files --progress --threads %d --outdir %s %s',
                      threads, fastq_dir, sra_id)
        result <- tryCatch({
          system(cmd, intern = TRUE)
        }, error = function(e) {
          warning(sprintf("Error converting %s: %s", sra_id, e$message))
          return(NULL)
        })
        
        # Check if conversion was successful
        if (!is.null(result)) {
          if (length(grep("error|failed", result, ignore.case = TRUE)) > 0) {
            warning(sprintf("Error converting %s: %s", sra_id, paste(result, collapse="\n")))
          } else {
            # Verify FASTQ files were created
            fastq_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq$"), full.names = TRUE)
            if (length(fastq_files) == 0) {
              warning(sprintf("No FASTQ files found after conversion for %s", sra_id))
            } else {
              cat(sprintf("Successfully created FASTQ files: %s\n", 
                         paste(basename(fastq_files), collapse = ", ")))
              
              # Compress the FASTQ files
              for (fastq_file in fastq_files) {
                cat(sprintf("Compressing %s\n", fastq_file))
                cmd <- sprintf("pigz -p %d %s", threads, fastq_file)
                system(cmd)
              }
            }
          }
        }
      }
    }
  }
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    # Try to get GSE ID from command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript download_data.R GSE223515")
    }
  }
  
  tryCatch({
    # Get base directory
    base_dir <- get_base_dir()
    cat(sprintf("Using base directory: %s\n", base_dir))
    
    # Get SRX IDs from GSM IDs
    cat(sprintf("Getting SRX IDs for %s...\n", gse_id))
    srx_ids <- get_srx_ids(gse_id)
    if (is.null(srx_ids)) {
      cat(sprintf("Failed to get SRX IDs for %s. This may be due to rate limiting or API issues.\n", gse_id))
      cat("Please try again later or check if the GSE ID is valid and accessible.\n")
      cat(sprintf("You can verify the GSE ID at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n", gse_id))
      return(NULL)
    }
    
    # Convert SRX to SRA IDs
    cat("Converting SRX to SRA IDs...\n")
    sra_ids <- convert_srx_to_sra(srx_ids)
    if (is.null(sra_ids)) {
      cat(sprintf("Failed to convert SRX to SRA IDs for %s. This may be due to rate limiting or API issues.\n", gse_id))
      cat("Please try again later.\n")
      return(NULL)
    }
    
    # Download SRA files
    cat("Downloading SRA files...\n")
    download_sra_files(gse_id, sra_ids)
    
    # Convert SRA to FASTQ
    cat("Converting SRA to FASTQ...\n")
    convert_sra_to_fastq(gse_id, sra_ids)
    
    # Verify files were created
    cat("Verifying files were created...\n")
    gse_dir <- file.path(base_dir, gse_id)
    if (dir.exists(gse_dir)) {
      cat(sprintf("GSE directory exists: %s\n", gse_dir))
      cat("Contents of GSE directory:\n")
      system(sprintf("ls -la %s", gse_dir))
      
      samples_dir <- file.path(gse_dir, "samples")
      if (dir.exists(samples_dir)) {
        cat(sprintf("Samples directory exists: %s\n", samples_dir))
        cat("Contents of samples directory:\n")
        system(sprintf("ls -la %s", samples_dir))
        
        # Check each sample directory
        sample_dirs <- list.dirs(samples_dir, recursive = FALSE)
        for (sample_dir in sample_dirs) {
          cat(sprintf("Contents of sample directory %s:\n", sample_dir))
          system(sprintf("ls -la %s", sample_dir))
          
          # Check SRA directory
          sra_dir <- file.path(sample_dir, "SRA")
          if (dir.exists(sra_dir)) {
            cat(sprintf("Contents of SRA directory %s:\n", sra_dir))
            system(sprintf("ls -la %s", sra_dir))
            
            # Check FASTQ directory
            fastq_dir <- file.path(sra_dir, "FASTQ")
            if (dir.exists(fastq_dir)) {
              cat(sprintf("Contents of FASTQ directory %s:\n", fastq_dir))
              system(sprintf("ls -la %s", fastq_dir))
            } else {
              cat(sprintf("FASTQ directory does not exist: %s\n", fastq_dir))
            }
          } else {
            cat(sprintf("SRA directory does not exist: %s\n", sra_dir))
          }
        }
      } else {
        cat(sprintf("Samples directory does not exist: %s\n", samples_dir))
      }
    } else {
      cat(sprintf("GSE directory does not exist: %s\n", gse_dir))
    }
    
    cat("Download and conversion complete!\n")
  }, error = function(e) {
    cat(sprintf("Error during download process: %s\n", e$message))
    cat("Please check if the GSE ID is valid and accessible.\n")
    cat(sprintf("You can verify the GSE ID at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n", gse_id))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 