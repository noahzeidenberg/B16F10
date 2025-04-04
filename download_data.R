# SRA Download and Conversion Pipeline
# This script downloads SRA files and converts them to FASTQ format

# Function to check if a command is available
check_command <- function(cmd) {
  tryCatch({
    system(sprintf("which %s", cmd), intern = TRUE, ignore.stderr = TRUE)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# Verify required commands are available
tryCatch({
  # Check for required commands
  if (!check_command("prefetch") || !check_command("fasterq-dump")) {
    stop("Required commands not found. Please ensure sra-toolkit is loaded in the SLURM script.")
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
  if (.Platform$OS.type == "windows") {
    df <- system(sprintf('wmic logicaldisk where "DeviceID=%s" get size', path), intern = TRUE)
    df <- strsplit(df, "\n")[[1]][2]
    return(as.numeric(gsub("[^0-9.]", "", df)))
  } else {
    df <- system(sprintf('df -B1G %s | awk \'NR==2 {print $4}\'', path), intern = TRUE)
    return(as.numeric(df))
  }
}

# Function to determine the appropriate base directory
get_base_dir <- function() {
  current_dir <- getwd()
  if (Sys.getenv("SLURM_TMPDIR") != "") {
    return(ifelse(Sys.getenv("USE_TMPDIR", "1") == "1", 
                 Sys.getenv("SLURM_TMPDIR"), 
                 Sys.getenv("SLURM_SUBMIT_DIR")))
  }
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
  
  cat(sprintf("Attempting to create sample directory: %s\n", sample_dir))
  
  # Check if parent directory exists and is writable
  parent_dir <- dirname(sample_dir)
  if (!dir.exists(parent_dir)) {
    cat(sprintf("Parent directory does not exist: %s\n", parent_dir))
    cat("Attempting to create parent directory...\n")
    if (!dir.create(parent_dir, showWarnings = TRUE, recursive = TRUE)) {
      cat(sprintf("Failed to create parent directory: %s\n", parent_dir))
      cat("Parent directory permissions:\n")
      system(sprintf("ls -ld %s", dirname(parent_dir)), intern = TRUE) |> cat(sep = "\n")
      stop(sprintf("Failed to create parent directory: %s", parent_dir))
    }
    cat("Successfully created parent directory\n")
  }
  
  # Check if parent directory is writable
  test_file <- file.path(parent_dir, ".test_write")
  if (!file.create(test_file)) {
    cat(sprintf("Parent directory is not writable: %s\n", parent_dir))
    cat("Parent directory permissions:\n")
    system(sprintf("ls -ld %s", parent_dir), intern = TRUE) |> cat(sep = "\n")
    stop(sprintf("Parent directory is not writable: %s", parent_dir))
  }
  file.remove(test_file)
  cat("Parent directory is writable\n")
  
  # Create sample directory with error checking
  if (!dir.exists(sample_dir)) {
    cat("Creating sample directory...\n")
    if (!dir.create(sample_dir, showWarnings = TRUE, recursive = TRUE)) {
      cat(sprintf("Failed to create sample directory: %s\n", sample_dir))
      cat("Current directory permissions:\n")
      system("pwd", intern = TRUE) |> cat(sep = "\n")
      system("ls -la", intern = TRUE) |> cat(sep = "\n")
      stop(sprintf("Failed to create sample directory: %s", sample_dir))
    }
    cat("Successfully created sample directory\n")
  } else {
    cat("Sample directory already exists\n")
  }
  
  # Create subdirectories
  dirs <- c(
    file.path(sample_dir, "SRA"),
    file.path(sample_dir, "SRA", "FASTQ")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      cat(sprintf("Creating directory: %s\n", dir))
      if (!dir.create(dir, showWarnings = TRUE, recursive = TRUE)) {
        cat(sprintf("Failed to create directory: %s\n", dir))
        cat("Parent directory permissions:\n")
        system(sprintf("ls -ld %s", dirname(dir)), intern = TRUE) |> cat(sep = "\n")
        stop(sprintf("Failed to create directory: %s", dir))
      }
      cat(sprintf("Successfully created directory: %s\n", dir))
    } else {
      cat(sprintf("Directory already exists: %s\n", dir))
    }
  }
  
  cat(sprintf("Created sample directory structure at: %s\n", sample_dir))
  return(sample_dir)
}

# Add rate limiting mechanism
last_request_time <- Sys.time()
request_count <- 0
MAX_REQUESTS_PER_MINUTE <- 50  # Increased from 10 to 50 requests per minute
MIN_REQUEST_INTERVAL <- 0.2  # Minimum time between requests in seconds

# Function to enforce rate limiting
enforce_rate_limit <- function() {
  current_time <- Sys.time()
  elapsed <- as.numeric(difftime(current_time, last_request_time, units = "secs"))
  
  # Always ensure minimum interval between requests
  if (elapsed < MIN_REQUEST_INTERVAL) {
    wait_time <- MIN_REQUEST_INTERVAL - elapsed
    cat(sprintf("Enforcing minimum request interval. Waiting %.1f seconds...\n", wait_time))
    Sys.sleep(wait_time)
    current_time <- Sys.time()
    elapsed <- MIN_REQUEST_INTERVAL
  }
  
  # Reset counter if we're in a new minute
  if (elapsed >= 60) {
    request_count <<- 0
    last_request_time <<- current_time
  } else if (request_count >= MAX_REQUESTS_PER_MINUTE) {
    # Wait until the next minute
    wait_time <- 60 - elapsed
    cat(sprintf("Rate limit reached (%d requests in %.1f seconds). Waiting %.1f seconds...\n", 
                request_count, elapsed, wait_time))
    Sys.sleep(wait_time)
    request_count <<- 0
    last_request_time <<- Sys.time()
  }
  
  request_count <<- request_count + 1
  cat(sprintf("Request count: %d/%d (%.1f seconds elapsed)\n", 
              request_count, MAX_REQUESTS_PER_MINUTE, elapsed))
}

# Function to make API request with rate limiting
make_api_request <- function(fn, ...) {
  max_retries <- 3
  retry_count <- 0
  
  while (retry_count < max_retries) {
    enforce_rate_limit()
    tryCatch({
      result <- fn(...)
      return(result)
    }, error = function(e) {
      retry_count <<- retry_count + 1
      
      if (grepl("HTTP 429|HTTP 403|Failed to perform HTTP request", e$message)) {
        # Rate limit hit, wait longer and retry
        wait_time <- 10 * retry_count  # Exponential backoff
        cat(sprintf("Rate limit hit (attempt %d of %d). Waiting %d seconds...\n", 
                   retry_count, max_retries, wait_time))
        Sys.sleep(wait_time)
      } else {
        cat(sprintf("Error in API request (attempt %d of %d): %s\n", 
                   retry_count, max_retries, e$message))
        if (retry_count >= max_retries) {
          stop(e)
        }
      }
    })
  }
  
  stop("Maximum retries reached for API request")
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
      
      # Get list of all entities
      all_entities <- GSMList(gse)
      
      # Filter out any non-GSM entries (GPL, GDS, etc.)
      gsm_list <- all_entities[grep("^GSM", names(all_entities))]
      total_gsm <- length(gsm_list)
      
      # Log only the GSM entities found
      cat(sprintf("Found %d GSM samples to process\n", total_gsm))
      
      if (total_gsm == 0) {
        cat("No GSM samples found. Skipping...\n")
        return(NULL)
      }
      
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
            # Extract SRX ID from the SRA link
            srx_id <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
            # Clean up the SRX ID
            srx_id <- gsub("\\[.*\\]", "", srx_id)  # Remove any brackets and their contents
            srx_id <- trimws(srx_id)  # Remove any whitespace
            cat(sprintf("Found SRX ID: %s\n", srx_id))
            return(srx_id)
          }
          cat("No SRX ID found for this GSM\n")
          return(NA)
        }, error = function(e) {
          cat(sprintf("Error processing GSM %s: %s\n", gsm_name, e$message))
          # Print additional debugging information
          cat("Current working directory:\n")
          system("pwd", intern = TRUE) |> cat(sep = "\n")
          cat("Directory contents:\n")
          system("ls -la", intern = TRUE) |> cat(sep = "\n")
          cat("Samples directory contents:\n")
          system(sprintf("ls -la %s/samples", gse_dir), intern = TRUE) |> cat(sep = "\n")
          return(NA)
        })
      })
      
      # Remove any NA values
      srx_ids <- srx_ids[!is.na(srx_ids)]
      
      if (length(srx_ids) == 0) {
        cat(sprintf("No SRX IDs found for GSE ID: %s. Skipping...\n", gse_id))
        return(NULL)
      }
      
      cat(sprintf("Successfully processed %d GSM samples and found %d SRX IDs\n", 
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
        
        # First, try to get the SRA ID directly using the SRX ID
        xml_result <- make_api_request(rentrez::entrez_search, 
                                     db = "sra", 
                                     term = paste0(srx_id, "[Accession]"),
                                     retmax = 1)
        
        if (length(xml_result$ids) > 0) {
          # If we found an ID, use it to fetch the details
          xml_result <- make_api_request(rentrez::entrez_fetch, 
                                       db = "sra", 
                                       id = xml_result$ids[1], 
                                       rettype = "xml")
        } else {
          # If direct search failed, try the old method
          xml_result <- make_api_request(rentrez::entrez_fetch, 
                                       db = "sra", 
                                       id = srx_id, 
                                       rettype = "xml")
        }
        
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

# Function to download SRA files using prefetch
download_sra_files <- function(gse_id, sra_ids) {
  if (!check_command("prefetch")) {
    stop("prefetch command not found. Please ensure SRA toolkit is installed and in PATH")
  }
  
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Check available space
  available_space <- check_disk_space(base_dir)
  if (!is.na(available_space) && available_space < 50) {
    stop(sprintf("Insufficient disk space. Only %.1fGB available. Need at least 50GB.", available_space))
  }
  
  # Download SRA files for each sample
  for (sra_id in sra_ids) {
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      gsm_id <- basename(gsm_dir)
      sra_dir <- file.path(gsm_dir, "SRA")
      sra_file <- file.path(sra_dir, sra_id, paste0(sra_id, ".sra"))
      
      # Check if this SRA ID belongs to this GSM
      gsm_rds <- file.path(gsm_dir, paste0(gsm_id, ".rds"))
      if (file.exists(gsm_rds)) {
        gsm_obj <- readRDS(gsm_rds)
        relations <- Meta(gsm_obj)$relation
        sra_link <- relations[grep("SRA:", relations)]
        if (length(sra_link) > 0) {
          srx_id <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
          xml_result <- make_api_request(rentrez::entrez_fetch, db = "sra", id = srx_id, rettype = "xml")
          xml_doc <- xml2::read_xml(xml_result)
          srr_ids <- xml2::xml_find_all(xml_doc, ".//RUN") |> xml2::xml_attr("accession")
          if (sra_id %in% srr_ids) {
            cat(sprintf("Downloading SRA %s for GSM %s\n", sra_id, gsm_id))
            dir.create(sra_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Download with prefetch
            cmd <- sprintf('prefetch --max-size 100G -O %s %s 2>&1', sra_dir, sra_id)
            result <- system(cmd, intern = TRUE)
            cat("prefetch output:\n", paste(result, collapse = "\n"), "\n")
            
            # Verify download
            if (file.exists(sra_file) && file.info(sra_file)$size > 0) {
              cat(sprintf("Successfully downloaded SRA file: %s (size: %.2f MB)\n", 
                         sra_file, file.info(sra_file)$size/1024/1024))
              break
            } else {
              warning(sprintf("Failed to download SRA file: %s", sra_id))
            }
          }
        }
      }
    }
  }
}

# Function to convert SRA to FASTQ
convert_sra_to_fastq <- function(gse_id, sra_ids, threads = 8) {
  if (!check_command("fasterq-dump")) {
    stop("fasterq-dump command not found. Please ensure SRA toolkit is installed and in PATH")
  }
  
  if (!check_command("pigz")) {
    stop("pigz command not found. Please ensure pigz is installed and in PATH")
  }
  
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Track which SRA IDs have been successfully converted
  converted_sra_ids <- c()
  
  for (sra_id in sra_ids) {
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_file <- file.path(gsm_dir, "SRA", sra_id, paste0(sra_id, ".sra"))
      fastq_dir <- file.path(gsm_dir, "SRA", "FASTQ")
      
      # Check if FASTQ files already exist for this SRA ID
      if (dir.exists(fastq_dir)) {
        existing_fastq_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq\\.gz$"), full.names = TRUE)
        if (length(existing_fastq_files) > 0) {
          cat(sprintf("FASTQ files already exist for SRA %s in %s. Skipping conversion.\n", 
                     sra_id, fastq_dir))
          converted_sra_ids <- c(converted_sra_ids, sra_id)
          break
        }
      }
      
      if (file.exists(sra_file)) {
        cat(sprintf("Converting SRA %s to FASTQ in %s\n", sra_id, fastq_dir))
        dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Convert with timeout and retry logic
        max_retries <- 3
        for (retry in 1:max_retries) {
          if (retry > 1) {
            cat(sprintf("Retry %d: Waiting %d seconds...\n", retry, 30 * retry))
            Sys.sleep(30 * retry)
          }
          
          # First use fasterq-dump to create uncompressed FASTQ files
          cmd <- sprintf('timeout 3600 fasterq-dump --split-files --progress --threads %d --outdir %s %s 2>&1',
                        threads, fastq_dir, sra_id)
          result <- system(cmd, intern = TRUE)
          cat("fasterq-dump output:\n", paste(result, collapse = "\n"), "\n")
          
          # Check for uncompressed FASTQ files
          fastq_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq$"), full.names = TRUE)
          if (length(fastq_files) > 0) {
            cat(sprintf("Successfully created FASTQ files, now compressing with pigz...\n"))
            
            # Compress each FASTQ file with pigz
            for (fastq_file in fastq_files) {
              cat(sprintf("Compressing %s with pigz...\n", basename(fastq_file)))
              compress_cmd <- sprintf('pigz -p %d %s', threads, fastq_file)
              compress_result <- system(compress_cmd, intern = TRUE)
              if (length(compress_result) > 0) {
                cat("pigz output:\n", paste(compress_result, collapse = "\n"), "\n")
              }
            }
            
            # Verify compressed files exist
            compressed_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq\\.gz$"), full.names = TRUE)
            if (length(compressed_files) > 0) {
              cat(sprintf("Successfully created compressed FASTQ files: %s\n", 
                         paste(basename(compressed_files), collapse = ", ")))
              
              # Clean up uncompressed files
              uncompressed_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq$"), full.names = TRUE)
              if (length(uncompressed_files) > 0) {
                cat("Cleaning up uncompressed FASTQ files...\n")
                for (file in uncompressed_files) {
                  if (file.remove(file)) {
                    cat(sprintf("Removed uncompressed file: %s\n", basename(file)))
                  } else {
                    warning(sprintf("Failed to remove uncompressed file: %s", basename(file)))
                  }
                }
              }
              
              converted_sra_ids <- c(converted_sra_ids, sra_id)
              break
            } else {
              warning(sprintf("Failed to compress FASTQ files for SRA %s", sra_id))
            }
          } else if (retry == max_retries) {
            warning(sprintf("Failed to convert SRA file after %d attempts: %s", max_retries, sra_id))
          }
        }
      }
    }
  }
  
  # Report on conversion status
  if (length(converted_sra_ids) == length(sra_ids)) {
    cat("All SRA files have been successfully converted to FASTQ.\n")
  } else {
    cat(sprintf("Converted %d out of %d SRA files to FASTQ.\n", 
                length(converted_sra_ids), length(sra_ids)))
    cat("The following SRA files were not converted:\n")
    for (sra_id in setdiff(sra_ids, converted_sra_ids)) {
      cat(sprintf("  - %s\n", sra_id))
    }
  }
}

# Function to check if FASTQ files exist for a GSE ID
check_fastq_files <- function(gse_id) {
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  if (!dir.exists(gse_dir)) {
    return(FALSE)
  }
  
  sample_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
  if (length(sample_dirs) == 0) {
    return(FALSE)
  }
  
  for (sample_dir in sample_dirs) {
    fastq_dir <- file.path(sample_dir, "SRA", "FASTQ")
    if (dir.exists(fastq_dir)) {
      fastq_files <- list.files(fastq_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
      if (length(fastq_files) > 0) {
        cat(sprintf("Found existing FASTQ files in %s\n", fastq_dir))
        return(TRUE)
      }
    }
  }
  
  return(FALSE)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript download_data.R GSE223515")
    }
  }
  
  tryCatch({
    base_dir <- get_base_dir()
    cat(sprintf("Using base directory: %s\n", base_dir))
    
    # Print environment information
    cat("=== Environment Information ===\n")
    cat("Current working directory:\n")
    system("pwd", intern = TRUE) |> cat(sep = "\n")
    cat("Directory contents:\n")
    system("ls -la", intern = TRUE) |> cat(sep = "\n")
    cat("User and permissions:\n")
    system("id", intern = TRUE) |> cat(sep = "\n")
    cat("=== End Environment Information ===\n")
    
    # Check if FASTQ files already exist
    if (check_fastq_files(gse_id)) {
      cat(sprintf("Some FASTQ files already exist for %s. Will check each SRA ID individually.\n", gse_id))
    }
    
    # Get SRX IDs
    cat(sprintf("Getting SRX IDs for %s...\n", gse_id))
    srx_ids <- get_srx_ids(gse_id)
    if (is.null(srx_ids)) {
      cat(sprintf("No valid GSM samples found for %s. This could be because:\n", gse_id))
      cat("  1. The dataset is not high-throughput sequencing data\n")
      cat("  2. The dataset contains only GPL (platform) entries\n")
      cat("  3. The dataset structure is not as expected\n")
      cat("Please check the dataset manually or try a different GSE ID.\n")
      return(invisible(NULL))
    }
    
    # Convert SRX to SRA IDs
    cat("Converting SRX to SRA IDs...\n")
    sra_ids <- convert_srx_to_sra(srx_ids)
    if (is.null(sra_ids)) {
      stop(sprintf("Failed to convert SRX to SRA IDs for %s", gse_id))
    }
    
    # Download and convert files
    cat("Downloading SRA files...\n")
    download_sra_files(gse_id, sra_ids)
    
    cat("Converting SRA to FASTQ...\n")
    convert_sra_to_fastq(gse_id, sra_ids)
    
    cat("Download and conversion complete!\n")
  }, error = function(e) {
    cat(sprintf("Error during download process: %s\n", e$message))
    cat("=== Error Details ===\n")
    cat("Current working directory:\n")
    system("pwd", intern = TRUE) |> cat(sep = "\n")
    cat("Directory contents:\n")
    system("ls -la", intern = TRUE) |> cat(sep = "\n")
    cat("User and permissions:\n")
    system("id", intern = TRUE) |> cat(sep = "\n")
    cat("=== End Error Details ===\n")
    stop(sprintf("Error during download process: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 