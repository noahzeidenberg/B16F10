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
get_srx_ids <- function(gse_id, gse_dir, api_key = NULL) {
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
      # Set getGPL=FALSE to prevent downloading and including GPL information
      gse <- make_api_request(getGEO, gse_id, GSEMatrix = FALSE, getGPL = FALSE)
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
      return(list(
        gsm_ids = names(gsm_list),
        srx_ids = srx_ids
      ))
      
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
convert_srx_to_sra <- function(srx_ids, api_key = NULL) {
  sra_ids <- list()
  srx_to_sra_mapping <- list()  # New: Create a mapping between SRX and SRA IDs
  
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
          # Store the mapping between SRX and SRA IDs
          srx_to_sra_mapping[[srx_id]] <- valid_srr_ids
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
  
  # Return both the list of all SRA IDs and the mapping
  return(list(
    sra_ids = all_sra_ids,
    srx_to_sra_mapping = srx_to_sra_mapping
  ))
}

# Function to download SRA files
download_sra_files <- function(sra_data, gse_dir, api_key = NULL) {
  if (!check_command("prefetch")) {
    stop("prefetch command not found. Please ensure SRA toolkit is installed and in PATH")
  }
  
  # Extract the SRA IDs and mapping
  sra_ids <- sra_data$sra_ids
  srx_to_sra_mapping <- sra_data$srx_to_sra_mapping
  
  # Debug: Check if GSE directory exists
  if (!dir.exists(gse_dir)) {
    stop(sprintf("GSE directory does not exist: %s", gse_dir))
  }
  
  # Debug: List contents of GSE directory
  cat(sprintf("Contents of GSE directory %s:\n", gse_dir))
  gse_contents <- list.dirs(gse_dir, recursive = FALSE, full.names = TRUE)
  for (item in gse_contents) {
    cat(sprintf("  - %s\n", item))
  }
  
  # Find all GSM directories
  gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
  cat(sprintf("Found %d GSM directories\n", length(gsm_dirs)))
  
  # Debug: List all GSM directories
  cat("GSM directories:\n")
  for (dir in gsm_dirs) {
    cat(sprintf("  - %s\n", dir))
  }
  
  # Create a mapping from GSM IDs to their directories
  gsm_to_dir_mapping <- list()
  for (gsm_dir in gsm_dirs) {
    gsm_id <- basename(gsm_dir)
    gsm_to_dir_mapping[[gsm_id]] <- gsm_dir
  }
  
  # Create a mapping from SRX IDs to GSM IDs by reading the GSM objects
  srx_to_gsm_mapping <- list()
  for (gsm_dir in gsm_dirs) {
    gsm_id <- basename(gsm_dir)
    rds_file <- file.path(gsm_dir, paste0(gsm_id, ".rds"))
    
    if (file.exists(rds_file)) {
      gsm <- readRDS(rds_file)
      relations <- Meta(gsm)$relation
      sra_link <- relations[grep("SRA:", relations)]
      
      if (length(sra_link) > 0) {
        # Extract SRX ID from the SRA link
        srx_id <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
        # Clean up the SRX ID
        srx_id <- gsub("\\[.*\\]", "", srx_id)  # Remove any brackets and their contents
        srx_id <- trimws(srx_id)  # Remove any whitespace
        
        # Store the mapping
        srx_to_gsm_mapping[[srx_id]] <- gsm_id
        cat(sprintf("Mapped SRX %s to GSM %s\n", srx_id, gsm_id))
      }
    }
  }
  
  # Download SRA files for each sample
  for (sra_id in sra_ids) {
    cat(sprintf("Downloading SRA %s\n", sra_id))
    
    # Find which SRX IDs are associated with this SRA ID
    associated_srx_ids <- names(srx_to_sra_mapping)[sapply(srx_to_sra_mapping, function(x) sra_id %in% x)]
    
    if (length(associated_srx_ids) == 0) {
      cat(sprintf("Warning: No SRX IDs found for SRA ID %s\n", sra_id))
      next
    }
    
    # Find which GSM IDs are associated with these SRX IDs
    associated_gsm_ids <- sapply(associated_srx_ids, function(srx_id) srx_to_gsm_mapping[[srx_id]])
    associated_gsm_ids <- associated_gsm_ids[!is.na(associated_gsm_ids)]
    
    if (length(associated_gsm_ids) == 0) {
      cat(sprintf("Warning: No GSM IDs found for SRA ID %s\n", sra_id))
      next
    }
    
    cat(sprintf("SRA %s is associated with GSM IDs: %s\n", 
                sra_id, paste(associated_gsm_ids, collapse = ", ")))
    
    # Download the SRA file to each associated GSM directory
    for (gsm_id in associated_gsm_ids) {
      gsm_dir <- gsm_to_dir_mapping[[gsm_id]]
      sra_dir <- file.path(gsm_dir, "SRA")
      
      # Debug: Check if SRA directory exists
      if (!dir.exists(sra_dir)) {
        cat(sprintf("SRA directory does not exist: %s\n", sra_dir))
        next
      }
      
      # Create SRA ID directory
      sra_id_dir <- file.path(sra_dir, sra_id)
      if (!dir.exists(sra_id_dir)) {
        cat(sprintf("Creating SRA ID directory: %s\n", sra_id_dir))
        dir.create(sra_id_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Check if SRA file already exists
      sra_file <- file.path(sra_id_dir, paste0(sra_id, ".sra"))
      if (file.exists(sra_file) && file.info(sra_file)$size > 0) {
        cat(sprintf("SRA file already exists: %s (size: %.2f MB)\n", 
                   sra_file, file.info(sra_file)$size/1024/1024))
        next
      }
      
      # Download SRA file
      cat(sprintf("Downloading SRA %s for GSM %s\n", sra_id, gsm_id))
      
      # Download with prefetch with retry logic
      max_retries <- 3
      retry_count <- 0
      success <- FALSE
      
      while (!success && retry_count < max_retries) {
        if (retry_count > 0) {
          delay_time <- 30 * retry_count
          cat(sprintf("Retry %d: Waiting %d seconds before trying again...\n", 
                     retry_count, delay_time))
          Sys.sleep(delay_time)
        }
        
        cmd <- sprintf('prefetch --max-size 100G --progress --verify yes --resume yes --heartbeat 1 --force no -O %s %s 2>&1', sra_id_dir, sra_id)
        result <- system(cmd, intern = TRUE)
        cat("prefetch output:\n", paste(result, collapse = "\n"), "\n")
        
        # Check for common error patterns in output
        if (any(grepl("error|failed|timeout|connection refused", tolower(result)))) {
          cat("Error detected in prefetch output, will retry...\n")
          retry_count <- retry_count + 1
          next
        }
        
        # Verify download
        if (file.exists(sra_file) && file.info(sra_file)$size > 0) {
          cat(sprintf("Successfully downloaded SRA file: %s (size: %.2f MB)\n", 
                     sra_file, file.info(sra_file)$size/1024/1024))
          success <- TRUE
        } else {
          warning(sprintf("Failed to download SRA file: %s", sra_id))
          retry_count <- retry_count + 1
        }
      }
      
      if (!success) {
        warning(sprintf("Failed to download SRA file after %d attempts: %s", max_retries, sra_id))
      }
    }
  }
}

# Function to convert SRA to FASTQ
convert_sra_to_fastq <- function(gse_id, sra_ids, threads = NULL) {
  if (!check_command("fasterq-dump")) {
    stop("fasterq-dump command not found. Please ensure SRA toolkit is installed and in PATH")
  }
  
  if (!check_command("pigz")) {
    stop("pigz command not found. Please ensure pigz is installed and in PATH")
  }
  
  # Determine optimal number of threads based on available resources
  if (is.null(threads)) {
    # Get number of CPUs from SLURM environment variable
    slurm_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
    # Use 75% of available CPUs to leave some for system processes
    threads <- max(1, floor(slurm_cpus * 0.75))
    cat(sprintf("Using %d threads for fasterq-dump (from %d available CPUs)\n", 
                threads, slurm_cpus))
  }
  
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Debug: Check if GSE directory exists
  if (!dir.exists(gse_dir)) {
    stop(sprintf("GSE directory does not exist: %s", gse_dir))
  }
  
  # Debug: List contents of GSE directory
  cat(sprintf("Contents of GSE directory %s:\n", gse_dir))
  gse_contents <- list.dirs(gse_dir, recursive = FALSE, full.names = TRUE)
  for (item in gse_contents) {
    cat(sprintf("  - %s\n", item))
  }
  
  # Create a temporary directory for faster I/O
  temp_dir <- file.path(base_dir, "temp_fastq")
  if (!dir.exists(temp_dir)) {
    cat(sprintf("Creating temporary directory: %s\n", temp_dir))
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Track which SRA IDs have been successfully converted
  converted_sra_ids <- c()
  
  # Process SRA files in parallel if possible
  if (length(sra_ids) > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cat(sprintf("Processing %d SRA files in parallel...\n", length(sra_ids)))
    
    # Determine number of parallel processes based on available resources
    # Use at most 4 parallel processes or 75% of available CPUs, whichever is smaller
    slurm_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
    max_parallel <- min(4, max(1, floor(slurm_cpus / 2)))
    num_cores <- min(max_parallel, parallel::detectCores() - 1)
    if (num_cores < 1) num_cores <- 1
    
    cat(sprintf("Using %d parallel processes for SRA conversion\n", num_cores))
    
    # Create a cluster
    cl <- parallel::makeCluster(num_cores)
    
    # Export necessary variables to the cluster
    parallel::clusterExport(cl, c("gse_dir", "temp_dir", "threads"), envir = environment())
    
    # Process SRA files in parallel with error handling
    results <- parallel::parSapply(cl, sra_ids, function(sra_id) {
      tryCatch({
        cat(sprintf("Starting conversion of SRA %s\n", sra_id))
        result <- process_single_sra(sra_id, gse_dir, temp_dir, threads)
        cat(sprintf("Completed conversion of SRA %s: %s\n", sra_id, ifelse(result, "SUCCESS", "FAILED")))
        return(result)
      }, error = function(e) {
        cat(sprintf("Error processing SRA %s: %s\n", sra_id, e$message))
        return(FALSE)
      })
    })
    
    # Stop the cluster
    parallel::stopCluster(cl)
    
    # Check results
    converted_sra_ids <- sra_ids[results]
    
  } else {
    # Process SRA files sequentially
    for (sra_id in sra_ids) {
      cat(sprintf("Processing SRA %s sequentially\n", sra_id))
      if (process_single_sra(sra_id, gse_dir, temp_dir, threads)) {
        converted_sra_ids <- c(converted_sra_ids, sra_id)
      }
    }
  }
  
  # Clean up temporary directory
  if (dir.exists(temp_dir)) {
    cat(sprintf("Cleaning up temporary directory: %s\n", temp_dir))
    unlink(temp_dir, recursive = TRUE)
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
  
  return(converted_sra_ids)
}

# Function to process a single SRA file
process_single_sra <- function(sra_id, gse_dir, temp_dir, threads) {
  cat(sprintf("Processing SRA %s in directory %s\n", sra_id, gse_dir))
  
  # Find all GSM directories
  gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
  cat(sprintf("Found %d GSM directories\n", length(gsm_dirs)))
  
  # Debug: List all GSM directories
  cat("GSM directories:\n")
  for (dir in gsm_dirs) {
    cat(sprintf("  - %s\n", dir))
  }
  
  # Create a mapping from GSM IDs to their directories
  gsm_to_dir_mapping <- list()
  for (gsm_dir in gsm_dirs) {
    gsm_id <- basename(gsm_dir)
    gsm_to_dir_mapping[[gsm_id]] <- gsm_dir
  }
  
  # Create a mapping from SRX IDs to GSM IDs by reading the GSM objects
  srx_to_gsm_mapping <- list()
  for (gsm_dir in gsm_dirs) {
    gsm_id <- basename(gsm_dir)
    rds_file <- file.path(gsm_dir, paste0(gsm_id, ".rds"))
    
    if (file.exists(rds_file)) {
      gsm <- readRDS(rds_file)
      relations <- Meta(gsm)$relation
      sra_link <- relations[grep("SRA:", relations)]
      
      if (length(sra_link) > 0) {
        # Extract SRX ID from the SRA link
        srx_id <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
        # Clean up the SRX ID
        srx_id <- gsub("\\[.*\\]", "", srx_id)  # Remove any brackets and their contents
        srx_id <- trimws(srx_id)  # Remove any whitespace
        
        # Store the mapping
        srx_to_gsm_mapping[[srx_id]] <- gsm_id
        cat(sprintf("Mapped SRX %s to GSM %s\n", srx_id, gsm_id))
      }
    }
  }
  
  # Find which GSM directories contain this SRA ID
  gsm_ids_with_sra <- c()
  for (gsm_dir in gsm_dirs) {
    gsm_id <- basename(gsm_dir)
    sra_dir <- file.path(gsm_dir, "SRA")
    sra_file <- file.path(sra_dir, sra_id, paste0(sra_id, ".sra"))
    
    if (file.exists(sra_file) && file.info(sra_file)$size > 0) {
      gsm_ids_with_sra <- c(gsm_ids_with_sra, gsm_id)
    }
  }
  
  if (length(gsm_ids_with_sra) == 0) {
    cat(sprintf("No GSM directories found with valid SRA file for SRA ID: %s\n", sra_id))
    return(FALSE)
  }
  
  cat(sprintf("Found SRA %s in %d GSM directories: %s\n", 
              sra_id, length(gsm_ids_with_sra), paste(gsm_ids_with_sra, collapse = ", ")))
  
  # Process the SRA file for each GSM directory
  success_count <- 0
  for (gsm_id in gsm_ids_with_sra) {
    gsm_dir <- gsm_to_dir_mapping[[gsm_id]]
    cat(sprintf("Processing SRA %s for GSM %s\n", sra_id, gsm_id))
    
    # Construct paths
    sra_dir <- file.path(gsm_dir, "SRA")
    sra_file <- file.path(sra_dir, sra_id, paste0(sra_id, ".sra"))
    fastq_dir <- file.path(sra_dir, "FASTQ")
    
    cat(sprintf("Looking for SRA file: %s\n", sra_file))
    
    # Debug: Check if SRA directory exists
    if (!dir.exists(sra_dir)) {
      cat(sprintf("SRA directory does not exist: %s\n", sra_dir))
      next
    }
    
    # Debug: List contents of SRA directory
    cat(sprintf("Contents of SRA directory %s:\n", sra_dir))
    sra_contents <- list.dirs(sra_dir, recursive = FALSE, full.names = TRUE)
    for (item in sra_contents) {
      cat(sprintf("  - %s\n", item))
    }
    
    # Check if FASTQ files already exist for this SRA ID
    if (dir.exists(fastq_dir)) {
      existing_fastq_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq\\.gz$"), full.names = TRUE)
      if (length(existing_fastq_files) > 0) {
        cat(sprintf("FASTQ files already exist for SRA %s in %s. Skipping conversion.\n", 
                   sra_id, fastq_dir))
        success_count <- success_count + 1
        next
      }
    }
    
    # Check if SRA file exists and has size > 0
    if (file.exists(sra_file)) {
      file_size <- file.info(sra_file)$size
      cat(sprintf("Found SRA file: %s (size: %.2f MB)\n", sra_file, file_size/1024/1024))
      
      if (file_size > 0) {
        cat(sprintf("Converting SRA %s to FASTQ in %s\n", sra_id, fastq_dir))
        
        # Create FASTQ directory if it doesn't exist
        if (!dir.exists(fastq_dir)) {
          cat(sprintf("Creating FASTQ directory: %s\n", fastq_dir))
          dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        # Convert with timeout and retry logic
        max_retries <- 3
        retry_count <- 0
        success <- FALSE
        
        while (!success && retry_count < max_retries) {
          if (retry_count > 0) {
            delay_time <- 30 * retry_count
            cat(sprintf("Retry %d: Waiting %d seconds before trying again...\n", 
                       retry_count, delay_time))
            Sys.sleep(delay_time)
          }
          
          # Create a temporary directory for this SRA
          sra_temp_dir <- file.path(temp_dir, sra_id)
          if (!dir.exists(sra_temp_dir)) {
            cat(sprintf("Creating temporary directory: %s\n", sra_temp_dir))
            dir.create(sra_temp_dir, recursive = TRUE, showWarnings = FALSE)
          }
          
          # First use fasterq-dump to create uncompressed FASTQ files with optimized settings
          # Use more aggressive memory settings for faster processing
          cmd <- sprintf('timeout 7200 fasterq-dump --split-files --split-spot --progress --threads %d --mem 2G --bufsize 16M --curcache 200M --disk-limit 200G --temp %s --outdir %s %s 2>&1',
                        threads, sra_temp_dir, fastq_dir, sra_id)
          
          cat(sprintf("Running command: %s\n", cmd))
          result <- system(cmd, intern = TRUE)
          cat("fasterq-dump output:\n", paste(result, collapse = "\n"), "\n")
          
          # Check for common error patterns in output
          if (any(grepl("error|failed|timeout|connection refused|disk space", tolower(result)))) {
            cat("Error detected in fasterq-dump output, will retry...\n")
            retry_count <- retry_count + 1
            next
          }
          
          # Check for uncompressed FASTQ files
          fastq_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq$"), full.names = TRUE)
          cat(sprintf("Found %d uncompressed FASTQ files\n", length(fastq_files)))
          
          if (length(fastq_files) > 0) {
            cat(sprintf("Successfully created FASTQ files, now compressing with pigz...\n"))
            
            # Compress each FASTQ file with pigz in parallel
            if (length(fastq_files) > 1 && requireNamespace("parallel", quietly = TRUE)) {
              # Use parallel compression for multiple files
              # Use all available threads for compression
              num_cores <- min(length(fastq_files), parallel::detectCores())
              if (num_cores < 1) num_cores <- 1
              
              cat(sprintf("Using %d cores for parallel compression\n", num_cores))
              
              cl <- parallel::makeCluster(num_cores)
              parallel::clusterExport(cl, c("threads"), envir = environment())
              
              parallel::parSapply(cl, fastq_files, function(fastq_file) {
                cat(sprintf("Compressing %s with pigz...\n", basename(fastq_file)))
                compress_cmd <- sprintf('pigz -p %d %s', threads, fastq_file)
                system(compress_cmd, intern = TRUE)
              })
              
              parallel::stopCluster(cl)
            } else {
              # Sequential compression
              for (fastq_file in fastq_files) {
                cat(sprintf("Compressing %s with pigz...\n", basename(fastq_file)))
                compress_cmd <- sprintf('pigz -p %d %s', threads, fastq_file)
                compress_result <- system(compress_cmd, intern = TRUE)
                if (length(compress_result) > 0) {
                  cat("pigz output:\n", paste(compress_result, collapse = "\n"), "\n")
                }
              }
            }
            
            # Verify compressed files exist
            compressed_files <- list.files(fastq_dir, pattern = paste0(sra_id, ".*\\.fastq\\.gz$"), full.names = TRUE)
            cat(sprintf("Found %d compressed FASTQ files\n", length(compressed_files)))
            
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
              
              # Clean up temporary directory
              if (dir.exists(sra_temp_dir)) {
                unlink(sra_temp_dir, recursive = TRUE)
              }
              
              success <- TRUE
              success_count <- success_count + 1
              break
            } else {
              warning(sprintf("Failed to compress FASTQ files for SRA %s", sra_id))
              retry_count <- retry_count + 1
            }
          } else {
            warning(sprintf("No FASTQ files created for SRA %s", sra_id))
            retry_count <- retry_count + 1
          }
        }
        
        if (!success) {
          stop(sprintf("Failed to convert SRA file after %d attempts: %s", max_retries, sra_id))
        }
      } else {
        cat(sprintf("SRA file exists but has zero size: %s\n", sra_file))
      }
    } else {
      cat(sprintf("SRA file not found: %s\n", sra_file))
      
      # Debug: Check if the SRA ID directory exists
      sra_id_dir <- file.path(sra_dir, sra_id)
      if (dir.exists(sra_id_dir)) {
        cat(sprintf("SRA ID directory exists: %s\n", sra_id_dir))
        cat("Contents of SRA ID directory:\n")
        sra_id_contents <- list.files(sra_id_dir, full.names = TRUE)
        for (item in sra_id_contents) {
          cat(sprintf("  - %s\n", item))
        }
      } else {
        cat(sprintf("SRA ID directory does not exist: %s\n", sra_id_dir))
      }
    }
  }
  
  if (success_count == 0) {
    cat(sprintf("Failed to process SRA ID %s for any GSM directory\n", sra_id))
    return(FALSE)
  }
  
  cat(sprintf("Successfully processed SRA ID %s for %d out of %d GSM directories\n", 
              sra_id, success_count, length(gsm_ids_with_sra)))
  return(TRUE)
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

# Main function to download and process SRA files
download_and_process_sra <- function(gse_id, api_key = NULL, threads = 4) {
  # Set up logging
  log_file <- file.path(getwd(), paste0(gse_id, "_download.log"))
  cat(sprintf("Starting download and processing for %s\n", gse_id), file = log_file)
  
  # Check for required commands
  check_required_commands()
  
  # Load required packages
  load_required_packages()
  
  # Validate API key
  if (!is.null(api_key)) {
    validate_api_key(api_key)
  }
  
  # Create directory structure
  gse_dir <- create_directory_structure(gse_id)
  
  # Get GSM IDs and SRX IDs
  cat("Getting GSM IDs and SRX IDs...\n", file = log_file, append = TRUE)
  gsm_data <- get_srx_ids(gse_id, gse_dir, api_key)
  
  if (length(gsm_data$gsm_ids) == 0) {
    cat("No GSM IDs found. Exiting.\n", file = log_file, append = TRUE)
    return(FALSE)
  }
  
  # Convert SRX IDs to SRA IDs
  cat("Converting SRX IDs to SRA IDs...\n", file = log_file, append = TRUE)
  sra_data <- convert_srx_to_sra(gsm_data$srx_ids, api_key)
  
  if (length(sra_data$sra_ids) == 0) {
    cat("No SRA IDs found. Exiting.\n", file = log_file, append = TRUE)
    return(FALSE)
  }
  
  # Download SRA files
  cat("Downloading SRA files...\n", file = log_file, append = TRUE)
  download_sra_files(sra_data, gse_dir, api_key)
  
  # Create temporary directory for processing
  temp_dir <- file.path(gse_dir, "temp")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Process SRA files
  cat("Processing SRA files...\n", file = log_file, append = TRUE)
  
  # Get all unique SRA IDs
  all_sra_ids <- unique(unlist(sra_data$srx_to_sra_mapping))
  cat(sprintf("Found %d unique SRA IDs to process\n", length(all_sra_ids)), file = log_file, append = TRUE)
  
  # Process each SRA ID
  success_count <- 0
  for (sra_id in all_sra_ids) {
    cat(sprintf("Processing SRA ID: %s\n", sra_id), file = log_file, append = TRUE)
    
    # Process the SRA file
    result <- process_single_sra(sra_id, gse_dir, temp_dir, threads)
    
    if (result) {
      success_count <- success_count + 1
      cat(sprintf("Successfully processed SRA ID: %s\n", sra_id), file = log_file, append = TRUE)
    } else {
      cat(sprintf("Failed to process SRA ID: %s\n", sra_id), file = log_file, append = TRUE)
    }
  }
  
  # Clean up temporary directory
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
  
  # Log final results
  cat(sprintf("Processed %d out of %d SRA IDs successfully\n", 
              success_count, length(all_sra_ids)), file = log_file, append = TRUE)
  
  return(success_count == length(all_sra_ids))
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
    sra_data <- convert_srx_to_sra(srx_ids)
    if (is.null(sra_data)) {
      stop(sprintf("Failed to convert SRX to SRA IDs for %s", gse_id))
    }
    
    # Download and convert files
    cat("Downloading SRA files...\n")
    download_sra_files(gse_id, sra_data)
    
    cat("Converting SRA to FASTQ...\n")
    convert_sra_to_fastq(gse_id, sra_data$sra_ids)
    
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