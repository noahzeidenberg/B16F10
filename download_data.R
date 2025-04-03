# SRA Download and Conversion Pipeline
# This script downloads SRA files and converts them to FASTQ format

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
  # First check if we're running in SLURM
  if (Sys.getenv("SLURM_TMPDIR") != "") {
    return(Sys.getenv("SLURM_TMPDIR"))
  }
  
  # If not in SLURM, try to use the current working directory
  current_dir <- getwd()
  
  # Check if we have write permissions in the current directory
  test_file <- file.path(current_dir, ".test_write")
  if (file.create(test_file)) {
    file.remove(test_file)
    return(current_dir)
  }
  
  # If we don't have write permissions, use a temporary directory
  temp_dir <- file.path(tempdir(), "B16F10")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  return(temp_dir)
}

# Function to create directory structure for a GSE ID
create_gse_structure <- function(gse_id) {
  # Create main GSE directory in the appropriate directory
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Create directory with error checking
  if (!dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)) {
    stop(sprintf("Failed to create GSE directory: %s", gse_dir))
  }
  
  # Create subdirectories
  dirs <- c(
    file.path(gse_dir, "samples"),
    file.path(gse_dir, "logs")
  )
  
  for (dir in dirs) {
    if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) {
      stop(sprintf("Failed to create directory: %s", dir))
    }
  }
  
  cat(sprintf("Created directory structure at: %s\n", gse_dir))
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

# Function to get SRX IDs from GSM IDs and save GSM objects
get_srx_ids <- function(gse_id) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Get the GEO object without matrix
  cat(sprintf("Fetching GSM information for %s...\n", gse_id))
  gse <- getGEO(gse_id, GSEMatrix = FALSE)
  
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
}

# Function to convert SRX IDs to SRA IDs using rentrez
convert_srx_to_sra <- function(srx_ids) {
  # Initialize list to store SRA IDs (using a list instead of vector to handle multiple SRRs)
  sra_ids <- list()
  
  # Process each SRX ID
  for (i in seq_along(srx_ids)) {
    srx_id <- srx_ids[i]
    cat(sprintf("Processing SRX ID %d of %d: %s\n", i, length(srx_ids), srx_id))
    
    tryCatch({
      # Use rentrez to get SRR ID(s)
      xml_result <- rentrez::entrez_fetch(db = "sra", id = srx_id, rettype = "xml")
      xml_doc <- xml2::read_xml(xml_result)
      srr_ids <- xml2::xml_find_all(xml_doc, ".//RUN") |> xml2::xml_attr("accession")
      
      # Filter for valid SRR IDs
      valid_srr_ids <- srr_ids[grepl("^SRR", srr_ids)]
      
      if (length(valid_srr_ids) > 0) {
        sra_ids[[srx_id]] <- valid_srr_ids
        cat(sprintf("Found %d SRA ID(s): %s\n", 
                   length(valid_srr_ids), 
                   paste(valid_srr_ids, collapse = ", ")))
      } else {
        cat(sprintf("No valid SRA IDs found for SRX ID: %s\n", srx_id))
      }
    }, error = function(e) {
      cat(sprintf("Error processing SRX ID %s: %s\n", srx_id, e$message))
    })
    
    # Add a small delay to avoid hitting rate limits
    Sys.sleep(0.5)
  }
  
  # Flatten the list of SRA IDs
  all_sra_ids <- unlist(sra_ids)
  
  if (length(all_sra_ids) == 0) {
    stop("No SRA IDs found for any of the SRX IDs")
  }
  
  return(all_sra_ids)
}

# Function to download SRA files using prefetch
download_sra_files <- function(gse_id, sra_ids) {
  # Use the existing GSE directory structure
  base_dir <- get_base_dir()
  gse_dir <- file.path(base_dir, gse_id)
  
  # Check available space in the directory
  available_space <- check_disk_space(base_dir)
  if (!is.na(available_space) && available_space < 50) {  # Require at least 50GB free
    stop(sprintf("Insufficient disk space in directory. Only %.1fGB available. Need at least 50GB.", 
                 available_space))
  }
  
  # Download SRA files for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_dir <- file.path(gsm_dir, "SRA")
      
      # Check space again before each download
      available_space <- check_disk_space(base_dir)
      if (!is.na(available_space) && available_space < 20) {  # Require at least 20GB free for each file
        stop(sprintf("Insufficient disk space for downloading %s. Only %.1fGB available.", 
                     sra_id, available_space))
      }
      
      cmd <- sprintf('prefetch --max-size 100G -O %s -p %s', 
                     sra_dir, sra_id)
      result <- tryCatch({
        system(cmd, intern = TRUE)
      }, error = function(e) {
        warning(sprintf("Error downloading %s: %s", sra_id, e$message))
        return(NULL)
      })
      
      # Check if download was successful
      if (!is.null(result) && length(grep("error|failed", result, ignore.case = TRUE)) > 0) {
        warning(sprintf("Error downloading %s: %s", sra_id, paste(result, collapse="\n")))
      }
    }
  }
}

# Function to convert SRA to FASTQ
convert_sra_to_fastq <- function(gse_id, sra_ids, threads = 8) {
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
        # Convert to FASTQ
        cmd <- sprintf('fasterq-dump --split-files --progress --threads %d --outdir %s %s',
                      threads, fastq_dir, sra_id)
        system(cmd, intern = TRUE)
        
        # Compress the FASTQ files
        fastq_files <- list.files(fastq_dir, pattern = "\\.fastq$", full.names = TRUE)
        for (fastq_file in fastq_files) {
          cmd <- sprintf("pigz -p %d %s", threads, fastq_file)
          system(cmd)
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
      stop("Failed to get SRX IDs")
    }
    
    # Convert SRX to SRA IDs
    cat("Converting SRX to SRA IDs...\n")
    sra_ids <- convert_srx_to_sra(srx_ids)
    if (is.null(sra_ids)) {
      stop("Failed to convert SRX to SRA IDs")
    }
    
    # Download SRA files
    cat("Downloading SRA files...\n")
    download_sra_files(gse_id, sra_ids)
    
    # Convert SRA to FASTQ
    cat("Converting SRA to FASTQ...\n")
    convert_sra_to_fastq(gse_id, sra_ids)
    
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