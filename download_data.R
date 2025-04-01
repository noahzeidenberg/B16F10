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

# Function to create directory structure for a GSE ID
create_gse_structure <- function(gse_id) {
  # Create main GSE directory
  gse_dir <- file.path(".", gse_id)
  dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories
  dirs <- c(
    file.path(gse_dir, "samples"),
    file.path(gse_dir, "logs")
  )
  
  for (dir in dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  return(gse_dir)
}

# Function to create sample directory structure
create_sample_structure <- function(gse_dir, gsm_id) {
  # Create sample directory
  sample_dir <- file.path(gse_dir, "samples", gsm_id)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories
  dirs <- c(
    file.path(sample_dir, "SRA"),
    file.path(sample_dir, "SRA", "FASTQ")
  )
  
  for (dir in dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  
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
  
  # Save each GSM object and extract SRX IDs
  srx_ids <- sapply(names(gsm_list), function(gsm_name) {
    gsm <- gsm_list[[gsm_name]]
    
    # Create sample directory structure
    sample_dir <- create_sample_structure(gse_dir, gsm_name)
    
    # Save GSM object
    saveRDS(gsm, file.path(sample_dir, paste0(gsm_name, ".rds")))
    
    # Extract SRX ID
    relations <- Meta(gsm)$relation
    sra_link <- relations[grep("SRA:", relations)]
    if (length(sra_link) > 0) {
      return(sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link))
    }
    return(NA)
  })
  
  # Remove any NA values
  srx_ids <- srx_ids[!is.na(srx_ids)]
  
  if (length(srx_ids) == 0) {
    cat(sprintf("No SRX IDs found for GSE ID: %s. Skipping...\n", gse_id))
    return(NULL)
  }
  
  return(srx_ids)
}

# Function to convert SRX IDs to SRA IDs using rentrez
convert_srx_to_sra <- function(srx_ids) {
  # Initialize vector to store SRA IDs
  sra_ids <- character(length(srx_ids))
  
  # Process each SRX ID
  for (i in seq_along(srx_ids)) {
    srx_id <- srx_ids[i]
    cat(sprintf("Processing SRX ID %d of %d: %s\n", i, length(srx_ids), srx_id))
    
    tryCatch({
      # Use rentrez to get SRR ID
      xml_result <- rentrez::entrez_fetch(db = "sra", id = srx_id, rettype = "xml")
      xml_doc <- xml2::read_xml(xml_result)
      srr_id <- xml2::xml_find_all(xml_doc, ".//RUN") |> xml2::xml_attr("accession")
      
      if (length(srr_id) > 0 && grepl("^SRR", srr_id)) {
        sra_ids[i] <- srr_id
        cat(sprintf("Found SRA ID: %s\n", srr_id))
      } else {
        cat(sprintf("No valid SRA ID found for SRX ID: %s\n", srx_id))
        sra_ids[i] <- NA
      }
    }, error = function(e) {
      cat(sprintf("Error processing SRX ID %s: %s\n", srx_id, e$message))
      sra_ids[i] <- NA
    })
    
    # Add a small delay to avoid hitting rate limits
    Sys.sleep(0.5)
  }
  
  # Remove any NA values
  sra_ids <- sra_ids[!is.na(sra_ids)]
  
  if (length(sra_ids) == 0) {
    stop("No SRA IDs found for any of the SRX IDs")
  }
  
  return(sra_ids)
}

# Function to download SRA files using prefetch
download_sra_files <- function(gse_id, sra_ids) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Download SRA files for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_dir <- file.path(gsm_dir, "SRA")
      cmd <- sprintf('prefetch --max-size 100G -O %s -p %s', 
                     sra_dir, sra_id)
      system(cmd, intern = TRUE)
    }
  }
}

# Function to convert SRA to FASTQ
convert_sra_to_fastq <- function(gse_id, sra_ids, threads = 8) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Convert SRA to FASTQ for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_file <- file.path(gsm_dir, "SRA", paste0(sra_id, ".sra"))
      fastq_dir <- file.path(gsm_dir, "SRA", "FASTQ")
      
      if (file.exists(sra_file)) {
        cmd <- sprintf('fasterq-dump --split-files --progress --threads %d --outdir %s %s',
                       threads, fastq_dir, sra_id)
        system(cmd, intern = TRUE)
      }
    }
  }
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    stop("Please provide a valid GSE ID. For example: main('GSE223515')")
  }
  
  tryCatch({
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
    cat("You can verify the GSE ID at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n", gse_id)
  })
}

# Example usage:
# main("GSE223515")  # Replace with your GSE ID 