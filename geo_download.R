# Load required libraries
library(GEOquery)
library(SRAdb)
library(dotenv)
library(rentrez)

# Load environment variables from .env file
load_dot_env()

# Initialize with API key
api_key <- Sys.getenv("NCBI_API_KEY_1")
if (api_key == "") {
  stop("No NCBI API key found in environment variables. Please check your .env file.")
}
set_entrez_key(api_key)

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
    
    # Search for the SRX
    srx_search <- rentrez::entrez_search(db = "sra", term = srx)
    if (length(srx_search$ids) == 0) return(NULL)
    
    # Get the summary
    srx_summary <- rentrez::entrez_summary(db = "sra", id = srx_search$ids)
    if (is.null(srx_summary$runs)) return(NULL)
    
    # Parse run accessions
    srr_ids <- unlist(strsplit(srx_summary$runs, ";"))
    srr_ids <- gsub("^.*?acc=([^,]+).*$", "\\1", srr_ids)
    
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
    gsm_data <- getGEO(gsm)
    
    # Extract SRA accession from relations
    relations <- gsm_data@header$relation
    sra_link <- relations[grep("SRA:", relations)]
    
    if (length(sra_link) == 0) return(NULL)
    
    sra_acc <- sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link)
    message(paste("  Found SRA accession:", sra_acc))
    return(sra_acc)
  }, error = function(e) {
    message(paste("Warning: Failed to get SRA info for", gsm, ":", e$message))
    return(NULL)
  })
}

# Main download function
download_geo_data <- function(accession, output_dir) {
  tryCatch({
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get GEO data
    message("\nDownloading GEO metadata...")
    gset <- getGEO(accession, GSEMatrix = TRUE)
    if (is.null(gset) || length(gset) == 0) {
      stop("Failed to retrieve GEO dataset")
    }
    
    # Save metadata
    metadata <- pData(gset[[1]])
    write.csv(metadata, file = file.path(output_dir, paste0(accession, "_metadata.csv")))
    
    # Get GSM accessions
    gsm_accessions <- metadata$geo_accession
    message(paste("\nFound", length(gsm_accessions), "GSM accessions"))
    
    # Create samples directory
    samples_dir <- file.path(output_dir, "samples")
    dir.create(samples_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Process each GSM
    for (gsm in gsm_accessions) {
      message(paste("\nProcessing", gsm))
      gsm_dir <- file.path(samples_dir, gsm)
      dir.create(gsm_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Get SRX accession
      srx_id <- get_srx_from_gsm(gsm)
      if (is.null(srx_id)) next
      
      # Get SRR accessions
      srr_ids <- get_srr_from_srx(srx_id)
      if (is.null(srr_ids)) next
      
      # Download FASTQ files
      fastq_dir <- file.path(gsm_dir, "FASTQ")
      dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
      
      for (srr in srr_ids) {
        message(paste("  Downloading", srr))
        cmd <- paste("module load sra-toolkit &&",
                    "prefetch", srr, "&&",
                    "fasterq-dump", srr,
                    "--outdir", fastq_dir,
                    "--split-files",
                    "--threads 8")
        
        if (system(cmd) != 0) {
          message(paste("  Warning: Failed to download", srr))
        } else {
          # Compress FASTQ files
          system(paste("gzip", file.path(fastq_dir, paste0(srr, "_*.fastq"))))
        }
      }
    }
    
    message("\nSuccessfully downloaded GEO data for", accession)
  }, error = function(e) {
    stop(paste("Failed to download GEO data:", e$message))
  })
}

# Execute download
download_geo_data(accession, output_dir) 