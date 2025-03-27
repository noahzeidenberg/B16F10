# Load required libraries
library(GEOquery)
library(SRAdb)
library(dotenv)
library(rentrez)

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

# Initialize with first API key and set it for rentrez
initial_key <- get_next_api_key()
Sys.setenv(ENTREZ_KEY = initial_key)
set_entrez_key(initial_key)

# Update the handle_rate_limit function to also set rentrez key when rotating
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
        set_entrez_key(new_key)  # Add this line to update rentrez key
        message(paste("Rotating to next API key:", substr(new_key, 1, 8), "..."))
        
        Sys.sleep(delay)
      } else {
        stop(e)
      }
    })
  }
  stop("Max retries reached")
}

# Function to get SRR accessions from an SRX ID
get_srr_from_srx <- function(srx) {
  tryCatch({
    message(paste("  Getting SRR accessions for", srx))
    
    # Use rentrez instead of command line tools
    srx_search <- rentrez::entrez_search(db = "sra", term = srx)
    if (length(srx_search$ids) == 0) {
      message("  No SRA records found")
      return(NULL)
    }
    
    # Get the SRA records
    srx_summary <- rentrez::entrez_summary(db = "sra", id = srx_search$ids)
    
    # Extract run information
    runs <- srx_summary$runs
    if (is.null(runs) || length(runs) == 0) {
      message("  No SRR accessions found")
      return(NULL)
    }
    
    # Parse the run accessions
    srr_ids <- unlist(strsplit(runs, ","))
    srr_ids <- gsub("^.*?acc=([^;]+).*$", "\\1", srr_ids)
    
    message(paste("  Found SRR accessions:", paste(srr_ids, collapse=", ")))
    return(srr_ids)
  }, error = function(e) {
    message(paste("Warning: Failed to get SRR info for", srx, ":", e$message))
    return(NULL)
  })
} 