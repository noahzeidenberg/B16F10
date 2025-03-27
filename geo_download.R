# Function to get SRR accessions from an SRX ID
get_srr_from_srx <- function(srx) {
  tryCatch({
    message(paste("  Getting SRR accessions for", srx))
    
    # Use rentrez to search for the SRX
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
    srr_ids <- unlist(strsplit(runs, ";"))
    srr_ids <- gsub("^.*?acc=([^,]+).*$", "\\1", srr_ids)
    
    message(paste("  Found SRR accessions:", paste(srr_ids, collapse=", ")))
    return(srr_ids)
  }, error = function(e) {
    message(paste("Warning: Failed to get SRR info for", srx, ":", e$message))
    return(NULL)
  })
}