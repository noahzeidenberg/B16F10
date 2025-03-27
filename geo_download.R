# Load required libraries
library(GEOquery)
library(SRAdb)
library(dotenv)
library(rentrez)

# Download GEO data
download_geo_data <- function(accession, output_dir) {
  tryCatch({
    # Get GEO data with rate limit handling and better error handling
    message("\nDownloading GEO metadata...")
    gset <- NULL
    
    # Try different methods to get GEO data
    for (method in c("auto", "curl", "wget", "libcurl")) {
      tryCatch({
        message(paste("Trying download method:", method))
        options('download.file.method.GEOquery' = method)
        
        gset <- handle_rate_limit(function() {
          getGEO(accession, GSEMatrix = TRUE, AnnotGPL = FALSE)
        })
        
        if (!is.null(gset)) {
          message(paste("Successfully downloaded using method:", method))
          break
        }
      }, error = function(e) {
        message(paste("Method", method, "failed:", e$message))
      })
    }
    
    if (is.null(gset)) {
      # Try one last time with direct URL construction
      direct_url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=full", accession)
      message(paste("Trying direct URL:", direct_url))
      
      gset <- handle_rate_limit(function() {
        getGEO(filename = direct_url, GSEMatrix = TRUE, AnnotGPL = FALSE)
      })
      
      if (is.null(gset)) {
        stop("Failed to download GEO data using all available methods")
      }
    }
    
    # Verify we have valid data
    if (length(gset) == 0 || is.null(gset[[1]])) {
      stop("Downloaded GEO data appears to be empty or invalid")
    }
    
    # Get platform information in multiple ways
    platform_info <- list()
    
    # Try getting platform from annotation
    if (!is.null(gset[[1]]@annotation) && gset[[1]]@annotation != "") {
      platform_info$annotation <- gset[[1]]@annotation
      message("Platform from annotation: ", platform_info$annotation)
    }
    
    # Try getting platform from platform_id
    if (!is.null(gset[[1]]$platform_id)) {
      platform_info$platform_id <- unique(gset[[1]]$platform_id)
      message("Platform from platform_id: ", paste(platform_info$platform_id, collapse=", "))
    }

    # Create output directory if it doesn't exist
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save platform information
    if (length(platform_info) > 0) {
      saveRDS(platform_info, file = file.path(output_dir, paste0(accession, "_platform_info.rds")))
      write.csv(data.frame(
        source = names(platform_info),
        platform = unlist(platform_info)
      ), file = file.path(output_dir, paste0(accession, "_platform_info.csv")))
    }

    # ... rest of existing function ...
  }, error = function(e) {
    stop(paste("Failed to download GEO data:", e$message))
  })
}

# Function to get SRR accessions from an SRX ID
get_srr_from_srx <- function(srx) {
  tryCatch({
    message(paste("  Getting SRR accessions for", srx))
    
    # Use rentrez to search and link
    srx_search <- handle_rate_limit(function() {
      entrez_search(db = "sra", term = srx)
    })
    
    if (srx_search$count == 0) {
      message("  No SRA records found")
      return(NULL)
    }
    
    # Get the linked SRA records
    sra_links <- handle_rate_limit(function() {
      entrez_link(dbfrom = "sra", db = "sra", id = srx_search$ids)
    })
    
    if (length(sra_links$links$sra_sra) == 0) {
      message("  No SRR accessions found")
      return(NULL)
    }
    
    # Fetch the SRA records
    sra_summary <- handle_rate_limit(function() {
      entrez_summary(db = "sra", id = sra_links$links$sra_sra)
    })
    
    # Extract SRR accessions
    srr_ids <- vapply(sra_summary, function(x) x$runs$run_acc, character(1))
    
    message(paste("  Found SRR accessions:", paste(srr_ids, collapse=", ")))
    return(srr_ids)
  }, error = function(e) {
    message(paste("Warning: Failed to get SRR info for", srx, ":", e$message))
    return(NULL)
  })
} 