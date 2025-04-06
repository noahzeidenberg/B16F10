# Design Matrix Generation Script
# This script extracts sample information from GEO and creates design matrices for differential expression analysis

# Load required packages
cat("=== Starting Design Matrix Generation ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "GEOquery",
  "rentrez",
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  cat(sprintf("Loading package: %s\n", pkg))
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

# Function to extract sample information from GEO
extract_sample_info <- function(gse_id) {
  cat(sprintf("Extracting sample information for %s...\n", gse_id))
  
  # Create output directory
  base_dir <- getwd()
  output_dir <- file.path(base_dir, "results", "design_matrices")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if design matrix already exists
  design_file <- file.path(output_dir, paste0(gse_id, "_design_matrix.txt"))
  if (file.exists(design_file)) {
    cat(sprintf("Design matrix already exists for %s. Skipping...\n", gse_id))
    return(TRUE)
  }
  
  # Get the GEO object
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
      gse <- getGEO(gse_id, GSEMatrix = FALSE, getGPL = FALSE)
      success <- TRUE
      
      # Get list of all entities
      all_entities <- GSMList(gse)
      
      # Filter out any non-GSM entries (GPL, GDS, etc.)
      gsm_list <- all_entities[grep("^GSM", names(all_entities))]
      total_gsm <- length(gsm_list)
      
      # Log only the GSM entities found
      cat(sprintf("Found %d GSM samples to process\n", total_gsm))
      
      if (total_gsm == 0) {
        cat("No GSM samples found. Skipping...\n")
        return(FALSE)
      }
      
      # Extract sample information
      sample_info <- data.frame(
        Sample_geo_accession = character(),
        Sample_title = character(),
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(names(gsm_list))) {
        gsm_name <- names(gsm_list)[i]
        gsm <- gsm_list[[gsm_name]]
        
        cat(sprintf("Processing GSM %d of %d: %s\n", i, total_gsm, gsm_name))
        
        # Extract sample title
        sample_title <- Meta(gsm)$title
        if (is.null(sample_title)) {
          sample_title <- gsm_name
        }
        
        # Add to data frame
        sample_info <- rbind(sample_info, data.frame(
          Sample_geo_accession = gsm_name,
          Sample_title = sample_title,
          stringsAsFactors = FALSE
        ))
      }
      
      # Save sample information
      sample_file <- file.path(output_dir, paste0(gse_id, "_sample_info.txt"))
      write.table(sample_info, sample_file, sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("Saved sample information to: %s\n", sample_file))
      
      # Generate design matrix
      # This is a simple example - in practice, you would need to parse the sample titles
      # to extract meaningful group information
      design_matrix <- data.frame(
        Sample_geo_accession = sample_info$Sample_geo_accession,
        Group = rep("Group1", nrow(sample_info)),  # Placeholder - replace with actual group information
        stringsAsFactors = FALSE
      )
      
      # Save design matrix
      write.table(design_matrix, design_file, sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("Saved design matrix to: %s\n", design_file))
      
      return(TRUE)
      
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
        return(FALSE)
      }
    })
  }
  
  return(FALSE)
}

# Function to parse sample titles to extract group information
parse_sample_titles <- function(sample_info) {
  # This function should be customized based on the naming conventions in your datasets
  # For example, if sample titles follow the pattern "Group_Replicate", you could extract the group
  
  # Example implementation:
  groups <- sapply(sample_info$Sample_title, function(title) {
    # Extract group information based on your naming convention
    # For example, if titles are like "Control_1", "Treatment_1", etc.
    parts <- strsplit(title, "_")[[1]]
    if (length(parts) > 1) {
      return(parts[1])
    } else {
      return("Unknown")
    }
  })
  
  return(groups)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript generate_design_matrices.R GSE223515")
    }
  }
  
  tryCatch({
    cat(sprintf("Generating design matrix for %s...\n", gse_id))
    success <- extract_sample_info(gse_id)
    
    if (success) {
      cat(sprintf("Successfully generated design matrix for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to generate design matrix for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during design matrix generation: %s\n", e$message))
    stop(sprintf("Error during design matrix generation: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 