#!/usr/bin/env Rscript

# Load required libraries
library(GEOquery)  # For platform info
library(affy)      # For Affymetrix arrays
library(oligo)     # For modern Affymetrix arrays
library(DESeq2)    # For RNA-seq data

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript process_methylation.R <input_dir> <output_dir>")
}

input_dir <- args[1]
output_dir <- args[2]

# Function to determine platform type
get_platform_type <- function(platform_id) {
  platform_info <- getGEO(platform_id)
  
  # Check platform technology
  tech <- platform_info@header$technology
  manufacturer <- platform_info@header$manufacturer
  
  message("Platform: ", platform_id)
  message("Technology: ", tech)
  message("Manufacturer: ", manufacturer)
  
  if (grepl("high-throughput sequencing", tolower(tech))) {
    return("RNA-seq")
  } else if (grepl("affymetrix", tolower(manufacturer))) {
    return("affymetrix")
  } else if (grepl("illumina", tolower(manufacturer))) {
    return("illumina")
  } else {
    return("unknown")
  }
}

# Process expression data
process_expression <- function(input_dir, output_dir, platform_id) {
  tryCatch({
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Determine platform type
    platform_type <- get_platform_type(platform_id)
    message("Detected platform type: ", platform_type)
    
    if (platform_type == "affymetrix") {
      # Process Affymetrix array data
      cel_files <- list.files(input_dir, pattern = "\\.cel$|\\.CEL$|\\.cel.gz$|\\.CEL.gz$", 
                            full.names = TRUE, recursive = TRUE)
      if (length(cel_files) > 0) {
        message("Processing Affymetrix array data...")
        raw_data <- read.celfiles(cel_files)
        eset <- rma(raw_data)
        expression_matrix <- exprs(eset)
      } else {
        stop("No CEL files found")
      }
      
    } else if (platform_type == "illumina") {
      # Process Illumina array data
      idat_files <- list.files(input_dir, pattern = "\\.idat$|\\.idat.gz$", 
                             full.names = TRUE, recursive = TRUE)
      if (length(idat_files) > 0) {
        message("Processing Illumina array data...")
        raw_data <- read.idat(idat_files)
        expression_matrix <- normalizeIllumina(raw_data)
      } else {
        stop("No IDAT files found")
      }
      
    } else if (platform_type == "RNA-seq") {
      # Process RNA-seq data
      count_files <- list.files(input_dir, pattern = "counts\\.txt$", 
                              full.names = TRUE, recursive = TRUE)
      if (length(count_files) > 0) {
        message("Processing RNA-seq data...")
        counts <- read.delim(count_files[1], row.names=1)
        dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = data.frame(row.names=colnames(counts)),
                                    design = ~ 1)
        vsd <- vst(dds)
        expression_matrix <- assay(vsd)
      } else {
        stop("No count files found")
      }
    } else {
      stop("Unsupported platform type: ", platform_type)
    }
    
    # Save results
    message("Saving processed data...")
    saveRDS(expression_matrix, file = file.path(output_dir, "expression_matrix.rds"))
    write.csv(expression_matrix, file = file.path(output_dir, "expression_matrix.csv"))
    
    # Generate QC plots
    message("Generating QC plots...")
    pdf(file.path(output_dir, "qc_plots.pdf"))
    boxplot(expression_matrix, main="Expression Distribution", las=2)
    hist(expression_matrix, main="Expression Histogram")
    dev.off()
    
    message("Successfully processed expression data")
    
  }, error = function(e) {
    stop("Failed to process expression data: ", e$message)
  })
}

# Execute processing
process_expression(input_dir, output_dir, platform_id) 