#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("GEOquery", "DESeq2", "clusterProfiler", "org.Mm.eg.db", "tidyverse")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Read the GDS table subset
gds_df <- read.csv("gds_table_subset.csv")

# Function to check if a GSE is RNA-seq
is_rnaseq <- function(gse_id) {
  tryCatch({
    gse <- getGEO(gse_id)
    # Check if it's RNA-seq by looking at the technology type and library strategy
    tech <- unique(unlist(lapply(gse, function(x) x@experimentData@other$type)))
    lib_strategy <- unique(unlist(lapply(gse, function(x) x@phenoData@data$library_strategy)))
    return(any(grepl("RNA-seq|RNA sequencing", tech, ignore.case = TRUE)) ||
           any(grepl("RNA-Seq", lib_strategy, ignore.case = TRUE)))
  }, error = function(e) {
    message(paste("Error checking RNA-seq status for", gse_id, ":", e$message))
    return(FALSE)
  })
}

# Function to process a single RNA-seq dataset
process_rnaseq <- function(gse_id) {
  message(paste("Processing", gse_id))
  
  tryCatch({
    # Download the data
    gse <- getGEO(gse_id)
    
    # Extract expression matrix and phenotype data
    expr_matrix <- exprs(gse[[1]])
    pheno_data <- pData(gse[[1]])
    
    # Get experiment metadata
    exp_metadata <- gse[[1]]@experimentData
    
    # Create results directory for this dataset
    results_dir <- file.path("results", gse_id)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save metadata
    metadata_list <- list(
      experiment = list(
        title = exp_metadata@title,
        abstract = exp_metadata@abstract,
        pubMedIds = exp_metadata@pubMedIds,
        type = exp_metadata@other$type,
        platform_id = exp_metadata@other$platform_id,
        overall_design = exp_metadata@other$overall_design
      ),
      samples = pheno_data
    )
    
    saveRDS(metadata_list, file = file.path(results_dir, paste0(gse_id, "_metadata.rds")))
    
    # Check for supplementary files with counts
    supp_files <- exp_metadata@other$supplementary_file
    if (!is.null(supp_files)) {
      message("Found supplementary files")
      message(paste("Files:", paste(supp_files, collapse="\n    ")))
      
      # Look for count files
      count_files <- supp_files[grep("counts|Counts", supp_files)]
      if (length(count_files) > 0) {
        message("Found count files")
        # Download and process count files
        for (file in count_files) {
          message(paste("Downloading:", file))
          tryCatch({
            download.file(file, destfile = file.path(results_dir, basename(file)), mode = "wb")
          }, error = function(e) {
            message(paste("Failed to download:", file, "-", e$message))
          })
        }
      }
    }
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = round(expr_matrix),
      colData = pheno_data,
      design = ~ 1  # You'll need to adjust this based on your experimental design
    )
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds)
    
    # Convert to data frame
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    # Perform KEGG pathway analysis
    gene_list <- res_df$gene[!is.na(res_df$padj) & res_df$padj < 0.05]
    
    if (length(gene_list) > 0) {
      kegg_result <- enrichKEGG(
        gene = gene_list,
        organism = 'mmu',
        keyType = 'kegg',
        pvalueCutoff = 0.05
      )
      
      # Save results
      write.csv(res_df, file = file.path(results_dir, paste0(gse_id, "_DESeq2_results.csv")))
      write.csv(as.data.frame(kegg_result), file = file.path(results_dir, paste0(gse_id, "_KEGG_results.csv")))
      
      # Save expression matrix and phenotype data
      write.csv(expr_matrix, file = file.path(results_dir, paste0(gse_id, "_expression_matrix.csv")))
      write.csv(pheno_data, file = file.path(results_dir, paste0(gse_id, "_phenotype_data.csv")))
      
      # Save DESeq2 object
      saveRDS(dds, file = file.path(results_dir, paste0(gse_id, "_DESeq2_object.rds")))
      
      message(paste("Successfully processed and saved results for", gse_id))
    } else {
      message(paste("No significant genes found for", gse_id))
    }
    
    return(list(deseq_results = res_df, kegg_results = kegg_result))
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
    return(NULL)
  })
}

# Create main results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Process each dataset
results <- list()
for (gse_id in gds_df$Accession) {
  tryCatch({
    if (is_rnaseq(gse_id)) {
      message(paste("\nProcessing RNA-seq dataset:", gse_id))
      results[[gse_id]] <- process_rnaseq(gse_id)
    } else {
      message(paste("\nSkipping non-RNA-seq dataset:", gse_id))
    }
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
  })
}

# Save session info
sink("results/session_info.txt")
sessionInfo()
sink() 