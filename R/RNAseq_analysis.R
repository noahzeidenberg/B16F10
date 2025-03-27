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

# Read the PMID dataframe
pmid_df <- read.csv("PMID_df_test.csv")

# Function to check if a GSE is RNA-seq
is_rnaseq <- function(gse_id) {
  tryCatch({
    gse <- getGEO(gse_id)
    # Check if it's RNA-seq by looking at the technology type
    tech <- unique(unlist(lapply(gse, function(x) x@experimentData@other$type)))
    return(any(grepl("RNA-seq|RNA sequencing", tech, ignore.case = TRUE)))
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
    # Convert gene IDs to ENTREZ IDs (you might need to adjust this based on your data)
    gene_list <- res_df$gene[!is.na(res_df$padj) & res_df$padj < 0.05]
    
    if (length(gene_list) > 0) {
      kegg_result <- enrichKEGG(
        gene = gene_list,
        organism = 'mmu',
        keyType = 'kegg',
        pvalueCutoff = 0.05
      )
      
      # Create results directory for this dataset
      results_dir <- file.path("results", gse_id)
      dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
      
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
for (gse_id in pmid_df$Accession_ID) {
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