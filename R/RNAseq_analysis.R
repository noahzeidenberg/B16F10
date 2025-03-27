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
pmid_df <- read.csv("PMID_df.csv")

# Function to check if a GSE is RNA-seq
is_rnaseq <- function(gse_id) {
  gse <- getGEO(gse_id)
  # Check if it's RNA-seq by looking at the technology type
  tech <- unique(unlist(lapply(gse, function(x) x@experimentData@other$type)))
  return(any(grepl("RNA-seq|RNA sequencing", tech, ignore.case = TRUE)))
}

# Function to process a single RNA-seq dataset
process_rnaseq <- function(gse_id) {
  message(paste("Processing", gse_id))
  
  # Download the data
  gse <- getGEO(gse_id)
  
  # Extract expression matrix and phenotype data
  expr_matrix <- exprs(gse[[1]])
  pheno_data <- pData(gse[[1]])
  
  # Create DESeq2 object
  # Note: This is a simplified example - you'll need to adjust based on your specific data structure
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
    
    # Save results
    write.csv(res_df, file = paste0("results/", gse_id, "_DESeq2_results.csv"))
    write.csv(as.data.frame(kegg_result), file = paste0("results/", gse_id, "_KEGG_results.csv"))
  }
  
  return(list(deseq_results = res_df, kegg_results = kegg_result))
}

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Process each dataset
results <- list()
for (gse_id in pmid_df$Accession_ID) {
  tryCatch({
    if (is_rnaseq(gse_id)) {
      results[[gse_id]] <- process_rnaseq(gse_id)
    }
  }, error = function(e) {
    message(paste("Error processing", gse_id, ":", e$message))
  })
}

# Save session info
sink("results/session_info.txt")
sessionInfo()
sink() 