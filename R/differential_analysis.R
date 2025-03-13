#!/usr/bin/env Rscript

# Load required libraries
library(DESeq2)
library(limma)
library(edgeR)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript differential_analysis.R <input_file> <metadata_file> <data_type> <output_dir>")
}

input_file <- args[1]
metadata_file <- args[2]
data_type <- args[3]
output_dir <- args[4]

# Differential analysis functions
run_deseq2_analysis <- function(counts_file, metadata_file, output_dir) {
  tryCatch({
    # Load data
    counts <- read.table(counts_file, header = TRUE, row.names = 1)
    metadata <- read.csv(metadata_file, row.names = 1)
    
    # Create DESeq dataset
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = ~ condition
    )
    
    # Filter low counts
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    # Run DESeq
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds)
    
    # Save results
    saveRDS(dds, file = file.path(output_dir, "deseq2_object.rds"))
    saveRDS(res, file = file.path(output_dir, "deseq2_results.rds"))
    write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))
    
    # Generate plots
    pdf(file.path(output_dir, "deseq2_plots.pdf"))
    plotMA(res)
    plotDispEsts(dds)
    dev.off()
    
    message("Successfully completed DESeq2 analysis")
  }, error = function(e) {
    stop("Failed to run DESeq2 analysis: ", e$message)
  })
}

run_limma_analysis <- function(beta_file, metadata_file, output_dir) {
  tryCatch({
    # Load data
    beta <- readRDS(beta_file)
    metadata <- read.csv(metadata_file)
    
    # Create design matrix
    design <- model.matrix(~ condition, data = metadata)
    
    # Fit linear model
    fit <- lmFit(beta, design)
    fit <- eBayes(fit)
    
    # Get results
    results <- topTable(fit, coef = 2, number = Inf)
    
    # Save results
    saveRDS(fit, file = file.path(output_dir, "limma_fit.rds"))
    saveRDS(results, file = file.path(output_dir, "limma_results.rds"))
    write.csv(results, file = file.path(output_dir, "limma_results.csv"))
    
    # Generate plots
    pdf(file.path(output_dir, "limma_plots.pdf"))
    volcanoplot(fit, coef = 2)
    dev.off()
    
    message("Successfully completed limma analysis")
  }, error = function(e) {
    stop("Failed to run limma analysis: ", e$message)
  })
}

# Run appropriate analysis based on data type
if (data_type == "rnaseq") {
  run_deseq2_analysis(input_file, metadata_file, output_dir)
} else if (data_type == "methylation") {
  run_limma_analysis(input_file, metadata_file, output_dir)
} else {
  stop("Invalid data type. Must be 'rnaseq' or 'methylation'")
} 