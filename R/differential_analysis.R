#!/usr/bin/env Rscript

# Load required libraries
library(edgeR)
library(limma)
library(sva)
library(tidyverse)

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
run_rnaseq_analysis <- function(counts_file, metadata_file, output_dir) {
  tryCatch({
    # Load data
    counts <- read.table(counts_file, header = TRUE, row.names = 1)
    metadata <- read.csv(metadata_file, row.names = 1)
    
    # Create DGEList object
    dge <- DGEList(counts = counts)
    
    # Filter low counts
    cpm <- cpm(dge)
    keepers <- rowSums(cpm > 1) >= min(table(metadata$condition))
    dge <- dge[keepers,]
    
    # Combat-Seq batch correction
    if ("batch" %in% colnames(metadata)) {
      counts_corrected <- ComBat_seq(counts = dge$counts, 
                                   batch = metadata$batch,
                                   group = metadata$condition)
      dge$counts <- counts_corrected
    }
    
    # GeTMM normalization
    dge <- calcNormFactors(dge, method = "TMM")
    
    # Create design matrix
    design <- model.matrix(~ condition, data = metadata)
    
    # Voom transformation
    v <- voom(dge, design, plot = TRUE)
    
    # Fit linear model
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    # Get results
    results <- topTable(fit, coef = 2, number = Inf)
    
    # Save results
    saveRDS(fit, file = file.path(output_dir, "edgeR_fit.rds"))
    saveRDS(results, file = file.path(output_dir, "edgeR_results.rds"))
    write.csv(results, file = file.path(output_dir, "edgeR_results.csv"))
    
    # Generate plots
    pdf(file.path(output_dir, "edgeR_plots.pdf"))
    plotMA(fit)
    plotSA(fit)
    dev.off()
    
    message("Successfully completed edgeR analysis")
  }, error = function(e) {
    stop("Failed to run edgeR analysis: ", e$message)
  })
}

run_methylation_analysis <- function(beta_file, metadata_file, output_dir) {
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
  run_rnaseq_analysis(input_file, metadata_file, output_dir)
} else if (data_type == "methylation") {
  run_methylation_analysis(input_file, metadata_file, output_dir)
} else {
  stop("Invalid data type. Must be 'rnaseq' or 'methylation'")
} 