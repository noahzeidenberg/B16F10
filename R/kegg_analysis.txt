#!/usr/bin/env Rscript

# Load required libraries
library(gage)
library(pathview)
library(org.Mm.eg.db)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript kegg_analysis.R <results_file> <data_type> <output_dir>")
}

results_file <- args[1]
data_type <- args[2]
output_dir <- args[3]

# KEGG pathway analysis
run_kegg_analysis <- function(results_file, data_type, output_dir) {
  tryCatch({
    # Load results
    results <- readRDS(results_file)
    
    # Prepare data for KEGG analysis
    if (data_type == "rnaseq") {
      fc <- results$log2FoldChange
      names(fc) <- rownames(results)
    } else {
      fc <- results$logFC
      names(fc) <- rownames(results)
    }
    
    # Convert gene symbols to Entrez IDs
    symbols <- names(fc)
    entrez_ids <- mapIds(org.Mm.eg.db, symbols, "ENTREZID", "SYMBOL")
    fc <- fc[!is.na(entrez_ids)]
    names(fc) <- entrez_ids[!is.na(entrez_ids)]
    
    # Load KEGG pathways
    data(kegg.gs)
    
    # Run KEGG pathway analysis
    kegg_res <- gage(fc, gsets = kegg.gs, ref = NULL, samp = NULL)
    
    # Generate pathway visualizations for top pathways
    top_pathways <- rownames(kegg_res$greater)[1:5]
    for (pathway in top_pathways) {
      pathview(
        gene.data = fc,
        pathway.id = substr(pathway, 1, 8),
        species = "mmu",
        out.suffix = "pathview",
        kegg.dir = output_dir
      )
    }
    
    # Save results
    saveRDS(kegg_res, file = file.path(output_dir, "kegg_results.rds"))
    write.csv(kegg_res$greater, file = file.path(output_dir, "kegg_upregulated.csv"))
    write.csv(kegg_res$less, file = file.path(output_dir, "kegg_downregulated.csv"))
    
    message("Successfully completed KEGG pathway analysis")
  }, error = function(e) {
    stop("Failed to run KEGG pathway analysis: ", e$message)
  })
}

# Execute analysis
run_kegg_analysis(results_file, data_type, output_dir) 