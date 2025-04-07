# Differential Expression Analysis Script
# This script performs differential expression analysis using edgeR

# Load required packages
cat("=== Starting Differential Expression Analysis ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "pheatmap",
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  cat(sprintf("Loading package: %s\n", pkg))
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Function to perform differential expression analysis
perform_de_analysis <- function(gse_id) {
  cat(sprintf("Performing differential expression analysis for %s...\n", gse_id))
  
  # Set up directories
  base_dir <- getwd()
  batch_correction_dir <- file.path(base_dir, "results", "batch_correction")
  design_dir <- file.path(base_dir, "results", "design_matrices")
  output_dir <- file.path(base_dir, "results", "differential_expression")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if analysis has already been performed
  results_file <- file.path(output_dir, paste0(gse_id, "_de_results.rds"))
  if (file.exists(results_file)) {
    cat(sprintf("Differential expression analysis already performed for %s. Skipping...\n", gse_id))
    return(TRUE)
  }
  
  # Load batch-corrected counts
  batch_corrected_file <- file.path(batch_correction_dir, "batch_corrected_counts.rds")
  if (!file.exists(batch_corrected_file)) {
    cat(sprintf("Batch-corrected counts not found. Cannot proceed.\n"))
    return(FALSE)
  }
  
  cat("Loading batch-corrected counts...\n")
  batch_corrected_data <- readRDS(batch_corrected_file)
  
  # Extract counts and batch info
  counts <- batch_corrected_data$corrected_counts
  batch_info <- batch_corrected_data$batch_info
  
  # Filter counts for the specific GSE ID
  gse_samples <- batch_info$sample[grep(paste0("^", gse_id, "_"), batch_info$sample)]
  if (length(gse_samples) == 0) {
    cat(sprintf("No samples found for %s in batch-corrected data. Cannot proceed.\n", gse_id))
    return(FALSE)
  }
  
  # Extract counts for this GSE
  gse_counts <- counts[, gse_samples, drop = FALSE]
  
  # Load design matrix
  design_file <- file.path(design_dir, paste0(gse_id, "_design_matrix.txt"))
  if (!file.exists(design_file)) {
    cat(sprintf("Design matrix not found for %s. Cannot proceed.\n", gse_id))
    return(FALSE)
  }
  
  cat("Loading design matrix...\n")
  design_matrix <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Check if we have the necessary data
  if (is.null(gse_counts) || nrow(gse_counts) == 0) {
    cat("No counts found. Cannot proceed.\n")
    return(FALSE)
  }
  
  # Create DGEList object
  cat("Creating DGEList object...\n")
  y <- DGEList(counts = gse_counts)
  
  # Filter low count genes
  cat("Filtering low count genes...\n")
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("Kept %d genes out of %d\n", sum(keep), length(keep)))
  
  # Create design matrix for edgeR
  cat("Creating design matrix for edgeR...\n")
  
  # Check if we have group information
  if ("Group" %in% colnames(design_matrix)) {
    # Create factor from group information
    group <- factor(design_matrix$Group)
    
    # Create design matrix
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    
    # Estimate dispersion
    cat("Estimating dispersion...\n")
    y <- estimateDisp(y, design)
    
    # Plot dispersion
    cat("Plotting dispersion...\n")
    pdf(file.path(output_dir, paste0(gse_id, "_dispersion_plot.pdf")))
    plotBCV(y)
    dev.off()
    
    # Fit model
    cat("Fitting model...\n")
    fit <- glmQLFit(y, design)
    
    # Create contrasts for all pairwise comparisons
    cat("Creating contrasts...\n")
    groups <- levels(group)
    contrasts <- list()
    
    if (length(groups) > 1) {
      for (i in 1:(length(groups) - 1)) {
        for (j in (i + 1):length(groups)) {
          contrast_name <- paste0(groups[j], "_vs_", groups[i])
          contrasts[[contrast_name]] <- makeContrasts(
            paste0(group, groups[j], " - ", group, groups[i]),
            levels = design
          )
        }
      }
    }
    
    # Perform differential expression analysis for each contrast
    results <- list()
    
    for (contrast_name in names(contrasts)) {
      cat(sprintf("Testing contrast: %s\n", contrast_name))
      
      # Test for differential expression
      qlf <- glmQLFTest(fit, contrast = contrasts[[contrast_name]])
      
      # Get results
      res <- topTags(qlf, n = Inf)
      
      # Add to results list
      results[[contrast_name]] <- res$table
      
      # Save results
      res_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_de_results.txt"))
      write.table(res$table, res_file, sep = "\t", quote = FALSE)
      cat(sprintf("Saved results to: %s\n", res_file))
      
      # Create MA plot
      cat(sprintf("Creating MA plot for %s...\n", contrast_name))
      pdf(file.path(output_dir, paste0(gse_id, "_", contrast_name, "_ma_plot.pdf")))
      plotMD(qlf)
      abline(h = c(-1, 1), col = "blue")
      dev.off()
      
      # Create volcano plot
      cat(sprintf("Creating volcano plot for %s...\n", contrast_name))
      pdf(file.path(output_dir, paste0(gse_id, "_", contrast_name, "_volcano_plot.pdf")))
      with(res$table, plot(logFC, -log10(FDR), pch = 20, main = paste("Volcano Plot:", contrast_name),
                          xlim = c(-5, 5), ylim = c(0, 10)))
      with(subset(res$table, FDR < 0.05 & abs(logFC) > 1), points(logFC, -log10(FDR), pch = 20, col = "red"))
      abline(h = -log10(0.05), col = "blue", lty = 2)
      abline(v = c(-1, 1), col = "blue", lty = 2)
      dev.off()
    }
    
    # Save all results
    saveRDS(results, results_file)
    cat(sprintf("Saved all results to: %s\n", results_file))
    
    # Create heatmap of differentially expressed genes
    if (length(results) > 0) {
      cat("Creating heatmap of differentially expressed genes...\n")
      
      # Get differentially expressed genes (FDR < 0.05 and |logFC| > 1)
      de_genes <- list()
      for (contrast_name in names(results)) {
        de_genes[[contrast_name]] <- rownames(results[[contrast_name]])[
          results[[contrast_name]]$FDR < 0.05 & abs(results[[contrast_name]]$logFC) > 1
        ]
      }
      
      # Get unique differentially expressed genes
      all_de_genes <- unique(unlist(de_genes))
      
      if (length(all_de_genes) > 0) {
        # Extract expression data for differentially expressed genes
        de_expr <- y$counts[all_de_genes, ]
        
        # Scale for better visualization
        de_expr_scaled <- t(scale(t(de_expr)))
        
        # Create heatmap
        pdf(file.path(output_dir, paste0(gse_id, "_de_heatmap.pdf")))
        pheatmap(de_expr_scaled, 
                 main = "Differentially Expressed Genes",
                 scale = "none",
                 clustering_method = "ward.D2",
                 show_rownames = FALSE)
        dev.off()
        cat(sprintf("Created heatmap with %d differentially expressed genes\n", length(all_de_genes)))
      } else {
        cat("No differentially expressed genes found for any contrast\n")
      }
    }
    
    return(TRUE)
  } else {
    cat("No group information found in design matrix. Cannot proceed with differential expression analysis.\n")
    return(FALSE)
  }
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript differential_expression.R GSE223515")
    }
  }
  
  tryCatch({
    cat(sprintf("Performing differential expression analysis for %s...\n", gse_id))
    success <- perform_de_analysis(gse_id)
    
    if (success) {
      cat(sprintf("Successfully performed differential expression analysis for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to perform differential expression analysis for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during differential expression analysis: %s\n", e$message))
    stop(sprintf("Error during differential expression analysis: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 