# Differential Gene Expression Analysis Pipeline for GSE223515
# This script performs a complete differential expression analysis workflow

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "DESeq2",
  "tidyverse",
  "pheatmap",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Mm.eg.db",
  "ggplot2",
  "dplyr",
  "tidyr",
  "Rsubread",
  "tximport",
  "GenomicFeatures",
  "AnnotationHub",
  "GEOquery",
  "rentrez",
  "xml2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Function to create directory structure for a GSE ID
create_gse_structure <- function(gse_id) {
  # Create main GSE directory
  gse_dir <- file.path(".", gse_id)
  dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories
  dirs <- c(
    file.path(gse_dir, "samples"),
    file.path(gse_dir, "logs"),
    file.path(gse_dir, "results")
  )
  
  for (dir in dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  return(gse_dir)
}

# Function to create sample directory structure
create_sample_structure <- function(gse_dir, gsm_id) {
  # Create sample directory
  sample_dir <- file.path(gse_dir, "samples", gsm_id)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories
  dirs <- c(
    file.path(sample_dir, "SRA"),
    file.path(sample_dir, "SRA", "FASTQ"),
    file.path(sample_dir, "alignment"),
    file.path(sample_dir, "quantification")
  )
  
  for (dir in dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  return(sample_dir)
}

# Function to safely save RDS file
safe_save_rds <- function(object, file) {
  if (file.exists(file)) {
    backup_file <- paste0(file, ".backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    cat(sprintf("File %s already exists. Creating backup at %s\n", file, backup_file))
    file.copy(file, backup_file)
  }
  saveRDS(object, file)
}

# Function to safely write table
safe_write_table <- function(x, file) {
  if (file.exists(file)) {
    backup_file <- paste0(file, ".backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    cat(sprintf("File %s already exists. Creating backup at %s\n", file, backup_file))
    file.copy(file, backup_file)
  }
  write.table(x, file = file, sep = "\t", quote = FALSE)
}

# Function to get expression matrix and sample information from GSE ID
get_expression_matrix <- function(gse_id) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Fetch the GEO object
  cat(sprintf("Fetching GEO object for %s...\n", gse_id))
  geo_object <- getGEO(gse_id, GSEMatrix = TRUE)
  
  # Save GEO object
  safe_save_rds(geo_object, file.path(gse_dir, paste0(gse_id, "_geo_object.rds")))
  
  # Extract expression matrix
  expr_matrix <- exprs(geo_object[[1]])
  
  # Extract phenotype data
  pheno_data <- pData(phenoData(geo_object[[1]]))
  
  return(list(expr_matrix = expr_matrix, pheno_data = pheno_data))
}

# Function to get SRX IDs from GSM IDs and save GSM objects
get_srx_ids <- function(gse_id) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Get the GEO object without matrix
  cat(sprintf("Fetching GSM information for %s...\n", gse_id))
  gse <- getGEO(gse_id, GSEMatrix = FALSE)
  
  # Check if this is HTS data. Currently I'm only using HTS data but could add array data later...
  data_type <- Meta(gse)$type
  if (is.null(data_type) || !grepl("high throughput sequencing", data_type, ignore.case = TRUE)) {
    cat(sprintf("Dataset %s is not High-Throughput Sequencing data. Found type: %s. Skipping...\n", 
                 gse_id, ifelse(is.null(data_type), "unknown", data_type)))
    return(NULL)
  }
  
  # Get list of GSM objects
  gsm_list <- GSMList(gse)
  
  # Save each GSM object and extract SRX IDs
  srx_ids <- sapply(names(gsm_list), function(gsm_name) {
    gsm <- gsm_list[[gsm_name]]
    
    # Create sample directory structure
    sample_dir <- create_sample_structure(gse_dir, gsm_name)
    
    # Save GSM object
    safe_save_rds(gsm, file.path(sample_dir, paste0(gsm_name, ".rds")))
    
    # Extract SRX ID
    relations <- Meta(gsm)$relation
    sra_link <- relations[grep("SRA:", relations)]
    if (length(sra_link) > 0) {
      return(sub("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "", sra_link))
    }
    return(NA)
  })
  
  # Remove any NA values
  srx_ids <- srx_ids[!is.na(srx_ids)]
  
  if (length(srx_ids) == 0) {
    cat(sprintf("No SRX IDs found for GSE ID: %s. Skipping...\n", gse_id))
    return(NULL)
  }
  
  return(srx_ids)
}

# Function to convert SRX IDs to SRA IDs using rentrez
convert_srx_to_sra <- function(srx_ids) {
  # Create output directory if it doesn't exist
  dir.create("sra_files", showWarnings = FALSE, recursive = TRUE)
  
  # Initialize vector to store SRA IDs
  sra_ids <- character(length(srx_ids))
  
  # Process each SRX ID
  for (i in seq_along(srx_ids)) {
    srx_id <- srx_ids[i]
    cat(sprintf("Processing SRX ID %d of %d: %s\n", i, length(srx_ids), srx_id))
    
    tryCatch({
      # Use rentrez to get SRR ID
      xml_result <- rentrez::entrez_fetch(db = "sra", id = srx_id, rettype = "xml")
      xml_doc <- xml2::read_xml(xml_result)
      srr_id <- xml2::xml_find_all(xml_doc, ".//RUN") |> xml2::xml_attr("accession")
      
      if (length(srr_id) > 0 && grepl("^SRR", srr_id)) {
        sra_ids[i] <- srr_id
        cat(sprintf("Found SRA ID: %s\n", srr_id))
      } else {
        cat(sprintf("No valid SRA ID found for SRX ID: %s\n", srx_id))
        sra_ids[i] <- NA
      }
    }, error = function(e) {
      cat(sprintf("Error processing SRX ID %s: %s\n", srx_id, e$message))
      sra_ids[i] <- NA
    })
    
    # Add a small delay to avoid hitting rate limits
    Sys.sleep(0.5)
  }
  
  # Remove any NA values
  sra_ids <- sra_ids[!is.na(sra_ids)]
  
  if (length(sra_ids) == 0) {
    stop("No SRA IDs found for any of the SRX IDs")
  }
  
  return(sra_ids)
}

# Function to download SRA files using prefetch
download_sra_files <- function(gse_id, sra_ids) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Download SRA files for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_dir <- file.path(gsm_dir, "SRA")
      cmd <- sprintf('wsl prefetch --max-size 500G -O %s -X 100G -f yes -p %s', 
                     sra_dir, sra_id)
      system(cmd, intern = TRUE)
    }
  }
}

# Function to convert SRA to FASTQ
convert_sra_to_fastq <- function(gse_id, sra_ids, threads = 8) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Convert SRA to FASTQ for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      sra_file <- file.path(gsm_dir, "SRA", paste0(sra_id, ".sra"))
      fastq_dir <- file.path(gsm_dir, "SRA", "FASTQ")
      
      if (file.exists(sra_file)) {
        cmd <- sprintf('fasterq-dump --split-files --progress --threads %d --outdir %s %s',
                       threads, fastq_dir, sra_id)
        system(cmd, intern = TRUE)
      }
    }
  }
}

# Function to align reads
align_reads <- function(gse_id, sra_ids, reference_index) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Align reads for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      fastq_dir <- file.path(gsm_dir, "SRA", "FASTQ")
      alignment_dir <- file.path(gsm_dir, "alignment")
      
      # Align reads using Rsubread
      align(index = reference_index,
            readfile1 = file.path(fastq_dir, paste0(sra_id, "_1.fastq.gz")),
            readfile2 = file.path(fastq_dir, paste0(sra_id, "_2.fastq.gz")),
            output_file = file.path(alignment_dir, paste0(sra_id, ".bam")),
            nthreads = 4)
    }
  }
}

# Function to quantify reads
quantify_reads <- function(gse_id, sra_ids, gtf_file) {
  gse_dir <- create_gse_structure(gse_id)
  
  # Get all BAM files
  bam_files <- list.files(file.path(gse_dir, "samples"), 
                         pattern = "*.bam", 
                         recursive = TRUE, 
                         full.names = TRUE)
  
  # Get gene-level counts
  counts <- featureCounts(files = bam_files,
                         annot.ext = gtf_file,
                         isGTFAnnotationFile = TRUE,
                         isPairedEnd = TRUE,
                         nthreads = 4)
  
  # Save counts for each sample
  for (sra_id in sra_ids) {
    # Find the corresponding GSM directory
    gsm_dirs <- list.dirs(file.path(gse_dir, "samples"), recursive = FALSE)
    for (gsm_dir in gsm_dirs) {
      quant_dir <- file.path(gsm_dir, "quantification")
      safe_write_table(counts$counts[, sra_id], 
                      file.path(quant_dir, paste0(sra_id, "_counts.txt")))
    }
  }
  
  return(counts$counts)
}

# Function to perform quality control
perform_qc <- function(expr_matrix, pheno_data) {
  # Create PCA plot
  tryCatch({
    # Scale the data
    scaled_data <- scale(t(expr_matrix))
    
    # Perform PCA
    pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    
    # Calculate variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    
    # Create PCA data frame
    pca_data <- data.frame(
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      Group = pheno_data$group
    )
    
    # Create PCA plot
    pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
      geom_point() +
      theme_minimal() +
      labs(title = "PCA Plot of Samples",
           x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"))
    
    ggsave("plots/pca_plot.pdf", pca_plot)
    
    # Create sample correlation heatmap
    cor_matrix <- cor(expr_matrix)
    pheatmap(cor_matrix, 
             filename = "plots/sample_correlation_heatmap.pdf",
             show_rownames = TRUE,
             show_colnames = TRUE)
  }, error = function(e) {
    cat("Error in PCA analysis:", e$message, "\n")
    cat("Attempting to continue with the rest of the analysis...\n")
  })
}

# Function to perform differential expression analysis
perform_differential_expression <- function(expr_matrix, pheno_data) {
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = round(expr_matrix),
    colData = pheno_data,
    design = ~ group
  )
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds)
  
  # Save results
  safe_save_rds(res, "differential_expression_results.rds")
  
  return(res)
}

# Function to create visualization plots
create_visualizations <- function(res) {
  # Volcano plot
  volcano_plot <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 3.0,
    labSize = 6.0)
  
  ggsave("plots/volcano_plot.pdf", volcano_plot)
  
  # MA plot
  plotMA(res, ylim = c(-5, 5))
  dev.copy(pdf, "plots/ma_plot.pdf")
  dev.off()
}

# Function to perform gene set enrichment analysis
perform_gsea <- function(res) {
  # Convert gene IDs to ENTREZ IDs
  gene_list <- res$log2FoldChange
  names(gene_list) <- rownames(res)
  
  # Perform GSEA
  gsea_result <- gseGO(geneList = gene_list,
                       ont = "BP",
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENSEMBL",
                       pvalueCutoff = 0.05)
  
  # Save results
  safe_save_rds(gsea_result, "gsea_results.rds")
  
  # Create GSEA plot
  gsea_plot <- dotplot(gsea_result, showCategory = 20)
  ggsave("plots/gsea_plot.pdf", gsea_plot)
}

# Function to prepare data for differential expression analysis
prepare_data <- function(gse_id = "GSE223515", reference_index = "mm10") {
  cat("Starting data preparation...\n")
  
  # Get SRX IDs from GSM IDs
  srx_ids <- get_srx_ids(gse_id)
  if (is.null(srx_ids)) {
    stop("Failed to get SRX IDs")
  }
  
  # Convert SRX to SRA IDs
  sra_ids <- convert_srx_to_sra(srx_ids)
  if (is.null(sra_ids)) {
    stop("Failed to convert SRX to SRA IDs")
  }
  
  # Download SRA files
  cat("Downloading SRA files...\n")
  download_sra_files(gse_id, sra_ids)
  
  # Convert SRA to FASTQ
  cat("Converting SRA to FASTQ...\n")
  convert_sra_to_fastq(gse_id, sra_ids)
  
  # Align reads
  cat("Aligning reads...\n")
  align_reads(gse_id, sra_ids, reference_index)
  
  # Quantify genes
  cat("Quantifying genes...\n")
  quantify_reads(gse_id, sra_ids, "mm10")
  
  # Get expression matrix and sample information
  cat("Getting expression matrix and sample information...\n")
  data <- get_expression_matrix(gse_id)
  
  return(data)
}

# Function to process multiple GSE IDs
process_gse_ids <- function(gse_ids) {
  # Create main directories
  dir.create("scripts", showWarnings = FALSE)
  dir.create("logs", showWarnings = FALSE)
  dir.create("results", showWarnings = FALSE)
  
  # Initialize list to store results
  results <- list()
  
  # Process each GSE ID
  for (gse_id in gse_ids) {
    cat(sprintf("\nProcessing GSE ID: %s\n", gse_id))
    
    # Get SRX IDs
    srx_ids <- get_srx_ids(gse_id)
    
    # If we got valid SRX IDs, convert to SRA IDs
    if (!is.null(srx_ids)) {
      tryCatch({
        sra_ids <- convert_srx_to_sra(srx_ids)
        results[[gse_id]] <- list(
          srx_ids = srx_ids,
          sra_ids = sra_ids
        )
        cat(sprintf("Successfully processed GSE ID: %s\n", gse_id))
      }, error = function(e) {
        cat(sprintf("Error processing GSE ID %s: %s\n", gse_id, e$message))
      })
    }
  }
  
  return(results)
}

# Main workflow
main <- function() {
  # Create necessary directories
  dir.create("plots", showWarnings = FALSE)
  
  # Prepare data
  cat("Preparing data...\n")
  data <- prepare_data()
  
  # Perform QC
  cat("Performing quality control...\n")
  perform_qc(data$expr_matrix, data$pheno_data)
  
  # Perform differential expression analysis
  cat("Performing differential expression analysis...\n")
  res <- perform_differential_expression(data$expr_matrix, data$pheno_data)
  
  # Create visualizations
  cat("Creating visualizations...\n")
  create_visualizations(res)
  
  # Perform GSEA
  cat("Performing gene set enrichment analysis...\n")
  perform_gsea(res)
  
  cat("Analysis complete! Check the plots directory for results.\n")
}

# Run the pipeline
main()
