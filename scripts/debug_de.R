# Debugging script for differential expression analysis
# This script is designed to be run line-by-line in R

# ===== Setup =====
# Set working directory to your B16F10 folder
setwd("~/scratch/B16F10")

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "pheatmap",
  "stringr",
  "data.table",
  "plotly",
  "dplyr",
  "htmlwidgets",
  "RColorBrewer"  # Added for color palettes
)

# Install and load packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ===== Parameters =====
# Set your GSE ID
gse_id <- "GSE287957"  # Change this to your GSE ID
base_dir <- path.expand("~/scratch/B16F10")
gse_dir <- file.path(base_dir, gse_id)

# Create output directory
output_dir <- file.path(gse_dir, "results", "differential_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===== Load Data =====
# Load normalized data
normalization_dir <- file.path(gse_dir, "results", "normalization")
normalized_file <- file.path(normalization_dir, "normalized_counts.rds")

# Check if file exists
if (!file.exists(normalized_file)) {
  stop("Normalized counts file not found!")
}

# Load the data
normalized_data <- readRDS(normalized_file)
print("Structure of normalized data:")
str(normalized_data)

# Extract counts
counts <- normalized_data$normalized_counts
print("Dimensions of counts matrix:")
dim(counts)
print("First few rows and columns of counts:")
head(counts[, 1:min(5, ncol(counts))])

# Check if counts are negative (log-transformed)
print("Summary of counts values:")
summary(as.vector(counts))

# Handle negative counts (log-transformed values)
if (any(counts < 0)) {
  print("Negative counts detected. Converting from log2(CPM/kb) to CPM/kb...")
  
  # First, check if we have raw counts available
  if ("raw_counts" %in% names(normalized_data)) {
    print("Using raw counts from normalized data...")
    counts <- normalized_data$raw_counts
    print("Raw counts summary:")
    print(summary(as.vector(counts)))
  } else {
    # Convert from log2(CPM/kb) to CPM/kb
    print("Converting from log2(CPM/kb) to CPM/kb...")
    cpm_kb <- 2^counts
    
    # Then, multiply by gene length (in kb) to get CPM
    if ("gene_lengths" %in% names(normalized_data)) {
      gene_lengths_kb <- normalized_data$gene_lengths$length / 1000
      
      # Ensure gene lengths match the counts matrix
      if (length(gene_lengths_kb) == nrow(counts)) {
        # Convert from CPM/kb to CPM
        counts <- sweep(cpm_kb, 1, gene_lengths_kb, "*")
        print("Converted to CPM values. Summary:")
        print(summary(as.vector(counts)))
      } else {
        stop(sprintf("Gene lengths (%d) do not match counts matrix dimensions (%d)", 
                    length(gene_lengths_kb), nrow(counts)))
      }
    } else {
      stop("Gene lengths not found in normalized data")
    }
  }
}

# ===== Load Design Matrix =====
# Function to find design matrix
find_design_matrix <- function(base_dir, gse_id) {
  possible_locations <- c(
    file.path(base_dir, "sample_design", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    file.path(base_dir, "results", "design_matrices", paste0(gse_id, "_design_matrix.txt")),
    file.path(base_dir, gse_id, "design_matrices", paste0(gse_id, "_design_matrix.txt"))
  )
  
  for (location in possible_locations) {
    if (file.exists(location)) {
      return(location)
    }
  }
  return(NULL)
}

# Find and load design matrix
design_file <- find_design_matrix(base_dir, gse_id)
if (is.null(design_file)) {
  stop("Design matrix file not found!")
}

design_matrix <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print("Structure of design matrix:")
str(design_matrix)

# ===== Sample Mapping =====
# Get sample names from counts matrix
samples <- colnames(counts)
print("Sample names:")
print(samples)

# Extract SRA IDs and GSM IDs
sra_ids <- character(length(samples))
gsm_ids <- character(length(samples))

for (i in 1:length(samples)) {
  sample_name <- samples[i]
  
  # Check if the sample name follows the new format (GSEID_GSMID_SRRID)
  if (grepl(paste0("^", gse_id, "_GSM[0-9]+_SRR"), sample_name)) {
    # Extract GSM ID and SRR ID from the new format
    parts <- strsplit(sample_name, "_")[[1]]
    gsm_ids[i] <- parts[2]  # GSM ID is the second part
    sra_ids[i] <- parts[3]  # SRR ID is the third part
  } else {
    # Fall back to the old format (GSEID_SRRID)
    sra_ids[i] <- gsub("_Aligned.sortedByCoord.out.bam", "", sample_name)
    sra_ids[i] <- gsub(paste0(gse_id, "_"), "", sra_ids[i])
    gsm_ids[i] <- NA  # We'll need to look up the GSM ID in the mapping file
  }
}

print("SRA IDs:")
print(sra_ids)
print("GSM IDs:")
print(gsm_ids)

# Try to find mapping file for missing GSM IDs
mapping_file <- file.path(base_dir, "sample_design", "sample_design", "sra_to_geo_mapping.txt")
if (file.exists(mapping_file) && any(is.na(gsm_ids))) {
  print("Found mapping file. Updating missing GSM IDs...")
  mapping <- read.table(mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (all(c("SRA_ID", "GEO_ID") %in% colnames(mapping))) {
    missing_gsm_indices <- which(is.na(gsm_ids))
    missing_sra_ids <- sra_ids[missing_gsm_indices]
    mapping <- mapping[mapping$SRA_ID %in% missing_sra_ids, ]
    
    if (nrow(mapping) > 0) {
      print(sprintf("Found %d mappings between SRA IDs and GEO IDs", nrow(mapping)))
      
      for (i in 1:nrow(mapping)) {
        sra_id <- mapping$SRA_ID[i]
        geo_id <- mapping$GEO_ID[i]
        sra_index <- which(sra_ids == sra_id)
        if (length(sra_index) > 0) {
          gsm_ids[sra_index] <- geo_id
        }
      }
    } else {
      print("No matching SRA IDs found in mapping file")
    }
  } else {
    print("Mapping file does not have required columns (SRA_ID, GEO_ID)")
  }
} else {
  if (!file.exists(mapping_file)) {
    print("Mapping file not found")
  }
  if (!any(is.na(gsm_ids))) {
    print("No missing GSM IDs to update")
  }
}

print("Updated GSM IDs:")
print(gsm_ids)

# ===== Create Group Factor =====
# Match samples to groups
sample_to_group <- data.frame(
  Sample_Index = integer(),
  Group = character(),
  stringsAsFactors = FALSE
)

for (i in 1:length(samples)) {
  gsm_id <- gsm_ids[i]
  if (!is.na(gsm_id)) {
    match_idx <- which(design_matrix$Sample_geo_accession == gsm_id)
    if (length(match_idx) > 0) {
      group <- design_matrix$Group[match_idx]
      sample_to_group <- rbind(sample_to_group, data.frame(
        Sample_Index = i,
        Group = group,
        stringsAsFactors = FALSE
      ))
    }
  }
}

print("Sample to group mapping:")
print(sample_to_group)

# Filter counts matrix
if (nrow(sample_to_group) > 0) {
  valid_indices <- sample_to_group$Sample_Index
  valid_samples <- samples[valid_indices]
  counts <- counts[, valid_samples]
  group_factor <- factor(sample_to_group$Group)
} else {
  stop("No valid samples found with matching GSM IDs in design matrix!")
}

print("Group factor:")
print(group_factor)
print("Group factor table:")
print(table(group_factor))

# ===== Create DGEList =====
# Create DGEList object
y <- DGEList(counts = counts)
print("DGEList structure:")
str(y)

# Check if we're using CPM values
if (all(y$counts >= 0) && max(y$counts) > 100) {
  print("Using CPM values. Setting library sizes to 1e6...")
  y$samples$lib.size <- rep(1e6, ncol(y))
}

# Filter low count genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
print("Number of genes kept after filtering:")
print(sum(keep))

# Create design matrix for edgeR
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
print("Design matrix:")
print(design)

# Estimate dispersion
y <- estimateDisp(y, design)
print("Dispersion estimation complete")

# Fit model
fit <- glmQLFit(y, design)
print("Model fitting complete")

# ===== Create Contrasts =====
# Get groups
groups <- colnames(design)
print("Groups for contrasts:")
print(groups)

# Find control and treatment groups
control_groups <- groups[grep("control", groups, ignore.case = TRUE)]
treatment_groups <- groups[grep("treatment", groups, ignore.case = TRUE)]

print("Control groups:")
print(control_groups)
print("Treatment groups:")
print(treatment_groups)

# Create contrasts
contrasts <- list()
if (length(control_groups) > 0 && length(treatment_groups) > 0) {
  control_group <- control_groups[1]
  for (treatment_group in treatment_groups) {
    contrast_name <- paste0(treatment_group, "_vs_", control_group)
    print(sprintf("Creating contrast: %s", contrast_name))
    contrasts[[contrast_name]] <- makeContrasts(
      contrasts = paste0(treatment_group, " - ", control_group),
      levels = design
    )
  }
} else {
  # If no control/treatment groups found, create contrasts for each group vs the first group
  control_group <- groups[1]
  for (i in 2:length(groups)) {
    contrast_name <- paste0(groups[i], "_vs_", control_group)
    print(sprintf("Creating contrast: %s", contrast_name))
    contrasts[[contrast_name]] <- makeContrasts(
      contrasts = paste0(groups[i], " - ", control_group),
      levels = design
    )
  }
}

print("Created contrasts:")
print(names(contrasts))

# ===== Test One Contrast =====
# Test the first contrast
if (length(contrasts) > 0) {
  contrast_name <- names(contrasts)[1]
  print(sprintf("Testing contrast: %s", contrast_name))
  
  # Test for differential expression
  qlf <- glmQLFTest(fit, contrast = contrasts[[contrast_name]])
  
  # Get results
  res <- topTags(qlf, n = Inf)
  print("Structure of results:")
  str(res)
  
  # Save results
  res_file <- file.path(output_dir, paste0(gse_id, "_", contrast_name, "_de_results.txt"))
  write.table(res$table, res_file, sep = "\t", quote = FALSE)
  print(sprintf("Saved results to: %s", res_file))
  
  # Create MA plot data
  print("Creating MA plot data")
  ma_data <- data.frame(
    logCPM = qlf$table$logCPM,
    logFC = qlf$table$logFC,
    FDR = qlf$table$FDR,
    genes = rownames(qlf$table)
  )
  print("MA plot data structure:")
  str(ma_data)
  
  # Create volcano plot data
  print("Creating volcano plot data")
  volcano_data <- data.frame(
    logFC = res$table$logFC,
    FDR = res$table$FDR,
    genes = rownames(res$table)
  )
  volcano_data$logP <- -log10(volcano_data$FDR)
  print("Volcano plot data structure:")
  str(volcano_data)
}

print("Debug script completed successfully!")

# ===== Load Results =====
# First, check if we have results from the previous analysis
res_file <- list.files(output_dir, pattern = "_de_results.txt$", full.names = TRUE)[1]
if (!file.exists(res_file)) {
  stop("No differential expression results found. Run differential_expression.R first!")
}

# Load the results
print("Loading differential expression results...")
res_table <- read.table(res_file, header = TRUE, sep = "\t", row.names = 1)
print("Structure of results table:")
str(res_table)

# ===== Prepare Plot Data =====
# Create MA plot data
print("Creating MA plot data...")
ma_data <- data.frame(
  logCPM = res_table$logCPM,
  logFC = res_table$logFC,
  FDR = res_table$FDR,
  genes = rownames(res_table)
)

# Add significance information for MA plot
ma_data$Significance <- dplyr::case_when(
  ma_data$FDR < 0.05 & ma_data$logFC > 1 ~ "Upregulated",
  ma_data$FDR < 0.05 & ma_data$logFC < -1 ~ "Downregulated",
  TRUE ~ "Not Significant"
)

print("MA plot data structure:")
str(ma_data)
print("Summary of MA plot data:")
summary(ma_data)
print("Table of significance:")
table(ma_data$Significance)

# Create volcano plot data
print("Creating volcano plot data...")
volcano_data <- data.frame(
  logFC = res_table$logFC,
  FDR = res_table$FDR,
  genes = rownames(res_table)
)
volcano_data$logP <- -log10(volcano_data$FDR)

# Add significance information for volcano plot
volcano_data$Significance <- dplyr::case_when(
  volcano_data$FDR < 0.05 & volcano_data$logFC > 1 ~ "Upregulated",
  volcano_data$FDR < 0.05 & volcano_data$logFC < -1 ~ "Downregulated",
  TRUE ~ "Not Significant"
)

print("Volcano plot data structure:")
str(volcano_data)
print("Summary of volcano plot data:")
summary(volcano_data)
print("Table of significance:")
table(volcano_data$Significance)

# ===== Plot Functions =====
# Function to create MA plot
create_ma_plot <- function(data, title = "MA Plot") {
  # Define colors for the points
  color_map <- c(
    "Upregulated" = "#c46666",
    "Downregulated" = "#1a80bb",
    "Not Significant" = "#b0b0b0"
  )
  
  # Create interactive MA plot
  p <- plot_ly(
    data = data,
    x = ~logCPM,
    y = ~logFC,
    text = ~paste(
      "Gene:", genes,
      "<br>logCPM:", round(logCPM, 2),
      "<br>logFC:", round(logFC, 2),
      "<br>FDR:", formatC(FDR, format = "e", digits = 2)
    ),
    hoverinfo = "text",
    mode = "markers",
    marker = list(size = 6)
  ) %>%
    add_markers(
      color = ~Significance,
      colors = color_map
    ) %>%
    layout(
      title = title,
      xaxis = list(title = "logCPM"),
      yaxis = list(title = "logFC"),
      shapes = list(
        list(type = "line", x0 = min(data$logCPM), x1 = max(data$logCPM),
             y0 = 1, y1 = 1, line = list(dash = "dash", color = "black"), opacity = 0.3),
        list(type = "line", x0 = min(data$logCPM), x1 = max(data$logCPM),
             y0 = -1, y1 = -1, line = list(dash = "dash", color = "black"), opacity = 0.3)
      )
    )
  
  return(p)
}

# Function to create volcano plot
create_volcano_plot <- function(data, title = "Volcano Plot") {
  # Define colors for the points
  color_map <- c(
    "Upregulated" = "#c46666",
    "Downregulated" = "#1a80bb",
    "Not Significant" = "#b0b0b0"
  )
  
  # Create interactive volcano plot
  p <- plot_ly(
    data = data,
    x = ~logFC,
    y = ~logP,
    text = ~paste(
      "Gene:", genes,
      "<br>logFC:", round(logFC, 2),
      "<br>FDR:", formatC(FDR, format = "e", digits = 2)
    ),
    hoverinfo = "text",
    mode = "markers",
    marker = list(size = 6)
  ) %>%
    add_markers(
      color = ~Significance,
      colors = color_map
    ) %>%
    layout(
      title = title,
      xaxis = list(title = "log2 Fold Change"),
      yaxis = list(title = "-log10(FDR)"),
      shapes = list(
        list(type = "line", x0 = -1, x1 = -1, y0 = 0, y1 = max(data$logP),
             line = list(dash = "dot"), opacity = 0.5),
        list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = max(data$logP),
             line = list(dash = "dot"), opacity = 0.5),
        list(type = "line", x0 = min(data$logFC), x1 = max(data$logFC),
             y0 = -log10(0.05), y1 = -log10(0.05),
             line = list(dash = "dash", color = "black"), opacity = 0.3)
      )
    )
  
  return(p)
}

# ===== Create and Save Plots =====
# Try creating MA plot
print("Creating MA plot...")
tryCatch({
  ma_plot <- create_ma_plot(ma_data)
  print("MA plot created successfully")
  print("Saving MA plot...")
  htmlwidgets::saveWidget(
    ma_plot,
    file.path(output_dir, paste0(gse_id, "_ma_plot.html")),
    selfcontained = TRUE
  )
  print("MA plot saved successfully")
}, error = function(e) {
  print(paste("Error creating MA plot:", e$message))
  print("MA plot data summary:")
  print(summary(ma_data))
})

# Try creating volcano plot
print("Creating volcano plot...")
tryCatch({
  volcano_plot <- create_volcano_plot(volcano_data)
  print("Volcano plot created successfully")
  print("Saving volcano plot...")
  htmlwidgets::saveWidget(
    volcano_plot,
    file.path(output_dir, paste0(gse_id, "_volcano_plot.html")),
    selfcontained = TRUE
  )
  print("Volcano plot saved successfully")
}, error = function(e) {
  print(paste("Error creating volcano plot:", e$message))
  print("Volcano plot data summary:")
  print(summary(volcano_data))
})

# ===== Create Static Plots (Fallback) =====
# Create static MA plot using ggplot2
print("Creating static MA plot...")
tryCatch({
  static_ma <- ggplot(ma_data, aes(x = logCPM, y = logFC, color = Significance)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(values = c(
      "Upregulated" = "#c46666",
      "Downregulated" = "#1a80bb",
      "Not Significant" = "#b0b0b0"
    )) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.3) +
    theme_bw() +
    labs(title = "MA Plot", x = "logCPM", y = "logFC")
  
  ggsave(
    file.path(output_dir, paste0(gse_id, "_ma_plot.pdf")),
    static_ma,
    width = 10,
    height = 8
  )
  print("Static MA plot saved successfully")
}, error = function(e) {
  print(paste("Error creating static MA plot:", e$message))
})

# Create static volcano plot using ggplot2
print("Creating static volcano plot...")
tryCatch({
  static_volcano <- ggplot(volcano_data, aes(x = logFC, y = logP, color = Significance)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(values = c(
      "Upregulated" = "#c46666",
      "Downregulated" = "#1a80bb",
      "Not Significant" = "#b0b0b0"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.3) +
    theme_bw() +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(FDR)")
  
  ggsave(
    file.path(output_dir, paste0(gse_id, "_volcano_plot.pdf")),
    static_volcano,
    width = 10,
    height = 8
  )
  print("Static volcano plot saved successfully")
}, error = function(e) {
  print(paste("Error creating static volcano plot:", e$message))
})

print("Debug script completed successfully!") 