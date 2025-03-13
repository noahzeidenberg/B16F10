# Methods

## Data Collection and Processing
RNA sequencing and methylation array data from B16F10 melanoma experiments were obtained from the Gene Expression Omnibus (GEO) database. A custom Python/R pipeline was developed to process and analyze these datasets in a standardized manner. The pipeline automatically detects data types and applies appropriate processing steps for each dataset type.

## RNA Sequencing Analysis
For RNA-seq datasets, raw reads were initially assessed for quality using FastQC (v0.11.9). Adapter sequences and low-quality bases were trimmed using Trimmomatic (PE mode) with parameters: LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36. Clean reads were aligned to the mouse reference genome (mm39) using STAR aligner with default parameters for splice-aware alignment. Gene-level quantification was performed using featureCounts (Subread package), counting only primary alignments and properly paired reads.

## Methylation Array Analysis
Methylation array data were processed using the minfi R package (v1.40.0). Raw IDAT files were imported and preprocessed using functional normalization (preprocessFunnorm) to account for technical variation. Quality control metrics were assessed using built-in minfi functions, and beta values were extracted for downstream analysis.

## Differential Analysis
Differential expression analysis for RNA-seq data was performed using DESeq2 (v1.34.0), incorporating experimental design factors and batch information where available. For methylation data, differential methylation analysis was conducted using limma (v3.50.0) with empirical Bayes moderation. In both cases, multiple testing correction was performed using the Benjamini-Hochberg method, with significance threshold set at FDR < 0.05.

## Pathway Analysis
Significantly differentially expressed genes or methylated regions were analyzed for pathway enrichment using the KEGG database through the gage R package. Pathway visualizations were generated using the pathview package, highlighting the directionality and magnitude of changes within significant pathways.

## Batch Effect Correction
To account for technical variation across multiple datasets, batch correction was performed separately for RNA-seq and methylation data using ComBat (sva package). The effectiveness of batch correction was assessed through principal component analysis before and after correction. For RNA-seq data, counts were normalized and log-transformed prior to batch correction.

## Implementation
The entire analysis pipeline was implemented in Python (v3.8+) and R (v4.1+), with parallel processing capabilities for handling multiple datasets simultaneously. All analysis steps were logged for reproducibility, and intermediate results were stored in a structured directory hierarchy. The pipeline source code is available at [repository URL]. 