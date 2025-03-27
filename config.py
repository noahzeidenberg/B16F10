#!/usr/bin/env python3

# Pipeline configuration

# System settings - derived from SLURM configuration
from b16f10_pipeline import get_available_cpus, get_available_memory

THREADS = get_available_cpus()  # Gets from SLURM_CPUS_PER_TASK or system
MEMORY = f"{get_available_memory()}G"  # Gets from SLURM_MEM_PER_CPU or defaults to 4G

# Reference genome settings
GENOME_REF = {
    "path": "/path/to/mouse_genome",
    "version": "mm10",
    "gtf": "/path/to/mouse_genome/genes.gtf",
    "star_index": "/path/to/mouse_genome/star_index"
}

# RNA-seq settings
RNASEQ_PARAMS = {
    "min_read_length": 36,
    "min_quality_score": 20, # Changed from 30 to 20
    "adapter_sequence": "TruSeq3-PE.fa",    # is this always the same?
    "star_params": {
        "outFilterMultimapNmax": 20,
        "alignSJoverhangMin": 8,
        "alignSJDBoverhangMin": 1,
        "outFilterMismatchNmax": 999,
        "outFilterMismatchNoverReadLmax": 0.04,
        "alignIntronMin": 20,
        "alignIntronMax": 1000000,
        "alignMatesGapMax": 1000000
    }
}

# Methylation array settings
METHYLATION_PARAMS = {
    "array_type": "450k",
    "normalization_method": "funnorm",
    "min_detection_p": 0.01,
    "min_beads": 3
}

# Differential expression settings
DE_PARAMS = {
    "min_count": 10,
    "min_samples": 2,
    "fdr_threshold": 0.05,
    "log2fc_threshold": 1.0
}

# KEGG pathway analysis settings
KEGG_PARAMS = {
    "organism": "mmu",
    "p_value_threshold": 0.05,
    "min_size": 10,
    "max_size": 500
}

# Batch correction settings
BATCH_PARAMS = {
    "method": "ComBat",
    "parametric": True,
    "prior_plots": True
}

# File paths
PATHS = {
    "metadata": "B16F10_RNASeq_NCBI_table.csv",
    "logs": "logs",
    "results": "results"
} 