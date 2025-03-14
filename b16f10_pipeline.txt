#!/usr/bin/env python3

import os
import sys
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import argparse
import json
from typing import Dict, List, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from functools import partial

# Import configuration
from config import *

# Create logs and checkpoint directories if they don't exist
Path("logs").mkdir(parents=True, exist_ok=True)
Path("checkpoints").mkdir(parents=True, exist_ok=True)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/pipeline.log'),
        logging.StreamHandler()
    ]
)

class PipelineCheckpoint:
    """Class to handle pipeline checkpointing"""
    def __init__(self, accession: str):
        self.accession = accession
        self.checkpoint_file = Path(f"checkpoints/{accession}_status.json")
        self.status = self._load_status()

    def _load_status(self) -> Dict:
        """Load checkpoint status from file"""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file) as f:
                return json.load(f)
        return {
            "download": False,
            "data_processing": False,
            "differential_analysis": False,
            "kegg_analysis": False,
            "completed": False
        }

    def save_status(self):
        """Save current status to checkpoint file"""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.status, f)

    def is_step_completed(self, step: str) -> bool:
        """Check if a step has been completed"""
        return self.status.get(step, False)

    def mark_step_completed(self, step: str):
        """Mark a step as completed"""
        self.status[step] = True
        self.save_status()

def get_available_cpus():
    """Get the number of available CPUs from SLURM or system"""
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        return int(os.environ['SLURM_CPUS_PER_TASK'])
    return multiprocessing.cpu_count()

def get_available_memory():
    """Get available memory from SLURM or system"""
    if 'SLURM_MEM_PER_CPU' in os.environ:
        return int(os.environ['SLURM_MEM_PER_CPU'])
    return 4  # Default to 4GB if not specified

def create_directory_structure(accession: str):
    """Create the directory structure for a given accession"""
    dirs = [
        f"{accession}",
        f"{accession}/samples",
        f"{accession}/alignment",
        f"{accession}/quantification",
        f"{accession}/results",
        "logs"
    ]
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
    logging.info(f"Created directory structure for {accession}")

def download_geo_data(accession: str, output_dir: str):
    """Download data from GEO using R script"""
    try:
        # Load required modules
        required_modules = ['r', 'sra-toolkit']
        for module in required_modules:
            if not load_module(module):
                logging.error(f"Required module {module} not found. Cannot proceed with download.")
                return False
            
        cmd = f"Rscript R/geo_download.R {accession} {output_dir}"
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f"Downloaded GEO data for {accession}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to download GEO data: {str(e)}")
        return False

def load_module(module_name: str) -> bool:
    """Load a module using the module system"""
    try:
        subprocess.run(f"module load {module_name}", shell=True, check=True)
        logging.info(f"Successfully loaded module: {module_name}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to load module {module_name}: {str(e)}")
        return False

def process_methylation(accession_dir: str, output_dir: str, max_workers: int = None):
    """Process methylation array data using minfi"""
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Run R script for methylation analysis
        cmd = f"Rscript R/process_methylation.R {accession_dir}/samples {output_dir}"
        logging.info(f"Running methylation analysis: {cmd}")
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        
        # Check for output files
        beta_file = os.path.join(output_dir, "beta_values.csv")
        if not os.path.exists(beta_file):
            raise FileNotFoundError("No beta values file was generated")
            
        logging.info("Successfully processed methylation data")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process methylation array data: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Failed to process methylation data: {str(e)}")
        return False

def process_rnaseq(accession_dir: str, output_dir: str, max_workers: int = None):
    """Process RNA-seq data using STAR and featureCounts"""
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Run R script for RNA-seq analysis
        cmd = f"Rscript R/process_rnaseq.R {accession_dir}/samples {output_dir}"
        logging.info(f"Running RNA-seq analysis: {cmd}")
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        
        # Check for output files
        counts_file = os.path.join(output_dir, "counts_matrix.csv")
        if not os.path.exists(counts_file):
            raise FileNotFoundError("No counts matrix was generated")
            
        logging.info("Successfully processed RNA-seq data")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process RNA-seq data: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Failed to process RNA-seq data: {str(e)}")
        return False

def process_methylation_array(accession: str):
    """Process methylation array data using R script"""
    try:
        cmd = f"Rscript R/process_methylation.R {accession}/samples {accession}/quantification"
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f"Completed methylation array processing for {accession}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process methylation array data: {str(e)}")
        return False

def differential_analysis(accession: str, data_type: str):
    """Perform differential analysis using R script"""
    try:
        if not load_module('r'):
            return False
            
        input_file = f"{accession}/quantification/{'counts.txt' if data_type == 'rnaseq' else 'beta_values.rds'}"
        metadata_file = f"{accession}/{accession}_metadata.csv"
        output_dir = f"{accession}/results"
        
        cmd = f"Rscript R/differential_analysis.R {input_file} {metadata_file} {data_type} {output_dir}"
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f"Completed differential analysis for {accession}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to perform differential analysis: {str(e)}")
        return False

def kegg_pathway_analysis(accession: str, data_type: str):
    """Perform KEGG pathway analysis using R script"""
    try:
        if not load_module('r'):
            return False
            
        results_file = f"{accession}/results/{'deseq2_results.rds' if data_type == 'rnaseq' else 'limma_results.rds'}"
        
        # Check if results file exists
        if not Path(results_file).exists():
            logging.error(f"Results file {results_file} not found. Ensure differential analysis completed successfully.")
            return False
            
        output_dir = f"{accession}/results/kegg"
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Check if R script exists
        r_script = "R/kegg_analysis.R"
        if not Path(r_script).exists():
            logging.error(f"KEGG analysis R script not found at {r_script}")
            return False
        
        cmd = f"Rscript {r_script} {results_file} {data_type} {output_dir}"
        subprocess.run(cmd, shell=True, check=True)
        logging.info(f"Completed KEGG pathway analysis for {accession}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to perform KEGG pathway analysis: {str(e)}")
        return False

def verify_batch_correction(data_type: str) -> bool:
    """Verify that batch correction completed successfully"""
    try:
        output_dir = Path("results/batch_correction")
        if data_type == "rnaseq":
            return (output_dir / "batch_corrected_counts.csv").exists()
        else:
            return (output_dir / "batch_corrected_beta_values.csv").exists()
    except Exception:
        return False

def batch_correction(accessions: List[str], data_type: str):
    """Perform batch correction using R script"""
    try:
        if not load_module('r'):
            return False
            
        # Create output directory
        output_dir = Path("results/batch_correction")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create temporary file with data file paths
        data_files = []
        batch_info = []
        
        # Verify all files exist before proceeding
        for i, acc in enumerate(accessions):
            if data_type == "rnaseq":
                file_path = Path(f"{acc}/quantification/counts.txt")
            else:
                file_path = Path(f"{acc}/quantification/beta_values.rds")
                
            if not file_path.exists():
                logging.error(f"Data file not found for {acc}: {file_path}")
                return False
                
            data_files.append(str(file_path))
            
            # Read metadata and extend batch info
            metadata_file = Path(f"{acc}/{acc}_metadata.csv")
            if not metadata_file.exists():
                logging.error(f"Metadata file not found for {acc}")
                return False
                
            try:
                sample_count = len(pd.read_csv(metadata_file))
                batch_info.extend([i] * sample_count)
            except Exception as e:
                logging.error(f"Error reading metadata for {acc}: {str(e)}")
                return False
        
        # Save data file paths to temporary file
        with open("data_files.txt", "w") as f:
            f.write("\n".join(data_files))
        
        # Save batch info
        pd.DataFrame({'batch': batch_info}).to_csv('batch_info.csv', index=False)
        
        # Run batch correction with file list
        cmd = f"Rscript R/batch_correction.R data_files.txt batch_info.csv {data_type} results/batch_correction"
        logging.info(f"Running batch correction command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        
        # Verify results
        if not verify_batch_correction(data_type):
            logging.error("Batch correction output files not found")
            return False
            
        # Clean up temporary files
        Path("data_files.txt").unlink()
        Path("batch_info.csv").unlink()
        
        logging.info("Completed batch correction successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to perform batch correction: {str(e)}")
        if hasattr(e, 'stderr'):
            logging.error(f"Stderr: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Failed to perform batch correction: {str(e)}")
        return False

def verify_download(accession: str, data_type: str) -> bool:
    """Verify that data was downloaded correctly"""
    try:
        samples_dir = Path(f"{accession}/samples")
        if not samples_dir.exists():
            return False
            
        if data_type == "rnaseq":
            # Check for either FASTQ files or SRA files
            fastq_files = list(samples_dir.glob("*/SRA/FASTQ/*.fastq.gz"))
            sra_files = list(samples_dir.glob("*/SRA/*.sra"))
            return len(fastq_files) > 0 or len(sra_files) > 0
        else:  # methylation
            # Check for either CEL files or IDAT files
            cel_files = list(samples_dir.glob("*/*.cel.gz"))
            idat_files = list(samples_dir.glob("*/*.idat"))
            chp_files = list(samples_dir.glob("*/*.chp.gz"))
            return len(cel_files) > 0 or len(idat_files) > 0 or len(chp_files) > 0
    except Exception as e:
        logging.error(f"Error during download verification: {str(e)}")
        return False

def verify_processing(accession: str, data_type: str) -> bool:
    """Verify that data processing completed successfully"""
    try:
        if data_type == "rnaseq":
            counts_file = Path(f"{accession}/quantification/counts_matrix.csv")
            return counts_file.exists()
        else:  # methylation
            beta_file = Path(f"{accession}/quantification/beta_values.csv")
            return beta_file.exists()
    except Exception:
        return False

def verify_differential_analysis(accession: str, data_type: str) -> bool:
    """Verify that differential analysis completed successfully"""
    try:
        results_file = Path(f"{accession}/results/{'deseq2_results.rds' if data_type == 'rnaseq' else 'limma_results.rds'}")
        return results_file.exists()
    except Exception:
        return False

def verify_kegg_analysis(accession: str) -> bool:
    """Verify that KEGG analysis completed successfully"""
    try:
        kegg_dir = Path(f"{accession}/results/kegg")
        return kegg_dir.exists() and len(list(kegg_dir.glob("*.pdf"))) > 0
    except Exception:
        return False

def process_single_accession(accession: str, max_workers: int = None) -> bool:
    """Process a single accession with resource management and checkpointing"""
    try:
        # Initialize checkpoint tracking
        checkpoint = PipelineCheckpoint(accession)
        
        # If already completed, skip
        if checkpoint.is_step_completed("completed"):
            logging.info(f"Accession {accession} already fully processed, skipping")
            return True
        
        # Set up resource limits
        if max_workers is None:
            max_workers = min(get_available_cpus(), 4)
        
        # Read metadata
        metadata = pd.read_csv(PATHS['metadata'])
        data_type = metadata[metadata['Accession'] == accession]['Type'].iloc[0]
        data_type = "rnaseq" if "sequencing" in data_type.lower() else "methylation"
        
        # Setup directories if not exist
        create_directory_structure(accession)
        
        # Download data if not already done
        if not checkpoint.is_step_completed("download"):
            if download_geo_data(accession, accession):
                if verify_download(accession, data_type):
                    checkpoint.mark_step_completed("download")
                else:
                    logging.error(f"Download verification failed for {accession}")
                    return False
            else:
                logging.error(f"Failed to download data for {accession}")
                return False
        
        # Process data if not already done
        if not checkpoint.is_step_completed("data_processing"):
            success = False
            if data_type == "rnaseq":
                success = process_rnaseq(accession, accession, max_workers)
            else:
                success = process_methylation(accession, accession, max_workers)
                
            if success and verify_processing(accession, data_type):
                checkpoint.mark_step_completed("data_processing")
            else:
                logging.error(f"Data processing failed for {accession}")
                return False
        
        # Perform differential analysis if not already done
        if not checkpoint.is_step_completed("differential_analysis"):
            if differential_analysis(accession, data_type):
                if verify_differential_analysis(accession, data_type):
                    checkpoint.mark_step_completed("differential_analysis")
                else:
                    logging.error(f"Differential analysis verification failed for {accession}")
                    return False
            else:
                logging.error(f"Differential analysis failed for {accession}")
                return False
        
        # KEGG pathway analysis if not already done
        if not checkpoint.is_step_completed("kegg_analysis"):
            if kegg_pathway_analysis(accession, data_type):
                if verify_kegg_analysis(accession):
                    checkpoint.mark_step_completed("kegg_analysis")
                else:
                    logging.error(f"KEGG analysis verification failed for {accession}")
                    return False
            else:
                logging.error(f"KEGG pathway analysis failed for {accession}")
                return False
        
        # Mark as fully completed
        checkpoint.mark_step_completed("completed")
        logging.info(f"Pipeline completed successfully for {accession}")
        return True
        
    except Exception as e:
        logging.error(f"Pipeline failed for {accession}: {str(e)}")
        return False

def process_all_accessions(max_parallel: int = None):
    """Process all accessions in parallel with resource management and checkpointing"""
    try:
        # Set up resource limits
        if max_parallel is None:
            max_parallel = min(get_available_cpus() // 4, 8)
        
        # Read metadata
        metadata = pd.read_csv(PATHS['metadata'])
        accessions = metadata['Accession'].unique()
        
        # Group accessions by data type and check existing progress
        rnaseq_accessions = []
        methylation_accessions = []
        completed_accessions = []
        
        for acc in accessions:
            checkpoint = PipelineCheckpoint(acc)
            if checkpoint.is_step_completed("completed"):
                completed_accessions.append(acc)
                data_type = metadata[metadata['Accession'] == acc]['Type'].iloc[0]
                if "sequencing" in data_type.lower():
                    rnaseq_accessions.append(acc)
                else:
                    methylation_accessions.append(acc)
                continue
                
            data_type = metadata[metadata['Accession'] == acc]['Type'].iloc[0]
            if "sequencing" in data_type.lower():
                rnaseq_accessions.append(acc)
            else:
                methylation_accessions.append(acc)
        
        # Skip already completed accessions
        accessions_to_process = [acc for acc in accessions if acc not in completed_accessions]
        if not accessions_to_process:
            logging.info("All accessions already processed")
            return True
            
        logging.info(f"Processing {len(accessions_to_process)} accessions, {len(completed_accessions)} already completed")
        
        # Process remaining accessions in parallel
        with ProcessPoolExecutor(max_workers=max_parallel) as executor:
            future_to_acc = {
                executor.submit(process_single_accession, acc, max_workers=4): acc 
                for acc in accessions_to_process
            }
            
            successful = []
            failed = []
            for future in as_completed(future_to_acc):
                acc = future_to_acc[future]
                try:
                    if future.result():
                        successful.append(acc)
                    else:
                        failed.append(acc)
                except Exception as e:
                    logging.error(f"Error processing {acc}: {str(e)}")
                    failed.append(acc)
        
        # Add successful accessions to the appropriate lists for batch correction
        for acc in successful:
            data_type = metadata[metadata['Accession'] == acc]['Type'].iloc[0]
            if "sequencing" in data_type.lower() and acc not in rnaseq_accessions:
                rnaseq_accessions.append(acc)
            elif acc not in methylation_accessions:
                methylation_accessions.append(acc)
        
        # Perform batch correction if needed
        if len(rnaseq_accessions) > 1:
            batch_correction(rnaseq_accessions, "rnaseq")
        if len(methylation_accessions) > 1:
            batch_correction(methylation_accessions, "methylation")
        
        # Log summary
        logging.info(f"Successfully processed: {len(successful)} datasets")
        logging.info(f"Previously completed: {len(completed_accessions)} datasets")
        logging.info(f"Failed to process: {len(failed)} datasets")
        if failed:
            logging.warning(f"Failed accessions: {', '.join(failed)}")
        
        return True
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        return False

def main():
    """Main entry point with command line argument parsing"""
    parser = argparse.ArgumentParser(description='B16F10 RNA-seq and Methylation Array Analysis Pipeline')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--accession', help='Process a single GSE accession')
    group.add_argument('--all', action='store_true', help='Process all accessions in the metadata file')
    parser.add_argument('-p', '--parallel', type=int, help='Maximum number of parallel processes when running all accessions')
    
    args = parser.parse_args()
    
    if args.all:
        process_all_accessions(args.parallel)
    else:
        process_single_accession(args.accession)

if __name__ == "__main__":
    main()
