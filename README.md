# B16F10 RNA-seq and Methylation Array Analysis Pipeline

This pipeline processes and analyzes RNA-seq and methylation array data from B16F10 melanoma cells. It handles data download, processing, differential analysis, and pathway analysis using the KEGG database.

## Requirements

### Python Dependencies
- pandas
- numpy
- rpy2
- pathlib
- logging

### R Dependencies
- GEOquery
- DESeq2
- limma
- edgeR
- minfi
- gage
- pathview
- sva
- biomaRt

### External Tools
- FastQC
- Trimmomatic
- STAR
- featureCounts

## Installation

1. Clone this repository:
```bash
git clone <repository_url>
cd B16F10
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Install R dependencies:
```R
install.packages(c("GEOquery", "DESeq2", "limma", "edgeR", "minfi", "gage", "pathview", "sva", "biomaRt"))
```

4. Install external tools (example for Ubuntu):
```bash
sudo apt-get install fastqc trimmomatic star subread
```

## Configuration

Edit `config.py` to set:
- System resources (threads, memory)
- Reference genome paths
- Analysis parameters
- File paths

## Usage

Run the pipeline for a single accession:
```bash
python b16f10_pipeline.py <accession_id>
```

Example:
```bash
python b16f10_pipeline.py GSE254073
```

## Pipeline Steps

1. **Setup**
   - Creates directory structure
   - Initializes R environment
   - Downloads data from GEO

2. **RNA-seq Processing**
   - Quality control (FastQC)
   - Adapter trimming (Trimmomatic)
   - Alignment (STAR)
   - Feature counting (featureCounts)

3. **Methylation Array Processing**
   - Data import and quality control
   - Normalization
   - Beta value calculation

4. **Differential Analysis**
   - RNA-seq: DESeq2
   - Methylation: limma

5. **KEGG Pathway Analysis**
   - Pathway enrichment
   - Visualization

6. **Batch Effect Correction**
   - ComBat for multiple dataset integration

## Output Structure

```
accession/
├── samples/
│   └── raw and processed sequencing data
├── alignment/
│   └── aligned BAM files
├── quantification/
│   └── count matrices or beta values
├── results/
│   ├── differential analysis results
│   ├── KEGG pathway analysis
│   └── quality control reports
└── logs/
    └── pipeline execution logs
```

## Notes

- The pipeline automatically detects data type (RNA-seq or methylation array)
- Batch effect correction is performed when analyzing multiple datasets
- All parameters can be customized in `config.py`
- Extensive logging is implemented for troubleshooting

## Citation

If you use this pipeline, please cite:
- DESeq2
- limma
- STAR
- featureCounts
- Other tools as appropriate

## License

[Add your license information here] 