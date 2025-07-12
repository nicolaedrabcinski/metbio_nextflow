# MetBio Nextflow Pipeline

Automated pipeline for viral strain quantification using kallisto/lr-kallisto for wastewater genomic surveillance.

## Overview

This pipeline implements the MetBio project's bioinformatics workflow for quantifying viral lineages in wastewater sequencing data. It uses kallisto's lr-kallisto algorithm optimized for long reads and high error rates, making it ideal for Oxford Nanopore sequencing data.

## Quick Start

```bash
# Auto-scan FASTQ directory
nextflow run main.nf --fastq_dir data/fastq --lineages_fasta lineages.fasta

# Use existing CSV
nextflow run main.nf --input samples.csv --lineages_fasta lineages.fasta

# Process only shuffled files (recommended for toy dataset)
nextflow run main.nf --fastq_dir data/reads_artic_small_overlaps --lineages_fasta data/lineages.fasta --filter_shuffled true
```

## Requirements

- Nextflow ≥ 24.04
- Python 3.6+ with pandas, matplotlib, seaborn, numpy
- kallisto ≥ 0.51 (auto-downloaded if missing)

## Input Formats

**Option 1: FASTQ Directory Auto-Detection**
```bash
--fastq_dir /path/to/fastq/files
```
Automatically finds all `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` files.

**Option 2: Sample Sheet CSV**
```csv
sample,fastq
mixture_1,/path/to/mixture_1_shuffled.fastq
mixture_2,/path/to/mixture_2_shuffled.fastq
```

**Reference Lineages**
- Multi-FASTA file containing viral lineage sequences
- Used to build kallisto index for quantification

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastq_dir` | - | Directory containing FASTQ files |
| `--input` | - | CSV file with sample definitions |
| `--lineages_fasta` | - | Reference lineages FASTA file |
| `--outdir` | `results` | Output directory name |
| `--filter_shuffled` | `false` | Only process files containing '_shuffled' |
| `--fragment_length` | `900` | Estimated fragment length for single-end reads |
| `--fragment_sd` | `30` | Standard deviation of fragment length |

## Pipeline Workflow

1. **Index Building**: Creates kallisto index from reference lineages
2. **Quantification**: Runs lr-kallisto on each sample using:
   - Single-end mode (`--single`)
   - Long read parameters optimized for Nanopore data
   - Fragment length estimation for wastewater samples
3. **Visualization**: Generates summary plots and abundance matrices

## Output Structure

```
results/[input_name]/
├── metadata/
│   └── auto_samples.csv              # Auto-generated sample sheet
├── quantification/
│   ├── mixture_1/
│   │   ├── abundance.tsv             # TPM abundance values
│   │   ├── abundance.h5              # HDF5 format results
│   │   └── run_info.json             # Kallisto run statistics
│   └── mixture_2/
│       └── ...
├── index/
│   └── lineages.idx                  # Kallisto index file
├── summary_percentages.csv           # Combined abundance matrix
├── top_lineages_barplot.png          # Top 10 lineages visualization
└── abundance_heatmap.png             # Sample vs lineage heatmap
```

## Key Features

- **lr-kallisto Integration**: Uses latest kallisto (v0.51+) with long-read support
- **Wastewater Optimized**: Parameters tuned for wastewater surveillance data
- **Flexible Input**: Auto-detection or manual sample specification
- **Comprehensive Output**: Quantification + visualization in one workflow
- **Reproducible**: Containerized with Nextflow for consistent results

## Usage Examples

```bash
# Basic run with toy dataset
nextflow run main.nf \
  --fastq_dir data/reads_artic_small_overlaps \
  --lineages_fasta data/lineages.fasta \
  --filter_shuffled true

# Production run with custom parameters
nextflow run main.nf \
  --input samples.csv \
  --lineages_fasta references/sars_cov2_lineages.fasta \
  --outdir analysis_2025 \
  --fragment_length 1200 \
  --fragment_sd 50

# Resume failed run
nextflow run main.nf -resume \
  --fastq_dir data/reads \
  --lineages_fasta lineages.fasta
```

## Performance Notes

- Processing time: ~20-30 seconds for 102 samples (toy dataset)
- Memory usage: Minimal for kallisto, moderate for visualization
- Parallel execution: Samples processed concurrently
- Scalability: Tested with 100+ samples

## Troubleshooting

**Common Issues:**
- Missing numpy: Ensure `numpy` is installed for visualization
- Kallisto version: Pipeline requires v0.51+ for lr-kallisto support
- File permissions: Ensure FASTQ files are readable
- Memory: Large datasets may require cluster execution

**Validation:**
The pipeline has been validated against MetBio project benchmarks showing high accuracy for both Illumina and Nanopore data, with lr-kallisto demonstrating optimal performance for long-read wastewater surveillance applications.

## Citation

If you use this pipeline, please cite:

> MetBio Project - "Metagenomics and Bioinformatics tools for Wastewater-based Genomic Surveillance of viral Pathogens for early prediction of public health risks"
> 
> Funded by EU NextGenerationEU Program, Romania's National Recovery and Resilience Plan
> Project number 760286/27.03.2024, code 167/31.07.2023

## Contributors

- Victor Gordeev (Pipeline Design & Validation)
- Nicolae Drabcinski (Nextflow Implementation)
- MetBio Consortium

## License

This project is part of the MetBio initiative funded under Romania's National Recovery and Resilience Plan through the EU NextGenerationEU program.