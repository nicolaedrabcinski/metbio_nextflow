# MetBio-WGSP Bioinformatics Package for Wastewater Surveillance

Automated pipeline for viral strain quantification using kallisto and Freyja for wastewater genomic surveillance.

## Overview

This pipeline implements the bioinformatics workflow developed as part of the MetBio-WGSP project for wastewater sequencing data analysis. In its current state, the pipeline is capable of detecting and quantifying pathogen lineages in wastewater sequencing data (using the tools [kallisto](https://github.com/pachterlab/kallisto) and [Freyja](https://github.com/andersen-lab/Freyja) as detection/quantification engines). While kallisto takes a user-provided list of candidate pathogen genomes, Freyja uses an embedded database of genome references for the lineages of 10 different pathogens -- the pathogen needs to be specified during analysis. Haplotype reconstruction modules and other functionality, such as validation of lineage presence (including for low-abundance lineages), will be implemented in subsequent stages.

## Quick Start

```bash
# Auto-scan FASTQ directory
nextflow run main.nf --fastq_dir data/fastq --lineages_fasta data/lineages.fasta

# Use existing CSV
nextflow run main.nf --input samples.csv --lineages_fasta data/lineages.fasta
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
| `--tool` | - | Specifies which of the bioinformatic instruments will be used. Currently available options are "kallisto" and "Freyja" |
| `--fastq_dir` | - | Directory with FASTQ files. Not needed if `--input` is specified |
| `--input` | - | CSV file with samples. Not needed if `--fastq_dir` is specified |
| `--lineages_fasta` | - | Multi-fasta file containing all candidate pathogen genome references |
| `--outdir` | `results` | Output directory |
| `--k-mer-size` | `31` | k-mer length (integer (odd), max value: 31). Lower k-mer lengths might enable higher sensitivity of lineage detection (this is true especially for shorter reads), while potentially increasing the rate of falsely assigned reads |
| `--single` | not included | Quantify single-end reads as opposed to paired-end reads |
| `--fragment-length` | - | Estimated average length of fragments in the sequencing library. Equals the average length of amplicons in case the amplicons were not additionally fragmented during sequencing library preparation |
| `--sd` | - | Estimated standard deviation of fragment length, or of the amplicons in case these were not additionally fragmented (default: -l, -s values are estimated from paired end data, but are required when using `--single`) |
| `--threads` | `1` | Number of threads to use for index construction and quantification (integer) |
=======
| `--fastq_dir` | - | Directory containing FASTQ files |
| `--input` | - | CSV file with sample definitions |
| `--lineages_fasta` | - | Reference lineages FASTA file |
| `--outdir` | `results` | Output directory name |
| `--filter_shuffled` | `false` | Only process files containing '_shuffled' |
| `--fragment_length` | `1000` | Estimated fragment length for single-end reads |
| `--fragment_sd` | `30` | Standard deviation of fragment length |
>>>>>>> e0c5f307b62d78f1c9692a0731006855049bae16

## Pipeline Workflow

kallisto:

1. **Index Building**: Creates kallisto index from reference lineages based on input genome references
2. **Pseudoalignment**: Pseudoalignment of the reads (through comparisons of k-mers statistics), followed by estimation of the relative abundances of each genomic sequence (via expectation maximization)
3. **Visualization**: Generation of summary plots and abundance matrices

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

## Usage Examples

```bash
# Basic run with toy dataset
nextflow run main.nf \
  --fastq_dir data/reads_artic_small_overlaps \
  --lineages_fasta data/lineages.fasta \
  --fragment-length 1000 \
  --sd 30
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
The pipeline has been validated against MetBio project benchmarks showing high accuracy for both Illumina and Nanopore data, with kallisto demonstrating optimal performance for wastewater surveillance applications.

## Contributors

- Victor Gordeev (Pipeline Design & Validation)
- Nicolae Drabcinski (Nextflow Implementation)

## Acknowledgment

The development of this software package is supported by the grant of the Ministry of Research, Innovation and Digitization, under Romania's National Recovery and Resilience Plan, Funded by the European Union NextGenerationEU program no. 760286/27.03.2024, code 167/31.07.2023, within Pillar III, Component C9, Investment 8.

The pipeline was developed as part of the project "Metagenomics and Bioinformatics tools for Wastewater-based Genomic Surveillance of viral Pathogens for early prediction of public health risks (MetBio-WGSP)"
