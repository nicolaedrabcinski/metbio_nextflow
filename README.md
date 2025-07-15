# MetBio-WGSP Bioinformatics Package for Wastewater Surveillance

A Nextflow pipeline for wastewater genomic surveillance intended for comprehensive and flexible analysis of targeted sequencing data and sensitive detection and quantification of known and novel pathogens.

## Overview

This pipeline implements the bioinformatics workflow developed as part of the MetBio-WGSP project for wastewater sequencing data analysis. In its current state, the pipeline is capable of detecting and quantifying complex mixture of pathogens including closely related pathogen lineages from wastewater sequencing data using the tool kallisto as detection/quantification engine. Apart from the preprocessed FASTQ files, Kallisto takes a user-provided list of candidate reference pathogen genomes. The pipeline also incorporates a toy dataset (in the data folder of this repository) with precisely controlled and known mixture compositions, which allow testing the current and future functionality of this package, including instruments for pathogen/lineage detection, quantification, and haplotype reconstruction. Other components of the quantification module, haplotype reconstruction modules and other functionality, such as validation of lineage presence (including for low-abundance lineages), will be implemented in subsequent stages of the project.

## Quick Start

```bash
# Auto-scan FASTQ directory
nextflow run main.nf --fastq_dir data/fastq --lineages_fasta data/lineages.fasta

# Use existing CSV
nextflow run main.nf --input samples.csv --lineages_fasta data/lineages.fasta
```

## Requirements

- Nextflow ≥24.04
- Python 3.6+ with pandas, matplotlib, seaborn, numpy
- kallisto ≥0.51 (auto-downloaded if missing)

## Input Data

### Required Files

1. **FASTQ files**: Raw sequencing data (single or paired-end)
   - For Nanopore: Use files with "shuffled" in the filename
   - Supported formats: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`

2. **Reference sequences**: FASTA file containing viral lineage genomes
   - Location: `data/lineages/` directory
   - Format: Multi-FASTA with lineage genomes

### Test Dataset

The pipeline includes a toy dataset for testing:
- Long-read sequences with different overlap lengths
- SARS-CoV-2 lineage references
- Expected outputs for validation

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
| `--fastq_dir` | - | Directory with FASTQ files. Not needed if `--input` is specified |
| `--input` | - | CSV file with samples. Not needed if `--fastq_dir` is specified |
| `--lineages_fasta` | - | Multi-fasta file containing all candidate pathogen genome references |
| `--outdir` | `results` | Output directory |
| `--k-mer-size` | `31` | k-mer length (integer (odd), max value: 31). Lower k-mer lengths might enable higher sensitivity of lineage detection (this is true especially for shorter reads), while potentially increasing the rate of falsely assigned reads |
| `--single` | not included | Quantify single-end reads as opposed to paired-end reads |
| `--fragment-length` | - | Estimated average length of fragments in the sequencing library. Equals the average length of amplicons in case the amplicons were not additionally fragmented during sequencing library preparation |
| `--sd` | - | Estimated standard deviation of fragment length, or of the amplicons in case these were not additionally fragmented (default: --fragment_length, --sd values are estimated from paired end data, but are required when using `--single`) |
| `--threads` | `1` | Number of threads to use for index construction and quantification (integer) |

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

## Description of the embedded ONT toy dataset

We simulated and embedded this dataset to allow testing the functionality of the pipeline (including detection, quantification and assembly of pathogen genomes). It contains 51 mixtures with precisely controlled mixture compositions and a known number of reads for each lineage of each mixture. Importantly, the dataset contains mixtures of various complexity, including samples containing a single lineage per mixture, the latter being intended as convenient controls for testing false positive results. The simulated FASTQ files, the exact composition of the ground truth mixtures, including the exact genome sequence for each lineage, but also the sequences of the primers used and the corresponding amplicons obtained for each lineage can be accessed from the "data" folder. The simulation procedure was as follows. Subsets of the ARTIC 5.3.2 primers covering the S gene, were chosen so as to achieve amplicon lengths of ~1000 bp and either small or big amplicon overlaps (~100 bp and ~400 bp respectively). All amplicons were successfully generated from the reference lineages included in the data folder using the publicly available in_silico_PCR tool (https://github.com/egonozer/in_silico_pcr) and the primer sequences. 3) These amplicons were used by the NanoSim simulator to generate long reads covering the full length of the amplicons (although this is a simplification compared to real ONT data since read lengths follow a distribution, this simplification is very useful to ensure an intuitive and known coverage level across the entire sequenced region). We trained the model for sequencing errors and base quality scores used by the NanoSim simulator on a sample obtained by sequencing synthetic RNA of the Wuhan strain of SARS-CoV-2 on a MinION device, R10.4.1 chemistry, followed by basecalling with guppy, in High Accuracy Mode. The public dataset containing the sample is accessible from https://zenodo.org/records/7786559). The authors report for this control sample an average PHRED score of 14.4 (with an accuracy of roughly 96.37%).

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
Kallisto has been extensively tested against roughly 20 other competing tools as part of the MetBio-WGSP project benchmarks, showing high accuracy for both Illumina and especially for long Nanopore reads, as well as high adaptability for monitoring pathogens other than SARS-CoV-2.

## Contributors

- Victor Gordeev (Pipeline Design & Validation)
- Nicolae Drabcinski (Nextflow Implementation)

## Acknowledgment

The development of this software package is supported by the grant of the Ministry of Research, Innovation and Digitization, under Romania's National Recovery and Resilience Plan, Funded by the European Union NextGenerationEU program no. 760286/27.03.2024, code 167/31.07.2023, within Pillar III, Component C9, Investment 8.

The pipeline was developed as part of the project "Metagenomics and Bioinformatics tools for Wastewater-based Genomic Surveillance of viral Pathogens for early prediction of public health risks (MetBio-WGSP)"
