# MetBio Nextflow Pipeline

Automated pipeline for viral strain quantification using kallisto/lr-kallisto.

## Quick Start

```bash
# Auto-scan FASTQ directory
nextflow run main.nf --fastq_dir data/fastq --lineages_fasta lineages.fasta

# Use existing CSV
nextflow run main.nf --input samples.csv --lineages_fasta lineages.fasta
```

## Requirements

- Nextflow ≥ 24.04
- Python 3.6+
- kallisto (auto-downloaded if missing)

## Input

**Option 1: FASTQ Directory**
```bash
--fastq_dir /path/to/fastq/files
```
Automatically finds all `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` files.

**Option 2: CSV File**
```csv
sample,fastq
sample1,/path/to/sample1.fastq
sample2,/path/to/sample2.fastq
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastq_dir` | - | Directory with FASTQ files |
| `--input` | - | CSV file with samples |
| `--lineages_fasta` | - | Reference lineages FASTA |
| `--outdir` | `results` | Output directory |
| `--filter_shuffled` | `false` | Only use files with '_shuffled' |

## Output

```
results/[input_name]/
├── metadata/
│   └── auto_samples.csv
├── quantification/
│   ├── sample1/
│   │   ├── abundance.tsv
│   │   └── run_info.json
│   └── sample2/
└── index/
    └── lineages.idx
```

## Examples

```bash
# Process all FASTQ files
nextflow run main.nf --fastq_dir data/reads --lineages_fasta lineages.fasta

# Custom output directory
nextflow run main.nf --fastq_dir data/reads --lineages_fasta lineages.fasta --outdir my_analysis
```

## License

MetBio Project - NextGenerationEU Program