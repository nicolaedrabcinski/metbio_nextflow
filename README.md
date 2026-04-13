# MetBio-WGSP

**Bioinformatics pipeline for wastewater genomic surveillance of pathogens.**

A Nextflow DSL2 pipeline for comprehensive analysis of wastewater sequencing data (Oxford Nanopore / Illumina), integrating tools for ONT read preprocessing, taxonomic profiling, lineage quantification, variant analysis, and haplotype reconstruction.

---

## Pipeline overview

```
Raw ONT FASTQ reads
        │
        ▼  [optional: --run_ont_preprocess]
┌─────────────────────────────────┐
│       ONT Preprocessing         │
│  NanoPlot (raw QC)              │
│      │                          │
│  Porechop_ABI (adapter trim)    │
│      │                          │
│  Filtlong (quality/length filter│
│      │                          │
│  [Host depletion - optional]    │
│      │                          │
│  NanoPlot (clean QC)            │
└─────────────────────────────────┘
        │
        ├──────────────────────────────────────────────────┐
        │                                                  │
        ▼  [if --run_kraken2]                              ▼  [if --reference_fasta]
┌────────────────────┐                         ┌──────────────────────────────────┐
│  Kraken2 + Bracken │                         │  Freyja                          │
│  Taxonomic         │                         │  minimap2 align → variants →     │
│  profiling         │                         │  demix → aggregate → plots       │
│  ↓                 │                         │       │            │              │
│  Krona tree        │                         │       ▼            ▼              │
│  Sankey/heatmap    │                         │  HaploDMF      VILOCA            │
│  plots             │                         │  haplotype     haplotype         │
│                    │                         │  reconstruction reconstruction   │
│                    │                         │       │                           │
└────────────────────┘                         │       ▼                           │
                                               │  Pangolin lineage classification │
        ▼  [if --lineages_fasta]               └──────────────────────────────────┘
┌────────────────────┐
│  lr-Kallisto       │
│  Lineage           │
│  quantification    │
└────────────────────┘
```

Tools are selected automatically based on provided parameters — only tools with the necessary inputs are run.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.04.0
- [Docker](https://www.docker.com/) (recommended)
- Python 3.8+ with `pandas`, `matplotlib`, `seaborn`, `plotly`, `numpy`

All bioinformatics tools run as Docker containers pulled automatically:

| Tool | Container | Purpose |
|------|-----------|---------|
| NanoPlot | `staphb/nanoplot:latest` | ONT read quality reports |
| Porechop_ABI | `quay.io/biocontainers/porechop_abi:0.5.1--py310h275bdba_0` | ONT adapter trimming |
| Filtlong | `staphb/filtlong:latest` | ONT quality/length filtering |
| minimap2 + samtools | `mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:...` | Read alignment |
| iVar | `staphb/ivar:latest` | Amplicon primer trimming |
| Freyja | `staphb/freyja:latest` | Variant calling, demixing, aggregation |
| Kallisto | `quay.io/biocontainers/kallisto:0.51.1--h2b92561_2` | Lineage quantification (lr-kallisto) |
| Kraken2 | `staphb/kraken2:latest` | Taxonomic classification |
| Bracken | `staphb/bracken:latest` | Abundance re-estimation |
| Krona | `nanozoo/krona:2.7.1--e7615f7` | Interactive taxonomy charts |
| VILOCA | `viloca:latest` (custom build) | Local haplotype reconstruction |
| HaploDMF | `haplodmf:latest` (custom build) | Deep learning haplotype reconstruction |
| Pangolin | `staphb/pangolin:latest` | SARS-CoV-2 lineage classification |

---

## Quick start

### 1. ONT preprocessing only

```bash
nextflow run main.nf \
  --fastq_dir data/reads \
  --run_ont_preprocess \
  --run_kraken2 \
  --kraken2_db databases/kraken2 \
  -profile docker
```

### 2. SARS-CoV-2 variant analysis (Freyja)

```bash
nextflow run main.nf \
  --fastq_dir data/reads \
  --reference_fasta data/NC_045512_2.fasta \
  -profile docker
```

### 3. Haplotype reconstruction (HaploDMF + Pangolin)

```bash
nextflow run main.nf \
  --fastq_dir data/reads \
  --reference_fasta data/NC_045512_2.fasta \
  --run_haplodmf \
  --run_pangolin \
  -profile docker
```

### 4. Full pipeline (ONT preprocessing + all tools)

```bash
nextflow run main.nf \
  --fastq_dir data/reads \
  --run_ont_preprocess \
  --reference_fasta data/NC_045512_2.fasta \
  --lineages_fasta data/lineages.fasta \
  --run_kraken2 --kraken2_db databases/kraken2 \
  --run_haplodmf \
  --run_viloca \
  --run_pangolin \
  -profile docker
```

### 5. Amplicon data with primer trimming (iVar)

```bash
nextflow run main.nf \
  --fastq_dir data/reads \
  --reference_fasta data/NC_045512_2.fasta \
  --run_ivar_trim \
  --primer_bed data/artic_v5.3.2.bed \
  --run_haplodmf \
  --run_pangolin \
  -profile docker
```

---

## Input data

### Supported formats

- `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`
- Oxford Nanopore Technologies (ONT): MinION, GridION, PromethION
- Illumina: MiSeq, NextSeq (use `--minimap2_preset sr` for alignment)

### Input options

**Option 1 — FASTQ directory (auto-detection):**
```bash
--fastq_dir /path/to/fastq/files
```
The pipeline auto-generates a sample sheet from all FASTQ files in the directory. Sample names are derived from filenames.

**Option 2 — Sample sheet CSV:**
```bash
--input samples.csv
```
```csv
sample,fastq
sample_01,/path/to/sample_01.fastq.gz
sample_02,/path/to/sample_02.fastq.gz
```

---

## Parameters

### Input / Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--fastq_dir` | Directory with FASTQ files | — |
| `--input` | CSV sample sheet | — |
| `--outdir` | Output directory | `results` |
| `--threads` | Number of threads | 32 |

### Tool selection

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reference_fasta` | Reference genome FASTA → enables Freyja | — |
| `--lineages_fasta` | Multi-FASTA lineages → enables Kallisto | — |
| `--run_kraken2` | Enable Kraken2 + Bracken | `false` |
| `--kraken2_db` | Path to Kraken2 database | — |
| `--run_ont_preprocess` | Enable ONT preprocessing | `false` |
| `--run_haplodmf` | Enable HaploDMF haplotype reconstruction | `false` |
| `--run_viloca` | Enable VILOCA haplotype reconstruction | `false` |
| `--run_pangolin` | Enable Pangolin (requires `--run_haplodmf`) | `true` |
| `--run_ivar_trim` | Enable iVar primer trimming (amplicon data) | `false` |
| `--primer_bed` | BED file with primer coordinates | — |

### ONT preprocessing parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--ont_min_quality` | Filtlong: minimum Phred quality score | `8` |
| `--ont_min_length` | Filtlong: minimum read length (bp) | `200` |
| `--ont_target_bases` | Filtlong: target total bases (e.g. `500000000`) | — |
| `--ont_skip_adapter_trim` | Skip Porechop_ABI adapter trimming | `false` |
| `--ont_skip_quality_filter` | Skip Filtlong filtering | `false` |
| `--ont_host_ref` | Host genome FASTA for depletion (e.g. hg38) | — |
| `--ont_skip_host_depletion` | Skip host read removal | `true` |

### Freyja parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--minimap2_preset` | Alignment preset (`map-ont` / `sr`) | `map-ont` |
| `--refname` | Reference sequence name | `NC_045512.2` |
| `--minq` | Minimum base quality for variant calling | `20` |
| `--freyja_eps` | Epsilon for demix (abundance threshold) | `0.001` |
| `--covcut` | Coverage cutoff | `10` |
| `--depthcutoff` | Depth cutoff | `0` |
| `--confirmedonly` | Use only confirmed lineages | `false` |
| `--pathogen` | Pathogen type | `SARS-CoV-2` |
| `--mincov` | Minimum coverage for aggregation | `60` |
| `--thresh` | Abundance threshold for plotting | `0` |

### Kraken2 / Bracken parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--kraken2_confidence` | Classification confidence threshold | `0.0` |
| `--kraken2_min_hit_groups` | Minimum hit groups | `2` |
| `--bracken_level` | Taxonomic level (D/P/C/O/F/G/S) | `S` |
| `--bracken_read_len` | Read length for Bracken model | `150` |
| `--bracken_threshold` | Minimum reads for re-estimation | `10` |

### HaploDMF parameters

HaploDMF reconstructs haplotypes from aligned reads using deep matrix factorisation. All parameters below are passed directly to the HaploDMF tool.

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| `--haplodmf_map_qual` | `-mq` | Minimum mapping quality | `0` |
| `--haplodmf_start_pos` | `-sp` | Starting genome position for reconstruction | — |
| `--haplodmf_end_pos` | `-ep` | Ending genome position for reconstruction | — |
| `--haplodmf_abundance` | `-a` | Filter out haplotypes below this relative abundance | `0` |
| `--haplodmf_error_rate` | `-e` | Expected sequencing error rate | `0.1` |
| `--haplodmf_signi_level` | `-s` | Significance level for binomial tests (SNV calling) | `0.05` |
| `--haplodmf_cond_pro` | `-cp` | Max conditional probability threshold | `0.65` |
| `--haplodmf_fre_snv` | `-f` | Dominant base frequency threshold for SNV filtering | `0.80` |
| `--haplodmf_smallest_snv` | `-ss` | Minimum SNV sites required for haplotype construction | `20` |
| `--haplodmf_weight` | `-w` | Weight for loss term 2 in DMF optimisation | `0.2` |
| `--haplodmf_learning_rate` | `-lr` | Learning rate for DMF | `0.001` |
| `--haplodmf_epoch` | `-epo` | Number of DMF training epochs | `20` |
| `--haplodmf_batch_size` | `-bs` | Batch size for DMF training | `128` |
| `--haplodmf_algorithm` | `-al` | Clustering algorithm: `ward` or `kmeans` | `ward` |
| `--haplodmf_depth` | `-d` | Minimum read depth for consensus sequence generation | `5` |
| `--haplodmf_cluster_thres1` | `-c1` | Edit distance similarity threshold 1 (haplotype merging) | `0.95` |
| `--haplodmf_cluster_thres2` | `-c2` | Edit distance similarity threshold 2 (haplotype merging) | `0.90` |
| `--haplodmf_largest_iter` | `-li` | Maximum clustering iterations | `20` |
| `--haplodmf_input_snv` | `-is` | Path to file with pre-defined SNV sites | — |

> **Note:** If `--haplodmf_smallest_snv` cannot be satisfied (too few SNV sites in the sample), HaploDMF is skipped for that sample and the pipeline continues without error (`errorStrategy = 'ignore'`).

### iVar primer trimming parameters

Applicable when `--run_ivar_trim true` is set (amplicon data only). iVar trims primer sequences from BAM alignments using a BED file of primer coordinates.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--primer_bed` | BED file with amplicon primer coordinates | — |
| `--ivar_min_length` | Minimum read length after trimming | `30` |
| `--ivar_min_quality` | Minimum base quality for soft-clipping | `20` |

---

## Output structure

```
results/<dataset_name>/
│
├── metadata/
│   └── auto_samples.csv              # Auto-generated sample sheet
│
├── qc/                               # ONT QC (--run_ont_preprocess)
│   ├── nanoplot_raw/<sample>/        # NanoPlot reports before preprocessing
│   ├── adapter_trimmed/<sample>/     # Porechop_ABI trimmed reads
│   ├── filtered/<sample>/            # Filtlong filtered reads
│   ├── nanoplot_clean/<sample>/      # NanoPlot reports after preprocessing
│   └── ont_qc_summary.tsv            # Summary table: reads/bases/quality before & after
│
├── kraken2/
│   ├── <sample>/
│   │   ├── <sample>.kraken2.report   # Kraken2 standard report
│   │   └── <sample>.kraken2.out      # Per-read classification
│   ├── bracken/<sample>/
│   │   ├── <sample>.bracken.report   # Bracken re-estimated report
│   │   └── <sample>.bracken.out      # Species abundance table
│   ├── combined_kraken2_reports/     # All reports (Pavian-compatible)
│   ├── krona/                        # Interactive Krona HTML charts
│   ├── plots/
│   │   ├── taxonomy_sankey_<sample>.png
│   │   ├── taxonomy_heatmap.png
│   │   ├── taxonomy_stacked_bar.png
│   │   ├── taxonomy_top_species.png
│   │   └── taxonomy_alpha_diversity.png
│   └── kraken2_summary.csv           # Classification rate per sample
│
├── alignments/                       # BAM files (minimap2)
├── variants/                         # Variant calls (Freyja)
├── demix/                            # Lineage demixing results (Freyja)
├── aggregated_freyja.tsv             # Aggregated Freyja results
│
├── plots/
│   ├── freyja_plot.pdf               # Freyja lineage abundance plot
│   ├── results_table.csv             # WHO variant summary table
│   └── who_plots/
│       ├── who_plot_stacked.pdf
│       └── who_plot_stacked.png
│
├── haplodmf/                         # HaploDMF haplotype reconstruction
│   ├── <sample>/
│   │   └── <sample>_haplodmf_outputs/
│   │       ├── <sample>_haplotypes.fasta   # Reconstructed haplotype sequences
│   │       ├── <sample>_acgt.txt           # Per-position nucleotide counts
│   │       └── <sample>_haplodmf.log       # HaploDMF run log
│   └── haplodmf_summary.csv              # Aggregated table: sample, haplotype count, haplotype info
│
├── viloca/                           # VILOCA haplotype reconstruction
│   ├── <sample>/
│   └── viloca_summary.csv
│
├── pangolin/                         # Pangolin lineage classification
│   ├── <sample>/
│   │   └── <sample>_lineage_report.csv    # Per-haplotype lineage assignments
│   └── pangolin_summary.csv              # Aggregated table across all samples
│
├── quantification/                   # Kallisto lr-kallisto results
│   └── <sample>/abundance.tsv
│
├── report.html                       # Nextflow execution report
├── timeline.html                     # Execution timeline
└── trace.txt                         # Per-process resource usage
```

---

## Databases

### Kraken2 database

A pre-built database is required. The pipeline was tested with the **PlusPF** database (standard genomes + protozoa + fungi):

```bash
# Download pre-built database (recommended, ~80 GB)
mkdir -p databases/kraken2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz
tar -xzf k2_pluspf_20240605.tar.gz -C databases/kraken2
```

### SARS-CoV-2 reference

```bash
# NC_045512.2 (Wuhan-Hu-1)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
gunzip GCF_009858895.2_ASM985889v3_genomic.fna.gz
```

### Freyja lineage database

Freyja's lineage database is updated automatically on each run. To update manually:

```bash
docker run --rm -v $PWD:/data staphb/freyja:latest freyja update --outdir /data/freyja_db
```

---

## Interpreting HaploDMF + Pangolin results

### HaploDMF output

Each reconstructed haplotype in `<sample>_haplotypes.fasta` has a header encoding key metadata:

```
>haplotype_0_length_29903_abundance_0.65_number_of_reads_15000_depth_1200.5
```

| Field | Meaning |
|-------|---------|
| `haplotype_N` | Haplotype index (0-based) |
| `length_N` | Consensus sequence length (bp) |
| `abundance_N` | Estimated relative abundance (0–1) |
| `number_of_reads_N` | Number of reads supporting this haplotype |
| `depth_N` | Mean read depth across the haplotype |

The `haplodmf_summary.csv` aggregates results per sample:

| Column | Meaning |
|--------|---------|
| `sample` | Sample name |
| `haplotype_count` | Number of haplotypes reconstructed |
| `haplotypes_info` | Semicolon-separated list of `haplotype;abundance;reads` per haplotype |

### Pangolin output

Each haplotype FASTA is classified by Pangolin. Key columns in `<sample>_lineage_report.csv`:

| Column | Meaning |
|--------|---------|
| `taxon` | Haplotype name (from FASTA header, includes abundance) |
| `lineage` | Assigned Pango lineage (e.g. `BA.2`, `XBB.1.5`) |
| `conflict` | Conflict score (lower = more confident) |
| `scorpio_call` | WHO variant name if applicable (e.g. `Omicron (BA.2)`) |
| `scorpio_support` | Fraction of defining mutations present |
| `qc_status` | `pass` or `fail` — sequences shorter than ~25 kb often fail |
| `qc_notes` | Reason for QC failure (e.g. `Ambiguous_content`) |

> **Tip:** Haplotypes with `qc_status = fail` and `Ambiguous_content` typically arise from low-depth regions. The abundance value in the haplotype name is still valid; the lineage call should be treated with caution. Cross-check with Freyja demix results for confirmation.

---

## Benchmark validation

### SARS-CoV-2 lineage mixture dataset (ijms_2023)

The pipeline was validated on **simulated SARS-CoV-2 lineage mixtures** from Kovacikova et al. 2023 ([doi:10.3390/ijms24241718](https://doi.org/10.3390/ijms/24/24/17184)), sourced from Zenodo ([record 7786559](https://zenodo.org/records/7786559)).

The dataset contains **20 wastewater mixtures** with precisely defined lineage compositions:

- **10 R9 mixtures** — simulated using NanoSim trained on MinION R9.4.1 reads
- **10 R10 mixtures** — simulated using NanoSim trained on MinION R10.4.1 reads
- ARTIC v5.3.2 amplicon scheme (with overlapping amplicons)
- Lineages include: Omicron (BA.1/BA.2), Delta (AY.*), Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1), and ancestral variants

**Pipeline run commands:**

```bash
# R9 chemistry
nextflow run main.nf \
  --fastq_dir data/ijms_2023/R9 \
  --reference_fasta data/NC_045512_2.fasta \
  --run_haplodmf --run_pangolin \
  --outdir results/R9 -profile docker

# R10 chemistry
nextflow run main.nf \
  --fastq_dir data/ijms_2023/R10 \
  --reference_fasta data/NC_045512_2.fasta \
  --run_haplodmf --run_pangolin \
  --outdir results/R10 -profile docker
```

**Execution times:** R9 — 8h 2m; R10 — 4h 37m (10 samples each, single-process mode)

Results are located in `results/R9/` and `results/R10/`. Ground-truth mixture compositions are in `data/mixture_compositions.xlsx`.

---

### ZymoBIOMICS Gut Microbiome Standard (Kraken2 validation)

The pipeline was validated on the **ZymoBIOMICS Gut Microbiome Standard** (SRR13128014, ~2M ONT reads):

| Step | Result |
|------|--------|
| Porechop_ABI | No adapters found (data already trimmed) |
| Filtlong | Removed low-quality/short reads |
| Kraken2 classification rate | **99.97%** (1,978,068 / 1,978,646 reads) |
| All 10 Zymo species detected | ✔ |

Top species detected by Bracken vs expected Zymo Gut Microbiome Standard composition:

| Species | Bracken | Expected |
|---------|---------|----------|
| *Faecalibacterium prausnitzii* | 21.1% | ~10.4% |
| *Bacteroides fragilis* | 19.3% | ~10.4% |
| *Veillonella rogosae* | 16.3% | ~10.4% |
| *Escherichia coli* | 16.1% | ~10.4% |
| *Prevotella corporis* | 7.6% | ~10.4% |
| *Roseburia hominis* | 5.7% | ~7.8% |
| *Fusobacterium nucleatum* | 3.7% | ~10.4% |
| *Clostridioides difficile* | 2.7% | ~10.4% |
| *Akkermansia muciniphila* | 2.0% | ~10.4% |

All expected species are detected. Proportion deviations are typical for Kraken2/Bracken with ONT long reads due to genome size differences and read length biases (no genome-length normalization applied).

---

## Pavian visualization

Kraken2/Bracken reports in `combined_kraken2_reports/` are compatible with [Pavian](https://github.com/fbreitwieser/pavian):

```r
# Install
if (!require(remotes)) install.packages("remotes")
remotes::install_github("fbreitwieser/pavian")

# Run
pavian::runApp(port=5000)
```

Upload `.kraken2.report` and `.bracken.report` files from the output directory.

---

## Execution profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run all tools in Docker containers (recommended) |
| `local` | Use locally installed tools |
| `mixed` | Mix of Docker and local tools |

---

## Custom Docker images

VILOCA and HaploDMF require custom Docker images. They are built automatically on first run when using `-profile docker`:

```bash
# Build manually if needed
docker build -t viloca:latest docker/viloca/
docker build -t haplodmf:latest docker/haplodmf/
```

---

## Test dataset

The repository includes a simulated ONT toy dataset for testing:

```bash
# Small overlaps (~100 bp)
nextflow run main.nf \
  --fastq_dir data/reads_artic_small_overlaps \
  --reference_fasta data/NC_045512_2.fasta \
  --run_haplodmf --run_pangolin \
  -profile docker

# Big overlaps (~400 bp)
nextflow run main.nf \
  --fastq_dir data/reads_artic_big_overlaps \
  --reference_fasta data/NC_045512_2.fasta \
  --run_haplodmf --run_pangolin \
  -profile docker
```

The test data contains 51 mixtures of SARS-CoV-2 lineages with precisely controlled compositions, simulated using NanoSim trained on MinION R10.4.1 data (ARTIC 5.3.2 primers).

---

## Troubleshooting

**Porechop_ABI out of memory on large files:**
Porechop_ABI loads the entire file into RAM for k-mer counting. For files >10 GB, increase memory in `nextflow.config`:
```groovy
withName: 'PORECHOP_ABI' {
    memory = '80 GB'
}
```
Do not use `--ab_initio` flag on large files — it requires significantly more RAM.

**Kraken2 slow on ONT data:**
Kraken2 processes each read independently and is not optimized for long reads. Expected throughput: ~50,000–100,000 reads/hour on a standard server. For large datasets (>1M reads), plan for 10–20 hours of runtime.

**Resume a failed pipeline:**
```bash
nextflow run main.nf [same parameters] -resume
```
Nextflow caches completed tasks. Only failed or new tasks will re-run.

---

## Methods

### Pipeline overview

Sequencing data analysis was performed using MetBio-WGSP, a Nextflow DSL2 pipeline (v24.04.4) developed for wastewater-based genomic surveillance of SARS-CoV-2. All bioinformatics tools were executed in isolated Docker containers to ensure reproducibility. The pipeline accepts Oxford Nanopore Technologies (ONT) long-read FASTQ files and processes them through sequential modules: optional read preprocessing, taxonomic profiling, variant-based lineage deconvolution, and haplotype-based lineage reconstruction.

### Read preprocessing (optional)

Raw ONT reads were optionally subjected to a preprocessing workflow. Adapter sequences were trimmed using Porechop_ABI v0.5.1, which performs automatic adapter detection via k-mer profiling. Quality and length filtering was performed with Filtlong v0.3.1 (minimum Phred quality score: 8, minimum read length: 200 bp). Read quality was assessed before and after preprocessing using NanoPlot v1.46.2. For amplicon-based data, primer sequences were trimmed from BAM alignments using iVar v1.4.4 with a BED file of primer coordinates.

### Taxonomic profiling

For metagenomic classification, reads were classified against the Kraken2 PlusPF database (standard genomes plus protozoa and fungi) using Kraken2 v2.17.1. Species-level abundance re-estimation was performed with Bracken v3.0.1. Taxonomic results were visualised as Sankey plots, stacked bar charts, heatmaps, and interactive Krona HTML reports.

### Read alignment and variant calling

Reads were aligned to the SARS-CoV-2 reference genome (NC_045512.2, Wuhan-Hu-1) using minimap2 v2.24 with the `map-ont` preset for ONT data. Alignment files were processed with samtools v1.14. Variant calling and lineage demixing were performed using Freyja v2.0.3. Briefly, Freyja calls variants and estimates per-sample lineage abundances by solving a constrained least-squares optimisation problem against a curated barcode database of lineage-defining mutations. Lineage barcode databases were updated prior to analysis using `freyja update`.

### Haplotype reconstruction

Viral haplotypes were reconstructed from aligned reads using HaploDMF v1.2.0, a deep matrix factorisation (DMF) approach that jointly learns haplotype sequences and their abundances from single-nucleotide variant (SNV) co-occurrence patterns. HaploDMF was run with default parameters (error rate: 0.1, minimum SNV sites: 20, clustering algorithm: Ward linkage). Samples with insufficient SNV sites were automatically skipped. As an alternative approach, haplotype reconstruction was also performed using VILOCA v1.2.0 (Variational Inference for LOngitudinal Cancer and viral Analysis), which employs a Bayesian statistical model for local haplotype reconstruction.

### Lineage classification

Reconstructed haplotype sequences from HaploDMF and VILOCA were assigned SARS-CoV-2 Pango lineages using Pangolin v4.4 (Phylogenetic Assignment of Named Global Outbreak LINeages). Pangolin was run in PUSHER mode, which uses a trained classifier for rapid assignment. Lineage assignments and confidence scores were aggregated across all haplotypes per sample into summary tables.

### Lineage quantification (lr-Kallisto)

As a complementary quantification approach, lineage-level abundances were estimated by pseudo-alignment of long reads against a multi-FASTA reference of lineage consensus sequences using lr-Kallisto (Kallisto v0.51.1 with long-read mode), which provides rapid probabilistic assignment without requiring full alignment.

### Software and reproducibility

The complete pipeline is implemented in Nextflow DSL2 (v24.04.4) and is publicly available at [https://github.com/nicolaedrabcinski/metbio_nextflow](https://github.com/nicolaedrabcinski/metbio_nextflow). All tools run as Docker containers. Tool versions: Nextflow 24.04.4, minimap2 2.24, samtools 1.14, Freyja 2.0.3, HaploDMF 1.2.0, VILOCA 1.2.0, Pangolin 4.4, Kallisto 0.51.1, Kraken2 2.17.1, Bracken 3.0.1, Porechop_ABI 0.5.1, Filtlong 0.3.1, NanoPlot 1.46.2, iVar 1.4.4.

---

## Contributors

- Victor Gordeev (Pipeline Design & Validation)
- Nicolae Drabcinski (Nextflow Implementation)

## Acknowledgment

The development of this software package is supported by the grant of the Ministry of Research, Innovation and Digitization, under Romania's National Recovery and Resilience Plan, funded by the European Union NextGenerationEU program no. 760286/27.03.2024, code 167/31.07.2023, within Pillar III, Component C9, Investment 8.

The pipeline was developed as part of the project "Metagenomics and Bioinformatics tools for Wastewater-based Genomic Surveillance of viral Pathogens for early prediction of public health risks (MetBio-WGSP)".
