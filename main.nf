#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters by default
params.tool = null
params.input = null
params.fastq_dir = null
params.lineages_fasta = null
params.reference_fasta = null
params.fasta_dir = null
params.index = null
params.outdir = 'results'
params.help = false
params.filter_shuffled = false
params.kmer_size = 31
params.fragment_length = null
params.sd = null
params.threads = 1
params.create_plots = true

// Создать папку результатов на основе имени входной папки
def getOutputDir() {
    if (params.fastq_dir) {
        def fastq_dir_name = new File(params.fastq_dir).getName()
        return "${params.outdir}/${fastq_dir_name}"
    } else if (params.input) {
        def input_name = new File(params.input).getName().replaceAll(/\.[^\.]+$/, '')
        return "${params.outdir}/${input_name}"
    } else {
        return params.outdir
    }
}

def final_outdir = getOutputDir()
params.final_outdir = final_outdir

// Help message
def helpMessage() {
    log.info """
    ========================================
    MetBio-WGSP Bioinformatics Pipeline - Usage
    ========================================
    
    Option 1 - Auto-generate CSV from FASTQ directory:
    nextflow run main.nf --tool kallisto --fastq_dir /path/to/fastq/files --lineages_fasta lineages.fasta
    nextflow run main.nf --tool freyja --fastq_dir /path/to/fastq/files --reference_fasta ref.fasta
    
    Option 2 - Use existing CSV:
    nextflow run main.nf --tool kallisto --input samples.csv --lineages_fasta lineages.fasta
    
    Parameters:
    --tool              Tool to use: "kallisto" or "freyja" (REQUIRED)
    --fastq_dir         Directory containing FASTQ files (auto-generates CSV)
    --input             CSV file with sample information (manual mode)  
    --lineages_fasta    FASTA file with reference lineages (for kallisto)
    --reference_fasta   Reference genome FASTA (for freyja)
    --fasta_dir         Directory with individual FASTA files (legacy)
    --index             Pre-built kallisto index file
    --outdir            Output directory (default: results)
    --threads           Number of threads (default: 1)
    
    Freyja-specific parameters:
    --refname           Reference name (default: NC_045512.2)
    --minq              Minimum base quality (default: 20)
    --freyja_eps        Epsilon for demix (default: 0.00000001)
    --covcut            Coverage cutoff (default: 10)
    --depthcutoff       Depth cutoff (default: 0)
    --confirmedonly     Use only confirmed lineages (default: false)
    --pathogen          Pathogen type (default: SC2)
    --mincov            Minimum coverage for aggregate (default: 60)
    --thresh            Threshold for plotting (default: 0.01)
    
    Examples:
    # Kallisto workflow
    nextflow run main.nf --tool kallisto --fastq_dir data/reads_artic_small_overlaps --lineages_fasta data/lineages.fasta
    
    # Freyja workflow  
    nextflow run main.nf --tool freyja --fastq_dir data/reads_artic_small_overlaps --reference_fasta data/NC_045512_2.fasta -profile docker
    
    # Freyja with custom parameters
    nextflow run main.nf --tool freyja --fastq_dir data/reads --reference_fasta ref.fasta --covcut 5 --mincov 50 --threads 8 -profile docker
    ========================================
    """.stripIndent()
}

// Validation
if (params.help || !params.tool) {
    helpMessage()
    exit 0
}

if (params.tool != "kallisto" && params.tool != "freyja") {
    log.error "Invalid tool '${params.tool}'. Must be 'kallisto' or 'freyja'"
    exit 1
}

if (!params.fastq_dir && !params.input) {
    log.error "You must specify either --fastq_dir or --input"
    exit 1
}

if (params.tool == "kallisto" && !params.lineages_fasta && !params.fasta_dir && !params.index) {
    log.error "Kallisto requires --lineages_fasta, --fasta_dir, or --index"
    exit 1
}

if (params.tool == "freyja" && !params.reference_fasta) {
    log.error "Freyja requires --reference_fasta"
    exit 1
}

// Validate file existence
if (params.reference_fasta && !file(params.reference_fasta).exists()) {
    log.error "Reference FASTA file not found: ${params.reference_fasta}"
    exit 1
}

if (params.lineages_fasta && !file(params.lineages_fasta).exists()) {
    log.error "Lineages FASTA file not found: ${params.lineages_fasta}"
    exit 1
}

// Logging
log.info """
========================================
MetBio-WGSP Bioinformatics Pipeline
========================================
tool           : ${params.tool}
fastq_dir      : ${params.fastq_dir ?: 'Not specified'}
input          : ${params.input ?: 'Not specified'}
lineages_fasta : ${params.lineages_fasta ?: 'Not specified'}
reference_fasta: ${params.reference_fasta ?: 'Not specified'}
outdir         : ${final_outdir}
threads        : ${params.threads}
========================================
"""

// Import modules
include { DOWNLOAD_KALLISTO } from './modules/download_kallisto.nf'
include { GENERATE_CSV_FROM_FASTQ_DIR } from './modules/generate_csv_kallisto.nf'
include { KALLISTO_INDEX } from './modules/kallisto_index.nf'
include { LR_KALLISTO } from './modules/lr_kallisto.nf'
include { FREYJA_ALIGN; FREYJA_VARIANTS; FREYJA_DEMIX; FREYJA_UPDATE; FREYJA_AGGREGATE; FREYJA_PLOT } from './modules/workflow_freyja.nf'
include { CREATE_PLOTS } from './modules/create_plots.nf'

// Helper process to pre-pull Docker images
process PULL_DOCKER_IMAGES {
    tag "Pre-pulling Docker images"
    
    output:
    path "docker_ready.flag", emit: flag
    
    when:
    workflow.profile.contains('docker')
    
    script:
    """
    echo "🐳 Pre-pulling required Docker images for reproducibility..."
    
    if [ "${params.tool}" == "freyja" ]; then
        echo "📦 Pulling Freyja Docker image..."
        docker pull staphb/freyja:latest || {
            echo "⚠️  Warning: Failed to pull staphb/freyja:latest"
            echo "Pipeline will attempt to use cached image or pull during execution"
        }
        
        echo "📦 Pulling minimap2 Docker image..."
        docker pull quay.io/biocontainers/minimap2:2.28--he4a0461_2 || {
            echo "⚠️  Warning: Failed to pull minimap2 image"
            echo "Pipeline will attempt to use cached image or pull during execution"
        }
    elif [ "${params.tool}" == "kallisto" ]; then
        echo "📦 Pulling Kallisto Docker image..."
        docker pull quay.io/biocontainers/kallisto:0.51.1--h6ad9e81_1 || {
            echo "⚠️  Warning: Failed to pull kallisto image"
            echo "Pipeline will attempt to use cached image or pull during execution"
        }
    fi
    
    touch docker_ready.flag
    echo "✅ Docker images ready"
    """
}

workflow {
    // Pre-pull Docker images if using docker profile (for reproducibility)
    if (workflow.profile.contains('docker')) {
        log.info "🐳 Docker profile detected - pre-pulling images for reproducibility..."
        PULL_DOCKER_IMAGES()
    }
    
    // Prepare input samples
    if (params.fastq_dir) {
        log.info "Auto-generating CSV from FASTQ directory: ${params.fastq_dir}"
        fastq_directory = Channel.fromPath(params.fastq_dir, type: 'dir')
        csv_file = GENERATE_CSV_FROM_FASTQ_DIR(fastq_directory)
        
        input_ch = csv_file.csv_file
            .splitCsv(header: true)
            .map { row -> [row.sample, file(row.fastq)] }
    } else if (params.input) {
        log.info "Using existing CSV file: ${params.input}"
        input_ch = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [row.sample, file(row.fastq)] }
    }

    // Branch workflow based on tool selection
    if (params.tool == "kallisto") {
        // Kallisto workflow
        if (file('./tools/kallisto/kallisto').exists()) {
            kallisto_path = Channel.value(file('./tools/kallisto/kallisto'))
            log.info "Using existing kallisto: ./tools/kallisto/kallisto"
        } else {
            DOWNLOAD_KALLISTO()
            kallisto_path = DOWNLOAD_KALLISTO.out.kallisto_path
            log.info "Downloaded kallisto"
        }
        
        // Create or use existing index
        if (params.index) {
            index_ch = Channel.value(file(params.index))
            log.info "Using existing index: ${params.index}"
        } else if (params.lineages_fasta) {
            log.info "Creating index from lineages.fasta: ${params.lineages_fasta}"
            lineages_file = Channel.value(file(params.lineages_fasta))
            KALLISTO_INDEX(lineages_file, kallisto_path)
            index_ch = KALLISTO_INDEX.out.index
        } else if (params.fasta_dir) {
            fasta_files = Channel
                .fromPath("${params.fasta_dir}/*.fasta")
                .collect()
            
            log.info "Creating index from FASTA files in: ${params.fasta_dir}"
            KALLISTO_INDEX(fasta_files, kallisto_path)
            index_ch = KALLISTO_INDEX.out.index
        }
        
        // Run lr-kallisto
        quantification_results = LR_KALLISTO(
            input_ch,
            index_ch,
            kallisto_path
        )
        
        results_ch = quantification_results.results
        
    } else if (params.tool == "freyja") {
        // Freyja workflow
        log.info "Running Freyja workflow"
        
        // Update Freyja database first (optional but recommended)
        // FREYJA_UPDATE()
        
        // Create file channel for reference
        reference_file = Channel.value(file(params.reference_fasta))
        
        // Alignment step (uses minimap2 from Docker)
        aligned_ch = FREYJA_ALIGN(
            input_ch,
            reference_file
        )
        
        // Variant calling (uses freyja from Docker)
        variants_ch = FREYJA_VARIANTS(
            aligned_ch.aligned,
            reference_file
        )
        
        // Demixing (uses freyja from Docker)
        demix_ch = FREYJA_DEMIX(
            variants_ch.variants_depths
        )
        
        // Aggregate results from all samples
        all_demix = demix_ch.results
            .map { sample, tsv -> tsv }
            .unique { it.name }
            .collect()
        
        aggregated_ch = FREYJA_AGGREGATE(all_demix)
        
        // Create Freyja plots
        FREYJA_PLOT(aggregated_ch.aggregated)
        
        results_ch = demix_ch.results
    }
    
    // Create plots (works for both workflows, but mainly for Kallisto)
    if (params.create_plots && params.tool == "kallisto") {
        results_with_unique_names = results_ch
            .map { sample, dir -> 
                def unique_name = dir.toString().replaceAll(/.*\/work\//, '').replaceAll(/\//, '_')
                return tuple(unique_name, dir)
            }
            .collectFile(name: 'results_paths.txt', newLine: true) { unique_name, dir ->
                "${unique_name}\t${dir}"
            }
        
        CREATE_PLOTS(results_with_unique_names)
    }
    
    // Output status.
    results_ch.subscribe { sample, result_dir ->
        log.info "✅ Sample ${sample} processed: ${result_dir}"
    }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline execution completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'ERROR'}
    Tool used: ${params.tool}
    Results: ${final_outdir}
    Execution time: ${workflow.duration}
    ========================================
    """
}