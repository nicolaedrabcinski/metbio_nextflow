#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters by default
params.tool = null                    // NEW: kallisto or freyja
params.input = null
params.fastq_dir = null
params.lineages_fasta = null
params.reference_fasta = null         // NEW: For Freyja
params.fasta_dir = null
params.index = null
params.outdir = 'results'
params.help = false
params.filter_shuffled = false
params.kmer_size = 31
params.fragment_length = null         // NEW: For kallisto
params.sd = null                      // NEW: For kallisto
params.threads = 1                    // NEW: Number of threads
params.freyja_eps = "0.00000001"      // NEW: Freyja epsilon
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
    --filter_shuffled   Only use files with '_shuffled' in name (default: false)
    --kmer_size         K-mer size for kallisto (default: 31)
    --fragment_length   Fragment length for kallisto single-end mode
    --sd                Fragment standard deviation for kallisto
    --threads           Number of threads (default: 1)
    --freyja_eps        Freyja epsilon parameter (default: 0.00000001)
    --create_plots      Create visualization plots (default: true)
    
    Examples:
    # Kallisto workflow
    nextflow run main.nf --tool kallisto --fastq_dir data/reads_artic_small_overlaps --lineages_fasta data/lineages.fasta
    
    # Freyja workflow  
    nextflow run main.nf --tool freyja --fastq_dir data/reads_artic_small_overlaps --reference_fasta data/NC_045512_2.fasta
    
    # Kallisto with custom parameters
    nextflow run main.nf --tool kallisto --fastq_dir data/reads --lineages_fasta lineages.fasta --kmer_size 21 --fragment_length 1000 --sd 30
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
kmer_size      : ${params.kmer_size}
fragment_length: ${params.fragment_length ?: 'Not specified'}
fragment_sd    : ${params.sd ?: 'Not specified'}
threads        : ${params.threads}
========================================
"""

// Import modules
include { DOWNLOAD_KALLISTO } from './modules/download_kallisto.nf'
include { GENERATE_CSV_FROM_FASTQ_DIR } from './modules/generate_csv_kallisto.nf'
include { KALLISTO_INDEX } from './modules/kallisto_index.nf'
include { LR_KALLISTO } from './modules/lr_kallisto.nf'
include { INSTALL_FREYJA } from './modules/install_freyja.nf'
include { INSTALL_MINIMAP2 } from './modules/install_minimap2.nf'
include { FREYJA_ALIGN; FREYJA_VARIANTS; FREYJA_DEMIX } from './modules/workflow_freyja.nf'
include { CREATE_PLOTS } from './modules/create_plots.nf'

workflow {
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
        // Kallisto workflow (existing working code)
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
        // Freyja workflow (new)
        log.info "Running Freyja workflow"
        
        INSTALL_FREYJA()
        INSTALL_MINIMAP2()
        
        // Alignment step
        FREYJA_ALIGN(
            input_ch,
            INSTALL_MINIMAP2.out.minimap2_bin,
            params.reference_fasta
        )
        
        // Variant calling
        FREYJA_VARIANTS(
            FREYJA_ALIGN.out.bam,
            INSTALL_FREYJA.out.freyja_bin,
            params.reference_fasta
        )
        
        // Demixing
        FREYJA_DEMIX(
            FREYJA_VARIANTS.out.variants_depths,
            INSTALL_FREYJA.out.freyja_bin
        )
        
        results_ch = FREYJA_DEMIX.out.results
    }
    
    // Create plots (works for both workflows)
    if (params.create_plots) {
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
    
    // Output status
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