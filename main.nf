#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Derive output directory: use --final_outdir if explicitly set, otherwise auto-derive
def getOutputDir() {
    if (params.final_outdir) {
        return params.final_outdir
    } else if (params.fastq_dir) {
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

// Determine which tools to run based on provided references
def run_kallisto = (params.lineages_fasta || params.fasta_dir || params.index) as boolean
def run_freyja = params.reference_fasta as boolean
def run_kraken2 = (params.run_kraken2 && params.kraken2_db) as boolean
def run_ont_preprocess = params.run_ont_preprocess as boolean
def tools_list = [
    run_ont_preprocess ? 'ont-preprocess' : null,
    run_kallisto       ? 'kallisto'        : null,
    run_freyja         ? 'freyja'          : null,
    run_kraken2        ? 'kraken2+bracken' : null,
].findAll().join(', ')

// Help message
def helpMessage() {
    log.info """
    ========================================
    MetBio-WGSP Bioinformatics Pipeline - Usage
    ========================================

    Tools are selected automatically based on provided reference files:
    - Provide --lineages_fasta (or --fasta_dir / --index) to run Kallisto
    - Provide --reference_fasta to run Freyja
    - Provide both to run Kallisto and Freyja in parallel

    Option 1 - Auto-generate CSV from FASTQ directory:
    nextflow run main.nf --fastq_dir /path/to/fastq/files --lineages_fasta lineages.fasta
    nextflow run main.nf --fastq_dir /path/to/fastq/files --reference_fasta ref.fasta -profile docker

    Option 2 - Use existing CSV:
    nextflow run main.nf --input samples.csv --lineages_fasta lineages.fasta

    Parameters:
    --fastq_dir         Directory containing FASTQ files (auto-generates CSV)
    --input             CSV file with sample information (manual mode)
    --lineages_fasta    FASTA file with reference lineages (enables Kallisto)
    --reference_fasta   Reference genome FASTA (enables Freyja)
    --fasta_dir         Directory with individual FASTA files (enables Kallisto)
    --index             Pre-built kallisto index file (enables Kallisto)
    --outdir            Output directory (default: results)
    --threads           Number of threads (default: 32)

    Freyja-specific parameters:
    --refname           Reference name (default: NC_045512.2)
    --minq              Minimum base quality (default: 20)
    --freyja_eps        Epsilon for demix (default: 0.001)
    --covcut            Coverage cutoff (default: 10)
    --depthcutoff       Depth cutoff (default: 0)
    --confirmedonly     Use only confirmed lineages (default: false)
    --pathogen          Pathogen type (default: SARS-CoV-2)
    --mincov            Minimum coverage for aggregate (default: 60)
    --thresh            Threshold for plotting (default: 0)
    --run_viloca        Enable VILOCA local haplotype reconstruction (default: false)
    --run_haplodmf      Enable HaploDMF haplotype reconstruction (default: false)
    --run_pangolin      Enable Pangolin lineage classification on haplotypes (default: true)
    --run_kraken2       Enable Kraken2+Bracken taxonomic profiling (default: false)
    --kraken2_db        Path to Kraken2 database (required if --run_kraken2)
    --bracken_level     Taxonomic level for Bracken: D/P/C/O/F/G/S (default: S)
    --bracken_read_len  Read length for Bracken (default: 150)

    ONT preprocessing (--run_ont_preprocess):
      Pipeline: NanoPlot(raw) → Porechop_ABI → Filtlong → [Host depletion] → NanoPlot(clean)
    --run_ont_preprocess        Enable ONT preprocessing (default: false)
    --ont_min_quality           Filtlong: min Phred score (default: 8)
    --ont_min_length            Filtlong: min read length bp (default: 200)
    --ont_target_bases          Filtlong: target total bases, e.g. 500000000 (optional)
    --ont_skip_adapter_trim     Skip Porechop_ABI adapter trimming
    --ont_skip_quality_filter   Skip Filtlong quality/length filtering
    --ont_host_ref              Host reference FASTA for depletion (e.g. hg38.fa)
    --ont_skip_host_depletion   Skip host read removal (default: true)
    --run_ivar_trim             Enable iVar primer trimming after alignment (amplicon)
    --primer_bed                BED file with primer coordinates (required if --run_ivar_trim)

    Examples:
    # Kallisto only
    nextflow run main.nf --fastq_dir data/reads_artic_small_overlaps --lineages_fasta data/lineages.fasta

    # Freyja only
    nextflow run main.nf --fastq_dir data/reads_artic_small_overlaps --reference_fasta data/NC_045512_2.fasta -profile docker

    # Both tools in parallel
    nextflow run main.nf --fastq_dir data/reads_artic_small_overlaps --lineages_fasta data/lineages.fasta --reference_fasta data/NC_045512_2.fasta -profile mixed

    # Freyja with VILOCA
    nextflow run main.nf --fastq_dir data/reads_artic_small_overlaps --reference_fasta data/NC_045512_2.fasta --run_viloca -profile mixed

    # Kraken2+Bracken taxonomic profiling
    nextflow run main.nf --fastq_dir data/reads --run_kraken2 --kraken2_db /path/to/kraken2_db -profile docker
    ========================================
    """.stripIndent()
}

// Validation
if (params.help || (!run_kallisto && !run_freyja && !run_kraken2)) {
    helpMessage()
    exit 0
}

if (!params.fastq_dir && !params.input) {
    log.error "You must specify either --fastq_dir or --input"
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

if (run_kraken2 && !file(params.kraken2_db).exists()) {
    log.error "Kraken2 database not found: ${params.kraken2_db}"
    exit 1
}

// Logging
log.info """
========================================
MetBio-WGSP Bioinformatics Pipeline
========================================
tools          : ${tools_list}
fastq_dir      : ${params.fastq_dir ?: 'Not specified'}
input          : ${params.input ?: 'Not specified'}
lineages_fasta : ${params.lineages_fasta ?: 'Not specified'}
reference_fasta: ${params.reference_fasta ?: 'Not specified'}
kraken2_db     : ${params.kraken2_db ?: 'Not specified'}
outdir         : ${final_outdir}
threads        : ${params.threads}
========================================
"""

// Import modules
include { ONT_PREPROCESS; IVAR_TRIM } from './modules/workflow_ont_preprocess'
include { DOWNLOAD_KALLISTO } from './modules/download_kallisto.nf'
include { GENERATE_CSV_FROM_FASTQ_DIR } from './modules/generate_csv_kallisto.nf'
include { KALLISTO_INDEX } from './modules/kallisto_index.nf'
include { LR_KALLISTO } from './modules/lr_kallisto.nf'
include {
    FREYJA_ALIGN;
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
    FREYJA_AGGREGATE;
    FREYJA_PLOT;
    COVERAGE_PLOT
} from './modules/workflow_freyja'
include { VILOCA_WORKFLOW } from './modules/workflow_viloca.nf'
include { HAPLODMF_WORKFLOW } from './modules/workflow_haplodmf.nf'
include { PANGOLIN_WORKFLOW } from './modules/workflow_pangolin.nf'
include { KRAKEN2_WORKFLOW } from './modules/workflow_kraken2.nf'
include { CREATE_PLOTS } from './modules/create_plots'
include { GENERATE_RESULT_TABLE_FREYJA } from './modules/generate_result_table_freyja'
include { CREATE_WHO_PLOTS } from './modules/create_who_plots'

// Docker image management: check availability, pull external or build custom
def ensureDockerImage(String image, String dockerfileDir = null) {
    // Check if image already exists locally
    def check = ['docker', 'image', 'inspect', image].execute()
    check.waitFor()
    if (check.exitValue() == 0) {
        log.info "  ${image} - OK"
        return
    }
    // Image not found: build from Dockerfile or pull from registry
    if (dockerfileDir && file("${dockerfileDir}/Dockerfile").exists()) {
        log.info "  Building ${image} from ${dockerfileDir}/Dockerfile..."
        def build = ['docker', 'build', '-t', image, dockerfileDir].execute()
        build.consumeProcessOutput(System.out, System.err)
        build.waitFor()
        if (build.exitValue() != 0) {
            log.error "  Failed to build Docker image: ${image}"
            exit 1
        }
        log.info "  Built ${image} successfully"
    } else {
        log.info "  Pulling ${image}..."
        def pull = ['docker', 'pull', image].execute()
        pull.consumeProcessOutput(System.out, System.err)
        pull.waitFor()
        if (pull.exitValue() != 0) {
            log.warn "  Failed to pull ${image} - will try to use cached version"
        }
    }
}

workflow {
    // Ensure Docker images are available (pull external / build custom from docker/)
    if (workflow.profile.contains('docker')) {
        log.info "Checking Docker images..."
        if (run_freyja) {
            ensureDockerImage('quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0')
            ensureDockerImage('staphb/freyja:latest')
            if (params.run_viloca) {
                ensureDockerImage(params.viloca_container, 'docker/viloca')
            }
            if (params.run_haplodmf) {
                ensureDockerImage(params.haplodmf_container, 'docker/haplodmf')
            }
            if (params.run_pangolin) {
                ensureDockerImage(params.pangolin_container)
            }
        }
        if (run_kallisto) {
            ensureDockerImage('quay.io/biocontainers/kallisto:0.51.1--h2b92561_2')
        }
        if (run_kraken2) {
            ensureDockerImage('staphb/kraken2:latest')
            ensureDockerImage('staphb/bracken:latest')
            ensureDockerImage('nanozoo/krona:2.7.1--e7615f7')
        }
        if (run_ont_preprocess) {
            ensureDockerImage('staphb/nanoplot:latest')
            if (!params.ont_skip_adapter_trim)   ensureDockerImage('staphb/porechop:latest')
            if (!params.ont_skip_quality_filter) ensureDockerImage('staphb/chopper:latest')
            if (params.run_ivar_trim)            ensureDockerImage('staphb/ivar:latest')
        }
        ensureDockerImage('jupyter/scipy-notebook:latest')
        log.info "All Docker images ready"
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

    // ONT preprocessing (optional, runs before all other tools)
    if (run_ont_preprocess) {
        log.info "Running ONT preprocessing (NanoPlot + Porechop_ABI + Chopper)"
        ONT_PREPROCESS(input_ch)
        input_ch = ONT_PREPROCESS.out.reads
    }

    // Collect results from all tools
    all_results_ch = Channel.empty()

    // Kallisto workflow
    if (run_kallisto) {
        log.info "Running Kallisto workflow"

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

        kallisto_results_ch = quantification_results.results
        all_results_ch = all_results_ch.mix(kallisto_results_ch)

        // Create Kallisto plots
        if (params.create_plots) {
            results_with_unique_names = kallisto_results_ch
                .map { sample, dir ->
                    def unique_name = dir.toString().replaceAll(/.*\/work\//, '').replaceAll(/\//, '_')
                    return tuple(unique_name, dir)
                }
                .collectFile(name: 'results_paths.txt', newLine: true) { unique_name, dir ->
                    "${unique_name}\t${dir}"
                }

            CREATE_PLOTS(results_with_unique_names)
        }
    }

    // Freyja workflow
    if (run_freyja) {
        log.info "Running Freyja workflow"

        // Create file channel for reference
        reference_file = Channel.value(file(params.reference_fasta))

        // Alignment step
        aligned_ch = FREYJA_ALIGN(
            input_ch,
            reference_file
        )

        // Optional iVar primer trimming (amplicon data)
        if (params.run_ivar_trim && params.primer_bed) {
            primer_bed_ch = Channel.value(file(params.primer_bed))
            bam_for_variants = IVAR_TRIM(
                aligned_ch.aligned,
                primer_bed_ch
            ).trimmed_bam
        } else {
            bam_for_variants = aligned_ch.aligned
        }

        // Variant calling
        variants_ch = FREYJA_VARIANTS(
            bam_for_variants,
            reference_file
        )

        // Demixing
        demix_ch = FREYJA_DEMIX(
            variants_ch.variants_depths
        )

        // Aggregate results from all samples
        all_demix = demix_ch.results
            .map { sample, tsv -> tsv }
            .unique { it.name }
            .collect()

        aggregated_ch = FREYJA_AGGREGATE(all_demix)

        // Generate Freyja results table CSV after aggregation
        GENERATE_RESULT_TABLE_FREYJA(aggregated_ch.aggregated)

        if (params.create_plots) {
            // Create Freyja plots
            FREYJA_PLOT(aggregated_ch.aggregated)

            // Generate per-sample coverage plots (Panel C, Pavian-style)
            COVERAGE_PLOT(aligned_ch.aligned)

            // Create WHO classification plots
            CREATE_WHO_PLOTS(GENERATE_RESULT_TABLE_FREYJA.out.results_table)
        }

        freyja_results_ch = demix_ch.results
        all_results_ch = all_results_ch.mix(freyja_results_ch)

        // Run VILOCA if enabled
        if (params.run_viloca) {
            log.info "Running VILOCA local haplotype reconstruction..."
            viloca_results = VILOCA_WORKFLOW(
                aligned_ch.aligned,
                reference_file
            )
            all_results_ch = all_results_ch.mix(viloca_results.results)
        }

        // Run HaploDMF if enabled
        if (params.run_haplodmf) {
            log.info "Running HaploDMF haplotype reconstruction..."
            haplodmf_results = HAPLODMF_WORKFLOW(
                aligned_ch.aligned,
                reference_file
            )
            all_results_ch = all_results_ch.mix(haplodmf_results.results)

            // Run Pangolin lineage classification on HaploDMF haplotypes
            if (params.run_pangolin) {
                log.info "Running Pangolin lineage classification on HaploDMF haplotypes..."
                haplodmf_haplotypes_ch = haplodmf_results.results
                    .flatMap { sample_id, outdir ->
                        def fastas = file("${outdir}/*_haplotypes.fasta")
                        fastas ? fastas.collect { f -> [sample_id, f] } : []
                    }

                pangolin_results = PANGOLIN_WORKFLOW(haplodmf_haplotypes_ch)
                all_results_ch = all_results_ch.mix(pangolin_results.reports)
            }
        }
    }

    // Kraken2 + Bracken taxonomic profiling workflow
    if (run_kraken2) {
        log.info "Running Kraken2 + Bracken taxonomic profiling"

        kraken2_db_ch = Channel.value(file(params.kraken2_db))

        KRAKEN2_WORKFLOW(
            input_ch,
            kraken2_db_ch
        )
    }

    // Output status
    all_results_ch.subscribe { sample, result_dir ->
        log.info "Sample ${sample} processed: ${result_dir}"
    }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline execution completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'ERROR'}
    Tools: ${tools_list}
    Results: ${final_outdir}
    Execution time: ${workflow.duration}
    ========================================
    """
}
