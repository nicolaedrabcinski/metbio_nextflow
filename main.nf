#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters by default
params.input = null
params.fastq_dir = null  // Новый параметр для автогенерации CSV
params.lineages_fasta = null
params.fasta_dir = null  // Для обратной совместимости
params.index = null
params.outdir = 'results'
params.help = false  // Исправлено: добавлена инициализация
params.filter_shuffled = false  // Изменено: по умолчанию использовать ВСЕ файлы
params.kmer_size = 31  // Новый параметр: размер k-mer для kallisto
params.create_plots = true  // Создавать ли графики результатов

// Создать папку результатов на основе имени входной папки
def getOutputDir() {
    if (params.fastq_dir) {
        // Извлечь имя папки из полного пути
        def fastq_dir_name = new File(params.fastq_dir).getName()
        return "${params.outdir}/${fastq_dir_name}"
    } else if (params.input) {
        // Если используется CSV, извлечь имя без расширения
        def input_name = new File(params.input).getName().replaceAll(/\.[^\.]+$/, '')
        return "${params.outdir}/${input_name}"
    } else {
        return params.outdir
    }
}

def final_outdir = getOutputDir()

// Передать final_outdir в params для использования в модулях
params.final_outdir = final_outdir

// Help message
def helpMessage() {
    log.info """
    ========================================
    MetBio Nextflow Pipeline - Usage
    ========================================
    
    Option 1 - Auto-generate CSV from FASTQ directory:
    nextflow run main.nf --fastq_dir /path/to/fastq/files --lineages_fasta lineages.fasta
    
    Option 2 - Use existing CSV:
    nextflow run main.nf --input samples.csv --lineages_fasta lineages.fasta
    
    Parameters:
    --fastq_dir         Directory containing FASTQ files (auto-generates CSV)
    --input             CSV file with sample information (manual mode)  
    --lineages_fasta    FASTA file with reference lineages
    --fasta_dir         Directory with individual FASTA files (legacy)
    --index             Pre-built kallisto index file
    --outdir            Output directory (default: results)
    --filter_shuffled   Only use files with '_shuffled' in name (default: false)
    --kmer_size         K-mer size for kallisto (default: 31)
    --create_plots      Create visualization plots (default: true)
    
    Note: Results will be saved in subdirectory named after input folder/file
    
    Examples:
    # Auto-scan ALL FASTQ files in directory
    nextflow run main.nf --fastq_dir data/reads_artic_big_overlaps --lineages_fasta lineages.fasta
    
    # Use custom k-mer size
    nextflow run main.nf --fastq_dir data/reads --lineages_fasta lineages.fasta --kmer_size 21
    
    # Process short overlaps data
    nextflow run main.nf --fastq_dir data/reads_artic_short_overlaps --lineages_fasta lineages.fasta
    
    # Use only _shuffled files
    nextflow run main.nf --fastq_dir data/reads --lineages_fasta lineages.fasta --filter_shuffled true
    
    # Use existing CSV
    nextflow run main.nf --input full.csv --lineages_fasta lineages.fasta  
    # Results in: results/full/
    ========================================
    """.stripIndent()
}

// Show help if no required parameters
if (params.help || (!params.fastq_dir && !params.input) || (!params.lineages_fasta && !params.fasta_dir && !params.index)) {
    helpMessage()
    exit 0
}

// Logging
log.info """
========================================
MetBio Nextflow Pipeline
========================================
fastq_dir      : ${params.fastq_dir}
input          : ${params.input}
lineages_fasta : ${params.lineages_fasta}
fasta_dir      : ${params.fasta_dir}
index          : ${params.index}
outdir         : ${final_outdir}
========================================
"""

// Импорт всех модулей
include { DOWNLOAD_KALLISTO } from './modules/download_kallisto.nf'
include { GENERATE_CSV_FROM_FASTQ_DIR } from './modules/generate_csv_kallisto.nf'
include { KALLISTO_INDEX } from './modules/kallisto_index.nf'
include { LR_KALLISTO } from './modules/lr_kallisto.nf'
include { CREATE_PLOTS } from './modules/create_plots.nf'

workflow {
    // Использовать существующий kallisto или скачать
    if (file('./tools/kallisto/kallisto').exists()) {
        kallisto_path = Channel.value(file('./tools/kallisto/kallisto'))
        log.info "Using existing kallisto: ./tools/kallisto/kallisto"
    } else {
        DOWNLOAD_KALLISTO()
        kallisto_path = DOWNLOAD_KALLISTO.out.kallisto_path
        log.info "Downloaded kallisto"
    }
    
    // Определить источник входных данных
    if (params.fastq_dir) {
        // Автогенерация CSV из папки с FASTQ
        log.info "Auto-generating CSV from FASTQ directory: ${params.fastq_dir}"
        fastq_directory = Channel.fromPath(params.fastq_dir, type: 'dir')
        csv_file = GENERATE_CSV_FROM_FASTQ_DIR(fastq_directory)
        
        // Читать сгенерированный CSV
        input_ch = csv_file.csv_file
            .splitCsv(header: true)
            .map { row -> 
                [row.sample, file(row.fastq)]
            }
    } else if (params.input) {
        // Использовать существующий CSV файл
        log.info "Using existing CSV file: ${params.input}"
        input_ch = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> 
                [row.sample, file(row.fastq)]
            }
    } else {
        error "You must specify either --fastq_dir or --input"
    }
    
    // Создать или использовать существующий индекс
    if (params.index) {
        // Использовать готовый индекс
        index_ch = Channel.value(file(params.index))
        log.info "Using existing index: ${params.index}"
    } else if (params.lineages_fasta) {
        // Использовать готовый lineages.fasta файл
        log.info "Creating index from lineages.fasta: ${params.lineages_fasta}"
        lineages_file = Channel.value(file(params.lineages_fasta))
        KALLISTO_INDEX(lineages_file, kallisto_path)
        index_ch = KALLISTO_INDEX.out.index
    } else if (params.fasta_dir) {
        // Создать новый индекс из отдельных FASTA файлов (старый способ)
        fasta_files = Channel
            .fromPath("${params.fasta_dir}/*.fasta")
            .collect()
        
        log.info "Creating index from FASTA files in: ${params.fasta_dir}"
        KALLISTO_INDEX(fasta_files, kallisto_path)
        index_ch = KALLISTO_INDEX.out.index
    } else {
        error "You must specify --lineages_fasta, --fasta_dir, or --index"
    }
    
    // Запустить lr-kallisto для каждого образца
    quantification_results = LR_KALLISTO(
        input_ch,
        index_ch,
        kallisto_path
    )
    
    // Создать графики и сводные данные
    if (params.create_plots) {
        // Собрать все результаты в один канал, используя уникальные имена
        results_with_unique_names = quantification_results.results
            .map { sample, dir -> 
                // Создать уникальное имя файла на основе пути к папке
                def unique_name = dir.toString().replaceAll(/.*\/work\//, '').replaceAll(/\//, '_')
                return tuple(unique_name, dir)
            }
            .collectFile(name: 'results_paths.txt', newLine: true) { unique_name, dir ->
                "${unique_name}\t${dir}"
            }
        
        CREATE_PLOTS(results_with_unique_names)
    }
    
    // Вывести статус завершения
    quantification_results.results
        .subscribe { sample, result_dir ->
            log.info "✅ Sample ${sample} processed: ${result_dir}"
        }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline is completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'ERROR'}
    Results: ${final_outdir}
    Execution time: ${workflow.duration}
    ========================================
    """
}