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
    
    Note: Results will be saved in subdirectory named after input folder/file
    
    Examples:
    # Auto-scan FASTQ directory
    nextflow run main.nf --fastq_dir data/reads_artic_big_overlaps --lineages_fasta lineages.fasta
    # Results in: results/reads_artic_big_overlaps/
    
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

// Импорт процессов
include { DOWNLOAD_KALLISTO } from './modules/download_kallisto.nf'
include { KALLISTO_INDEX } from './modules/kallisto_index.nf'
include { LR_KALLISTO } from './modules/lr_kallisto.nf'

// Процесс для автогенерации CSV из папки с FASTQ
process GENERATE_CSV_FROM_FASTQ_DIR {
    tag "Auto-generate CSV"
    publishDir "${final_outdir}/metadata", mode: 'copy'
    
    input:
    path fastq_dir
    
    output:
    path "auto_samples.csv"
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    import csv
    from pathlib import Path
    
    fastq_files = []
    
    # Поиск FASTQ файлов с различными расширениями
    extensions = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
    
    # Рекурсивный поиск во всех подпапках
    for ext in extensions:
        pattern = os.path.join('${fastq_dir}', '**', ext)
        fastq_files.extend(glob.glob(pattern, recursive=True))
    
    # Сортировка для стабильного порядка
    fastq_files.sort()
    
    print(f"Found {len(fastq_files)} FASTQ files:")
    for f in fastq_files:
        print(f"  - {f}")
    
    # Создание CSV файла
    with open('auto_samples.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'fastq'])  # Заголовки как в оригинальном CSV
        
        for fastq_file in fastq_files:
            # Извлечение имени образца из имени файла
            sample_name = Path(fastq_file).stem
            
            # Удаление распространенных суффиксов
            suffixes_to_remove = ['.fastq', '.fq', '.gz', '_long', '_short']
            for suffix in suffixes_to_remove:
                if sample_name.endswith(suffix):
                    sample_name = sample_name[:-len(suffix)]
            
            # Использование абсолютного пути к файлу
            abs_path = os.path.abspath(fastq_file)
            writer.writerow([sample_name, abs_path])
    
    print(f"Generated CSV file with {len(fastq_files)} samples")
    """
}

// Локальные процессы с правильными путями
process KALLISTO_INDEX_LOCAL {
    tag "Building kallisto index"
    publishDir "${final_outdir}/index", mode: 'copy'
    
    input:
    path lineages_fasta
    path kallisto_path
    
    output:
    path "lineages.idx"
    
    script:
    """
    ${kallisto_path} index \\
        -i lineages.idx \\
        -k 31 \\
        ${lineages_fasta}
    """
}

process LR_KALLISTO_LOCAL {
    tag "${sample}"
    publishDir "${final_outdir}/quantification", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq_file)
    path index_file
    path kallisto_path
    
    output:
    tuple val(sample), path("${sample}")
    
    script:
    """
    mkdir -p ${sample}
    
    ${kallisto_path} quant \\
        -i ${index_file} \\
        -o ${sample} \\
        --single \\
        -l 1000 \\
        -s 30 \\
        ${fastq_file}
    """
}

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
        input_ch = csv_file
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
        KALLISTO_INDEX_LOCAL(lineages_file, kallisto_path)
        index_ch = KALLISTO_INDEX_LOCAL.out
    } else if (params.fasta_dir) {
        // Создать новый индекс из отдельных FASTA файлов (старый способ)
        fasta_files = Channel
            .fromPath("${params.fasta_dir}/*.fasta")
            .collect()
        
        log.info "Creating index from FASTA files in: ${params.fasta_dir}"
        KALLISTO_INDEX_LOCAL(fasta_files, kallisto_path)
        index_ch = KALLISTO_INDEX_LOCAL.out
    } else {
        error "You must specify --lineages_fasta, --fasta_dir, or --index"
    }
    
    // Запустить lr-kallisto для каждого образца
    LR_KALLISTO_LOCAL(
        input_ch,
        index_ch,
        kallisto_path
    )
    
    // Вывести статус завершения
    LR_KALLISTO_LOCAL.out
        .subscribe { sample, result_dir ->
            log.info "Sample ${sample} processed: ${result_dir}"
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