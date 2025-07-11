#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Параметры по умолчанию
params.input = 'samples.csv'
params.fasta_dir = null
params.index = null
params.outdir = 'results'

// Логирование параметров
log.info """
========================================
MetBio Nextflow Pipeline
========================================
input        : ${params.input}
fasta_dir    : ${params.fasta_dir}
index        : ${params.index}
outdir       : ${params.outdir}
========================================
"""

// Импорт процессов
include { DOWNLOAD_KALLISTO } from './download_kallisto.nf'
include { KALLISTO_INDEX } from './kallisto_index.nf'
include { LR_KALLISTO } from './lr_kallisto.nf'

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
    
    // Читать входной CSV файл
    input_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            [row.sample, file(row.fastq)]
        }
    
    // Создать или использовать существующий индекс
    if (params.index) {
        // Использовать готовый индекс
        index_ch = Channel.value(file(params.index))
        log.info "Using existing index: ${params.index}"
    } else {
        // Создать новый индекс из FASTA файлов
        if (!params.fasta_dir) {
            error "Необходимо указать --fasta_dir или --index"
        }
        
        fasta_files = Channel
            .fromPath("${params.fasta_dir}/*.fasta")
            .collect()
        
        log.info "Creating index from FASTA files in: ${params.fasta_dir}"
        KALLISTO_INDEX(fasta_files, kallisto_path)
        index_ch = KALLISTO_INDEX.out.index
    }
    
    // Запустить lr-kallisto для каждого образца
    LR_KALLISTO(
        input_ch,
        index_ch,
        kallisto_path
    )
    
    // Вывести статус завершения
    LR_KALLISTO.out.abundance
        .subscribe { sample, abundance ->
            log.info "✅ Образец ${sample} обработан: ${abundance}"
        }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline завершен!
    ========================================
    Статус: ${workflow.success ? 'УСПЕШНО' : 'ОШИБКА'}
    Результаты: ${params.outdir}
    Время выполнения: ${workflow.duration}
    ========================================
    """
}