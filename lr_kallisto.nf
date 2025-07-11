#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process LR_KALLISTO {
    // Конфигурация процесса
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq_file)
    path index_file
    path kallisto_path
    
    output:
    tuple val(sample_id), path("${sample_id}/abundance.tsv"), emit: abundance
    tuple val(sample_id), path("${sample_id}/run_info.json"), emit: run_info
    
    script:
    """
    echo "Processing long reads: ${fastq_file}"
    
    # Создать выходную директорию
    mkdir -p ${sample_id}
    
    # Запустить lr-kallisto quantification
    ${kallisto_path} quant \\
        -i ${index_file} \\
        -o ${sample_id} \\
        --single \\
        -l 900 \\
        -s 30 \\
        -t ${task.cpus} \\
        ${fastq_file}
    
    echo "Completed processing: ${sample_id}"
    """
}