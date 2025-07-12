#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process LR_KALLISTO {
    tag "${sample}"
    publishDir "${params.final_outdir}/quantification", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq_file)
    path index_file
    path kallisto_path
    
    output:
    tuple val(sample), path("${sample}"), emit: results
    
    script:
    """
    mkdir -p ${sample}
    
    # Проверить, существует ли файл
    if [ ! -f "${fastq_file}" ]; then
        echo "❌ Error: FASTQ file ${fastq_file} not found"
        exit 1
    fi
    
    echo "🔬 Processing sample: ${sample}"
    echo "📄 Input file: ${fastq_file}"
    echo "📊 Index file: ${index_file}"
    
    ${kallisto_path} quant \\
        -i ${index_file} \\
        -o ${sample} \\
        --single \\
        -l 1000 \\
        -s 30 \\
        ${fastq_file}
    
    # Проверить, что результат создан
    if [ -f "${sample}/abundance.tsv" ]; then
        echo "✅ Sample ${sample} processed successfully"
        echo "📈 Results saved in ${sample}/"
    else
        echo "❌ Error: abundance.tsv not created for ${sample}"
        exit 1
    fi
    """
}