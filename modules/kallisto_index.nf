#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process KALLISTO_INDEX {
    publishDir "${params.outdir}/index", mode: 'copy'
    
    input:
    path lineages_fasta
    path kallisto_bin
    
    output:
    path "transcripts.kidx", emit: index
    
    script:
    """
    echo "=== Создание kallisto индекса ==="
    echo "Используем готовый файл: ${lineages_fasta}"
    echo ""
    
    echo "Проверка содержимого lineages.fasta:"
    echo "Размер файла: \$(wc -c < ${lineages_fasta}) байт"
    echo "Общее количество заголовков: \$(grep -c '^>' ${lineages_fasta})"
    echo ""
    echo "Все заголовки в файле:"
    grep '^>' ${lineages_fasta}
    echo ""
    
    # Создать индекс напрямую из lineages.fasta
    echo "Создание индекса kallisto..."
    ${kallisto_bin} index \\
        -i transcripts.kidx \\
        ${lineages_fasta}
    
    echo ""
    echo "=== Проверка созданного индекса ==="
    ls -la transcripts.kidx
    echo "Индекс создан успешно!"
    """
}