// kallisto_index.nf - Создание индекса kallisto

process KALLISTO_INDEX {
    publishDir "${params.outdir}/index", mode: 'copy'
    
    input:
    path fasta_files
    path kallisto_bin
    
    output:
    path "transcripts.kidx", emit: index
    
    script:
    """
    echo "=== Создание kallisto индекса ==="
    echo "FASTA файлы: $fasta_files"
    
    # Объединяем все FASTA файлы в один multi-fasta как просил Виктор
    cat $fasta_files > combined_transcripts.fasta
    
    # Создаем индекс как в команде Виктора
    ${kallisto_bin} index \\
        --index transcripts.kidx \\
        combined_transcripts.fasta
    
    echo "Индекс создан: transcripts.kidx"
    """
}