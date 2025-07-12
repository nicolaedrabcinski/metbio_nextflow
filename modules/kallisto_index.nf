#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process KALLISTO_INDEX {
    tag "Building kallisto index"
    publishDir "${params.final_outdir}/index", mode: 'copy'
    
    input:
    path lineages_fasta
    path kallisto_path
    
    output:
    path "lineages.idx", emit: index
    
    script:
    """
    ${kallisto_path} index \\
        -i lineages.idx \\
        -k ${params.kmer_size} \\
        ${lineages_fasta}
    
    echo "✅ Kallisto index created with k-mer size: ${params.kmer_size}"
    """
}