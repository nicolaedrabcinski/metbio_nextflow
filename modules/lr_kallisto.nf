// modules/lr_kallisto.nf

process LR_KALLISTO {
    tag "${sample}"
    publishDir "${params.final_outdir}/quantification", mode: 'copy'
    
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