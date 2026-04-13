process VILOCA_RUN {
    tag "VILOCA ${sample_id}"
    publishDir "${params.final_outdir}/viloca/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_fasta
    
    // Emit a single, unique directory per sample with all outputs to avoid name collisions
    output:
    tuple val(sample_id), path("${sample_id}_viloca_outputs"), emit: results
    tuple val(sample_id), path("${sample_id}_viloca_outputs/${sample_id}_viloca.log"), emit: logs
    
    script:
    def threads_arg = params.threads > 0 ? "-t ${params.threads}" : "-t 0"
    """
    echo "Running VILOCA for ${sample_id}..."
    echo "BAM file: ${bam}"
    echo "Reference: ${reference_fasta}"
    
    # Check if viloca is available
    if ! command -v viloca &> /dev/null; then
        echo "ERROR: viloca not found"
        exit 1
    fi
    
    # Check input files
    if [ ! -f "${bam}" ]; then
        echo "ERROR: BAM file ${bam} not found!"
        exit 1
    fi
    
    if [ ! -f "${reference_fasta}" ]; then
        echo "ERROR: Reference file ${reference_fasta} not found!"
        exit 1
    fi
    
    # Run VILOCA
    echo "Running VILOCA local haplotype reconstruction..."
    viloca run \\
        -b ${bam} \\
        -f ${reference_fasta} \\
        --mode use_quality_scores \\
        -of csv \\
        ${threads_arg} \\
        > ${sample_id}_viloca.log 2>&1
    rc=\$?

    # Ensure output structure exists and collect outputs into a single, unique folder
    mkdir -p snv
    outdir=${sample_id}_viloca_outputs
    mkdir -p "\${outdir}"

    # Move CSVs if present
    if compgen -G "*.csv" > /dev/null; then
        mv -- *.csv "\${outdir}/" 2>/dev/null || true
    fi

    # Move snv directory if created
    if [ -d "snv" ]; then
        mv snv "\${outdir}/" 2>/dev/null || true
    fi

    # Move the process log into the outputs folder for collection (exists even on failure)
    mv ${sample_id}_viloca.log "\${outdir}/" 2>/dev/null || true

    # Check results
    if [ \$rc -eq 0 ]; then
        echo "VILOCA completed successfully for ${sample_id}"
        csv_count=\$(find "\${outdir}" -name "*.csv" -type f | wc -l)
        echo "Generated CSV files: \$csv_count"
        if [ -d "\${outdir}/snv" ]; then
            snv_count=\$(find "\${outdir}/snv/" -type f | wc -l)
            echo "SNV files: \$snv_count"
        fi
    else
        echo "ERROR: VILOCA failed for ${sample_id} (rc=\$rc)"
        # print limited log for debugging
        head -n 200 "\${outdir}/${sample_id}_viloca.log" || true
        exit \$rc
    fi
    """
}

process VILOCA_AGGREGATE {
    tag "Aggregating VILOCA results"
    publishDir "${params.final_outdir}/viloca", mode: 'copy'

    input:
    path result_dirs  // collected *_viloca_outputs directories from VILOCA_RUN

    output:
    path("viloca_summary.csv"), emit: summary
    path("viloca_aggregated.log"), emit: log

    script:
    """
    echo "Aggregating VILOCA results..." > viloca_aggregated.log
    echo "sample,haplotype_count,snv_count,total_coverage" > viloca_summary.csv

    for result_dir in *_viloca_outputs; do
        if [ -d "\$result_dir" ]; then
            sample=\$(basename "\$result_dir" | sed 's/_viloca_outputs\$//')
            snv_count=0
            haplotype_count=0

            if [ -f "\$result_dir/snv/SNVs_0.010000_final.csv" ]; then
                snv_count=\$(tail -n +2 "\$result_dir/snv/SNVs_0.010000_final.csv" | wc -l)
            fi

            haplotype_count=\$(find "\$result_dir" -name "w-*.csv" -type f | wc -l)

            echo "\$sample,\$haplotype_count,\$snv_count,0" >> viloca_summary.csv
            echo "Processed \$sample: haplotypes=\$haplotype_count, SNVs=\$snv_count" >> viloca_aggregated.log
        fi
    done

    echo "Aggregation complete!" >> viloca_aggregated.log
    """
}

workflow VILOCA_WORKFLOW {
    take:
    aligned_ch    // Channel of [sample_id, bam, bai]
    reference_file

    main:
    viloca_results = VILOCA_RUN(aligned_ch, reference_file)

    // Collect all output directories through the channel — no filesystem race
    collected_dirs = viloca_results.results
        .map { sample_id, outdir -> outdir }
        .collect()

    aggregated = VILOCA_AGGREGATE(collected_dirs)

    emit:
    results = viloca_results.results
    logs    = viloca_results.logs
    summary = aggregated.summary
}
