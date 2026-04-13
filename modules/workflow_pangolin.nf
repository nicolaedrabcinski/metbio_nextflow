process PANGOLIN_CLASSIFY {
    tag "Pangolin ${sample_id}"
    publishDir "${params.final_outdir}/pangolin/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(haplotypes_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_lineage_report.csv"), emit: report

    script:
    """
    echo "Running Pangolin lineage classification for ${sample_id}..."

    pangolin ${haplotypes_fasta} \
        --outfile ${sample_id}_lineage_report.csv \
        --threads ${task.cpus} \
        --verbose

    echo "Pangolin completed for ${sample_id}"
    """
}

process PANGOLIN_AGGREGATE {
    tag "Aggregating Pangolin results"
    publishDir "${params.final_outdir}/pangolin", mode: 'copy'

    input:
    path reports  // collected *_lineage_report.csv files from PANGOLIN_CLASSIFY

    output:
    path("pangolin_summary.csv"), emit: summary

    script:
    """
    header_written=false

    for report in *_lineage_report.csv; do
        if [ -f "\$report" ]; then
            if [ "\$header_written" = false ]; then
                head -1 "\$report" > pangolin_summary.csv
                header_written=true
            fi
            tail -n +2 "\$report" >> pangolin_summary.csv
        fi
    done

    if [ "\$header_written" = false ]; then
        echo "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note" > pangolin_summary.csv
    fi
    """
}

process NEXTCLADE_CLASSIFY {
    tag "Nextclade ${sample_id}"
    publishDir "${params.final_outdir}/nextclade/${sample_id}", mode: 'copy'
    container 'nextstrain/nextclade:latest'

    input:
    tuple val(sample_id), path(haplotypes_fasta)
    path nextclade_dataset

    output:
    tuple val(sample_id), path("${sample_id}_nextclade.tsv"), emit: report

    script:
    """
    echo "Running Nextclade classification for ${sample_id}..."

    nextclade run \
        -D ${nextclade_dataset} \
        -O . \
        ${haplotypes_fasta} 2>&1

    mv nextclade.tsv ${sample_id}_nextclade.tsv

    echo "Nextclade completed for ${sample_id}"
    """
}

process NEXTCLADE_AGGREGATE {
    tag "Aggregating Nextclade results"
    publishDir "${params.final_outdir}/nextclade", mode: 'copy'

    input:
    tuple val(sample_ids), path(reports)  // collected [sample_id, *.tsv] tuples

    output:
    path("nextclade_summary.tsv"), emit: summary

    script:
    """
    header_written=false

    for report in *_nextclade.tsv; do
        if [ -f "\$report" ]; then
            sample=\$(basename "\$report" _nextclade.tsv)
            if [ "\$header_written" = false ]; then
                echo -e "sample\\t\$(head -1 \$report)" > nextclade_summary.tsv
                header_written=true
            fi
            tail -n +2 "\$report" | awk -v s="\$sample" 'BEGIN{OFS="\\t"} {print s, \$0}' >> nextclade_summary.tsv
        fi
    done

    if [ "\$header_written" = false ]; then
        echo -e "sample\\tindex\\tseqName\\tclade\\tclade_display\\tclade_who\\tclade_nextstrain\\tpartiallyAliased" > nextclade_summary.tsv
    fi
    """
}

workflow PANGOLIN_WORKFLOW {
    take:
    haplotypes_ch  // Channel of [sample_id, haplotypes.fasta]

    main:
    classified = PANGOLIN_CLASSIFY(haplotypes_ch)

    // Collect all report files through the channel — no filesystem race
    collected_reports = classified.report
        .map { sample_id, report -> report }
        .collect()

    aggregated = PANGOLIN_AGGREGATE(collected_reports)

    // Nextclade classification (if dataset provided)
    if (params.nextclade_dataset && file(params.nextclade_dataset).exists()) {
        nextclade_db = Channel.value(file(params.nextclade_dataset))
        nc_classified = NEXTCLADE_CLASSIFY(haplotypes_ch, nextclade_db)

        collected_nc = nc_classified.report.collect()
        NEXTCLADE_AGGREGATE(collected_nc)
    }

    emit:
    reports = classified.report
    summary = aggregated.summary
}
