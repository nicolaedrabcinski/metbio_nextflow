process HAPLODMF_RUN {
    tag "HaploDMF ${sample_id}"
    publishDir "${params.final_outdir}/haplodmf/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_fasta

    output:
    tuple val(sample_id), path("${sample_id}_haplodmf_outputs"), emit: results
    tuple val(sample_id), path("${sample_id}_haplodmf_outputs/${sample_id}_haplodmf.log"), emit: logs

    script:
    // Build optional arguments
    def opt_args = ""
    opt_args += " -t ${task.cpus}"
    opt_args += " -mq ${params.haplodmf_map_qual}"
    if (params.haplodmf_start_pos) opt_args += " -sp ${params.haplodmf_start_pos}"
    if (params.haplodmf_end_pos)   opt_args += " -ep ${params.haplodmf_end_pos}"
    opt_args += " -a ${params.haplodmf_abundance}"
    opt_args += " -e ${params.haplodmf_error_rate}"
    opt_args += " -s ${params.haplodmf_signi_level}"
    opt_args += " -c ${params.haplodmf_cond_pro}"
    opt_args += " -f ${params.haplodmf_fre_snv}"
    opt_args += " -ss ${params.haplodmf_smallest_snv}"
    opt_args += " -w ${params.haplodmf_weight}"
    opt_args += " -lr ${params.haplodmf_learning_rate}"
    opt_args += " -epo ${params.haplodmf_epoch}"
    opt_args += " -bs ${params.haplodmf_batch_size}"
    opt_args += " -al ${params.haplodmf_algorithm}"
    opt_args += " -d ${params.haplodmf_depth}"
    opt_args += " -c1 ${params.haplodmf_cluster_thres1}"
    opt_args += " -c2 ${params.haplodmf_cluster_thres2}"
    opt_args += " -li ${params.haplodmf_largest_iter}"
    if (params.haplodmf_input_snv) opt_args += " -is ${params.haplodmf_input_snv}"
    """
    echo "Running HaploDMF for ${sample_id}..."
    echo "BAM file: ${bam}"
    echo "Reference: ${reference_fasta}"

    # Downsample BAM to max 50000 reads to prevent OOM in HaploDMF matrix factorization
    TOTAL_READS=\$(samtools view -c -F 4 ${bam})
    if [ "\${TOTAL_READS}" -gt "50000" ]; then
        FRAC=\$(python3 -c "print(round(50000/\${TOTAL_READS}, 4))")
        samtools view -h -s "\${FRAC}" -o ${sample_id}.sam ${bam}
        echo "Downsampled \${TOTAL_READS} -> ~50000 reads (fraction \${FRAC})"
    else
        samtools view -h -o ${sample_id}.sam ${bam}
        echo "Using all \${TOTAL_READS} reads (below threshold)"
    fi

    # Save absolute paths before changing directory (haplodmf.sh uses relative ./src/ paths)
    WORK_DIR=\$(pwd)
    SAM_ABS=\$(realpath ${sample_id}.sam)
    REF_ABS=\$(realpath ${reference_fasta})

    # Run HaploDMF from its install directory so ./src/ paths resolve
    # Note: haplodmf.sh may return non-zero even on success, so we check output files instead
    cd /opt/haplodmf
    bash haplodmf.sh \
        -i \${SAM_ABS} \
        -r \${REF_ABS} \
        -o \${WORK_DIR}/${sample_id}_result \
        -p ${sample_id} \
        ${opt_args} \
        > \${WORK_DIR}/${sample_id}_haplodmf.log 2>&1 || true
    cd \${WORK_DIR}

    # Collect outputs into a single directory
    outdir=${sample_id}_haplodmf_outputs
    mkdir -p "\${outdir}"

    # Move haplotype FASTA
    if compgen -G "${sample_id}_result/*_haplotypes.fasta" > /dev/null; then
        cp ${sample_id}_result/*_haplotypes.fasta "\${outdir}/"
    fi

    # Move SNV file
    if compgen -G "${sample_id}_result/*_SNV.txt" > /dev/null; then
        cp ${sample_id}_result/*_SNV.txt "\${outdir}/"
    fi

    # Move ACGT counts
    if compgen -G "${sample_id}_result/*_acgt.txt" > /dev/null; then
        cp ${sample_id}_result/*_acgt.txt "\${outdir}/"
    fi

    # Move log
    mv ${sample_id}_haplodmf.log "\${outdir}/" 2>/dev/null || true

    # Clean up intermediate SAM
    rm -f ${sample_id}.sam

    # Check success by output files (haplodmf.sh exit code is unreliable)
    if compgen -G "\${outdir}/*_haplotypes.fasta" > /dev/null; then
        hap_count=\$(grep -c "^>" "\${outdir}"/*_haplotypes.fasta 2>/dev/null || echo 0)
        echo "HaploDMF completed successfully for ${sample_id}"
        echo "Reconstructed haplotypes: \$hap_count"
    else
        echo "ERROR: HaploDMF failed for ${sample_id} - no haplotypes produced"
        head -n 200 "\${outdir}/${sample_id}_haplodmf.log" 2>/dev/null || true
        exit 1
    fi
    """
}

process HAPLODMF_AGGREGATE {
    tag "Aggregating HaploDMF results"
    publishDir "${params.final_outdir}/haplodmf", mode: 'copy'

    input:
    path result_dirs  // collected *_haplodmf_outputs directories from HAPLODMF_RUN

    output:
    path("haplodmf_summary.csv"), emit: summary

    script:
    """
    echo "sample,haplotype_count,haplotypes_info" > haplodmf_summary.csv

    for result_dir in *_haplodmf_outputs; do
        if [ -d "\$result_dir" ]; then
            sample=\$(basename "\$result_dir" | sed 's/_haplodmf_outputs\$//')
            hap_count=0
            hap_info=""

            hap_file=\$(find "\$result_dir" -name "*_haplotypes.fasta" -type f | head -1)
            if [ -n "\$hap_file" ] && [ -f "\$hap_file" ]; then
                hap_count=\$(grep -c "^>" "\$hap_file" 2>/dev/null || echo 0)
                hap_info=\$(grep "^>" "\$hap_file" | tr '\\n' ';' | sed 's/;\$//')
            fi

            echo "\$sample,\$hap_count,\"\$hap_info\"" >> haplodmf_summary.csv
        fi
    done
    """
}

workflow HAPLODMF_WORKFLOW {
    take:
    aligned_ch      // Channel of [sample_id, bam, bai]
    reference_file

    main:
    haplodmf_results = HAPLODMF_RUN(aligned_ch, reference_file)

    // Collect all output directories through the channel — no filesystem race
    collected_dirs = haplodmf_results.results
        .map { sample_id, outdir -> outdir }
        .collect()

    aggregated = HAPLODMF_AGGREGATE(collected_dirs)

    emit:
    results = haplodmf_results.results
    logs    = haplodmf_results.logs
    summary = aggregated.summary
}
