process FREYJA_ALIGN {
    tag "Aligning ${sample_id}"
    publishDir "${params.final_outdir}/alignments", mode: 'copy'
    // НЕ указываем container здесь, потому что он указан в nextflow.config
    
    input:
    tuple val(sample_id), path(fastq)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: aligned
    
    script:
    """
    echo "Aligning ${sample_id} to reference..."
    echo "FASTQ file: ${fastq}"
    echo "Reference: ${reference_fasta}"
    
    # УБРАЛИ установку через conda - инструменты должны быть в Docker-образе!
    
    # Check if minimap2 is available
    if ! command -v minimap2 &> /dev/null; then
        echo "ERROR: minimap2 not found in container"
        exit 1
    fi
    
    if ! command -v samtools &> /dev/null; then
        echo "ERROR: samtools not found in container"
        exit 1
    fi
    
    # Check input files
    if [ ! -f "${fastq}" ]; then
        echo "ERROR: FASTQ file ${fastq} not found!"
        exit 1
    fi
    
    if [ ! -f "${reference_fasta}" ]; then
        echo "ERROR: Reference file ${reference_fasta} not found!"
        exit 1
    fi
    
    # Run minimap2 alignment
    echo "Running minimap2 alignment..."
    minimap2 -ax sr ${reference_fasta} ${fastq} | \\
        samtools view -bS - | \\
        samtools sort -o ${sample_id}.sorted.bam -
    
    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index ${sample_id}.sorted.bam
    
    # Check results
    if [ ! -f "${sample_id}.sorted.bam" ]; then
        echo "ERROR: BAM file was not created!"
        exit 1
    fi
    
    if [ ! -f "${sample_id}.sorted.bam.bai" ]; then
        echo "ERROR: BAM index was not created!"
        exit 1
    fi
    
    # Show alignment statistics
    echo "Alignment statistics:"
    samtools flagstat ${sample_id}.sorted.bam
    
    echo "Alignment completed for ${sample_id}"
    """
}

process FREYJA_VARIANTS {
    tag "Calling variants ${sample_id}"
    publishDir "${params.final_outdir}/variants", mode: 'copy'
    container 'staphb/freyja:latest'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_variants.tsv"), path("${sample_id}_depths"), emit: variants_depths
    
    script:
    def refname = params.refname ?: 'NC_045512.2'
    def minq = params.minq ?: 20
    """
    echo "Calling variants for ${sample_id}..."
    echo "BAM file: ${bam}"
    echo "Reference: ${reference_fasta}"
    
    # Check file existence
    if [ ! -f "${bam}" ]; then
        echo "ERROR: BAM file ${bam} not found!"
        exit 1
    fi
    
    if [ ! -f "${reference_fasta}" ]; then
        echo "ERROR: Reference file ${reference_fasta} not found!"
        exit 1
    fi
    
    # Check BAM file integrity
    echo "Validating BAM file..."
    samtools quickcheck ${bam}
    if [ \$? -ne 0 ]; then
        echo "ERROR: BAM file validation failed!"
        exit 1
    fi
    
    # Check that BAM file contains data
    echo "Checking BAM file contents..."
    READ_COUNT=\$(samtools view -c ${bam})
    echo "Total reads in BAM: \$READ_COUNT"
    
    if [ "\$READ_COUNT" -eq 0 ]; then
        echo "WARNING: BAM file contains no reads!"
        touch "${sample_id}_variants.tsv"
        touch "${sample_id}_depths"
        exit 0
    fi
    
    # Call variants with Freyja
    echo "Running Freyja variants..."
    freyja variants \\
        ${bam} \\
        --variants ${sample_id}_variants.tsv \\
        --depths ${sample_id}_depths \\
        --ref ${reference_fasta} \\
        --refname ${refname} \\
        --minq ${minq}
    
    # Check that files were created
    if [ ! -f "${sample_id}_variants.tsv" ]; then
        echo "ERROR: Variants file was not created!"
        exit 1
    fi
    
    if [ ! -f "${sample_id}_depths" ]; then
        echo "ERROR: Depths file was not created!"
        exit 1
    fi
    
    echo "Variant calling completed for ${sample_id}"
    echo "Variants file size: \$(wc -l < ${sample_id}_variants.tsv) lines"
    echo "Depths file size: \$(wc -l < ${sample_id}_depths) lines"
    """
}

process FREYJA_DEMIX {
    tag "Demixing ${sample_id}"
    publishDir "${params.final_outdir}/demix", mode: 'copy'
    container 'staphb/freyja:latest'
    
    input:
    tuple val(sample_id), path(variants), path(depths)
    
    output:
    tuple val(sample_id), path("${sample_id}_demixed.tsv"), emit: results
    
    script:
    def eps = params.freyja_eps ?: "0.00000001"
    def covcut = params.covcut ?: 10
    def depthcutoff = params.depthcutoff ?: 0
    def confirmedonly = params.confirmedonly ? '--confirmedonly' : ''
    def pathogen = params.pathogen ?: 'SARS-CoV-2'
    def barcodes = params.barcodes ? "--barcodes ${params.barcodes}" : ''
    """
    echo "Running Freyja demix for ${sample_id}..."
    
    # Check if variants file has content
    if [ ! -s "${variants}" ]; then
        echo "WARNING: Variants file is empty, creating empty demix result"
        echo -e "summarized\tlineages\tabundances\tresid\tcoverage" > ${sample_id}_demixed.tsv
        echo -e "0.0\tundetermined\t1.0\t0.0\t0.0" >> ${sample_id}_demixed.tsv
        exit 0
    fi
    
    # Run demix
    echo "Running demix analysis..."
    freyja demix \\
        ${variants} \\
        ${depths} \\
        --output ${sample_id}_demixed.tsv \\
        --eps ${eps} \\
        --covcut ${covcut} \\
        --depthcutoff ${depthcutoff} \\
        --pathogen ${pathogen} \\
        ${confirmedonly} \\
        ${barcodes}
    
    # Verify output
    if [ ! -f "${sample_id}_demixed.tsv" ]; then
        echo "ERROR: Demix output file was not created!"
        exit 1
    fi
    
    echo "Demix completed for ${sample_id}"
    echo "Output file size: \$(wc -l < ${sample_id}_demixed.tsv) lines"
    
    # Show sample of results
    echo "Sample of demix results:"
    head -5 ${sample_id}_demixed.tsv
    """
}

process FREYJA_UPDATE {
    tag "Updating Freyja database"
    publishDir "${params.final_outdir}/freyja_db", mode: 'copy'
    container 'staphb/freyja:latest'
    
    output:
    path "freyja_updated.flag", emit: flag
    path "usher_barcodes.csv", emit: barcodes, optional: true
    
    script:
    def outdir = params.freyja_outdir ?: './'
    def pathogen = params.pathogen ?: 'SC2'
    """
    echo "Updating Freyja lineage database..."
    freyja update \\
        --outdir ${outdir} \\
        --pathogen ${pathogen}
    
    # Create a flag file to indicate completion
    touch freyja_updated.flag
    
    echo "Freyja database update completed"
    echo "Database files:"
    ls -lh ${outdir}
    """
}

process FREYJA_AGGREGATE {
    tag "Aggregating Freyja results"
    publishDir "${params.final_outdir}", mode: 'copy'
    container 'staphb/freyja:latest'
    
    input:
    path demix_results
    
    output:
    path "aggregated_freyja.tsv", emit: aggregated
    
    script:
    """
    echo "Aggregating Freyja demix results..."
    
    # Create a directory with all demix files
    mkdir -p demix_files
    
    # Copy all demix TSV files
    for file in ${demix_results}; do
        if [[ \$file == *.tsv ]]; then
            cp "\$file" demix_files/
        fi
    done
    
    # Check if we have files to aggregate
    file_count=\$(ls demix_files/*.tsv 2>/dev/null | wc -l)
    if [ "\$file_count" -eq 0 ]; then
        echo "ERROR: No demix TSV files found to aggregate!"
        exit 1
    fi
    
    echo "Found \$file_count demix files to aggregate"
    
    # Run aggregate (согласно freyja 2.0.0 help)
    freyja aggregate \
        demix_files/ \
        --output aggregated_freyja.tsv
    
    echo "Aggregation completed"
    echo "Aggregated results:"
    head -10 aggregated_freyja.tsv
    """
}

process FREYJA_PLOT {
    tag "Creating Freyja plots"
    publishDir "${params.final_outdir}/plots", mode: 'copy'
    container 'staphb/freyja:latest'
    
    input:
    path aggregated_results
    
    output:
    path "freyja_plot.pdf", optional: true, emit: plot
    path "freyja_plot.png", optional: true, emit: plot_png
    
    script:
    // def mincov = params.mincov ?: 60
    def thresh = params.thresh ?: 0.01
    def pathogen = params.pathogen ?: 'SC2'
    def lineage = params.lineage ? "--lineage ${params.lineage}" : ''
    """
    echo "Creating Freyja visualization plots..."
    
    # Create plot
    freyja plot \\
        ${aggregated_results} \\
        --output freyja_plot.pdf \\
        --mincov 1 \\
        --thresh ${thresh} \\
        --pathogen ${pathogen} \\
        ${lineage} || {
        echo "⚠️  Warning: Freyja plot generation failed"
        echo "Possible reasons:"
        echo "  - Too many lineages detected (try --thresh 0.10 or higher)"
        echo "  - Insufficient coverage (samples below --mincov threshold)"
        echo ""
        echo "Pipeline continues - other results are available in: ${params.final_outdir}"
        exit 0
        }
    
    # Try to convert to PNG if possible
    if command -v convert &> /dev/null; then
        convert -density 300 freyja_plot.pdf freyja_plot.png || true
    fi
    
    echo "Plot creation completed"
    """
}