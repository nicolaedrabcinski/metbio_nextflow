process FREYJA_ALIGN {
    tag "Aligning ${sample_id}"
    publishDir "${params.final_outdir}/alignments", mode: 'copy'
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'

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
    minimap2 -ax ${params.minimap2_preset} ${reference_fasta} ${fastq} | \\
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
    tag "Variant calling for ${sample_id}"
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
    echo "Running variant calling for ${sample_id}..."
    echo "BAM file: ${bam}"
    echo "Reference: ${reference_fasta}"
    
    # Check if freyja is available
    if ! command -v freyja &> /dev/null; then
        echo "ERROR: freyja not found in container"
        exit 1
    fi
    
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
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(variants), path(depths)
    
    output:
    tuple val(sample_id), path("${sample_id}_demixed.tsv"), optional: true, emit: results
    
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
        echo "WARNING: Variants file is empty for ${sample_id}, skipping demix"
        exit 0
    fi

    # Check minimum variant count for solver
    variant_count=\$(tail -n +2 "${variants}" | wc -l)
    if [ "\$variant_count" -lt 2 ]; then
        echo "WARNING: Too few variants (\$variant_count) for ${sample_id}, skipping demix"
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
    
    # Run aggregate
    freyja aggregate \\
        demix_files/ \\
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
    def thresh = params.thresh ?: 0.01
    def pathogen = params.pathogen ?: 'SC2'
    def lineage = params.lineage ? "--lineage ${params.lineage}" : ''
    """
    echo "Creating Freyja visualization plots..."
    
    # Create plot
    freyja plot \\
        ${aggregated_results} \\
        --output freyja_plot \\
        --thresh ${thresh} \\
        --pathogen ${pathogen} \\
        ${lineage}
    
    echo "Plot creation completed"
    """
}

process COVERAGE_PLOT {
    tag "Coverage ${sample_id}"
    publishDir "${params.final_outdir}/coverage", mode: 'copy'
    container 'staphb/samtools:latest'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.png"), emit: plot
    tuple val(sample_id), path("${sample_id}.depth.tsv"),   emit: depth

    script:
    """
    echo "Calculating coverage depth for ${sample_id}..."

    # Get per-base depth (forward + reverse separately)
    samtools depth -a -d 0 ${bam} > ${sample_id}.depth.tsv

    # Get forward strand depth
    samtools view -F 16 -b ${bam} | samtools depth -a -d 0 - > ${sample_id}.fwd.tsv
    # Get reverse strand depth
    samtools view -f 16 -b ${bam} | samtools depth -a -d 0 - > ${sample_id}.rev.tsv

    python3 - << 'PYEOF'
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sample_id = "${sample_id}"

def load_depth(path):
    pos, depth = [], []
    try:
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) >= 3:
                    pos.append(int(parts[1]))
                    depth.append(int(parts[2]))
    except:
        pass
    return np.array(pos), np.array(depth)

pos_all, dep_all = load_depth("${sample_id}.depth.tsv")
pos_fwd, dep_fwd = load_depth("${sample_id}.fwd.tsv")
pos_rev, dep_rev = load_depth("${sample_id}.rev.tsv")

if len(pos_all) == 0:
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.text(0.5, 0.5, 'No coverage data', ha='center', va='center', fontsize=14)
    ax.axis('off')
    plt.savefig(f"{sample_id}.coverage.png", dpi=150)
    plt.close()
    sys.exit(0)

# Genome length and coverage stats
genome_len    = pos_all[-1] if len(pos_all) > 0 else 1
mean_cov      = float(np.mean(dep_all)) if len(dep_all) > 0 else 0
pct_covered   = float(np.mean(dep_all > 0) * 100) if len(dep_all) > 0 else 0
n_reads_fwd   = int(np.sum(dep_fwd)) if len(dep_fwd) > 0 else 0
n_reads_rev   = int(np.sum(dep_rev)) if len(dep_rev) > 0 else 0

fig, ax = plt.subplots(figsize=(14, 5))

# Fill forward strand (cyan, like Pavian)
if len(pos_fwd) > 0:
    ax.fill_between(pos_fwd, dep_fwd, color='#00BCD4', alpha=0.85,
                    linewidth=0, label='Forward strand')

# Fill reverse strand (pink, like Pavian) — plotted as negative for visual separation
# But shown as positive overlay (Pavian shows both stacked)
if len(pos_rev) > 0:
    ax.fill_between(pos_rev, dep_rev, color='#F06292', alpha=0.55,
                    linewidth=0, label='Reverse strand')

# Total depth line
ax.plot(pos_all, dep_all, color='#1565C0', linewidth=0.4, alpha=0.6)

# Annotations
ax.set_xlabel('Position (bp)', fontsize=11)
ax.set_ylabel('Coverage depth', fontsize=11)
ax.set_title(
    f'Coverage — {sample_id}   '
    f'[{pos_all[0]:,}–{genome_len:,} bp covered]   '
    f'avg cov: {mean_cov:.0f}x',
    fontsize=11, fontweight='bold'
)
ax.set_xlim(pos_all[0], genome_len)
ax.set_ylim(0)

fwd_patch = mpatches.Patch(color='#00BCD4', alpha=0.85, label=f'Forward ({n_reads_fwd:,} reads)')
rev_patch = mpatches.Patch(color='#F06292', alpha=0.55, label=f'Reverse ({n_reads_rev:,} reads)')
ax.legend(handles=[fwd_patch, rev_patch], loc='upper right', fontsize=9)

ax.grid(axis='y', alpha=0.25, linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add stats box (like Pavian)
stats_text = (f"{pos_all[0]:,} / {genome_len:,} bp covered\\n"
              f"avg cov: {mean_cov:.0f}x\\n"
              f"{pct_covered:.1f}% covered")
ax.text(0.98, 0.95, stats_text, transform=ax.transAxes,
        ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#cccccc', alpha=0.9))

plt.tight_layout()
plt.savefig(f"{sample_id}.coverage.png", dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Coverage plot saved: {sample_id}.coverage.png")
PYEOF

    echo "Coverage plot completed for ${sample_id}"
    """
}