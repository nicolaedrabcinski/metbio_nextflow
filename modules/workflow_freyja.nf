process FREYJA_ALIGN {
    tag "Aligning ${sample_id}"
    publishDir "${params.outdir}/${params.fastq_dir ? file(params.fastq_dir).name : 'analysis'}/alignment", mode: 'copy', pattern: "*.bam*"
    
    input:
    tuple val(sample_id), path(fastq)
    path minimap2_bin
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}.bam.bai", emit: bai
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "🧬 Aligning ${sample_id} with minimap2..."
    
    # Align reads to reference
    ${minimap2_bin}/minimap2 -ax map-ont ${reference_fasta} ${fastq} | \\
        samtools view -bS - | \\
        samtools sort -o ${sample_id}.bam -
    
    # Index BAM file
    samtools index ${sample_id}.bam
    
    echo "✅ Alignment completed for ${sample_id}"
    """
}

process FREYJA_VARIANTS {
    tag "Calling variants ${sample_id}"
    publishDir "${params.outdir}/${params.fastq_dir ? file(params.fastq_dir).name : 'analysis'}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path freyja_bin
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_variants.tsv"), path("${sample_id}_depths"), emit: variants_depths
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "🧬 Calling variants for ${sample_id}..."
    
    # Call variants using Freyja
    ${freyja_bin}/freyja variants \\
        --variants ${sample_id}_variants.tsv \\
        --depths ${sample_id}_depths \\
        --ref ${reference_fasta} \\
        ${bam}
    
    echo "✅ Variant calling completed for ${sample_id}"
    """
}

process FREYJA_DEMIX {
    tag "Demixing ${sample_id}"
    publishDir "${params.outdir}/${params.fastq_dir ? file(params.fastq_dir).name : 'analysis'}/quantification/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(variants), path(depths)
    path freyja_bin
    
    output:
    tuple val(sample_id), path("${sample_id}_demixed.tsv"), emit: results
    path "${sample_id}_demixed.tsv", emit: demixed_results
    
    script:
    def eps = params.freyja_eps ?: "0.00000001"
    """
    #!/bin/bash
    set -e
    
    echo "🧬 Demixing lineages for ${sample_id}..."
    
    # Perform lineage demixing
    ${freyja_bin}/freyja demix \\
        --eps ${eps} \\
        --output ${sample_id}_demixed.tsv \\
        ${variants} \\
        ${depths}
    
    # Convert to abundance format for compatibility
    echo -e "target_id\\test_counts\\teff_length\\ttpm" > abundance.tsv
    
    # Parse Freyja output and convert to kallisto-like format
    if [ -f "${sample_id}_demixed.tsv" ]; then
        python3 << 'EOF'
import pandas as pd
import sys

try:
    # Read Freyja output
    with open("${sample_id}_demixed.tsv", 'r') as f:
        lines = f.readlines()
    
    # Find summarized section
    lineages = []
    abundances = []
    
    for line in lines:
        if line.strip().startswith("summarized"):
            parts = line.strip().split("\\t")
            if len(parts) >= 2:
                lineage_data = parts[1].split()
                abundance_data = parts[2].split() if len(parts) > 2 else []
                
                for i, lineage in enumerate(lineage_data):
                    if i < len(abundance_data):
                        abundance = float(abundance_data[i])
                        if abundance > 0:
                            lineages.append(lineage)
                            abundances.append(abundance)
    
    # Write abundance file
    with open("abundance.tsv", 'w') as f:
        f.write("target_id\\test_counts\\teff_length\\ttpm\\n")
        for lineage, abundance in zip(lineages, abundances):
            # Convert abundance (0-1) to TPM-like values
            tpm = abundance * 1000000
            f.write(f"{lineage}\\t{abundance*1000}\\t1000\\t{tpm}\\n")
            
    print(f"✅ Converted {len(lineages)} lineages to abundance format")
    
except Exception as e:
    print(f"⚠️  Error converting Freyja output: {e}")
    # Create empty abundance file
    with open("abundance.tsv", 'w') as f:
        f.write("target_id\\test_counts\\teff_length\\ttpm\\n")
EOF
    else
        echo "Freyja demix output not found, creating empty abundance.tsv"
    fi
    
    echo "Demixing completed for ${sample_id}"
    """
}