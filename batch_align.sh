#!/bin/bash

REF="data/NC_045512_2.fasta"
OUTDIR="custom_results/alignments"
mkdir -p "$OUTDIR"

for folder in data/reads_artic_big_overlaps data/reads_artic_small_overlaps; do
  for i in {10..50}; do
    FASTQ="$folder/mixture_${i}.fastq"
    BAM="$OUTDIR/mixture_${i}_$(basename $folder).bam"
    CRAM="$OUTDIR/mixture_${i}_$(basename $folder).cram"
    if [ -f "$FASTQ" ]; then
      echo "Processing $FASTQ..."
      minimap2 -a "$REF" "$FASTQ" | samtools view -bS - | samtools sort -o "$BAM"
      samtools index "$BAM"
      samtools view -C -T "$REF" "$BAM" -o "$CRAM"
      samtools index "$CRAM"
    fi
  done
done
