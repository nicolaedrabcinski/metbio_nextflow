#!/bin/bash
# convert_cram_to_bam.sh

REF="data/NC_045512_2.fasta"
CRAM_DIR="custom_results/alignments"
BAM_DIR="custom_results/alignments_bam"

mkdir -p "$BAM_DIR"

echo "Конвертация CRAM в BAM..."
echo "---"

for cram in "$CRAM_DIR"/*.cram; do
    SAMPLE=$(basename "$cram" .cram)
    BAM="$BAM_DIR/${SAMPLE}.bam"
    
    echo "Конвертация $SAMPLE..."
    samtools view -b -T "$REF" -o "$BAM" "$cram"
    
    if [ $? -eq 0 ]; then
        echo "✓ $SAMPLE.bam создан"
        # Создаем индекс
        samtools index "$BAM"
        echo "✓ Индекс создан"
    else
        echo "✗ Ошибка при конвертации $SAMPLE"
    fi
    echo "---"
done

echo "Готово! BAM файлы в $BAM_DIR"