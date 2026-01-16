#!/bin/bash
# run_viloca_batch.sh

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
REF="$BASEDIR/data/NC_045512_2.fasta"
ALIGNDIR="$BASEDIR/custom_results/alignments_bam"  # Используем BAM вместо CRAM
OUTDIR="$BASEDIR/custom_results/viloca_results"
THREADS=32

mkdir -p "$OUTDIR"

if [ ! -f "$REF" ]; then
    echo "Error: Reference file $REF not found"
    exit 1
fi

if [ ! -d "$ALIGNDIR" ]; then
    echo "Error: Alignment directory $ALIGNDIR not found"
    exit 1
fi

echo "Используем потоков: $THREADS (0 = все доступные)"
echo "Обрабатываем семплы 11-50 (пропускаем 1-10)"
echo "---"

for bam in "$ALIGNDIR"/*.bam; do
    SAMPLE=$(basename "$bam" .bam)
    
    # Извлекаем номер mixture
    MIXTURE_NUM=$(echo "$SAMPLE" | grep -oP 'mixture_\K\d+')
    
    # Пропускаем семплы 1-10
    if [ ! -z "$MIXTURE_NUM" ] && [ "$MIXTURE_NUM" -le 10 ]; then
        echo "⊘ Пропускаем $SAMPLE (mixture $MIXTURE_NUM)"
        continue
    fi
    
    echo "Запуск VILOCA для $SAMPLE (mixture $MIXTURE_NUM)..."
    
    SAMPLE_DIR="$OUTDIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"
    
    ABS_BAM=$(realpath "$bam")
    ABS_REF=$(realpath "$REF")
    
    (cd "$SAMPLE_DIR" && viloca run \
        -b "$ABS_BAM" \
        -f "$ABS_REF" \
        --mode use_quality_scores \
        -of csv \
        -t "$THREADS" \
        > "${SAMPLE}_viloca.log" 2>&1)
    
    if [ $? -eq 0 ]; then
        echo "✓ $SAMPLE завершен успешно"
    else
        echo "✗ $SAMPLE завершился с ошибкой (см. $SAMPLE_DIR/${SAMPLE}_viloca.log)"
    fi
    echo "---"
done

echo "Готово! Результаты в $OUTDIR"