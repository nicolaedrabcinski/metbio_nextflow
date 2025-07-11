#!/bin/bash

echo "🚀 Быстрый тест исправленного pipeline"
echo "====================================="

# Создать CSV с тремя образцами
echo "sample,fastq" > quick_test_fixed.csv
echo "mixture_1_shuffled,data/fastq/mixture_1_shuffled.fastq" >> quick_test_fixed.csv
echo "mixture_2_shuffled,data/fastq/mixture_2_shuffled.fastq" >> quick_test_fixed.csv
echo "mixture_3_shuffled,data/fastq/mixture_3_shuffled.fastq" >> quick_test_fixed.csv

echo "📋 Тестируемые образцы:"
cat quick_test_fixed.csv
echo ""

# Проверить что файлы существуют
echo "🔍 Проверка файлов:"
echo "✅ Kallisto: $(ls -la tools/kallisto/kallisto 2>/dev/null || echo '❌ НЕ НАЙДЕН')"
echo "✅ FASTA файлы: $(ls data/lineages/*.fasta | wc -l) файлов"
echo "✅ FASTQ файлы для теста:"
for i in 1 2 3; do
    file="data/fastq/mixture_${i}_shuffled.fastq"
    if [ -f "$file" ]; then
        echo "   ✅ mixture_${i}_shuffled.fastq"
    else
        echo "   ❌ mixture_${i}_shuffled.fastq"
    fi
done

echo ""
echo "🧹 Очистка..."
rm -rf work/ .nextflow*

echo ""
echo "🚀 Запуск..."
nextflow run main.nf \
  --input quick_test_fixed.csv \
  --fasta_dir data/lineages \
  -profile local

echo ""
echo "📊 Результаты:"
success=0
for sample in mixture_1_shuffled mixture_2_shuffled mixture_3_shuffled; do
    if [ -f "results/${sample}/${sample}/abundance.tsv" ]; then
        lines=$(wc -l < "results/${sample}/${sample}/abundance.tsv")
        echo "✅ ${sample}: ${lines} строк"
        success=$((success + 1))
    else
        echo "❌ ${sample}: не найден"
    fi
done

echo ""
if [ $success -eq 3 ]; then
    echo "🎉 ВСЕ 3 ОБРАЗЦА УСПЕШНО ОБРАБОТАНЫ!"
    echo ""
    echo "🚀 Готов запускать полный анализ всех 51 образца:"
    echo "nextflow run main.nf --input victor_full.csv --fasta_dir data/lineages -profile local"
else
    echo "⚠️ Обработано только $success из 3 образцов"
    echo "🔍 Проверьте логи: nextflow log"
fi