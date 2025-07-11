#!/bin/bash

echo "🔍 Отладочный запуск - только 3 образца"
echo "======================================"

# Создать небольшой тестовый файл
echo "sample,fastq" > debug_test.csv
echo "mixture_1_shuffled,data/fastq/mixture_1_shuffled.fastq" >> debug_test.csv
echo "mixture_2_shuffled,data/fastq/mixture_2_shuffled.fastq" >> debug_test.csv
echo "mixture_3_shuffled,data/fastq/mixture_3_shuffled.fastq" >> debug_test.csv

echo "📋 Тестируемые образцы:"
cat debug_test.csv
echo ""

# Очистить предыдущие результаты работы
echo "🧹 Очистка work директории..."
rm -rf work/
rm -rf .nextflow*

echo "🚀 Запуск с отладкой..."
echo ""

# Запуск с verbose логированием
nextflow run main.nf \
  --input debug_test.csv \
  --fasta_dir data/lineages \
  -profile local \
  -with-timeline timeline.html \
  -with-report report.html \
  -with-trace \
  -with-dag dag.png

echo ""
echo "📊 Проверка результатов:"
for sample in mixture_1_shuffled mixture_2_shuffled mixture_3_shuffled; do
    if [ -f "results/${sample}/${sample}/abundance.tsv" ]; then
        lines=$(wc -l < "results/${sample}/${sample}/abundance.tsv")
        echo "✅ ${sample}: ${lines} строк"
    else
        echo "❌ ${sample}: не найден"
    fi
done

echo ""
echo "📁 Созданные файлы отчетов:"
echo "- timeline.html (временная диаграмма)"
echo "- report.html (подробный отчет)"
echo "- trace.txt (трассировка процессов)"
echo "- dag.png (граф выполнения)"

echo ""
echo "🔍 Если есть ошибки, проверьте:"
echo "1. nextflow log - показать логи последнего запуска"
echo "2. ls work/ - проверить рабочие директории"
echo "3. cat .nextflow.log - подробные логи"