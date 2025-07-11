#!/bin/bash

echo "🧬 MetBio Pipeline - Полный анализ"
echo "=================================="
echo "📂 Файл входных данных: victor_full.csv"
echo "🔬 Количество образцов: $(tail -n +2 victor_full.csv | wc -l)"
echo "🧪 Референсные геномы: $(ls data/lineages/*.fasta | wc -l) линий SARS-CoV-2"
echo "⚙️  Инструмент: lr-kallisto (для long reads)"
echo ""

# Проверка файлов
echo "🔍 Проверка входных данных..."
if [ ! -f "victor_full.csv" ]; then
    echo "❌ Файл victor_full.csv не найден!"
    exit 1
fi

if [ ! -d "data/lineages" ]; then
    echo "❌ Папка data/lineages не найдена!"
    exit 1
fi

fastq_count=$(ls data/fastq/*.fastq 2>/dev/null | wc -l)
if [ $fastq_count -eq 0 ]; then
    echo "❌ FASTQ файлы не найдены в data/fastq/"
    exit 1
fi

echo "✅ Все входные данные на месте"
echo ""

# Показать первые несколько образцов
echo "📋 Первые 5 образцов для анализа:"
head -6 victor_full.csv | tail -5
echo "   ... и еще $(expr $(tail -n +2 victor_full.csv | wc -l) - 5) образцов"
echo ""

# Запустить анализ
echo "🚀 Запуск полного анализа..."
echo "Команда: nextflow run main.nf --input victor_full.csv --fasta_dir data/lineages -profile local"
echo ""

# Время начала
start_time=$(date)
echo "⏰ Начало: $start_time"
echo ""

# Запуск
nextflow run main.nf \
  --input victor_full.csv \
  --fasta_dir data/lineages \
  -profile local

# Результат
echo ""
echo "📊 Анализ завершен!"
echo "⏰ Время завершения: $(date)"
echo ""

# Проверка результатов
echo "🔍 Проверка результатов..."
success_count=0
total_samples=$(tail -n +2 victor_full.csv | wc -l)

while IFS=, read -r sample fastq; do
    if [ "$sample" != "sample" ]; then  # Пропустить заголовок
        if [ -f "results/${sample}/${sample}/abundance.tsv" ]; then
            lines=$(wc -l < "results/${sample}/${sample}/abundance.tsv")
            echo "✅ ${sample}: ${lines} строк в abundance.tsv"
            success_count=$((success_count + 1))
        else
            echo "❌ ${sample}: результаты не найдены"
        fi
    fi
done < victor_full.csv

echo ""
echo "📈 Итоговая статистика:"
echo "✅ Успешно обработано: ${success_count}/${total_samples} образцов"
echo "📁 Результаты сохранены в: results/"
echo ""

if [ $success_count -eq $total_samples ]; then
    echo "🎉 ВСЕ ОБРАЗЦЫ УСПЕШНО ОБРАБОТАНЫ!"
else
    echo "⚠️  Некоторые образцы не обработались. Проверьте логи."
fi