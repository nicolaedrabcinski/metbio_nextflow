#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters by default
params.input = null
params.fastq_dir = null  // Новый параметр для автогенерации CSV
params.lineages_fasta = null
params.fasta_dir = null  // Для обратной совместимости
params.index = null
params.outdir = 'results'
params.help = false  // Исправлено: добавлена инициализация
params.filter_shuffled = false  // Изменено: по умолчанию использовать ВСЕ файлы
params.kmer_size = 31  // Новый параметр: размер k-mer для kallisto
params.create_plots = true  // Создавать ли графики результатов

// Создать папку результатов на основе имени входной папки
def getOutputDir() {
    if (params.fastq_dir) {
        // Извлечь имя папки из полного пути
        def fastq_dir_name = new File(params.fastq_dir).getName()
        return "${params.outdir}/${fastq_dir_name}"
    } else if (params.input) {
        // Если используется CSV, извлечь имя без расширения
        def input_name = new File(params.input).getName().replaceAll(/\.[^\.]+$/, '')
        return "${params.outdir}/${input_name}"
    } else {
        return params.outdir
    }
}

def final_outdir = getOutputDir()

// Help message
def helpMessage() {
    log.info """
    ========================================
    MetBio Nextflow Pipeline - Usage
    ========================================
    
    Option 1 - Auto-generate CSV from FASTQ directory:
    nextflow run main.nf --fastq_dir /path/to/fastq/files --lineages_fasta lineages.fasta
    
    Option 2 - Use existing CSV:
    nextflow run main.nf --input samples.csv --lineages_fasta lineages.fasta
    
    Parameters:
    --fastq_dir         Directory containing FASTQ files (auto-generates CSV)
    --input             CSV file with sample information (manual mode)  
    --lineages_fasta    FASTA file with reference lineages
    --fasta_dir         Directory with individual FASTA files (legacy)
    --index             Pre-built kallisto index file
    --outdir            Output directory (default: results)
    --filter_shuffled   Only use files with '_shuffled' in name (default: false)
    --kmer_size         K-mer size for kallisto (default: 31)
    --create_plots      Create visualization plots (default: true)
    
    Note: Results will be saved in subdirectory named after input folder/file
    
    Examples:
    # Auto-scan ALL FASTQ files in directory
    nextflow run main.nf --fastq_dir data/reads_artic_big_overlaps --lineages_fasta lineages.fasta
    
    # Use custom k-mer size
    nextflow run main.nf --fastq_dir data/reads --lineages_fasta lineages.fasta --kmer_size 21
    
    # Process short overlaps data
    nextflow run main.nf --fastq_dir data/reads_artic_short_overlaps --lineages_fasta lineages.fasta
    
    # Use only _shuffled files
    nextflow run main.nf --fastq_dir data/reads --lineages_fasta lineages.fasta --filter_shuffled true
    
    # Use existing CSV
    nextflow run main.nf --input full.csv --lineages_fasta lineages.fasta  
    # Results in: results/full/
    ========================================
    """.stripIndent()
}

// Show help if no required parameters
if (params.help || (!params.fastq_dir && !params.input) || (!params.lineages_fasta && !params.fasta_dir && !params.index)) {
    helpMessage()
    exit 0
}

// Logging
log.info """
========================================
MetBio Nextflow Pipeline
========================================
fastq_dir      : ${params.fastq_dir}
input          : ${params.input}
lineages_fasta : ${params.lineages_fasta}
fasta_dir      : ${params.fasta_dir}
index          : ${params.index}
outdir         : ${final_outdir}
========================================
"""

// Импорт процессов
include { DOWNLOAD_KALLISTO } from './modules/download_kallisto.nf'
include { KALLISTO_INDEX } from './modules/kallisto_index.nf'
include { LR_KALLISTO } from './modules/lr_kallisto.nf'

// Процесс для автогенерации CSV из папки с FASTQ
process GENERATE_CSV_FROM_FASTQ_DIR {
    tag "Auto-generate CSV"
    publishDir "${final_outdir}/metadata", mode: 'copy'
    
    input:
    path fastq_dir
    
    output:
    path "auto_samples.csv"
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    import csv
    from pathlib import Path
    
    fastq_files = []
    
    # Поиск FASTQ файлов с различными расширениями
    extensions = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
    
    # Рекурсивный поиск во всех подпапках
    for ext in extensions:
        pattern = os.path.join('${fastq_dir}', '**', ext)
        found_files = glob.glob(pattern, recursive=True)
        
        # Фильтровать только файлы с _shuffled в названии (согласно указаниям Виктора)
        filter_shuffled = ${params.filter_shuffled ? 'True' : 'False'}
        if filter_shuffled:
            shuffled_files = [f for f in found_files if '_shuffled' in os.path.basename(f)]
            fastq_files.extend(shuffled_files)
        else:
            fastq_files.extend(found_files)
    
    # Если shuffled файлы не найдены и фильтр включен, взять все файлы
    if not fastq_files and filter_shuffled:
        print("No '_shuffled' files found, using all FASTQ files...")
        for ext in extensions:
            pattern = os.path.join('${fastq_dir}', '**', ext)
            fastq_files.extend(glob.glob(pattern, recursive=True))
    
    # Сортировка для стабильного порядка
    fastq_files.sort()
    
    print(f"Found {len(fastq_files)} FASTQ files:")
    for f in fastq_files[:10]:  # Показать первые 10 для проверки
        print(f"  - {f}")
    if len(fastq_files) > 10:
        print(f"  ... and {len(fastq_files) - 10} more files")
    
    # Создание CSV файла
    with open('auto_samples.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'fastq'])  # Заголовки как в оригинальном CSV
        
        for fastq_file in fastq_files:
            # Извлечение имени образца из имени файла
            sample_name = Path(fastq_file).stem
            
            # Удаление распространенных суффиксов
            suffixes_to_remove = ['.fastq', '.fq', '.gz', '_shuffled', '_long', '_short']
            for suffix in suffixes_to_remove:
                if sample_name.endswith(suffix):
                    sample_name = sample_name[:-len(suffix)]
            
            # Использование абсолютного пути к файлу
            abs_path = os.path.abspath(fastq_file)
            writer.writerow([sample_name, abs_path])
    
    print(f"Generated CSV file with {len(fastq_files)} samples")
    print("Sample names preview:")
    with open('auto_samples.csv', 'r') as f:
        lines = f.readlines()
        for line in lines[1:6]:  # Показать первые 5 образцов
            print(f"  {line.strip()}")
        if len(lines) > 6:
            print(f"  ... and {len(lines) - 6} more samples")
    """
}

// Локальные процессы с правильными путями
process KALLISTO_INDEX_LOCAL {
    tag "Building kallisto index"
    publishDir "${final_outdir}/index", mode: 'copy'
    
    input:
    path lineages_fasta
    path kallisto_path
    
    output:
    path "lineages.idx"
    
    script:
    """
    ${kallisto_path} index \\
        -i lineages.idx \\
        -k ${params.kmer_size} \\
        ${lineages_fasta}
    """
}

process LR_KALLISTO_LOCAL {
    tag "${sample}"
    publishDir "${final_outdir}/quantification", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq_file)
    path index_file
    path kallisto_path
    
    output:
    tuple val(sample), path("${sample}")
    
    script:
    """
    mkdir -p ${sample}
    
    ${kallisto_path} quant \\
        -i ${index_file} \\
        -o ${sample} \\
        --single \\
        -l 1000 \\
        -s 30 \\
        ${fastq_file}
    """
}

// Процесс для создания графиков и конвертации TPM в проценты
process CREATE_PLOTS {
    tag "Creating visualization plots"
    publishDir "${final_outdir}/plots", mode: 'copy'
    
    input:
    path results_file
    
    output:
    path "*.png", optional: true
    path "summary_percentages.csv"
    
    when:
    params.create_plots
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    import glob
    from pathlib import Path
    
    # Читать список путей к результатам
    results_paths = []
    with open('${results_file}', 'r') as f:
        for line in f:
            if line.strip():
                unique_name, dir_path = line.strip().split('\\t')
                results_paths.append((unique_name, dir_path))
    
    print(f"Found {len(results_paths)} result directories")
    
    # Собрать все результаты
    all_results = []
    processed_samples = set()
    
    for unique_name, dir_path in results_paths:
        # Извлечь оригинальное имя образца из пути к папке
        sample_base_name = os.path.basename(dir_path)
        
        # Если образец уже обработан, добавить суффикс
        original_name = sample_base_name
        counter = 1
        while sample_base_name in processed_samples:
            sample_base_name = f"{original_name}_v{counter}"
            counter += 1
        
        processed_samples.add(sample_base_name)
        
        abundance_file = os.path.join(dir_path, 'abundance.tsv')
        
        if os.path.exists(abundance_file):
            try:
                df = pd.read_csv(abundance_file, sep='\\t')
                df['sample'] = sample_base_name  # Использовать уникальное имя
                # Конвертировать TPM в проценты: % = TPM / 1,000,000 * 100
                df['percentage'] = df['tpm'] / 1000000 * 100
                all_results.append(df)
                print(f"Processed {sample_base_name}: {len(df)} targets")
            except Exception as e:
                print(f"Error processing {dir_path}: {e}")
        else:
            print(f"abundance.tsv not found in {dir_path}")
    
    if all_results:
        # Объединить все результаты
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Сохранить сводную таблицу с процентами
        summary_df = combined_df.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
        summary_df.to_csv('summary_percentages.csv')
        
        print(f"Created summary with {len(summary_df)} targets across {len(summary_df.columns)} samples")
        
        # Создать barplot для топ-10 targets по среднему проценту
        mean_percentages = combined_df.groupby('target_id')['percentage'].mean().sort_values(ascending=False)
        top_targets = mean_percentages.head(10).index.tolist()
        
        if len(top_targets) > 0:
            plot_data = combined_df[combined_df['target_id'].isin(top_targets)]
            
            # Создать barplot
            plt.figure(figsize=(12, 8))
            avg_data = plot_data.groupby('target_id')['percentage'].mean().sort_values(ascending=False)
            plt.bar(range(len(avg_data)), avg_data.values)
            plt.title('Top 10 Lineages by Average Percentage')
            plt.xlabel('Lineage')
            plt.ylabel('Percentage (%)')
            plt.xticks(range(len(avg_data)), avg_data.index, rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig('top_lineages_barplot.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Создать heatmap для всех образцов и топ targets (только если много образцов)
            if len(plot_data['sample'].unique()) > 1:
                heatmap_data = plot_data.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
                
                # Ограничить количество колонок для читаемости
                if len(heatmap_data.columns) > 20:
                    # Взять образцы с наибольшим разнообразием
                    sample_variance = heatmap_data.var(axis=0).sort_values(ascending=False)
                    top_samples = sample_variance.head(20).index
                    heatmap_data = heatmap_data[top_samples]
                
                plt.figure(figsize=(max(10, len(heatmap_data.columns) * 0.4), 
                                  max(6, len(heatmap_data.index) * 0.3)))
                sns.heatmap(heatmap_data, annot=False, fmt='.1f', cmap='viridis', 
                           cbar_kws={'label': 'Percentage (%)'})
                plt.title('Lineage Abundance Across Samples (Top 10 Lineages)')
                plt.xlabel('Sample')
                plt.ylabel('Lineage')
                plt.tight_layout()
                plt.savefig('abundance_heatmap.png', dpi=300, bbox_inches='tight')
                plt.close()
        
        print(f"Processed {len(all_results)} samples total")
        print(f"Found {len(combined_df['target_id'].unique())} unique targets")
        print(f"Created visualization plots and summary CSV")
    else:
        print("No abundance files found")
        # Создать пустую сводную таблицу
        pd.DataFrame().to_csv('summary_percentages.csv')
    """
}

workflow {
    // Использовать существующий kallisto или скачать
    if (file('./tools/kallisto/kallisto').exists()) {
        kallisto_path = Channel.value(file('./tools/kallisto/kallisto'))
        log.info "Using existing kallisto: ./tools/kallisto/kallisto"
    } else {
        DOWNLOAD_KALLISTO()
        kallisto_path = DOWNLOAD_KALLISTO.out.kallisto_path
        log.info "Downloaded kallisto"
    }
    
    // Определить источник входных данных
    if (params.fastq_dir) {
        // Автогенерация CSV из папки с FASTQ
        log.info "Auto-generating CSV from FASTQ directory: ${params.fastq_dir}"
        fastq_directory = Channel.fromPath(params.fastq_dir, type: 'dir')
        csv_file = GENERATE_CSV_FROM_FASTQ_DIR(fastq_directory)
        
        // Читать сгенерированный CSV
        input_ch = csv_file
            .splitCsv(header: true)
            .map { row -> 
                [row.sample, file(row.fastq)]
            }
    } else if (params.input) {
        // Использовать существующий CSV файл
        log.info "Using existing CSV file: ${params.input}"
        input_ch = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> 
                [row.sample, file(row.fastq)]
            }
    } else {
        error "You must specify either --fastq_dir or --input"
    }
    
    // Создать или использовать существующий индекс
    if (params.index) {
        // Использовать готовый индекс
        index_ch = Channel.value(file(params.index))
        log.info "Using existing index: ${params.index}"
    } else if (params.lineages_fasta) {
        // Использовать готовый lineages.fasta файл
        log.info "Creating index from lineages.fasta: ${params.lineages_fasta}"
        lineages_file = Channel.value(file(params.lineages_fasta))
        KALLISTO_INDEX_LOCAL(lineages_file, kallisto_path)
        index_ch = KALLISTO_INDEX_LOCAL.out
    } else if (params.fasta_dir) {
        // Создать новый индекс из отдельных FASTA файлов (старый способ)
        fasta_files = Channel
            .fromPath("${params.fasta_dir}/*.fasta")
            .collect()
        
        log.info "Creating index from FASTA files in: ${params.fasta_dir}"
        KALLISTO_INDEX_LOCAL(fasta_files, kallisto_path)
        index_ch = KALLISTO_INDEX_LOCAL.out
    } else {
        error "You must specify --lineages_fasta, --fasta_dir, or --index"
    }
    
    // Запустить lr-kallisto для каждого образца
    quantification_results = LR_KALLISTO_LOCAL(
        input_ch,
        index_ch,
        kallisto_path
    )
    
    // Создать графики и сводные данные
    if (params.create_plots) {
        // Собрать все результаты в один канал, используя уникальные имена
        results_with_unique_names = quantification_results
            .map { sample, dir -> 
                // Создать уникальное имя файла на основе пути к папке
                def unique_name = dir.toString().replaceAll(/.*\/work\//, '').replaceAll(/\//, '_')
                return tuple(unique_name, dir)
            }
            .collectFile(name: 'results_paths.txt', newLine: true) { unique_name, dir ->
                "${unique_name}\t${dir}"
            }
        
        CREATE_PLOTS(results_with_unique_names)
    }
    
    // Вывести статус завершения
    quantification_results
        .subscribe { sample, result_dir ->
            log.info "✅ Sample ${sample} processed: ${result_dir}"
        }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline is completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'ERROR'}
    Results: ${final_outdir}
    Execution time: ${workflow.duration}
    ========================================
    """
}