#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GENERATE_CSV_FROM_FASTQ_DIR {
    tag "Auto-generate CSV"
    publishDir "${params.final_outdir}/metadata", mode: 'copy'
    
    input:
    path fastq_dir
    
    output:
    path "auto_samples.csv", emit: csv_file
    
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
        
        # Фильтровать только файлы с _shuffled в названии
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