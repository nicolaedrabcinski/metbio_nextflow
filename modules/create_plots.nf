#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_PLOTS {
    tag "Creating visualization plots"
    publishDir "${params.final_outdir}/plots", mode: 'copy'
    
    input:
    path results_file
    
    output:
    path "*.png", optional: true, emit: plots
    path "summary_percentages.csv", emit: summary
    
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
    
    print("🎨 Starting visualization creation...")
    
    # Читать список путей к результатам
    results_paths = []
    with open('${results_file}', 'r') as f:
        for line in f:
            if line.strip():
                unique_name, dir_path = line.strip().split('\\t')
                results_paths.append((unique_name, dir_path))
    
    print(f"📁 Found {len(results_paths)} result directories")
    
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
                df['sample'] = sample_base_name
                # Конвертировать TPM в проценты: % = TPM / 1,000,000 * 100
                df['percentage'] = df['tpm'] / 1000000 * 100
                all_results.append(df)
                print(f"✅ Processed {sample_base_name}: {len(df)} targets")
            except Exception as e:
                print(f"❌ Error processing {dir_path}: {e}")
        else:
            print(f"⚠️  abundance.tsv not found in {dir_path}")
    
    if all_results:
        # Объединить все результаты
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Сохранить сводную таблицу с процентами
        summary_df = combined_df.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
        summary_df.to_csv('summary_percentages.csv')
        
        print(f"📊 Created summary with {len(summary_df)} targets across {len(summary_df.columns)} samples")
        
        # Создать barplot для топ-10 targets по среднему проценту
        mean_percentages = combined_df.groupby('target_id')['percentage'].mean().sort_values(ascending=False)
        top_targets = mean_percentages.head(10).index.tolist()
        
        if len(top_targets) > 0:
            plot_data = combined_df[combined_df['target_id'].isin(top_targets)]
            
            # Создать barplot
            plt.figure(figsize=(12, 8))
            avg_data = plot_data.groupby('target_id')['percentage'].mean().sort_values(ascending=False)
            
            colors = plt.cm.viridis(np.linspace(0, 1, len(avg_data)))
            bars = plt.bar(range(len(avg_data)), avg_data.values, color=colors)
            
            plt.title('Top 10 Lineages by Average Percentage', fontsize=16, fontweight='bold')
            plt.xlabel('Lineage', fontsize=12)
            plt.ylabel('Percentage (%)', fontsize=12)
            plt.xticks(range(len(avg_data)), avg_data.index, rotation=45, ha='right')
            
            # Добавить значения на столбцы
            for i, (bar, value) in enumerate(zip(bars, avg_data.values)):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{value:.2f}%', ha='center', va='bottom', fontsize=10)
            
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig('top_lineages_barplot.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            print("🎨 Created barplot: top_lineages_barplot.png")
            
            # Создать heatmap для всех образцов и топ targets (только если много образцов)
            if len(plot_data['sample'].unique()) > 1:
                heatmap_data = plot_data.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
                
                # Ограничить количество колонок для читаемости
                if len(heatmap_data.columns) > 20:
                    # Взять образцы с наибольшим разнообразием
                    sample_variance = heatmap_data.var(axis=0).sort_values(ascending=False)
                    top_samples = sample_variance.head(20).index
                    heatmap_data = heatmap_data[top_samples]
                    print(f"📊 Limited heatmap to top 20 most variable samples")
                
                plt.figure(figsize=(max(10, len(heatmap_data.columns) * 0.4), 
                                  max(6, len(heatmap_data.index) * 0.3)))
                sns.heatmap(heatmap_data, annot=False, fmt='.1f', cmap='viridis', 
                           cbar_kws={'label': 'Percentage (%)'})
                plt.title('Lineage Abundance Across Samples (Top 10 Lineages)', 
                         fontsize=14, fontweight='bold')
                plt.xlabel('Sample', fontsize=12)
                plt.ylabel('Lineage', fontsize=12)
                plt.tight_layout()
                plt.savefig('abundance_heatmap.png', dpi=300, bbox_inches='tight')
                plt.close()
                
                print("🎨 Created heatmap: abundance_heatmap.png")
        
        print(f"✅ Processed {len(all_results)} samples total")
        print(f"🔬 Found {len(combined_df['target_id'].unique())} unique targets")
        print("📊 Created visualization plots and summary CSV")
        
        # Показать топ-5 lineages
        print("\\n🏆 Top 5 Lineages by Average Percentage:")
        top_5 = mean_percentages.head(5)
        for lineage, pct in top_5.items():
            print(f"   {lineage}: {pct:.3f}%")
            
    else:
        print("❌ No abundance files found")
        # Создать пустую сводную таблицу
        pd.DataFrame().to_csv('summary_percentages.csv')
        print("📄 Created empty summary_percentages.csv")
    """
}