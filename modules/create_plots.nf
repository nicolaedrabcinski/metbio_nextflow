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
    path "results_table.csv", optional: true, emit: table
    
    when:
    params.create_plots
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    import glob
    import re
    import ast
    from pathlib import Path
    
    print("Starting visualization creation...")
    
    # WHO variant names mapping for COVID variants
    who_names_mapping = {
        'B.1.1.7': 'Alpha',
        'B.1.351': 'Beta', 
        'P.1': 'Gamma',
        'B.1.617.2': 'Delta',
        'B.1.1.529': 'Omicron',
        'BA.1': 'Omicron BA.1',
        'BA.2': 'Omicron BA.2',
        'BA.3': 'Omicron BA.3',
        'BA.4': 'Omicron BA.4',
        'BA.5': 'Omicron BA.5',
        'BQ.1': 'Omicron BQ.1',
        'XBB': 'Omicron XBB',
        'XBB.1.5': 'Omicron XBB.1.5',
        'EG.5': 'Omicron EG.5',
        'JN.1': 'Omicron JN.1',
        'KP.2': 'Omicron KP.2',
        'KP.3': 'Omicron KP.3',
        'B.1.429': 'Epsilon',
        'B.1.427': 'Epsilon',
        'P.2': 'Zeta',
        'B.1.525': 'Eta',
        'B.1.526': 'Iota',
        'B.1.617.1': 'Kappa',
        'C.37': 'Lambda',
        'B.1.621': 'Mu',
        'P.3': 'Theta',
        
    }
    
    def get_who_name(target_name):
        if target_name in who_names_mapping:
            return who_names_mapping[target_name]
        
        # Check for partial matches for sublineages
        for lineage, who_name in who_names_mapping.items():
            if target_name.startswith(lineage + '.'):
                return f"{who_name} ({target_name})"
        
        return target_name
    
    def sort_sample_names(sample_names):
        # Sort sample names numerically when possible, with robust fallback
        def extract_number(sample_name):
            # Try to extract numbers from sample names for sorting
            import re
            numbers = re.findall(r'\\d+', str(sample_name))
            if numbers:
                # Use the first number found, or concatenate if multiple
                return int(numbers[0]) if len(numbers) == 1 else int(''.join(numbers[:2]))
            return float('inf')  # Non-numeric samples go to the end
        
        try:
            # Sort by extracted numbers, then alphabetically for ties
            return sorted(sample_names, key=lambda x: (extract_number(x), str(x)))
        except:
            # Fallback to alphabetical sorting if numerical sorting fails
            return sorted(sample_names)
    
    # Read list of result paths
    results_paths = []
    with open('${results_file}', 'r') as f:
        for line in f:
            if line.strip():
                unique_name, dir_path = line.strip().split('\\t')
                results_paths.append((unique_name, dir_path))
    
    print(f"Found {len(results_paths)} result directories")
    
    # Collect all results from Kallisto files
    all_results = []
    processed_samples = set()
    
    for unique_name, dir_path in results_paths:
        # Extract original sample name from directory path
        sample_base_name = os.path.basename(dir_path)
        
        # If sample already processed, add suffix
        original_name = sample_base_name
        counter = 1
        while sample_base_name in processed_samples:
            sample_base_name = f"{original_name}_v{counter}"
            counter += 1
        
        processed_samples.add(sample_base_name)
        
        # Look for Kallisto abundance.tsv files
        abundance_file = os.path.join(dir_path, 'abundance.tsv')
        
        found_data = False
        
        # Process Kallisto abundance.tsv file
        if os.path.exists(abundance_file):
            try:
                df = pd.read_csv(abundance_file, sep='\\t')
                df['sample'] = sample_base_name
                # Convert TPM to percentages: % = TPM / 1,000,000 * 100
                # df['percentage'] = df['tpm'] / 1000000 * 100
                df['percentage'] = df['tpm'] / df['tpm'].sum() * 100
                # Apply WHO names mapping for COVID variants, keep original names for genes
                df['target_id'] = df['target_id'].apply(get_who_name)
                all_results.extend(df[['sample', 'target_id', 'percentage']].to_dict('records'))
                print(f"Processed {sample_base_name} from abundance.tsv: {len(df)} targets")
                found_data = True
            except Exception as e:
                print(f"Error processing {dir_path}: {e}")
        
        if not found_data:
            print(f"No abundance.tsv file found in {dir_path}")
    
    if all_results:
        # Convert to DataFrame
        combined_df = pd.DataFrame(all_results)
        
        # Save summary table with percentages
        summary_df = combined_df.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
        
        # Sort sample columns numerically in summary CSV
        sorted_summary_samples = sort_sample_names(summary_df.columns.tolist())
        summary_df = summary_df[sorted_summary_samples]
        summary_df.to_csv('summary_percentages.csv')
        
        print(f"Created summary with {len(summary_df)} targets across {len(summary_df.columns)} samples")
        
        # Create barplot for top targets by average percentage
        mean_percentages = combined_df.groupby('target_id')['percentage'].mean().sort_values(ascending=False)
        top_targets = mean_percentages.head(10).index.tolist()
        
        if len(top_targets) > 0:
            plot_data = combined_df[combined_df['target_id'].isin(top_targets)]
            
            # Create barplot
            plt.figure(figsize=(14, 8))
            avg_data = plot_data.groupby('target_id')['percentage'].mean().sort_values(ascending=False)
            
            # Use more contrasting color scheme
            colors = plt.cm.plasma(np.linspace(0, 1, len(avg_data)))
            bars = plt.bar(range(len(avg_data)), avg_data.values, color=colors, edgecolor='black', linewidth=0.5)
            
            plt.title('Top Targets by Average Percentage', fontsize=16, fontweight='bold')
            plt.xlabel('Target Name', fontsize=12)
            plt.ylabel('Percentage (%)', fontsize=12)
            plt.xticks(range(len(avg_data)), avg_data.index, rotation=45, ha='right')
            
            # Add values on bars
            for i, (bar, value) in enumerate(zip(bars, avg_data.values)):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(avg_data.values) * 0.01,
                        f'{value:.2f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig('lineages_barplot.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            print("Created barplot: lineages_barplot.png")
            
            # Create heatmap for all samples and top targets (only if multiple samples)
            if len(plot_data['sample'].unique()) > 1:
                heatmap_data = plot_data.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
                
                # Limit columns for readability
                if len(heatmap_data.columns) > 20:
                    # Take samples with highest variance
                    sample_variance = heatmap_data.var(axis=0).sort_values(ascending=False)
                    top_samples = sample_variance.head(20).index
                    heatmap_data = heatmap_data[top_samples]
                    print(f"Limited heatmap to 20 most variable samples")
                
                # Sort sample columns numerically
                sorted_samples = sort_sample_names(heatmap_data.columns.tolist())
                heatmap_data = heatmap_data[sorted_samples]
                
                plt.figure(figsize=(max(12, len(heatmap_data.columns) * 0.6), 
                                  max(8, len(heatmap_data.index) * 0.5)))
                
                # Use more contrasting color scheme and add values in cells
                sns.heatmap(heatmap_data, 
                           annot=True,  # Show values in cells
                           fmt='.1f',   # Number format (1 decimal place)
                           cmap='RdYlBu_r',  # More contrasting color scheme
                           cbar_kws={'label': 'Percentage (%)'},
                           linewidths=0.5,  # Add lines between cells
                           linecolor='white',
                           annot_kws={'size': 8, 'weight': 'bold'},  # Text settings for cells
                           square=False,  # Don't make cells square
                           robust=True)   # Use robust color normalization
                
                plt.title('Target Abundance Across Samples', 
                         fontsize=14, fontweight='bold', pad=20)
                plt.xlabel('Sample', fontsize=12)
                plt.ylabel('Target Name', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig('abundance_heatmap.png', dpi=300, bbox_inches='tight')
                plt.close()
                
                print("Created heatmap: abundance_heatmap.png")
                
                # Create additional heatmap with all variants (if not too many)
                if len(combined_df['target_id'].unique()) <= 15:
                    full_heatmap_data = combined_df.pivot(index='target_id', columns='sample', values='percentage').fillna(0)
                    
                    # Sort sample columns numerically
                    sorted_samples_full = sort_sample_names(full_heatmap_data.columns.tolist())
                    full_heatmap_data = full_heatmap_data[sorted_samples_full]
                    
                    plt.figure(figsize=(max(12, len(full_heatmap_data.columns) * 0.6), 
                                      max(10, len(full_heatmap_data.index) * 0.4)))
                    
                    sns.heatmap(full_heatmap_data, 
                               annot=True, 
                               fmt='.2f', 
                               cmap='viridis',  # Alternative contrasting scheme
                               cbar_kws={'label': 'Percentage (%)'},
                               linewidths=0.3,
                               linecolor='gray',
                               annot_kws={'size': 7, 'weight': 'bold'},
                               square=False)
                    
                    plt.title('Complete Target Abundance Matrix', 
                             fontsize=14, fontweight='bold', pad=20)
                    plt.xlabel('Sample', fontsize=12)
                    plt.ylabel('Target Name', fontsize=12)
                    plt.xticks(rotation=45, ha='right')
                    plt.yticks(rotation=0)
                    plt.tight_layout()
                    plt.savefig('complete_abundance_heatmap.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    print("Created complete heatmap: complete_abundance_heatmap.png")
        
        print(f"Processed {len(all_results)} samples total")
        print(f"Found {len(combined_df['target_id'].unique())} unique targets")
        print("Created visualization plots and summary CSV")
        
        # Create and display formatted results table
        print("\\n" + "="*80)
        print("TARGET ABUNDANCE RESULTS TABLE")
        print("="*80)
        
        # Create results table with targets as columns and samples as rows
        results_table = summary_df.T  # Transpose to have samples as rows, targets as columns
        
        # Filter to show only targets with >0.1% in at least one sample
        significant_targets = []
        for col in results_table.columns:
            if results_table[col].max() > 0.1:  # Only targets with >0.1% abundance
                significant_targets.append(col)
        
        if significant_targets:
            display_table = results_table[significant_targets]
            
            # Save formatted table as CSV
            display_table.to_csv('results_table.csv', float_format='%.3f')
            
            # Display table in console with nice formatting
            print(f"\\nResults for {len(display_table)} samples and {len(significant_targets)} significant targets:")
            print("-" * 80)
            
            # Create header
            header = "Sample".ljust(10)
            for target in display_table.columns:
                header += target.ljust(12)
            print(header)
            print("-" * len(header))
            
            # Display data rows
            for sample in display_table.index:
                row = sample.ljust(10)
                for target in display_table.columns:
                    value = display_table.loc[sample, target]
                    row += f"{value:.3f}%".ljust(12)
                print(row)
            
            print("-" * len(header))
            print(f"Table saved as: results_table.csv")
        else:
            print("\\nNo significant targets found (>0.1% abundance)")
        
        # Show top 5 targets summary
        print("\\nTop 5 Targets by Average Percentage:")
        top_5 = mean_percentages.head(5)
        for target, pct in top_5.items():
            print(f"   {target}: {pct:.3f}%")
            
    else:
        print("No suitable data files found")
        # Create empty summary tables
        pd.DataFrame().to_csv('summary_percentages.csv')
        pd.DataFrame().to_csv('results_table.csv')
        print("Created empty summary files")
    """
}