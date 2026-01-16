process CREATE_WHO_PLOTS {
    tag "Creating WHO classification plots"
    publishDir "${params.final_outdir}/who_plots", mode: 'copy'
    container 'jupyter/scipy-notebook:latest'
    containerOptions '--user root'
    beforeScript '''
        export MPLCONFIGDIR=/tmp/matplotlib
        export TMPDIR=/tmp
    '''
    errorStrategy 'retry'
    maxRetries 2

    input:
    path results_table

    output:
    path "who_plot_stacked.pdf", emit: plot
    path "who_plot_stacked.png", emit: plot_png

    script:
    """
#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Set style
plt.style.use('default')  # Use default style instead of seaborn
# Set figure style manually
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

try:
    # Read the results table
    df = pd.read_csv("${results_table}")
    print("File contents:")
    print(df.head())
    print("Columns:", df.columns.tolist())
    
    # Set up the colors for WHO variants using standard, high-contrast colors
    colors = {
        'Alpha': '#E41A1C',    # Red
        'Beta': '#377EB8',     # Blue
        'Gamma': '#4DAF4A',    # Green
        'Delta': '#FF7F00',    # Orange
        'Epsilon': '#984EA3',  # Purple
        'Iota': '#F781BF',     # Pink
        'Kappa': '#A65628',    # Brown
        'Lambda': '#FFD700',   # Gold
        'Zeta': '#000000'      # Black
    }
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Get variant columns (exclude 'sample' column)
    variant_cols = [col for col in df.columns if col != 'sample']
    
    # Create the stacked bar plot
    bottom = np.zeros(len(df))
    bars = []
    labels = []
    
    # Plot each variant as a stacked component
    for variant in variant_cols:
        if df[variant].sum() > 0:  # Only plot variants that are present
            mask = df[variant] > 0  # Create mask for non-zero values
            bar = plt.bar(range(len(df)), df[variant], bottom=bottom, 
                         label=variant, color=colors.get(variant, '#808080'))
            bottom += df[variant]
            bars.append(bar)
            labels.append(variant)
    
    # Customize plot
    plt.title('SARS-CoV-2 WHO Variant Distribution by Sample', fontsize=14, pad=20)
    plt.xlabel('Sample', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    
    # Set sample names as x-axis labels
    plt.xticks(range(len(df)), df['sample'], rotation=45, ha='right')
    
    # Add legend
    plt.legend(title='WHO Variants', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Removed percentage labels for cleaner visualization
    
    # Adjust layout to prevent label cutoff
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save plots
    plt.savefig('who_plot_stacked.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('who_plot_stacked.png', dpi=300, bbox_inches='tight')
    print("Successfully created stacked bar plot")

except Exception as e:
    print(f"Error: {str(e)}")
    sys.exit(1)
"""
}