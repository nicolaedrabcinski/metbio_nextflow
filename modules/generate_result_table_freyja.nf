nextflow.enable.dsl = 2

process GENERATE_RESULT_TABLE_FREYJA {
    tag "Generate Freyja results table"
    publishDir "${params.final_outdir}/plots", mode: 'copy'

    input:
    path freyja_agg_tsv

    output:
    path "results_table.csv", emit: results_table

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    input_file = 'aggregated_freyja.tsv'
    print(f"Generating results_table_freyja.csv from: {input_file}")
    freyja_df = pd.read_csv(input_file, sep='\t', header=0)
    sample_names = []
    all_lineages = set()
    sample_abundances = {}
    for idx, row in freyja_df.iterrows():
        sample = row[0].replace('_variants.tsv','').replace('sample_','')
        lineages = str(row['lineages']).split()
        abundances = [float(x) for x in str(row['abundances']).split()]
        sample_names.append(sample)
        for lin in lineages:
            all_lineages.add(lin)
        sample_abundances[sample] = dict(zip(lineages, abundances))
    all_lineages = sorted(list(all_lineages))
    # Map lineages to WHO names
    who_mapping = {
        'B.1.1.7': 'Alpha',
        'B.1.351': 'Beta',
        'P.1': 'Gamma',
        'B.1.617.2': 'Delta',
        'B.1.429': 'Epsilon',
        'B.1.526': 'Iota',
        'B.1.617.1': 'Kappa',
        'C.37': 'Lambda',
        'P.2': 'Zeta'
    }

    # Initialize WHO table
    who_columns = sorted(set(who_mapping.values()))
    who_table = pd.DataFrame(0, index=sample_names, columns=who_columns)
    
    # Convert lineages to WHO classifications
    for sample in sample_names:
        sample_data = {}
        # Initialize counts for each WHO variant
        for who_name in who_columns:
            sample_data[who_name] = 0
            
        # Sum abundances for each WHO variant
        for lineage, abundance in sample_abundances[sample].items():
            if lineage in who_mapping:
                who_name = who_mapping[lineage]
                sample_data[who_name] += abundance * 100  # Convert to percentage
                
        # Update the table
        for who_name, total_abundance in sample_data.items():
            who_table.loc[sample, who_name] = total_abundance
            
    # Save the table with WHO classifications
    who_table.index.name = 'sample'  # Use 'sample' as in the example
    who_table.to_csv('results_table.csv', float_format='%.3f')
    print(f"WHO variant table saved as: results_table.csv")
    """
}
