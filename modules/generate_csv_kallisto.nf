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

    seen = set()
    fastq_files = []

    # Search for FASTQ files with various extensions
    extensions = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']

    # Recursive search in all subdirectories
    filter_shuffled = ${params.filter_shuffled ? 'True' : 'False'}
    for ext in extensions:
        pattern = os.path.join('${fastq_dir}', '**', ext)
        found_files = glob.glob(pattern, recursive=True)

        # Filter only files with _shuffled in the name if requested
        if filter_shuffled:
            found_files = [f for f in found_files if '_shuffled' in os.path.basename(f)]

        # Deduplicate: skip files already added by a previous extension pattern
        for f in found_files:
            real = os.path.realpath(f)
            if real not in seen:
                seen.add(real)
                fastq_files.append(f)

    # If no shuffled files found and filter is enabled, use all files
    if not fastq_files and filter_shuffled:
        print("No '_shuffled' files found, using all FASTQ files...")
        for ext in extensions:
            pattern = os.path.join('${fastq_dir}', '**', ext)
            for f in glob.glob(pattern, recursive=True):
                real = os.path.realpath(f)
                if real not in seen:
                    seen.add(real)
                    fastq_files.append(f)

    # Sort for stable ordering
    fastq_files.sort()

    print(f"Found {len(fastq_files)} FASTQ files:")
    for f in fastq_files[:10]:
        print(f"  - {f}")
    if len(fastq_files) > 10:
        print(f"  ... and {len(fastq_files) - 10} more files")

    # Create CSV file
    with open('auto_samples.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'fastq'])

        for fastq_file in fastq_files:
            # Extract sample name from filename
            sample_name = Path(fastq_file).stem

            # Remove common suffixes
            suffixes_to_remove = ['.fastq', '.fq', '.gz', '_shuffled', '_long', '_short']
            for suffix in suffixes_to_remove:
                if sample_name.endswith(suffix):
                    sample_name = sample_name[:-len(suffix)]

            # Use absolute path to file
            abs_path = os.path.abspath(fastq_file)
            writer.writerow([sample_name, abs_path])

    print(f"Generated CSV file with {len(fastq_files)} samples")
    print("Sample names preview:")
    with open('auto_samples.csv', 'r') as f:
        lines = f.readlines()
        for line in lines[1:6]:
            print(f"  {line.strip()}")
        if len(lines) > 6:
            print(f"  ... and {len(lines) - 6} more samples")
    """
}
