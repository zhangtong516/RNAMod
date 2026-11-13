#!/usr/bin/env python3
"""
Python version of STREME analysis script with modification type filtering.

This script processes RNA modification data, filters by modification type,
and runs STREME motif analysis with FIMO occurrence detection.
"""

import argparse
import csv
import os
import subprocess
import sys
import shutil
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run STREME analysis with modification type filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_streme.py --prefix batch2 --mod-type m6A
  python run_streme.py --prefix batch2 --mod-type m7G --window-size 150
  python run_streme.py --prefix batch2 --mod-type m6A --config custom_config.csv
        """
    )
    
    parser.add_argument(
        '--prefix',
        required=True,
        help='Output prefix for naming files and directories'
    )
    
    parser.add_argument(
        '--mod-type',
        required=True,
        choices=['m6A', 'm7G'],
        help='Modification type to filter (m6A or m7G)'
    )
    
    parser.add_argument(
        '--config',
        default='../batch2_config.csv',
        help='Path to configuration CSV file (default: ../batch2_config.csv)'
    )
    
    parser.add_argument(
        '--window-size',
        type=int,
        default=100,
        help='Window size around peak summit (default: 100)'
    )
    
    parser.add_argument(
        '--sample-rate',
        type=float,
        default=0.25,
        help='Sampling rate for sequences if too many (default: 0.25)'
    )
    
    parser.add_argument(
        '--fimo-threshold',
        type=float,
        default=0.01,
        help='FIMO p-value threshold (default: 0.01)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without executing'
    )
    
    return parser.parse_args()

def read_and_filter_config(config_file, mod_type):
    """Read CSV config file and filter by modification type."""
    logger.info(f"Reading config file: {config_file}")
    logger.info(f"Filtering for modification type: {mod_type}")
    
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")
    
    filtered_samples = []
    
    with open(config_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_name = row['sampleName']
            if mod_type in sample_name:
                filtered_samples.append({
                    'sample_name': sample_name,
                    'peak_file': row['peak_file'],
                    'liftover_file': row['liftover_file'],
                    'fasta_file': row['fasta_file']
                })
    
    logger.info(f"Found {len(filtered_samples)} samples with {mod_type} modification")
    for sample in filtered_samples:
        logger.info(f"  - {sample['sample_name']}")
    
    return filtered_samples

def run_command(cmd, dry_run=False, check=True):
    """Run a shell command with optional dry-run mode."""
    if dry_run:
        logger.info(f"[DRY RUN] Would execute: {cmd}")
        return None
    
    logger.info(f"Executing: {cmd}")
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=check,
            capture_output=True,
            text=True
        )
        if result.stdout:
            logger.debug(f"STDOUT: {result.stdout}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {cmd}")
        logger.error(f"Error: {e.stderr}")
        raise

def create_bed_files(samples, window_size, dry_run=False):
    """Create BED files from peak files."""
    logger.info("Creating BED files from peak files...")
    
    for sample in samples:
        sample_name = sample['sample_name']
        peak_file = sample['peak_file']
        bed_file = f"{sample_name}_positive_regions.bed"
        
        # Remove existing file if it exists
        if os.path.exists(bed_file) and not dry_run:
            os.remove(bed_file)
        
        # AWK command to create BED file
        awk_cmd = f"""awk -v window={window_size} 'BEGIN{{OFS="\\t"}} {{
            summit = $2 + $10;
            start = summit - window;
            end = summit + window;
            if(start < 0) start = 0;
            print $1, start, end, $4"_plus", $5, "+";
            print $1, start, end, $4"_minus", $5, "-";
        }}' {peak_file} > {bed_file}"""
        
        run_command(awk_cmd, dry_run)

def extract_sequences(samples, dry_run=False):
    """Extract sequences using bedtools."""
    logger.info("Extracting sequences using bedtools...")
    
    for sample in samples:
        sample_name = sample['sample_name']
        fasta_file = sample['fasta_file']
        bed_file = f"{sample_name}_positive_regions.bed"
        seq_file = f"{sample_name}_positive_sequences.fa"
        
        # Remove existing sequence file if it exists
        if os.path.exists(seq_file) and not dry_run:
            os.remove(seq_file)
        
        # Bedtools command
        bedtools_cmd = f"bedtools getfasta -s -fi {fasta_file} -bed {bed_file} -fo {seq_file}"
        run_command(bedtools_cmd, dry_run)

def combine_and_sample_sequences(prefix, mod_type, sample_rate, dry_run=False):
    """Combine all sequences and sample if necessary."""
    logger.info("Combining sequences...")
    
    combined_file = f"{prefix}_{mod_type}_sequences_combined.fa"
    subsample_file = f"{prefix}_{mod_type}_sequences_combined_subsample.fa"
    
    # Remove existing files
    for file in [combined_file, subsample_file]:
        if os.path.exists(file) and not dry_run:
            os.remove(file)
    
    # Combine all positive sequence files
    combine_cmd = f"cat *_positive_sequences.fa > {combined_file}"
    run_command(combine_cmd, dry_run)
    
    # Sample sequences using seqkit
    logger.info(f"Sampling sequences with rate: {sample_rate}")
    sample_cmd = f"seqkit sample -p {sample_rate} {combined_file} > {subsample_file}"
    run_command(sample_cmd, dry_run)
    
    return subsample_file

def run_streme_analysis(input_file, prefix, mod_type, dry_run=False):
    """Run STREME motif analysis."""
    logger.info("Running STREME analysis...")
    
    output_dir = f"streme_{prefix}_{mod_type}"
    
    # Remove existing output directory
    if os.path.exists(output_dir) and not dry_run:
        shutil.rmtree(output_dir)
    
    # Load singularity module and run STREME
    meme_image = "/home/users/astar/gis/zhangt/scratch/tools/meme_v5.5.8.sif"
    
    streme_cmd = f"""module load singularity && \\
singularity run {meme_image} streme \\
    --p {input_file} \\
    --rna \\
    --minw 5 \\
    --maxw 6 \\
    --nmotifs 15 \\
    --thresh 0.01 \\
    --niter 100 \\
    --order 2 \\
    --seed 123 \\
    --oc {output_dir}"""
    
    run_command(streme_cmd, dry_run)
    
    return output_dir

def run_fimo_analysis(streme_dir, prefix, mod_type, fimo_threshold, dry_run=False):
    """Run FIMO motif occurrence detection."""
    logger.info("Running FIMO analysis...")
    
    streme_xml = os.path.join(streme_dir, "streme.xml")
    if not os.path.exists(streme_xml) and not dry_run:
        logger.warning(f"STREME XML not found: {streme_xml}")
        return
    
    # Create FIMO output directory
    fimo_output_dir = f"{prefix}_{mod_type}_fimo_results"
    if os.path.exists(fimo_output_dir) and not dry_run:
        shutil.rmtree(fimo_output_dir)
    
    if not dry_run:
        os.makedirs(fimo_output_dir, exist_ok=True)
    
    # Find all positive sequence files for this mod type
    pos_files = []
    for file in os.listdir('.'):
        if file.endswith('_positive_sequences.fa'):
            # Check if this file contains the mod_type
            sample_name = file.replace('_positive_sequences.fa', '')
            if mod_type in sample_name:
                pos_files.append(file)
    
    logger.info(f"Found {len(pos_files)} sequence files for FIMO analysis")
    
    meme_image = "/home/users/astar/gis/zhangt/scratch/tools/meme_v5.5.8.sif"
    
    # Process each sequence file with FIMO
    for seq_file in pos_files:
        sample_name = seq_file.replace('_positive_sequences.fa', '')
        sample_output = os.path.join(fimo_output_dir, sample_name)
        
        if not dry_run:
            os.makedirs(sample_output, exist_ok=True)
        
        fimo_cmd = f"""module load singularity && \\
singularity run {meme_image} fimo \\
    --thresh {fimo_threshold} \\
    --oc {sample_output} \\
    --parse-genomic-coord \\
    {streme_xml} \\
    {seq_file}"""
        
        run_command(fimo_cmd, dry_run, check=False)  # Don't fail if no motifs found
    
    # Aggregate FIMO results
    aggregate_fimo_results(fimo_output_dir, prefix, mod_type, pos_files, dry_run)

def aggregate_fimo_results(fimo_output_dir, prefix, mod_type, pos_files, dry_run=False):
    """Aggregate FIMO results across all samples."""
    logger.info("Aggregating FIMO results...")
    
    aggregated_file = f"{prefix}_{mod_type}_fimo_occurrences_aggregated.tsv"
    
    if dry_run:
        logger.info(f"[DRY RUN] Would create aggregated file: {aggregated_file}")
        return
    
    # Create header
    with open(aggregated_file, 'w') as f:
        f.write("sample_name\tmotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n")
    
    total_occurrences = 0
    
    # Combine all results
    for seq_file in pos_files:
        sample_name = seq_file.replace('_positive_sequences.fa', '')
        sample_output = os.path.join(fimo_output_dir, sample_name)
        fimo_tsv = os.path.join(sample_output, "fimo.tsv")
        
        if os.path.exists(fimo_tsv):
            # Read FIMO results and add sample name
            with open(fimo_tsv, 'r') as f_in:
                lines = f_in.readlines()
                if len(lines) > 1:  # Has data beyond header
                    with open(aggregated_file, 'a') as f_out:
                        for line in lines[1:]:  # Skip header
                            f_out.write(f"{sample_name}\t{line}")
                            total_occurrences += 1
    
    logger.info(f"Total motif occurrences found: {total_occurrences}")
    
    if total_occurrences > 0:
        # Create summary files
        create_summary_statistics(aggregated_file, prefix, mod_type)

def create_summary_statistics(aggregated_file, prefix, mod_type):
    """Create summary statistics showing motif occurrence ratios per sample."""
    logger.info("Creating motif occurrence ratio statistics...")
    
    try:
        import pandas as pd
        import numpy as np
        use_pandas = True
        logger.info("Using pandas for motif ratio analysis")
    except ImportError:
        logger.warning("pandas not available, falling back to awk-based analysis")
        use_pandas = False
    
    if use_pandas:
        create_motif_ratio_statistics_pandas(aggregated_file, prefix, mod_type)
    else:
        create_motif_ratio_statistics_awk(aggregated_file, prefix, mod_type)

def create_motif_ratio_statistics_pandas(aggregated_file, prefix, mod_type):
    """Create motif ratio statistics using pandas."""
    import pandas as pd
    import numpy as np
    
    # Read the aggregated FIMO results
    logger.info(f"Reading aggregated results from: {aggregated_file}")
    
    try:
        df = pd.read_csv(aggregated_file, sep='\t')
        logger.info(f"Loaded {len(df)} motif occurrences")
    except Exception as e:
        logger.error(f"Failed to read aggregated file: {e}")
        return
    
    if len(df) == 0:
        logger.warning("No data found in aggregated file")
        return
    
    # Get total number of peaks (unique sequences) per sample
    logger.info("Calculating total peaks per sample...")
    
    # Count unique sequences per sample (these represent the peaks)
    total_peaks_per_sample = df.groupby('sample_name')['sequence_name'].nunique().to_dict()
    
    logger.info("Total peaks per sample:")
    for sample, count in total_peaks_per_sample.items():
        logger.info(f"  {sample}: {count} peaks")
    
    # Calculate motif occurrence ratios
    logger.info("Calculating motif occurrence ratios...")
    
    # For each motif-sample combination, count unique sequences with that motif
    motif_sample_peaks = df.groupby(['motif_id', 'sample_name'])['sequence_name'].nunique().reset_index()
    motif_sample_peaks.columns = ['motif_id', 'sample_name', 'peaks_with_motif']
    
    # Add total peaks for each sample
    motif_sample_peaks['total_peaks'] = motif_sample_peaks['sample_name'].map(total_peaks_per_sample)
    
    # Calculate ratio
    motif_sample_peaks['motif_ratio'] = motif_sample_peaks['peaks_with_motif'] / motif_sample_peaks['total_peaks']
    motif_sample_peaks['motif_percentage'] = motif_sample_peaks['motif_ratio'] * 100
    
    # Create pivot table for easy viewing: motifs as rows, samples as columns
    ratio_matrix = motif_sample_peaks.pivot(index='motif_id', columns='sample_name', values='motif_ratio').fillna(0)
    percentage_matrix = motif_sample_peaks.pivot(index='motif_id', columns='sample_name', values='motif_percentage').fillna(0)
    peaks_matrix = motif_sample_peaks.pivot(index='motif_id', columns='sample_name', values='peaks_with_motif').fillna(0)
    
    # Add summary statistics to matrices
    ratio_matrix['mean_ratio'] = ratio_matrix.mean(axis=1)
    ratio_matrix['std_ratio'] = ratio_matrix.std(axis=1)
    ratio_matrix['max_ratio'] = ratio_matrix.iloc[:, :-2].max(axis=1)
    ratio_matrix['samples_with_motif'] = (ratio_matrix.iloc[:, :-3] > 0).sum(axis=1)
    
    percentage_matrix['mean_percentage'] = percentage_matrix.mean(axis=1)
    percentage_matrix['std_percentage'] = percentage_matrix.std(axis=1)
    percentage_matrix['max_percentage'] = percentage_matrix.iloc[:, :-2].max(axis=1)
    
    peaks_matrix['total_peaks_with_motif'] = peaks_matrix.sum(axis=1)
    
    # Save results
    ratio_file = f"{prefix}_{mod_type}_motif_ratio_matrix.tsv"
    percentage_file = f"{prefix}_{mod_type}_motif_percentage_matrix.tsv"
    peaks_file = f"{prefix}_{mod_type}_motif_peaks_matrix.tsv"
    detailed_file = f"{prefix}_{mod_type}_motif_ratio_detailed.tsv"
    
    # Round to reasonable precision
    ratio_matrix.round(4).to_csv(ratio_file, sep='\t')
    percentage_matrix.round(2).to_csv(percentage_file, sep='\t')
    peaks_matrix.astype(int).to_csv(peaks_file, sep='\t')
    
    # Create detailed table with all information
    detailed_stats = motif_sample_peaks.copy()
    detailed_stats = detailed_stats.round({'motif_ratio': 4, 'motif_percentage': 2})
    detailed_stats = detailed_stats.sort_values(['motif_id', 'sample_name'])
    detailed_stats.to_csv(detailed_file, sep='\t', index=False)
    
    logger.info(f"Motif ratio analysis files created:")
    logger.info(f"  - Ratio matrix: {ratio_file}")
    logger.info(f"  - Percentage matrix: {percentage_file}")
    logger.info(f"  - Peaks count matrix: {peaks_file}")
    logger.info(f"  - Detailed statistics: {detailed_file}")
    
    # Create summary report
    summary_file = f"{prefix}_{mod_type}_motif_ratio_summary.tsv"
    
    # Calculate overall statistics per motif
    motif_summary = []
    for motif in ratio_matrix.index:
        motif_data = detailed_stats[detailed_stats['motif_id'] == motif]
        
        total_peaks_with_motif = motif_data['peaks_with_motif'].sum()
        total_peaks_all_samples = motif_data['total_peaks'].sum()
        overall_ratio = total_peaks_with_motif / total_peaks_all_samples if total_peaks_all_samples > 0 else 0
        
        samples_with_motif = (motif_data['peaks_with_motif'] > 0).sum()
        total_samples = len(motif_data)
        
        mean_ratio = motif_data['motif_ratio'].mean()
        std_ratio = motif_data['motif_ratio'].std()
        max_ratio = motif_data['motif_ratio'].max()
        min_ratio = motif_data['motif_ratio'].min()
        
        # Find sample with highest ratio
        max_sample = motif_data.loc[motif_data['motif_ratio'].idxmax(), 'sample_name']
        max_ratio_value = motif_data['motif_ratio'].max()
        
        motif_summary.append({
            'motif_id': motif,
            'total_peaks_with_motif': total_peaks_with_motif,
            'total_peaks_all_samples': total_peaks_all_samples,
            'overall_ratio': overall_ratio,
            'overall_percentage': overall_ratio * 100,
            'samples_with_motif': samples_with_motif,
            'total_samples': total_samples,
            'sample_coverage': samples_with_motif / total_samples,
            'mean_ratio_across_samples': mean_ratio,
            'std_ratio_across_samples': std_ratio,
            'min_ratio': min_ratio,
            'max_ratio': max_ratio,
            'best_sample': max_sample,
            'best_sample_ratio': max_ratio_value,
            'best_sample_percentage': max_ratio_value * 100
        })
    
    summary_df = pd.DataFrame(motif_summary)
    summary_df = summary_df.round({
        'overall_ratio': 4, 'overall_percentage': 2,
        'sample_coverage': 3, 'mean_ratio_across_samples': 4,
        'std_ratio_across_samples': 4, 'min_ratio': 4, 'max_ratio': 4,
        'best_sample_ratio': 4, 'best_sample_percentage': 2
    })
    
    # Sort by overall ratio (most prevalent motifs first)
    summary_df = summary_df.sort_values('overall_ratio', ascending=False)
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    logger.info(f"  - Summary report: {summary_file}")
    
    # Log key findings
    logger.info("\n=== Motif Ratio Analysis Summary ===")
    logger.info(f"Analyzed {len(summary_df)} motifs across {len(total_peaks_per_sample)} samples")
    
    if len(summary_df) > 0:
        logger.info(f"\nTop 5 most prevalent motifs (by overall ratio):")
        top_motifs = summary_df.head(5)
        for _, row in top_motifs.iterrows():
            logger.info(f"  {row['motif_id']}: {row['overall_percentage']:.2f}% of peaks ({row['total_peaks_with_motif']}/{row['total_peaks_all_samples']})")
            logger.info(f"    Best in {row['best_sample']}: {row['best_sample_percentage']:.2f}%")
        
        logger.info(f"\nMotifs with highest sample coverage:")
        high_coverage = summary_df.nlargest(5, 'sample_coverage')
        for _, row in high_coverage.iterrows():
            logger.info(f"  {row['motif_id']}: found in {row['samples_with_motif']}/{row['total_samples']} samples ({row['sample_coverage']:.1%})")

def create_motif_ratio_statistics_awk(aggregated_file, prefix, mod_type):
    """Create motif ratio statistics using awk (fallback)."""
    logger.info("Creating motif ratio statistics with awk...")
    
    ratio_file = f"{prefix}_{mod_type}_motif_ratio_summary.tsv"
    
    # Use awk to calculate ratios
    awk_cmd = f"""awk -F'\\t' '
    NR > 1 {{
        motif = $2
        sample = $1
        sequence = $4
        
        # Count unique sequences per motif-sample combination
        motif_sample_sequences[motif][sample][sequence] = 1
        # Count total sequences per sample
        sample_sequences[sample][sequence] = 1
    }}
    END {{
        print "motif_id\\tsample_name\\tpeaks_with_motif\\ttotal_peaks\\tmotif_ratio\\tmotif_percentage"
        
        for (motif in motif_sample_sequences) {{
            for (sample in motif_sample_sequences[motif]) {{
                # Count sequences with this motif in this sample
                peaks_with_motif = 0
                for (seq in motif_sample_sequences[motif][sample]) {{
                    peaks_with_motif++
                }}
                
                # Count total sequences in this sample
                total_peaks = 0
                for (seq in sample_sequences[sample]) {{
                    total_peaks++
                }}
                
                # Calculate ratio
                ratio = (total_peaks > 0) ? peaks_with_motif / total_peaks : 0
                percentage = ratio * 100
                
                print motif "\\t" sample "\\t" peaks_with_motif "\\t" total_peaks "\\t" ratio "\\t" percentage
            }}
        }}
    }}' {aggregated_file} > {ratio_file}"""
    
    run_command(awk_cmd)
    
    logger.info(f"Motif ratio statistics created: {ratio_file}")

def main():
    """Main function."""
    args = parse_arguments()
    
    logger.info("=== STREME Analysis with Python ===")
    logger.info(f"Prefix: {args.prefix}")
    logger.info(f"Modification type: {args.mod_type}")
    logger.info(f"Config file: {args.config}")
    logger.info(f"Window size: {args.window_size}")
    logger.info(f"Sample rate: {args.sample_rate}")
    logger.info(f"FIMO threshold: {args.fimo_threshold}")
    logger.info(f"Dry run: {args.dry_run}")
    logger.info("")
    
    try:
        # Step 1: Read and filter configuration
        samples = read_and_filter_config(args.config, args.mod_type)
        
        if not samples:
            logger.error(f"No samples found for modification type: {args.mod_type}")
            sys.exit(1)
        
        # Step 2: Create BED files
        create_bed_files(samples, args.window_size, args.dry_run)
        
        # Step 3: Extract sequences
        extract_sequences(samples, args.dry_run)
        
        # Step 4: Combine and sample sequences
        input_file = combine_and_sample_sequences(
            args.prefix, args.mod_type, args.sample_rate, args.dry_run
        )
        
        # Step 5: Run STREME analysis
        streme_dir = run_streme_analysis(
            input_file, args.prefix, args.mod_type, args.dry_run
        )
        
        # Step 6: Run FIMO analysis
        run_fimo_analysis(
            streme_dir, args.prefix, args.mod_type, args.fimo_threshold, args.dry_run
        )
        
        logger.info("=== Analysis Complete ===")
        
        if not args.dry_run:
            logger.info("Output files:")
            logger.info(f"  - STREME results: {streme_dir}/")
            logger.info(f"  - FIMO results: {args.prefix}_{args.mod_type}_fimo_results/")
            logger.info(f"  - Aggregated occurrences: {args.prefix}_{args.mod_type}_fimo_occurrences_aggregated.tsv")
            logger.info(f"  - Motif summary: {args.prefix}_{args.mod_type}_fimo_summary_by_motif.tsv")
            logger.info(f"  - Sample summary: {args.prefix}_{args.mod_type}_fimo_summary_by_sample.tsv")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()