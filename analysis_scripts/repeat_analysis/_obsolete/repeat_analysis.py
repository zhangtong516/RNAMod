#!/usr/bin/env python3
"""
M6A Signal Analysis on Repeat Elements
Analyze m6A modification patterns across different repeat element families
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
import subprocess
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class M6ARepeatAnalyzer:
    def __init__(self, species_data, repeat_annotations, bin_size=50, extend_region=2000):
        """
        Initialize M6A repeat analyzer
        
        Args:
            species_data: dict {species: {'peaks': peak_file, 'bigwig': bigwig_file}}
            repeat_annotations: dict {species: repeat_annotation_file}
            bin_size: Size of bins for signal calculation (bp)
            extend_region: Extension around repeat elements (bp)
        """
        self.species_data = species_data
        self.repeat_annotations = repeat_annotations
        self.bin_size = bin_size
        self.extend_region = extend_region
        self.species_names = list(species_data.keys())
        
        # Define repeat family classifications
        self.repeat_families = {
            'LINE': ['LINE', 'L1', 'L2', 'CR1', 'RTE', 'Penelope'],
            'SINE': ['SINE', 'Alu', 'MIR', 'B1', 'B2', 'B4', 'ID'],
            'LTR': ['LTR', 'ERV', 'ERVL', 'ERVK', 'Gypsy', 'Copia'],
            'DNA': ['DNA', 'hAT', 'TcMar', 'Merlin', 'Transib', 'PiggyBac'],
            'Satellite': ['Satellite', 'centr', 'macro'],
            'Simple_repeat': ['Simple_repeat', 'Low_complexity'],
            'Other': ['Unknown', 'Other', 'Unclassified']
        }
        
    def load_repeat_annotations(self, species):
        """Load and classify repeat annotations for a species"""
        repeat_file = self.repeat_annotations[species]
        
        # Try different formats (RepeatMasker out, BED, GTF)
        try:
            if repeat_file.endswith('.out'):
                # RepeatMasker .out format
                repeats = self.parse_repeatmasker_out(repeat_file)
            elif repeat_file.endswith('.gtf') or repeat_file.endswith('.gff'):
                # GTF/GFF format
                repeats = self.parse_gtf_repeats(repeat_file)
            else:
                # Assume BED format
                repeats = pd.read_csv(repeat_file, sep='\t', header=None,
                                    names=['chr', 'start', 'end', 'name', 'score', 
                                          'strand', 'family', 'class'])
        except Exception as e:
            print(f"Error loading repeats for {species}: {e}")
            return pd.DataFrame()
        
        # Classify repeat families
        repeats['family_class'] = repeats.apply(self.classify_repeat_family, axis=1)
        
        return repeats
    
    def parse_repeatmasker_out(self, out_file):
        """Parse RepeatMasker .out file"""
        repeats = []
        with open(out_file, 'r') as f:
            # Skip header lines
            for _ in range(3):
                next(f)
            
            for line in f:
                if line.strip():
                    fields = line.split()
                    if len(fields) >= 14:
                        repeats.append({
                            'chr': fields[4],
                            'start': int(fields[5]) - 1,  # Convert to 0-based
                            'end': int(fields[6]),
                            'name': fields[9],
                            'family': fields[10],
                            'class': fields[10].split('/')[0] if '/' in fields[10] else fields[10],
                            'strand': '+' if fields[8] == '+' else '-',
                            'score': float(fields[0]) if fields[0].replace('.', '').isdigit() else 0
                        })
        
        return pd.DataFrame(repeats)
    
    def parse_gtf_repeats(self, gtf_file):
        """Parse GTF format repeat annotations"""
        repeats = []
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] in ['repeat_region', 'transposable_element']:
                    # Parse attributes
                    attributes = {}
                    for attr in fields[8].split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attributes[key.strip()] = value.strip().strip('"')
                    
                    repeats.append({
                        'chr': fields[0],
                        'start': int(fields[3]) - 1,  # Convert to 0-based
                        'end': int(fields[4]),
                        'name': attributes.get('Name', attributes.get('ID', 'unknown')),
                        'family': attributes.get('repeat_family', attributes.get('family', 'unknown')),
                        'class': attributes.get('repeat_class', attributes.get('class', 'unknown')),
                        'strand': fields[6],
                        'score': float(fields[5]) if fields[5] != '.' else 0
                    })
        
        return pd.DataFrame(repeats)
    
    def classify_repeat_family(self, row):
        """Classify repeat into major families"""
        family = row['family'].upper()
        class_name = row['class'].upper() if pd.notna(row['class']) else ''
        
        for family_class, keywords in self.repeat_families.items():
            for keyword in keywords:
                if keyword.upper() in family or keyword.upper() in class_name:
                    return family_class
        
        return 'Other'
    
    def create_repeat_regions(self, repeats_df):
        """Create extended regions around repeat elements"""
        extended_repeats = repeats_df.copy()
        
        # Extend regions
        extended_repeats['extended_start'] = np.maximum(0, 
            extended_repeats['start'] - self.extend_region)
        extended_repeats['extended_end'] = (extended_repeats['end'] + 
                                          self.extend_region)
        
        # Create BED format for extended regions
        bed_regions = extended_repeats[['chr', 'extended_start', 'extended_end', 
                                       'name', 'score', 'strand', 'family_class']].copy()
        bed_regions.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'family_class']
        
        return bed_regions
    
    def calculate_signal_matrix(self, species, repeats_df, bigwig_file):
        """Calculate m6A signal matrix for repeat regions"""
        # Create temporary BED file
        temp_bed = f"temp_{species}_repeats.bed"
        bed_regions = self.create_repeat_regions(repeats_df)
        
        # Save to BED file
        bed_regions[['chr', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            temp_bed, sep='\t', header=False, index=False
        )
        
        # Calculate signal using computeMatrix (deepTools)
        matrix_file = f"temp_{species}_matrix.gz"
        
        # Number of bins for the extended region
        n_bins = (2 * self.extend_region + bed_regions['end'] - bed_regions['start']).iloc[0] // self.bin_size
        
        cmd = [
            'computeMatrix', 'scale-regions',
            '-S', bigwig_file,
            '-R', temp_bed,
            '-o', matrix_file,
            '--beforeRegionStartLength', str(self.extend_region),
            '--afterRegionStartLength', str(self.extend_region),
            '--binSize', str(self.bin_size),
            '--skipZeros'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Read matrix
            signal_matrix = self.read_deeptools_matrix(matrix_file)
            
            # Clean up temp files
            os.remove(temp_bed)
            os.remove(matrix_file)
            
            return signal_matrix
            
        except subprocess.CalledProcessError as e:
            print(f"Error calculating matrix for {species}: {e}")
            return None
        except FileNotFoundError:
            print("deepTools not found. Using alternative method...")
            return self.calculate_signal_matrix_alternative(species, repeats_df, bigwig_file)
    
    def calculate_signal_matrix_alternative(self, species, repeats_df, bigwig_file):
        """Alternative method using pyBigWig if deepTools is not available"""
        try:
            import pyBigWig
            
            # Open bigwig file
            bw = pyBigWig.open(bigwig_file)

            bed_regions = self.create_repeat_regions(repeats_df)
            region_starts = bed_regions['start'].values
            region_ends = bed_regions['end'].values
            region_chrs = bed_regions['chr'].values

            region_length = region_ends[0] - region_starts[0]
            n_bins = region_length // self.bin_size

            # Precompute bin offsets
            bin_offsets = np.arange(n_bins) * self.bin_size

            signal_matrix = np.zeros((len(bed_regions), n_bins), dtype=float)

            for idx, (chr_, start, end) in enumerate(zip(region_chrs, region_starts, region_ends)):
                bin_starts = start + bin_offsets
                bin_ends = np.minimum(end, bin_starts + self.bin_size)

                # pyBigWig.stats can take lists for batch queries
                try:
                    signals = bw.stats(chr_, list(bin_starts), list(bin_ends), type="mean")
                    # Replace None with 0
                    signal_matrix[idx, :] = [s if s is not None else 0 for s in signals]
                except Exception:
                    signal_matrix[idx, :] = 0

            bw.close()
            return signal_matrix
            
        except ImportError:
            print("pyBigWig not available. Please install deepTools or pyBigWig.")
            return None
    
    def read_deeptools_matrix(self, matrix_file):
        """Read deepTools computeMatrix output"""
        import gzip
        
        with gzip.open(matrix_file, 'rt') as f:
            # Skip header line
            header = f.readline()
            
            matrix_data = []
            for line in f:
                if line.startswith('@'):
                    continue
                values = line.strip().split('\t')[6:]  # Skip first 6 columns (metadata)
                matrix_data.append([float(x) if x != 'nan' else 0 for x in values])
        
        return np.array(matrix_data)
    
    def analyze_all_species(self):
        """Analyze m6A signal on repeats for all species"""
        results = {}
        
        for species in self.species_names:
            print(f"Analyzing {species}...")
            
            # Load repeat annotations
            repeats = self.load_repeat_annotations(species)
            if repeats.empty:
                print(f"No repeats loaded for {species}")
                continue
            
            # Get bigwig file
            bigwig_file = self.species_data[species].get('bigwig')
            if not bigwig_file or not os.path.exists(bigwig_file):
                print(f"BigWig file not found for {species}")
                continue
            
            # Calculate signal matrix
            signal_matrix = self.calculate_signal_matrix_alternative(species, repeats, bigwig_file)
            if signal_matrix is None:
                continue
            
            # Group by repeat family
            family_signals = {}
            family_counts = repeats['family_class'].value_counts()
            
            for family in self.repeat_families.keys():
                family_mask = repeats['family_class'] == family
                if family_mask.sum() > 0:
                    family_signals[family] = signal_matrix[family_mask]
            
            results[species] = {
                'repeats_df': repeats,
                'signal_matrix': signal_matrix,
                'family_signals': family_signals,
                'family_counts': family_counts
            }
        
        return results
    
    def plot_signal_heatmaps(self, results, output_dir='m6a_repeat_plots'):
        """Create signal heatmaps for each species and repeat family"""
        os.makedirs(output_dir, exist_ok=True)
        
        for species, data in results.items():
            fig, axes = plt.subplots(2, 4, figsize=(20, 10))
            axes = axes.flatten()
            
            family_signals = data['family_signals']
            family_counts = data['family_counts']
            
            for i, (family, signals) in enumerate(family_signals.items()):
                if i >= len(axes):
                    break
                
                if len(signals) == 0:
                    axes[i].text(0.5, 0.5, f'No {family} repeats', 
                               ha='center', va='center', transform=axes[i].transAxes)
                    axes[i].set_title(f'{family} (n=0)')
                    continue
                
                # Calculate mean signal for sorting
                mean_signals = np.mean(signals, axis=1)
                sorted_indices = np.argsort(mean_signals)[::-1]  # Descending order
                
                # Plot heatmap (show top 1000 repeats max for visibility)
                n_show = min(1000, len(signals))
                plot_signals = signals[sorted_indices[:n_show]]
                
                im = axes[i].imshow(plot_signals, aspect='auto', cmap='Reds', 
                                  interpolation='nearest')
                
                # Set labels
                n_bins = plot_signals.shape[1]
                x_ticks = [0, n_bins//4, n_bins//2, 3*n_bins//4, n_bins-1]
                x_labels = [f'-{self.extend_region//1000}kb', f'-{self.extend_region//2000}kb', 
                          'Center', f'+{self.extend_region//2000}kb', f'+{self.extend_region//1000}kb']
                
                axes[i].set_xticks(x_ticks)
                axes[i].set_xticklabels(x_labels)
                axes[i].set_ylabel('Repeat Elements (sorted by signal)')
                axes[i].set_title(f'{family} (n={family_counts.get(family, 0)})')
                
                # Add colorbar
                plt.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
            
            # Remove empty subplots
            for i in range(len(family_signals), len(axes)):
                axes[i].remove()
            
            plt.suptitle(f'M6A Signal on Repeat Elements - {species.title()}', 
                        fontsize=16, fontweight='bold')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{species}_repeat_heatmaps.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def plot_average_profiles(self, results, output_dir='m6a_repeat_plots'):
        """Plot average m6A signal profiles for each repeat family"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Prepare data for plotting
        families = list(self.repeat_families.keys())
        
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        axes = axes.flatten()
        
        for i, family in enumerate(families):
            if i >= len(axes):
                break
            
            for species in self.species_names:
                if species not in results:
                    continue
                
                family_signals = results[species]['family_signals'].get(family)
                if family_signals is None or len(family_signals) == 0:
                    continue
                
                # Calculate average profile
                avg_profile = np.mean(family_signals, axis=0)
                sem_profile = stats.sem(family_signals, axis=0)
                
                # X-axis (genomic position)
                n_bins = len(avg_profile)
                x_positions = np.linspace(-self.extend_region, self.extend_region, n_bins)
                
                # Plot profile
                axes[i].plot(x_positions, avg_profile, label=species.title(), linewidth=2)
                axes[i].fill_between(x_positions, 
                                   avg_profile - sem_profile, 
                                   avg_profile + sem_profile, 
                                   alpha=0.2)
            
            axes[i].set_xlabel('Distance from Repeat Center (bp)')
            axes[i].set_ylabel('Average M6A Signal')
            axes[i].set_title(f'{family} Repeats')
            axes[i].legend()
            axes[i].grid(True, alpha=0.3)
            
            # Add vertical line at center
            axes[i].axvline(x=0, color='black', linestyle='--', alpha=0.5)
        
        # Remove empty subplots
        for i in range(len(families), len(axes)):
            axes[i].remove()
        
        plt.suptitle('Average M6A Signal Profiles by Repeat Family', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/average_repeat_profiles.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_comparative_heatmap(self, results, output_file='m6a_repeat_plots/comparative_heatmap.png'):
        """Create comparative heatmap showing signal intensity across species and families"""
        # Calculate average signal for each species-family combination
        comparison_data = []
        
        for species in self.species_names:
            if species not in results:
                continue
            
            species_data = []
            for family in self.repeat_families.keys():
                family_signals = results[species]['family_signals'].get(family)
                if family_signals is not None and len(family_signals) > 0:
                    avg_signal = np.mean(family_signals)
                    species_data.append(avg_signal)
                else:
                    species_data.append(0)
            
            comparison_data.append(species_data)
        
        # Create DataFrame for heatmap
        comparison_df = pd.DataFrame(comparison_data, 
                                   index=[s.title() for s in self.species_names if s in results],
                                   columns=list(self.repeat_families.keys()))
        
        # Create heatmap
        plt.figure(figsize=(12, 8))
        mask = comparison_df == 0
        sns.heatmap(comparison_df, annot=True, fmt='.3f', cmap='Reds', 
                   mask=mask, cbar_kws={'label': 'Average M6A Signal'})
        
        plt.title('M6A Signal Intensity Across Species and Repeat Families', 
                 fontsize=14, fontweight='bold')
        plt.xlabel('Repeat Family')
        plt.ylabel('Species')
        plt.tight_layout()
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return comparison_df
    
    def statistical_analysis(self, results):
        """Perform statistical analysis of repeat family differences"""
        stats_results = {}
        
        for species in self.species_names:
            if species not in results:
                continue
            
            family_signals = results[species]['family_signals']
            species_stats = {}
            
            # Calculate statistics for each family
            for family, signals in family_signals.items():
                if len(signals) > 0:
                    avg_signals = np.mean(signals, axis=1)  # Average per repeat
                    species_stats[family] = {
                        'mean': np.mean(avg_signals),
                        'std': np.std(avg_signals),
                        'median': np.median(avg_signals),
                        'n_elements': len(signals),
                        'total_signal': np.sum(avg_signals)
                    }
            
            # Pairwise comparisons between families
            family_names = list(species_stats.keys())
            pairwise_comparisons = {}
            
            for i in range(len(family_names)):
                for j in range(i+1, len(family_names)):
                    fam1, fam2 = family_names[i], family_names[j]
                    
                    signals1 = np.mean(family_signals[fam1], axis=1)
                    signals2 = np.mean(family_signals[fam2], axis=1)
                    
                    # Mann-Whitney U test
                    statistic, p_value = stats.mannwhitneyu(signals1, signals2, 
                                                          alternative='two-sided')
                    
                    pairwise_comparisons[f'{fam1}_vs_{fam2}'] = {
                        'statistic': statistic,
                        'p_value': p_value,
                        'significant': p_value < 0.05,
                        'effect_size': (np.mean(signals1) - np.mean(signals2)) / 
                                     np.sqrt((np.std(signals1)**2 + np.std(signals2)**2) / 2)
                    }
            
            stats_results[species] = {
                'family_statistics': species_stats,
                'pairwise_comparisons': pairwise_comparisons
            }
        
        return stats_results
    
    def generate_comprehensive_report(self, results, stats_results, output_file='m6a_repeat_analysis_report.md'):
        """Generate comprehensive analysis report"""
        with open(output_file, 'w') as f:
            f.write("# M6A Signal Analysis on Repeat Elements\n\n")
            
            f.write("## Summary\n")
            f.write(f"- Species analyzed: {', '.join(self.species_names)}\n")
            f.write(f"- Repeat families: {', '.join(self.repeat_families.keys())}\n")
            f.write(f"- Analysis parameters: {self.bin_size}bp bins, ±{self.extend_region}bp extension\n\n")
            
            # Species-specific results
            for species in self.species_names:
                if species not in results:
                    continue
                
                f.write(f"## {species.title()} Results\n\n")
                
                # Repeat counts
                family_counts = results[species]['family_counts']
                f.write("### Repeat Element Counts\n")
                for family, count in family_counts.items():
                    f.write(f"- {family}: {count:,} elements\n")
                f.write("\n")
                
                # Statistical summary
                if species in stats_results:
                    species_stats = stats_results[species]['family_statistics']
                    f.write("### M6A Signal Statistics by Family\n")
                    f.write("| Family | Mean Signal | Std | Median | Elements |\n")
                    f.write("|--------|-------------|-----|--------|----------|\n")
                    
                    for family, stats_data in species_stats.items():
                        f.write(f"| {family} | {stats_data['mean']:.4f} | "
                               f"{stats_data['std']:.4f} | {stats_data['median']:.4f} | "
                               f"{stats_data['n_elements']:,} |\n")
                    f.write("\n")
                    
                    # Significant comparisons
                    f.write("### Significant Pairwise Comparisons (p < 0.05)\n")
                    pairwise = stats_results[species]['pairwise_comparisons']
                    significant_pairs = {k: v for k, v in pairwise.items() if v['significant']}
                    
                    if significant_pairs:
                        for comparison, data in significant_pairs.items():
                            f.write(f"- **{comparison.replace('_vs_', ' vs ')}**: "
                                   f"p = {data['p_value']:.2e}, effect size = {data['effect_size']:.3f}\n")
                    else:
                        f.write("- No significant differences found between repeat families\n")
                    f.write("\n")


# Example usage and shell script integration
# def create_analysis_script():
    """Create a shell script for running the analysis"""
    script_content = """#!/bin/bash
# M6A Repeat Analysis Pipeline

# Set up directories
mkdir -p data/repeats
mkdir -p data/bigwig
mkdir -p results

# Example data structure:
# data/
# ├── repeats/
# │   ├── human_repeats.out
# │   ├── monkey_repeats.out
# │   └── mouse_repeats.out
# └── bigwig/
#     ├── human_m6a.bw
#     ├── monkey_m6a.bw
#     └── mouse_m6a.bw

echo "Starting M6A repeat analysis..."

# Run Python analysis
python3 << 'EOF'
from m6a_repeat_analysis import M6ARepeatAnalyzer
import json

# Define species data
species_data = {
    'human': {
        'peaks': 'data/peaks/human_peaks.bed',
        'bigwig': 'data/bigwig/human_m6a.bw'
    },
    'monkey': {
        'peaks': 'data/peaks/monkey_peaks.bed', 
        'bigwig': 'data/bigwig/monkey_m6a.bw'
    },
    'mouse': {
        'peaks': 'data/peaks/mouse_peaks.bed',
        'bigwig': 'data/bigwig/mouse_m6a.bw'
    }
}

# Define repeat annotations
repeat_annotations = {
    'human': 'data/repeats/human_repeats.out',
    'monkey': 'data/repeats/monkey_repeats.out', 
    'mouse': 'data/repeats/mouse_repeats.out'
}

# Create analyzer
analyzer = M6ARepeatAnalyzer(species_data, repeat_annotations, 
                            bin_size=50, extend_region=2000)

# Run analysis
print("Analyzing m6A signal on repeats...")
results = analyzer.analyze_all_species()

# Generate plots
print("Creating visualizations...")
analyzer.plot_signal_heatmaps(results, 'results/heatmaps')
analyzer.plot_average_profiles(results, 'results/profiles') 
comparison_df = analyzer.plot_comparative_heatmap(results, 'results/comparative_analysis.png')

# Statistical analysis
print("Performing statistical analysis...")
stats_results = analyzer.statistical_analysis(results)

# Generate report
analyzer.generate_comprehensive_report(results, stats_results, 'results/m6a_repeat_report.md')

# Save results
with open('results/analysis_results.json', 'w') as f:
    # Convert numpy arrays to lists for JSON serialization
    json_results = {}
    for species, data in results.items():
        json_results[species] = {
            'family_counts': data['family_counts'].to_dict(),
            'repeats_summary': {
                'total_repeats': len(data['repeats_df']),
                'families_present': list(data['family_signals'].keys())
            }
        }
    json.dump(json_results, f, indent=2)

print("Analysis completed! Results saved in 'results/' directory")
EOF

echo "M6A repeat analysis completed!"
"""
    
    with open('run_m6a_repeat_analysis.sh', 'w') as f:
        f.write(script_content)
    
    # Make executable
    import stat
    st = os.stat('run_m6a_repeat_analysis.sh')
    os.chmod('run_m6a_repeat_analysis.sh', st.st_mode | stat.S_IEXEC)
    
    print("Analysis script created: run_m6a_repeat_analysis.sh")


if __name__ == "__main__":
    # Create the analysis script
    # create_analysis_script()
    
    print("M6A Repeat Analysis Tool initialized!")
    print("Usage:")
    print("1. Prepare your data in the expected directory structure")
    print("2. Run: ./run_m6a_repeat_analysis.sh")
    print("3. Or use the M6ARepeatAnalyzer class directly in Python")