#!/usr/bin/env python3
# window_based_analysis.py

import pandas as pd
import pybedtools
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

class WindowBasedAnalyzer:
    def __init__(self, peak_files, window_sizes=[100, 200, 500, 1000]):
        self.peak_files = peak_files
        self.window_sizes = window_sizes
        self.species_names = list(peak_files.keys())
        
    def create_peak_summits(self, bed_file, output_file):
        """Create summit coordinates from peak files"""
        peaks = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0, 1, 2, 3, 4],
                           names=['chr', 'start', 'end', 'name', 'score'])
        
        # Calculate summit as midpoint
        peaks['summit'] = ((peaks['start'] + peaks['end']) // 2).astype(int)
        
        # Create summit BED format
        summits = peaks[['chr', 'summit', 'name', 'score']].copy()
        summits['summit_end'] = summits['summit'] + 1
        summits = summits[['chr', 'summit', 'summit_end', 'name', 'score']]
        
        summits.to_csv(output_file, sep='\t', header=False, index=False)
        return output_file
    
    def create_windows(self, summit_file, window_size, output_file):
        """Create windows around summits"""
        summits = pybedtools.BedTool(summit_file)
        windows = summits.slop(b=window_size//2, g='"hg38.chrom."sizes')
        windows.saveas(output_file)
        return output_file
    
    def window_overlap_analysis(self):
        """Analyze overlap across different window sizes"""
        results = {}
        
        # Create summit files
        summit_files = {}
        for species, peak_file in self.peak_files.items():
            summit_file = f"summits/{species}_summits.bed"
            summit_files[species] = self.create_peak_summits(peak_file, summit_file)
        
        # Analyze for each window size
        for window_size in self.window_sizes:
            window_files = {}
            
            # Create windows for each species
            for species, summit_file in summit_files.items():
                window_file = f"windows/{species}_windows_{window_size}bp.bed"
                window_files[species] = self.create_windows(
                    summit_file, window_size, window_file
                )
            
            # Calculate overlaps
            window_analyzer = PeakOverlapAnalyzer(window_files)
            overlap_results = window_analyzer.pairwise_overlap_analysis()
            
            results[window_size] = overlap_results
        
        return results
    
    def plot_window_analysis(self, results, output_file):
        """Plot window analysis results"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Prepare data for plotting
        window_sizes = list(results.keys())
        species_pairs = list(results[window_sizes[0]].keys())
        
        # Jaccard index heatmap
        jaccard_data = []
        for window_size in window_sizes:
            row = []
            for pair in species_pairs:
                row.append(results[window_size][pair]['jaccard_index'])
            jaccard_data.append(row)
        
        jaccard_df = pd.DataFrame(jaccard_data, 
                                 index=window_sizes, 
                                 columns=species_pairs)
        
        sns.heatmap(jaccard_df, annot=True, fmt='.3f', 
                   cmap='viridis', ax=axes[0,0])
        axes[0,0].set_title('Jaccard Index by Window Size')
        axes[0,0].set_xlabel('Species Pairs')
        axes[0,0].set_ylabel('Window Size (bp)')
        
        # Line plot for each species pair
        for pair in species_pairs:
            jaccard_values = [results[ws][pair]['jaccard_index'] 
                            for ws in window_sizes]
            axes[0,1].plot(window_sizes, jaccard_values, 
                          marker='o', label=pair)
        
        axes[0,1].set_xlabel('Window Size (bp)')
        axes[0,1].set_ylabel('Jaccard Index')
        axes[0,1].set_title('Conservation vs Window Size')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

# Usage

peak_files = {
    "H9_G": "data/lifted_peaks/H9_G__m6A_hg38ToHg38.lifted.bed",
    "H9_I": "data/lifted_peaks/H9_I__m6A_hg38ToHg38.lifted.bed",
    "LG_I": "data/lifted_peaks/LG_I__m6A_nomLeu3ToHg38.lifted.bed",
    "LG_I": "data/lifted_peaks/LG_I__m7G_nomLeu3ToHg38.lifted.bed",
    "OU_G": "data/lifted_peaks/OU_G__m6A_nomLeu3ToHg38.lifted.bed",
    "OU_G": "data/lifted_peaks/OU_G__m7G_nomLeu3ToHg38.lifted.bed",
    "OU_I": "data/lifted_peaks/OU_I__m6A_nomLeu3ToHg38.lifted.bed",
    "OU_I": "data/lifted_peaks/OU_I__m7G_nomLeu3ToHg38.lifted.bed",
    "PM_G": "data/lifted_peaks/PM_G__m6A_nasLar1ToHg38.lifted.bed",
    "PM_I": "data/lifted_peaks/PM_I__m6A_nasLar1ToHg38.lifted.bed"

}
window_analyzer = WindowBasedAnalyzer(peak_files)
window_results = window_analyzer.window_overlap_analysis()
window_analyzer.plot_window_analysis(window_results, 'window_analysis.png')