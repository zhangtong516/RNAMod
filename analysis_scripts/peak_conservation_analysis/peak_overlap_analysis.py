#!/usr/bin/env python3
# peak_overlap_analysis.py

import pandas as pd
import numpy as np
from scipy import stats
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import json

class PeakOverlapAnalyzer:
    def __init__(self, peak_files):
        self.peak_files = peak_files
        self.species_names = list(peak_files.keys())
        self.bed_objects = {}
        self.load_peaks()
    
    def load_peaks(self):
        """Load peak files as pybedtools BedTool objects"""
        for species, file_path in self.peak_files.items():
            self.bed_objects[species] = pybedtools.BedTool(file_path)
    
    def calculate_jaccard_index(self, species1, species2, min_overlap=0.5):
        """Calculate Jaccard index between two species"""
        bed1 = self.bed_objects[species1]
        bed2 = self.bed_objects[species2]
        
        # Find intersections with minimum overlap
        intersection = bed1.intersect(bed2, f=min_overlap, r=True)
        
        # Calculate sizes
        intersection_size = len(intersection)
        union_size = len(bed1) + len(bed2) - intersection_size
        
        jaccard = intersection_size / union_size if union_size > 0 else 0
        return jaccard, intersection_size, union_size
    
    def pairwise_overlap_analysis(self, min_overlap=0.5):
        """Perform pairwise overlap analysis for all species pairs"""
        results = {}
        
        for species1, species2 in combinations(self.species_names, 2):
            jaccard, intersection, union = self.calculate_jaccard_index(
                species1, species2, min_overlap
            )
            
            results[f"{species1}_vs_{species2}"] = {
                'jaccard_index': jaccard,
                'intersection_size': intersection,
                'union_size': union,
                'species1': species1,
                'species2': species2
            }
        
        return results
    
    def hypergeometric_test(self, species1, species2, genome_size=3000000000):
        """Perform hypergeometric test for overlap significance"""
        bed1 = self.bed_objects[species1]
        bed2 = self.bed_objects[species2]
        
        # Calculate overlap
        overlap = bed1.intersect(bed2, u=True)
        
        # Parameters for hypergeometric test
        k = len(overlap)  # observed overlap
        n = len(bed1)     # peaks in species1
        K = len(bed2)     # peaks in species2
        N = genome_size   # total possible positions (approximation)
        
        # Hypergeometric test
        p_value = stats.hypergeom.sf(k-1, N, K, n)
        
        return p_value, k, n, K
    
    def significance_testing(self, genome_size=3000000000):
        """Test significance of overlaps for all species pairs"""
        results = {}
        
        for species1, species2 in combinations(self.species_names, 2):
            p_value, overlap, n1, n2 = self.hypergeometric_test(
                species1, species2, genome_size
            )
            
            results[f"{species1}_vs_{species2}"] = {
                'p_value': p_value,
                'overlap_count': overlap,
                'peaks_species1': n1,
                'peaks_species2': n2,
                'significant': p_value < 0.05
            }
        
        return results
    
    def multiple_testing_correction(self, p_values, method='fdr_bh'):
        """Apply multiple testing correction"""
        try:
            from statsmodels.stats.multitest import multipletests
            
            p_vals = list(p_values.values())
            species_pairs = list(p_values.keys())
            
            if method == 'fdr_bh':
                reject, p_adj, _, _ = multipletests(p_vals, method='fdr_bh')
            elif method == 'bonferroni':
                reject, p_adj, _, _ = multipletests(p_vals, method='bonferroni')
        
        except ImportError:
            # Fallback implementation without statsmodels
            print("Warning: statsmodels not available. Using simple Bonferroni correction.")
            p_vals = list(p_values.values())
            species_pairs = list(p_values.keys())
            
            if method == 'bonferroni' or method == 'fdr_bh':
                # Simple Bonferroni correction
                n_tests = len(p_vals)
                p_adj = [min(p * n_tests, 1.0) for p in p_vals]
                reject = [p < 0.05 for p in p_adj]
            
        corrected_results = {}
        for i, pair in enumerate(species_pairs):
            corrected_results[pair] = {
                'p_value': p_vals[i],
                'p_adjusted': p_adj[i],
                'significant': reject[i]
            }
        
        return corrected_results

    
# Usage example
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

analyzer = PeakOverlapAnalyzer(peak_files)
overlap_results = analyzer.pairwise_overlap_analysis()
significance_results = analyzer.significance_testing()

# Extract p-values for correction
p_values = {k: v['p_value'] for k, v in significance_results.items()}
corrected_results = analyzer.multiple_testing_correction(p_values)

# Save results
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

with open('overlap_analysis_results.json', 'w') as f:
    json.dump({
        'overlap_results': overlap_results,
        'significance_results': significance_results,
        'corrected_results': corrected_results
    }, f,  cls=NumpyEncoder, indent=2)

