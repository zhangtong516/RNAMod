#!/usr/bin/env python3
# phylogenetic_analysis.py

import pandas as pd
import numpy as np
import dendropy
import io
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns

class PhylogeneticAnalyzer:
    def __init__(self, peak_files):
        self.peak_files = peak_files
        self.species_names = list(peak_files.keys())
        self.distance_matrix = None
        
    def create_binary_matrix(self, window_size=200):
        """Create binary presence/absence matrix for all peaks"""
        # Get all unique genomic regions
        all_regions = set()
        species_peaks = {}
        
        for species, file_path in self.peak_files.items():
            peaks = pd.read_csv(file_path, sep='\t', header=None,
                               names=['chr', 'start', 'end', 'name', 'score'])
            
            # Create windows around peaks
            regions = []
            for _, peak in peaks.iterrows():
                center = (peak['start'] + peak['end']) // 2
                region = (peak['chr'], 
                         center - window_size//2, 
                         center + window_size//2)
                regions.append(region)
                all_regions.add(region)
            
            species_peaks[species] = set(regions)
        
        # Create binary matrix
        all_regions = sorted(list(all_regions))
        binary_matrix = np.zeros((len(self.species_names), len(all_regions)))
        
        for i, species in enumerate(self.species_names):
            for j, region in enumerate(all_regions):
                if region in species_peaks[species]:
                    binary_matrix[i, j] = 1
        
        return binary_matrix, all_regions
    
    def calculate_distance_matrix(self, binary_matrix, method='jaccard'):
        """Calculate distance matrix between species"""
        if method == 'jaccard':
            # Jaccard distance = 1 - Jaccard similarity
            distances = pdist(binary_matrix, metric='jaccard')
        elif method == 'hamming':
            distances = pdist(binary_matrix, metric='hamming')
        elif method == 'euclidean':
            distances = pdist(binary_matrix, metric='euclidean')
        
        self.distance_matrix = squareform(distances)
        return self.distance_matrix
    
    def construct_phylogenetic_tree(self, method='nj'):
        
        # Create taxon namespace
        taxon_namespace = dendropy.TaxonNamespace()
        
        # Add taxa to the namespace
        taxa = []
        for i in range(len(self.species_names)):
            taxon = dendropy.Taxon(label=self.species_names[i])
            taxon_namespace.add_taxon(taxon)
            taxa.append(taxon)
        
        # Create CSV format from distance matrix
        csv_content = io.StringIO()
        
        # Write header row with taxon names
        csv_content.write("," + ",".join(self.species_names) + "\n")
        
        # Write distance matrix rows
        for i, taxon_name in enumerate(self.species_names):
            row = [taxon_name] + [str(self.distance_matrix[i][j]) for j in range(len(self.species_names))]
            csv_content.write(",".join(row) + "\n")
        
        # Reset StringIO position to beginning
        csv_content.seek(0)
        
        # Create PhylogeneticDistanceMatrix from CSV
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=csv_content,
            taxon_namespace=taxon_namespace,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            delimiter=","
        )
        
        # Construct tree based on method
        if method == 'nj':
            tree = pdm.nj_tree()
        elif method == 'upgma':
            tree = pdm.upgma_tree()
        else:
            raise ValueError("Method must be 'nj' or 'upgma'")
        
        return tree
    def ancestral_state_reconstruction(self, binary_matrix, tree):
        """Perform ancestral state reconstruction"""
        # Simplified parsimony-based reconstruction
        # This is a basic implementation - more sophisticated methods available
        
        results = {}
        n_regions = binary_matrix.shape[1]
        
        # For each genomic region, reconstruct ancestral states
        for region_idx in range(n_regions):
            region_states = binary_matrix[:, region_idx]
            
            # Simple parsimony: ancestral state = most common state
            ancestral_state = 1 if np.sum(region_states) > len(region_states)/2 else 0
            
            results[region_idx] = {
                'ancestral_state': ancestral_state,
                'current_states': {species: int(region_states[i]) 
                                 for i, species in enumerate(self.species_names)},
                'conservation_score': np.mean(region_states == ancestral_state)
            }
        
        return results
    
    def identify_lineage_specific_changes(self, ancestral_results):
        """Identify lineage-specific gains and losses"""
        lineage_changes = {species: {'gains': [], 'losses': []} 
                          for species in self.species_names}
        
        for region_idx, result in ancestral_results.items():
            ancestral_state = result['ancestral_state']
            current_states = result['current_states']
            
            for species, current_state in current_states.items():
                if ancestral_state == 0 and current_state == 1:
                    lineage_changes[species]['gains'].append(region_idx)
                elif ancestral_state == 1 and current_state == 0:
                    lineage_changes[species]['losses'].append(region_idx)
        
        return lineage_changes
    
    def plot_phylogenetic_analysis(self, output_prefix):
        """Plot phylogenetic analysis results"""
        # Distance matrix heatmap
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Heatmap
        mask = np.triu(np.ones_like(self.distance_matrix, dtype=bool))
        sns.heatmap(self.distance_matrix, 
                   annot=True, 
                   fmt='.3f',
                   mask=mask,
                   xticklabels=self.species_names,
                   yticklabels=self.species_names,
                   cmap='viridis',
                   ax=axes[0])
        axes[0].set_title('Distance Matrix')
        
        # Dendrogram
        linkage_matrix = linkage(pdist(self.distance_matrix), method='average')
        dendrogram(linkage_matrix, 
                  labels=self.species_names, 
                  ax=axes[1])
        axes[1].set_title('Hierarchical Clustering')
        axes[1].set_xlabel('Species')
        axes[1].set_ylabel('Distance')
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_phylogenetic.png', dpi=300, bbox_inches='tight')
        plt.close()

# Usage

peak_files = {
    "H9_G": "data/lifted_peaks/H9_G__m6A_hg38ToHg38.lifted.bed",
    "H9_I": "data/lifted_peaks/H9_I__m6A_hg38ToHg38.lifted.bed",
    "LG_I": "data/lifted_peaks/LG_I__m6A_nomLeu3ToHg38.lifted.bed",
    "OU_G": "data/lifted_peaks/OU_G__m6A_nomLeu3ToHg38.lifted.bed",
    "OU_I": "data/lifted_peaks/OU_I__m6A_nomLeu3ToHg38.lifted.bed",
    "PM_G": "data/lifted_peaks/PM_G__m6A_nasLar1ToHg38.lifted.bed",
    "PM_I": "data/lifted_peaks/PM_I__m6A_nasLar1ToHg38.lifted.bed"

}

phylo_analyzer = PhylogeneticAnalyzer(peak_files)
binary_matrix, regions = phylo_analyzer.create_binary_matrix()
distance_matrix = phylo_analyzer.calculate_distance_matrix(binary_matrix)
tree = phylo_analyzer.construct_phylogenetic_tree()
ancestral_results = phylo_analyzer.ancestral_state_reconstruction(binary_matrix, tree)
lineage_changes = phylo_analyzer.identify_lineage_specific_changes(ancestral_results)

phylo_analyzer.plot_phylogenetic_analysis('phylogenetic_analysis')

## m7G

peak_files = {
    "LG_I": "data/lifted_peaks/LG_I__m7G_nomLeu3ToHg38.lifted.bed",
    "OU_G": "data/lifted_peaks/OU_G__m7G_nomLeu3ToHg38.lifted.bed",
    "OU_I": "data/lifted_peaks/OU_I__m7G_nomLeu3ToHg38.lifted.bed"
}
