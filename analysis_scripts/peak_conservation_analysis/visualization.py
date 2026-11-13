#!/usr/bin/env python3
# visualization_dashboard.py

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

class ConservationDashboard:
    def __init__(self, overlap_results, significance_results, 
                 window_results, phylo_results):
        self.overlap_results = overlap_results
        self.significance_results = significance_results
        self.window_results = window_results
        self.phylo_results = phylo_results
        
    def create_comprehensive_dashboard(self, output_file):
        """Create comprehensive conservation analysis dashboard"""
        fig = plt.figure(figsize=(20, 16))
        
        # Create grid layout
        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
        
        # 1. Jaccard Index Heatmap
        ax1 = fig.add_subplot(gs[0, :2])
        self.plot_jaccard_heatmap(ax1)
        
        # 2. Significance Test Results
        ax2 = fig.add_subplot(gs[0, 2:])
        self.plot_significance_results(ax2)
        
        # 3. Window Size Analysis
        ax3 = fig.add_subplot(gs[1, :2])
        self.plot_window_analysis(ax3)
        
        # 4. Conservation Distribution
        ax4 = fig.add_subplot(gs[1, 2:])
        self.plot_conservation_distribution(ax4)
        
        # 5. Phylogenetic Tree
        ax5 = fig.add_subplot(gs[2, :2])
        self.plot_phylogenetic_tree(ax5)
        
        # 6. Lineage-specific Changes
        ax6 = fig.add_subplot(gs[2, 2:])
        self.plot_lineage_changes(ax6)
        
        # 7. Summary Statistics
        ax7 = fig.add_subplot(gs[3, :])
        self.plot_summary_statistics(ax7)
        
        plt.suptitle('M6A Peak Conservation Analysis Dashboard', 
                    fontsize=16, fontweight='bold')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_jaccard_heatmap(self, ax):
        """Plot Jaccard index heatmap"""
        species_pairs = list(self.overlap_results.keys())
        jaccard_values = [self.overlap_results[pair]['jaccard_index'] 
                         for pair in species_pairs]
        
        # Create matrix for heatmap
        n_species = 3  # human, monkey, mouse
        matrix = np.zeros((n_species, n_species))
        
        # Fill matrix (this is simplified - adjust based on your species order)
        pair_indices = {
            'human_vs_monkey': (0, 1),
            'human_vs_mouse': (0, 2),
            'monkey_vs_mouse': (1, 2)
        }
        
        for pair, jaccard in zip(species_pairs, jaccard_values):
            if pair in pair_indices:
                i, j = pair_indices[pair]
                matrix[i, j] = jaccard
                matrix[j, i] = jaccard
        
        sns.heatmap(matrix, annot=True, fmt='.3f', cmap='viridis',
                   xticklabels=['Human', 'Monkey', 'Mouse'],
                   yticklabels=['Human', 'Monkey', 'Mouse'],
                   ax=ax)
        ax.set_title('Jaccard Index Between Species')
    
    def plot_significance_results(self, ax):
        """Plot significance test results"""
        pairs = list(self.significance_results.keys())
        p_values = [self.significance_results[pair]['p_value'] for pair in pairs]
        
        # Create bar plot
        colors = ['red' if p < 0.05 else 'blue' for p in p_values]
        bars = ax.bar(range(len(pairs)), -np.log10(p_values), color=colors)
        
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', 
                  label='p=0.05 threshold')
        ax.set_xticks(range(len(pairs)))
        ax.set_xticklabels([p.replace('_vs_', ' vs ') for p in pairs], 
                          rotation=45)
        ax.set_ylabel('-log10(p-value)')
        ax.set_title('Significance of Peak Overlaps')
        ax.legend()
    
    def plot_window_analysis(self, ax):
        """Plot window size analysis"""
        window_sizes = list(self.window_results.keys())
        species_pairs = list(self.window_results[window_sizes[0]].keys())
        
        for pair in species_pairs:
            jaccard_values = [self.window_results[ws][pair]['jaccard_index'] 
                            for ws in window_sizes]
            ax.plot(window_sizes, jaccard_values, marker='o', label=pair)
        
        ax.set_xlabel('Window Size (bp)')
        ax.set_ylabel('Jaccard Index')
        ax.set_title('Conservation vs Window Size')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def plot_conservation_distribution(self, ax):
        """Plot conservation score distribution"""
        # This would use results from ancestral state reconstruction
        # Simplified example
        conservation_scores = np.random.beta(2, 2, 1000)  # Example data
        
        ax.hist(conservation_scores, bins=30, alpha=0.7, color='skyblue')
        ax.set_xlabel('Conservation Score')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Conservation Scores')
        ax.axvline(x=np.mean(conservation_scores), color='red', 
                  linestyle='--', label=f'Mean: {np.mean(conservation_scores):.3f}')
        ax.legend()
    
    def plot_phylogenetic_tree(self, ax):
        """Plot simplified phylogenetic tree"""
        # Simplified tree visualization
        ax.text(0.1, 0.5, 'Human', fontsize=12, ha='left', va='center')
        ax.text(0.1, 0.3, 'Monkey', fontsize=12, ha='left', va='center')
        ax.text(0.1, 0.1, 'Mouse', fontsize=12, ha='left', va='center')
        
        # Draw tree branches (simplified)
        ax.plot([0.05, 0.08], [0.5, 0.5], 'k-', linewidth=2)
        ax.plot([0.05, 0.08], [0.3, 0.3], 'k-', linewidth=2)
        ax.plot([0.05, 0.08], [0.1, 0.1], 'k-', linewidth=2)
        
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 0.6)
        ax.set_title('Phylogenetic Relationship')
        ax.axis('off')
    
    def plot_lineage_changes(self, ax):
        """Plot lineage-specific changes"""
        # Example data
        species = ['Human', 'Monkey', 'Mouse']
        gains = [45, 67, 123]
        losses = [23, 34, 89]
        
        x = np.arange(len(species))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, gains, width, label='Gains', color='green')
        bars2 = ax.bar(x + width/2, losses, width, label='Losses', color='red')
        
        ax.set_xlabel('Species')
        ax.set_ylabel('Number of Changes')
        ax.set_title('Lineage-Specific m6A Changes')
        ax.set_xticks(x)
        ax.set_xticklabels(species)
        ax.legend()
    
    def plot_summary_statistics(self, ax):
        """Plot summary statistics table"""
        # Create summary table
        summary_data = {
            'Species': ['Human', 'Monkey', 'Mouse'],
            'Total Peaks': [1234, 1567, 2345],
            'Conserved Peaks': [567, 678, 789],
            'Conservation Rate': [0.46, 0.43, 0.34],
            'Unique Peaks': [234, 345, 567]
        }
        
        df = pd.DataFrame(summary_data)
        
        # Create table
        ax.axis('tight')
        ax.axis('off')
        table = ax.table(cellText=df.values,
                        colLabels=df.columns,
                        cellLoc='center',
                        loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        
        ax.set_title('Summary Statistics')

# Usage
# dashboard = ConservationDashboard(
#     overlap_results, significance_results, 
#     window_results, phylo_results
# )
# dashboard.create_comprehensive_dashboard('conservation_dashboard.png')