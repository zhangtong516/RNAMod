#!/usr/bin/env python3
import json
import pandas as pd
import os
import re
from collections import defaultdict
import sys

# Initialize results dictionary
results = defaultdict(dict)

# Process fastp JSON reports
fastp_json_reports = sys.argv[1]
star_logs = sys.argv[2]
summary_tsv = sys.argv[3]
summary_html = sys.argv[4]

for json_file in fastp_json_reports.split(','):
    sample_id = os.path.basename(json_file).replace('_fastp.json', '')
    with open(json_file, 'r') as f:
        data = json.load(f)
        
        # Extract read counts and filtering stats
        results[sample_id]['total_reads'] = data['summary']['before_filtering']['total_reads'] / 2 
        results[sample_id]['filtered_reads'] = data['summary']['after_filtering']['total_reads'] / 2
        results[sample_id]['filtering_rate'] = round((1 - data['summary']['after_filtering']['total_reads'] / data['summary']['before_filtering']['total_reads']) * 100, 2)
        results[sample_id]['duplication_rate'] = round(data['duplication']['rate'] * 100, 2) if 'duplication' in data else 'NA'

# Process STAR alignment logs
for log_file in star_logs.split(','):
    sample_id = os.path.basename(log_file).replace('_Log.final.out', '')
    with open(log_file, 'r') as f:
        content = f.read()
        
        # Extract mapping rates
        input_reads_match = re.search(r'Number of input reads \|\s+(\d+)', content)
        if input_reads_match:
            results[sample_id]['input_reads'] = int(input_reads_match.group(1))
            
        uniquely_mapped_match = re.search(r'Uniquely mapped reads % \|\s+([\d\.]+)%', content)
        if uniquely_mapped_match:
            results[sample_id]['uniquely_mapped_rate'] = float(uniquely_mapped_match.group(1))
            
        multi_mapped_match = re.search(r'% of reads mapped to multiple loci \|\s+([\d\.]+)%', content)
        if multi_mapped_match:
            results[sample_id]['multi_mapped_rate'] = float(multi_mapped_match.group(1))
            
        total_mapped_match = re.search(r'% of reads mapped to too many loci \|\s+([\d\.]+)%', content)
        if total_mapped_match:
            results[sample_id]['too_many_loci_rate'] = float(total_mapped_match.group(1))
            too_many_loci = float(total_mapped_match.group(1))
            # Calculate total mapping rate only if the required rates are available
            uniquely_mapped = results[sample_id].get('uniquely_mapped_rate', 0)
            multi_mapped = results[sample_id].get('multi_mapped_rate', 0)
            results[sample_id]['total_mapping_rate'] = uniquely_mapped + multi_mapped + too_many_loci

# Convert to DataFrame and save as TSV
df = pd.DataFrame.from_dict(results, orient='index')
df.index.name = 'Sample'
df.reset_index(inplace=True)

# Define all possible columns
all_columns = ['Sample', 'total_reads', 'filtered_reads', 'filtering_rate', 'duplication_rate', 
          'input_reads', 'uniquely_mapped_rate', 'multi_mapped_rate', 'too_many_loci_rate']

# Filter to only include columns that exist in the DataFrame
columns = ['Sample'] + [col for col in all_columns[1:] if col in df.columns]
df = df[columns]

# Save to TSV
df.to_csv(summary_tsv, sep='\t', index=False)

# Generate HTML report
html = """<!DOCTYPE html>
<html>
<head>
    <title>RNA Methylation Pipeline QC Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .header { background-color: #4CAF50; color: white; padding: 10px; margin-bottom: 20px; }
    </style>
</head>
<body>
    <div class="header">
        <h1>RNA Methylation Pipeline QC Summary</h1>
    </div>
    <table>
        <tr>
"""

# Add table headers
html += "<th>" + "</th><th>".join(columns) + "</th></tr>"

# Add data rows
for _, row in df.iterrows():
    html += "<tr>"
    for col in columns:
        # Format numeric values with comma separators
        if col != 'Sample' and isinstance(row[col], (int, float)) and not pd.isna(row[col]):
            if isinstance(row[col], int):
                # Format integers with comma separators
                formatted_value = f"{row[col]:,}"
            elif col.endswith('_rate'):
                # Keep percentage values with 2 decimal places
                formatted_value = f"{row[col]:.2f}"
            else:
                # Format floats with comma separators and 2 decimal places
                formatted_value = f"{row[col]:,.2f}"
            html += f"<td>{formatted_value}</td>"
        else:
            # Keep non-numeric values as is
            html += f"<td>{row[col]}</td>"
    html += "</tr>"

html += """        </table>
</body>
</html>
"""

with open(summary_html, 'w') as f:
    f.write(html)