process COLLECT_QC_METRICS {
    storeDir "${params.outdir}/qc", mode: 'copy'
    
    container 'https://depot.galaxyproject.org/singularity/python:3.9--1'

    input:
    path fastp_json_reports
    path star_logs

    output:
    path "qc_summary.tsv", emit: qc_summary
    path "qc_summary.html", emit: qc_html

    script:
    """
    #!/usr/bin/env python3
    import json
    import pandas as pd
    import os
    import re
    from collections import defaultdict

    # Initialize results dictionary
    results = defaultdict(dict)

    # Process fastp JSON reports
    for json_file in "${fastp_json_reports}".split():
        sample_id = os.path.basename(json_file).replace('_fastp.json', '')
        with open(json_file, 'r') as f:
            data = json.load(f)
            
            # Extract read counts and filtering stats
            results[sample_id]['total_reads'] = data['summary']['before_filtering']['total_reads']
            results[sample_id]['filtered_reads'] = data['summary']['after_filtering']['total_reads']
            results[sample_id]['filtering_rate'] = round((1 - data['summary']['after_filtering']['total_reads'] / data['summary']['before_filtering']['total_reads']) * 100, 2)
            results[sample_id]['duplication_rate'] = round(data['duplication']['rate'] * 100, 2) if 'duplication' in data else 'NA'

    # Process STAR alignment logs
    for log_file in "${star_logs}".split():
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
                too_many_loci = float(total_mapped_match.group(1))
                results[sample_id]['total_mapping_rate'] = results[sample_id]['uniquely_mapped_rate'] + results[sample_id]['multi_mapped_rate'] + too_many_loci

    # Convert to DataFrame and save as TSV
    df = pd.DataFrame.from_dict(results, orient='index')
    df.index.name = 'Sample'
    df.reset_index(inplace=True)
    
    # Reorder columns
    columns = ['Sample', 'total_reads', 'filtered_reads', 'filtering_rate', 'duplication_rate', 
              'input_reads', 'uniquely_mapped_rate', 'multi_mapped_rate', 'total_mapping_rate']
    df = df[columns]
    
    # Save to TSV
    df.to_csv('qc_summary.tsv', sep='\t', index=False)
    
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
            html += f"<td>{row[col]}</td>"
        html += "</tr>"
    
    html += """        </table>
    </body>
    </html>
    """
    
    with open('qc_summary.html', 'w') as f:
        f.write(html)
    """
}

workflow QC {
    take:
    fastp_reports
    star_logs

    main:
    COLLECT_QC_METRICS(fastp_reports.collect(), star_logs.collect())

    emit:
    qc_summary = COLLECT_QC_METRICS.out.qc_summary
    qc_html = COLLECT_QC_METRICS.out.qc_html
}