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
    python ${baseDir}/bin/collect_qc_metrics.py ${fastp_json_reports} ${star_logs}
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