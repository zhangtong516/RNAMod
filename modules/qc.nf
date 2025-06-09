process COLLECT_QC_METRICS {
    storeDir "${params.outdir}/qc"

    input:
    path fastp_json_reports
    path star_logs

    output:
    path "${params.batchName}_qc_summary.tsv", emit: qc_summary
    path "${params.batchName}_qc_summary.html", emit: qc_html

    script:
    def fastp_json_files = fastp_json_reports.join(',')
    def star_log_files = star_logs.join(',')
    def summary_tsv = "${params.batchName}_qc_summary.tsv"
    def summary_html = "${params.batchName}_qc_summary.html"

    """
    python ${baseDir}/bin/collect_qc_metrics.py ${fastp_json_files} ${star_log_files} ${summary_tsv} ${summary_html} 
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