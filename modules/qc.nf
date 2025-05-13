process COLLECT_QC_METRICS {
    storeDir "${params.outdir}/qc"

    input:
    path fastp_json_reports
    path star_logs

    output:
    path "qc_summary.tsv", emit: qc_summary
    path "qc_summary.html", emit: qc_html

    script:

    def fastp_json_files = fastp_json_reports.collect().map{ items -> items.join(',') }
    def star_log_files = star_logs.collect().map{ items -> items.join(',') }
    """
    python ${baseDir}/bin/collect_qc_metrics.py ${fastp_json_files} ${star_log_files}
    """
}

workflow QC {
    take:
    fastp_reports
    star_logs

    main:
    COLLECT_QC_METRICS(fastp_reports.collect().join(','), star_logs.collect().join(','))

    emit:
    qc_summary = COLLECT_QC_METRICS.out.qc_summary
    qc_html = COLLECT_QC_METRICS.out.qc_html
}