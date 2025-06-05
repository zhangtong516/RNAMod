process PEAK_ANNOTATION {
    storeDir "${params.outdir}/annotation"

    input:
    tuple val(sample_id), path(peaks), val(gtf_path)

    output:
    tuple val(sample_id), path("${sample_id}_annotated_peaks.csv"), emit: annotated_peaks

    script:
    """
    Rscript ${baseDir}/bin/annotate_peaks.R \
        ${peaks} \
        ${gtf_path} \
        ${sample_id}
    """
}