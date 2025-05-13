process PEAK_ANNOTATION {
    storeDir "${params.outdir}/annotation"

    input:
    tuple val(sample_id), path(peaks)
    path(gtf)

    output:
    tuple val(sample_id), path("${sample_id}_annotated_peaks.csv"), emit: annotated_peaks
    // tuple val(sample_id), path("${sample_id}_peak_distribution.pdf"), emit: peak_distribution

    script:
    """
    Rscript ${baseDir}/bin/annotate_peaks.R \
        ${peaks} \
        ${gtf} \
        ${sample_id}
    """
}