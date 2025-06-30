process INSERT_SIZE_METRICS {
    storeDir "${params.outdir}/metrics"

    input:
    tuple val(sample_id), path(bam), val(genome_prefix) 

    output:
    path "${sample_id}_insert_size_metrics.txt"
    path "${sample_id}_insert_size_histogram.pdf"

    script:
    """
    picard CollectInsertSizeMetrics \
        I=${bam} \
        O=${sample_id}_insert_size_metrics.txt \
        H=${sample_id}_insert_size_histogram.pdf
    """
}