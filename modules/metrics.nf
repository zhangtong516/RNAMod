process INSERT_SIZE_METRICS {
    storeDir "${params.outdir}/metrics"
    
    container 'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0'

    input:
    tuple val(sample_id), path(bam)

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