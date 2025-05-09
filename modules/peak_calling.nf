process MACS2_PEAK_CALLING {
    storeDir "${params.outdir}/peaks", mode: 'copy'
    
    container 'https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py39hf95cd2a_4'

    input:
    tuple val(sample_id), path(clip_bams), path(input_bams)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), emit: peaks
    tuple val(sample_id), path("${sample_id}_peaks.xls"), emit: peak_xls
    tuple val(sample_id), path("${sample_id}_summits.bed"), emit: summits
    tuple val(sample_id), path("${sample_id}_treat_pileup.bdg"), path("${sample_id}_control_lambda.bdg"), emit: bedgraphs

    script:
    def m6a_bam_files = m6a_bams.join(' ')
    def input_bam_files = input_bams.join(' ')
    
    """
    macs2 callpeak \
        -f BAM \
        -B \
        -q 0.01 \
        --nomodel \
        --extsize 100 \
        --keep-dup all \
        -t ${m6a_bam_files} \
        -c ${input_bam_files} \
        -n ${sample_id}
    
    # Filter high-confidence peaks (log-q-value < 10)
    awk '$9 > 10' ${sample_id}_peaks.narrowPeak > ${sample_id}_filtered_peaks.narrowPeak
    
    # Rename filtered peaks to standard output name for downstream processing
    mv ${sample_id}_filtered_peaks.narrowPeak ${sample_id}_peaks.narrowPeak
    """
}