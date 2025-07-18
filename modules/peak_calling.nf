process MACS2_PEAK_CALLING {
    storeDir "${params.outdir}/peaks"
    
    input:
    tuple val(sample_id), path(clip_bams), path(input_bams), val(genome_prefix)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), val(genome_prefix), emit: peaks
    tuple val(sample_id), path("${sample_id}_peaks.xls"), emit: peak_xls
    tuple val(sample_id), path("${sample_id}_summits.bed"), emit: summits
    tuple val(sample_id), path("${sample_id}_treat_pileup.bdg"), path("${sample_id}_control_lambda.bdg"), emit: bedgraphs

    script:
    def clip_bam_files = clip_bams.join(' ')
    def input_bam_files = input_bams.join(' ')
    
    """
    macs2 callpeak \
        -f BAM \
        -B \
        -q 0.1 \
        --nomodel \
        --extsize 100 \
        --keep-dup all \
        -t ${clip_bam_files} \
        -c ${input_bam_files} \
        -n ${sample_id}
    
    ### Will perform the filtering later
    # Filter high-confidence peaks (log-q-value < 10) 
    # awk '\$9 > 10' ${sample_id}_peaks.narrowPeak > ${sample_id}_filtered_peaks.narrowPeak
    
    # Rename filtered peaks to standard output name for downstream processing
    # mv ${sample_id}_filtered_peaks.narrowPeak ${sample_id}_peaks.narrowPeak
    """
}

process MACS2_PEAK_CALLING_SINGLE {
    storeDir "${params.outdir}/peaks_single_replicate"
    
    input:
    tuple val(sample_id), path(clip_bam), path(input_bam), val(genome_prefix) 

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), val(genome_prefix), emit: peaks
    tuple val(sample_id), path("${sample_id}_peaks.xls"), emit: peak_xls
    tuple val(sample_id), path("${sample_id}_summits.bed"), emit: summits
    tuple val(sample_id), path("${sample_id}_treat_pileup.bdg"), path("${sample_id}_control_lambda.bdg"), emit: bedgraphs

    script:
    """
    macs2 callpeak \
        -f BAM \
        -B \
        -q 0.1 \
        --nomodel \
        --extsize 100 \
        --keep-dup all \
        -t ${clip_bam} \
        -c ${input_bam} \
        -n ${sample_id}
    
    ### Will perform the filtering later
    # Filter high-confidence peaks (log-q-value < 10) 
    # awk '\$9 > 10' ${sample_id}_peaks.narrowPeak > ${sample_id}_filtered_peaks.narrowPeak
    
    # Rename filtered peaks to standard output name for downstream processing
    # mv ${sample_id}_filtered_peaks.narrowPeak ${sample_id}_peaks.narrowPeak
    """
}