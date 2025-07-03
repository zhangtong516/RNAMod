process PEAK_ANNOTATION {
    storeDir "${params.outdir}/annotation"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(peaks), val(genome_prefix)

    output:
    tuple val(sample_id), path("${sample_id}_annotated_peaks.csv"), emit: annotated_peaks
    tuple val(sample_id), path("${sample_id}_annotated_peaks.pdf"), emit: annotated_peaks_plots
    
    script:
    def gtf_file
    if (params.liftover) {
        gtf_file = params.genomes['hg38'].gtf
    } else {
        gtf_file = params.genomes[genome_prefix].gtf
    }
    
    """ 
    Rscript ${baseDir}/bin/annotate_peaks.R \
        ${peaks} \
        ${gtf_file} \
        ${sample_id}
    """
}