process PEAK_ANNOTATION {
    storeDir "${params.outdir}/annotation"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(peaks), val(geonome_dir), val(genome_prefix) 

    output:
    tuple val(sample_id), path("${sample_id}_annotated_peaks.csv"), emit: annotated_peaks
    tuple val(sample_id), path("${sample_id}_annotated_peaks.pdf"), emit: annotated_peaks_plots
    
    script:
    def gtf_file = "${params.reference_dir}/hg38/hg38.ncbiRefSeq.gtf" // since all peaks are liftovered to hg38 
    """ 
    Rscript ${baseDir}/bin/annotate_peaks.R \
        ${peaks} \
        ${gtf_file} \
        ${sample_id}
    """
}