process PEAK_ANNOTATION {
    storeDir "${params.outdir}/annotation"

    input:
    tuple val(sample_id), path(peaks), val(genome_dir), val(gtf_path), val(genome_fasta) 

    output:
    tuple val(sample_id), path("${sample_id}_annotated_peaks.csv"), emit: annotated_peaks
    tuple val(sample_id), path("${sample_id}_annotated_peaks.pdf"), emit: annotated_peaks_plots
    
    script:
    """
    Rscript ${baseDir}/bin/annotate_peaks.R \
        ${peaks} \
        ${gtf_path} \
        ${sample_id}
    """
}