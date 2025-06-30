process MOTIFS {
    storeDir "${params.outdir}/motifs"

    input:
    tuple val(sample_id), path(narrow_peaks), val(genome_prefix) 

    output:
    file("${sample_id}_motifs_logoPlot.pdf")
    file("${sample_id}_motif_summary.csv")

    script:    
    def genome_fasta = "${params.reference_dir}/${genome_prefix}/${genome_prefix}.fa"
    
    """
    bedtools getfasta -fi ${genome_fasta} -bed ${narrow_peaks} -fo ${sample_id}_narrowpeak.fa
    bedtools shuffle -i ${narrow_peaks} -g ${genome_fasta}.fai -excl ${narrow_peaks} > ${sample_id}_background.bed
    bedtools getfasta -fi ${genome_fasta} -bed ${sample_id}_background.bed -fo ${sample_id}_background.fa
    # singularity exec -B ./:/workspace -W /workspace \
    #      ${params.cacheDir}/meme.sif fasta-shuffle-letters ${sample_id}_narrowpeak.fa > shuffled_background.fa

    Rscript ${baseDir}/bin/find_motifs.R ${sample_id}_narrowpeak.fa ${sample_id}_background.fa ${sample_id}

    """
}