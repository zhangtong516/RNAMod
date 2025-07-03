process STAR_ALIGN {
    storeDir "${params.outdir}/aligned"

    input:
    tuple val(sample_id), path(reads), val(genome_prefix)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), val(genome_prefix), emit: aligned_bam
    path "${sample_id}_Log.final.out", emit: log_final
    
    script:
    def genome_dir = "${params.reference_dir}/${genome_prefix}/STAR_v2.7.10a/" 
    """
    module load star
    STAR \
        --genomeDir ${genome_dir} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --outFilterMatchNminOverLread 0 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNmin 30 \
        --outFilterMismatchNmax 3 \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --outFileNamePrefix ${sample_id}_
    """
}