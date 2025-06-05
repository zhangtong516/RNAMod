process FEATURE_COUNTS {
    storeDir "${params.outdir}/counts"

    input:
    tuple val(sample_id), path(bam), val(genome_dir), val(gtf_path)

    output:
    tuple val(sample_id), file("${sample_id}_counts.txt"), emit: count_files

    script:
    """
    featureCounts \
        -p \
        -C \
        -P \
        -B \
        -d 51 \
        -D 200 \
        -s 0 \
        -F GTF \
        -t exon \
        -T ${task.cpus} \
        -a ${gtf_path} \
        -o ${sample_id}_counts.txt \
        ${bam}
    """
}