process FEATURE_COUNTS {
    storeDir "${params.outdir}/counts", mode: 'copy'
    
    container 'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0'

    input:
    tuple val(sample_id), path(bam)

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
        -a ${params.gtf} \
        -o ${sample_id}_counts.txt \
        ${bam}
    """
}