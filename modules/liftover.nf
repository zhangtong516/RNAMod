process LIFTOVER {
    tag "$sample_id"
    label 'process_low'
    storeDir "${params.outdir}/liftover"


    input:
    tuple val(sample_id), path(peaks), val(genome_prefix)
    path chain_file

    output:
    tuple val(sample_id), path("*.lifted.bed"), val(genome_prefix), emit: lifted_peaks

    script:
    def chain_file = "${params.reference_dir}/${genome_prefix}/${genome_prefix}ToHg38.chain.gz"
    def out_file = "${sample_id}_${genome_prefix}ToHg38.lifted.bed"
    if (genome_prefix == 'hg38') {
        """
        cp ${peaks} ${out_file}
        """
    } else {
        """
        liftOver \
            ${peaks} \
            ${chain_file} \
            ${out_file} \
            ${sample_id}_${genome_prefix}.unmapped.bed
        """
    }
}