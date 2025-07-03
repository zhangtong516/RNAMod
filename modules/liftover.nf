process LIFTOVER {
    storeDir "${params.outdir}/liftover"

    input:
    tuple val(sample_id), path(peaks), val(genome_prefix)

    output:
    tuple val(sample_id), path("${sample_id}_${genome_prefix}ToHg38.lifted.bed"), val(genome_prefix), emit: lifted_peaks

    script:
    def chain_file = "${params.reference_dir}/${genome_prefix}/${genome_prefix}ToHg38.over.chain.gz"
    def out_file = "${sample_id}_${genome_prefix}ToHg38.lifted.bed"
    if (genome_prefix == 'hg38') {
        """
        cp ${peaks} ${out_file}
        """
    } else {
        """
        cut -f1-6  ${peaks} > tmp.bed 
        liftOver \
            tmp.bed \
            ${chain_file} \
            tmp_out.bed \
            tmp_out.unmapped
        Rscript ${baseDir}/bin/complete_lifted_bed.R tmp_out.bed ${peaks} tmp_out.unsorted.bed 
        bedtools sort -i tmp_out.unsorted.bed > ${out_file} 
        """
    }
}