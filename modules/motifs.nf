process MOTIFS {
    storeDir "${params.outdir}/motifs"

    input:
    tuple val(sample_id), path(narrow_peaks), val(genome_dir), val(gtf_file), val(genome_fasta)

    output:
    path "${sample_id}_motifs"
    path "${sample_id}_plot_motifs.pdf"
    path "${sample_id}_streme_summary_table.txt"
    path "${sample_id}_streme_motif_similarity_matrix.txt" 

    script:    
    """
    module load singularity  
    awk 'OFS="\t" {
        summit_pos = \$2 + \$10
        start = summit_pos - 100
        end = summit_pos + 100
        if (start < 0) start = 0
        print \$1, start, end, \$4, \$5, "."
    }' ${narrow_peaks} > ${sample_id}_summit_200bp.bed 

    bedtools getfasta -fi ${genome_fasta} -bed ${sample_id}_summit_200bp.bed -fo ${sample_id}_summit_200bp.fa
    
    singularity exec -B ./:/workspace -W /workspace \
        ${params.cacheDir}/meme.sif fasta-shuffle-letters  ${sample_id}_summit_200bp.fa > shuffled_background.fa

    singularity exec -B ./:/workspace -W /workspace \
        ${params.cacheDir}/meme.sif streme --p  ${sample_id}_summit_200bp.fa \
            --n shuffled_background.fa --rna --minw 5 --maxw 8 \
            --nmotifs 5 --thresh 0.0001 \
            --oc ${sample_id}_motifs
    Rscript ${baseDir}/bin/parse_streme.R ${sample_id}_motifs/streme.xml ${sample_id}

    """
}