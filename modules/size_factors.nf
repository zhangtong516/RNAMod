process CALCULATE_SIZE_FACTORS {
    storeDir "${params.outdir}/size_factors"
    
    
    
    input:
    path m6a_counts
    path input_counts
    
    output:
    path "size_factors.txt"
    
    script:
    """
    Rscript ${baseDir}/bin/calculate_size_factors.R \
        ${m6a_counts} \
        ${input_counts} \
        size_factors.txt
    """
}