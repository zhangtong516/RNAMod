process BAM_COVERAGE {
    storeDir "${params.outdir}/bigwig", mode: 'copy'
    
    container 'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0'

    input:
    tuple val(sample_id), path(bam), path(input_bam)
    path size_factors

    output:
    path "${sample_id}.bw"

    script:
    def genome_size = params.genome == 'mm9' ? 2620345972 : 2864785220
    def size_factor_map = size_factors.text.split('\n')
        .findAll { it.trim() && !it.startsWith('sample') }
        .collectEntries { line ->
            def (sample, factor) = line.split('\t')
            [(sample): factor]
        }
    def mod_factor = size_factor_map[sample_id]
    def input_sample = sample_id.replaceAll(/(m6A|m7G)/, 'input')
    def input_factor = size_factor_map[input_sample]
    """
    bamCompare \
        --bamfile1 ${bam} \
        --bamfile2 ${input_bam} \
        --outFileName ${sample_id}.bw \
        --scaleFactors ${mod_factor}:${input_factor} \
        --effectiveGenomeSize ${genome_size} \
        --numberOfProcessors ${task.cpus} \
        --operation ratio
    """
}