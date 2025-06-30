process FASTP {
    storeDir "${params.outdir}/trimmed"

    input:
    tuple val(sample_id), path(read1), path(read2), val(genome_prefix)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz"), val(genome_prefix), emit: trimmed_reads
    path "${sample_id}_fastp.json", emit: json_report
    path "${sample_id}_fastp.html", emit: html_report

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_R1.trimmed.fastq.gz \
        --out2 ${sample_id}_R2.trimmed.fastq.gz \
        --json ${sample_id}_fastp.json \
        --html ${sample_id}_fastp.html \
        --dedup \
        --detect_adapter_for_pe \
        --n_base_limit 5 \
        --length_required 50 \
        --thread ${task.cpus} \
        --qualified_quality_phred 30 \
        --thread ${task.cpus} 
    """
}
