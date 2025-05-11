#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// Default parameters
params.samplesheet = "$projectDir/samplesheet.csv"  // Add this line
params.outdir = "$projectDir/results"
params.threads = Runtime.runtime.availableProcessors()
params.help = false

// Import modules
include { FASTP } from './modules/trimming'
include { STAR_ALIGN } from './modules/alignment'
include { INSERT_SIZE_METRICS } from './modules/metrics'
include { FEATURE_COUNTS } from './modules/quantification'
include { BAM_COVERAGE } from './modules/visualization'
include { CALCULATE_SIZE_FACTORS } from './modules/size_factors'
include { QC } from './modules/qc'
include { MACS2_PEAK_CALLING } from './modules/peak_calling'
include { PEAK_ANNOTATION } from './modules/annotation'

// Print help message
if (params.help) {
    log.info"""
    ===================================================
    RNA methylation Analysis Pipeline
    ===================================================
    
    Usage:
    nextflow run main.nf [options]
    
    Options:
      --samplesheet     CSV file containing sample information (default: $params.samplesheet)
      --outdir          Output directory (default: $params.outdir)
      --batchName       Batch name for processing (default: $params.batchName)
      --help            Show this message
    """
    exit 0
}

// Create channel for input reads

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sampleName +"__" +row.libType +"__" +row.treatment +"__" +row.replicate,
                        file(row.r1), file(row.r2)) }
    .set { input_reads }

// Main workflow
workflow {
    // Trim adapters
    FASTP(input_reads)

    // Align reads
    STAR_ALIGN(FASTP.out.trimmed_reads)

    // Analyze insert size
    INSERT_SIZE_METRICS(STAR_ALIGN.out.aligned_bam)

    // Count features and calculate size factors
    FEATURE_COUNTS(STAR_ALIGN.out.aligned_bam)
    

    // CALCULATE_SIZE_FACTORS(
    //     FEATURE_COUNTS.out.count_files.collect{ it[1] }.filter{ it.name =~ /m6A|m7G/ },
    //     FEATURE_COUNTS.out.count_files.collect{ it[1] }.filter{ it.name =~ /input/ }
    // )
    
    // Generate coverage tracks
    STAR_ALIGN.out.aligned_bam
        .map{ key, value ->
            // Extract the "grouping" key: first three parts separated by '_'
            def groupKey = key.split('__')[0..2].join('__')
            def type = key.split('__')[2]  // this determines if it's an x or y type
            [groupKey, type, value]
        }
        .groupTuple(by: 0) // Use the 'by' parameter to specify grouping by the first element
        .map { groupKey, types, values ->
            // Split into x and y lists
            def xs = []
            def ys = []
            
            for (int i = 0; i < types.size(); i++) {
                if (types[i] == 'm6A|m7G') {
                    xs.add(values[i])
                } else if (types[i] == 'Input') {
                    ys.add(values[i])
                }
            }
            
            [groupKey, xs, ys]
        }
        .set { chanel_for_peak_calling }
    
    // Call peaks using MACS2
    MACS2_PEAK_CALLING(chanel_for_peak_calling)
    
    // Annotate peaks using ChIPseeker and visualize with Guitar
    PEAK_ANNOTATION(MACS2_PEAK_CALLING.out.peaks, params.gtf)
    
    // Run QC module to collect and summarize metrics
    ch_for_qc = 
    QC(
        FASTP.out.json_report,
        STAR_ALIGN.out.log_final
    )
}

