#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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


// Default parameters
params.samplesheet = "$projectDir/samplesheet.csv"  // Add this line
params.outdir = "$projectDir/results"
params.threads = Runtime.runtime.availableProcessors()
params.help = false


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
input_reads= process_sample_name(params.samplesheet, params.batchName, "reads")

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
        }.groupTuple( it[0])
        .map { groupKey, list ->
            // Split into x and y lists
            def xs = list.findAll { it[1] == 'm6A|m7G' }.collect { it[2] }
            def ys = list.findAll { it[1] == 'Input' }.collect { it[2] }
            [groupKey, xs, ys]
        }.set { chanel_for_peak_calling}

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

// extra functions to play with the sample name
// Sample is named in the format: [sample_name]_[libtype]_[replicate]
// e.g.  ES_DC_m7G_2 and ES_DC_input_2 
def process_sample_name(condition, mode="merge_replicate", sep="__") {
    def p = condition.toString().split(sep)
    if(mode == "sample") { result = p[0] + sep + p[1] }
    if(mode == "sample_lib") { result = p[0] + sep + p[1] + sep + p[2] }
    if(mode == "sample_lib_rep") { result = p[0] + sep + p[1] + sep + p[2] + sep + p[3] } 
    return(result)
}

def process_sampleinfo(tsvFile, batchName, mode) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def sample          = row.sampleName // Sample01
            def batch           = row.compBatch //example: batch1
            def libType         = row.libType // m6A or m7G etc. 
            def treatment       = row.treatment // Input or treated etc. 
            def replicate       = row.replicate // 1 or 2 etc. 
            def read1           = row.r1 // file with full path
            def read2           = row.r2 // file with full path
            def condition       = sample + "__" + libType + "__" + treatment 
            def unique_name     = condition + "__" + replicate

            if (mode == "reads" && batch == batchName) return [ unique_name, file(read1), file(read2) ]
        }
        .unique()
    } 