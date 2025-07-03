#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.samplesheet = "$projectDir/samplesheet.csv"  // Add this line
params.outdir = "$projectDir/results/rnaseq"
params.threads = Runtime.runtime.availableProcessors()
params.help = false

// Import modules
include { FASTP } from './modules/trimming'
include { STAR_ALIGN } from './modules/alignment'
include { FEATURE_COUNTS } from './modules/quantification'
include { QC } from './modules/qc'


// Print help message
if (params.help) {
    log.info"""
    ===================================================
    RNA seq Analysis Pipeline
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

// Create channel for input reads with genome and GTF paths
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> 
        tuple(row.sampleName +"__" +row.libType +"__" +row.treatment +"__" +row.replicate,
              file(row.r1), file(row.r2), row.genome_prefix)
    }
    .set { input_reads }

// Main workflow
workflow {
    // Trim adapters
    FASTP(input_reads)

    // Align reads with sample-specific genome
    STAR_ALIGN(FASTP.out.trimmed_reads)

    // Count features and calculate size factors
    FEATURE_COUNTS(STAR_ALIGN.out.aligned_bam)
    
    // Run QC module to collect and summarize metrics
    QC(
        FASTP.out.json_report,
        STAR_ALIGN.out.log_final
    )
}

