#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.samplesheet = "$projectDir/samplesheet.csv"  // Add this line
params.outdir = "$projectDir/results"
params.threads = Runtime.runtime.availableProcessors()
params.help = false
params.liftover = true // Add liftover option

// Import modules
include { FASTP } from './modules/trimming'
include { STAR_ALIGN } from './modules/alignment'
include { INSERT_SIZE_METRICS } from './modules/metrics'
include { FEATURE_COUNTS } from './modules/quantification'
// include { BAM_COVERAGE } from './modules/visualization'
include { CALCULATE_SIZE_FACTORS } from './modules/size_factors'
include { QC } from './modules/qc'
include { MACS2_PEAK_CALLING } from './modules/peak_calling'
include { MACS2_PEAK_CALLING_SINGLE } from './modules/peak_calling'
include { PEAK_ANNOTATION } from './modules/annotation' 
include { MOTIFS } from './modules/motifs'
include { LIFTOVER } from './modules/liftover' // Import the new module

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
      --liftover        Set to true to perform liftover (default: $params.liftover)
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

    // Analyze insert size
    INSERT_SIZE_METRICS(STAR_ALIGN.out.aligned_bam)

    // Count features and calculate size factors
    FEATURE_COUNTS(STAR_ALIGN.out.aligned_bam)
    

    // prepare for diffBind: non-merged single replicates peak files. 
    STAR_ALIGN.out.aligned_bam
    .map { item ->
        def parts = item[0].split('__')
        def prefix = parts[0] + '__' + parts[1] + "__" + parts[3]
        def group = parts[2]  // This is 'a' or 'b'
        
        // Return a tuple of [prefix, group, value]
        return [prefix, group, item[1], item[2]]
    }.branch { 
        inputs: it[1] =~ /Input/
        treated: true
    }.set{ branched_bam }

    ch_bam_inputs = branched_bam.inputs.map{ item -> return [item[0], item[2], item[3]] }
    ch_bam_treated = branched_bam.treated.map{ item -> return[item[0], item[2]] } 
    chanel_for_peak_calling_single = ch_bam_treated.join(ch_bam_inputs)
    // Call peaks using MACS2
    MACS2_PEAK_CALLING_SINGLE(chanel_for_peak_calling_single)  

    // Generate coverage tracks
    STAR_ALIGN.out.aligned_bam
    .map { item ->
        def parts = item[0].split('__')
        def prefix = parts[0] + '__' + parts[1]
        def group = parts[2]  // This is 'a' or 'b'

        // Return a tuple of [prefix, group, value]
        return [prefix, group, item[1], item[2]]
    }
    .groupTuple(by: [0, 1, 3])
    // .view { v -> "1st: ${v}" }
    .map { prefix, group, values, genome_prefix  ->
        // At this point, we have entries like: ['a_a', 'a', [1, 2]] and ['a_a', 'b', [3, 4]]
        return [prefix, group, genome_prefix, values]
    }
    .groupTuple(by: [0,1, 2] )  // Group by just the prefix now
    // .view{ v -> "2last: ${v}"}
    .map { prefix, group, genome_prefix, valuesList ->
        // Now we have: ['a_a', ['a', 'b'], [[1, 2], [3, 4]]]
        return [prefix, group, genome_prefix,  valuesList.flatten()]
    }
    .branch {
        inputs: it[1] =~ /Input/
        treated: true
    }.set{ branched_bam_list }

    ch_bam_list_inputs = branched_bam_list.inputs.map{ item -> return [item[0], item[3], item[2]] }
    ch_bam_list_treated = branched_bam_list.treated.map{ item -> return[item[0], item[3]] }
    chanel_for_peak_calling = ch_bam_list_treated.join(ch_bam_list_inputs)
    // Call peaks using MACS2
    MACS2_PEAK_CALLING(chanel_for_peak_calling)

    // Call motifs using STREME with genome_dir and gtf_path
    MOTIFS(MACS2_PEAK_CALLING.out.peaks)

    // Liftover non-hg38 peaks if requested
    if (params.liftover) {
        LIFTOVER(MACS2_PEAK_CALLING.out.peaks)
        PEAK_ANNOTATION(LIFTOVER.out.lifted_peaks)
    } else {
        PEAK_ANNOTATION(MACS2_PEAK_CALLING.out.peaks)
    }
    
    // Run QC module to collect and summarize metrics
    QC(
        FASTP.out.json_report,
        STAR_ALIGN.out.log_final
    )
}

