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

// Create channel for input reads with genome and GTF paths
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> 
        // Get the genome prefix or use default
        def prefix = row.genome_prefix ?: params.default_genome_prefix
        
        // Look up genome and GTF paths from the config
        def genome_dir = params.genomes.containsKey(prefix) ? params.genomes[prefix].genome : params.genomes[params.default_genome_prefix].genome
        def gtf_path = params.genomes.containsKey(prefix) ? params.genomes[prefix].gtf : params.genomes[params.default_genome_prefix].gtf
        
        tuple(row.sampleName +"__" +row.libType +"__" +row.treatment +"__" +row.replicate,
              file(row.r1), file(row.r2),
              genome_dir, gtf_path)
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
    

    // CALCULATE_SIZE_FACTORS(
    //     FEATURE_COUNTS.out.count_files.collect{ it[1] }.filter{ it.name =~ /m6A|m7G/ },
    //     FEATURE_COUNTS.out.count_files.collect{ it[1] }.filter{ it.name =~ /input/ }
    // )
    
    // Generate coverage tracks
    STAR_ALIGN.out.aligned_bam
    .map { item ->
        def parts = item[0].split('__')
        def prefix = parts[0] + '__' + parts[1]
        def group = parts[2]  // This is 'a' or 'b'
        
        // Return a tuple of [prefix, group, value]
        return [prefix, group, item[1]]
    }
    .groupTuple(by: [0, 1])  // Group by both prefix and a/b group
    .map { prefix, group, values ->
        // At this point, we have entries like: ['a_a', 'a', [1, 2]] and ['a_a', 'b', [3, 4]]
        return [prefix, group, values]
    }.groupTuple(by: 0)  // Group by just the prefix now
    .map { prefix, groups, valuesList ->
        // Now we have: ['a_a', ['a', 'b'], [[1, 2], [3, 4]]]
        return [prefix, valuesList.flatten()]
    }
    .map { prefix, values ->
        // Group values into pairs
        def aValues = []
        def bValues = []
        
        // Find index where 'b' values start
        def half = values.size() / 2
        aValues = values[0..<half]
        bValues = values[half..<values.size()]
        
        return [prefix, aValues, bValues]
    }
    .map { prefix, aValues, bValues ->
        // Format to match the desired output
        return [prefix, aValues, bValues]
    }.set { chanel_for_peak_calling } 
    // Call peaks using MACS2
    MACS2_PEAK_CALLING(chanel_for_peak_calling)
    
    // In the workflow section, replace the PEAK_ANNOTATION line with:
    
    // Annotate peaks using ChIPseeker and visualize with Guitar
    PEAK_ANNOTATION(
        MACS2_PEAK_CALLING.out.peaks.join(
            STAR_ALIGN.out.aligned_bam.map { sample_id, bam, genome_dir, gtf_path -> 
                [sample_id, gtf_path] 
            }
        )
    )
    
    // Run QC module to collect and summarize metrics
    QC(
        FASTP.out.json_report,
        STAR_ALIGN.out.log_final
    )
}

