# Peak Calling and Annotation

## Overview

This module implements peak calling for m6A-CT data using MACS2, followed by peak annotation with ChIPseeker and visualization with Guitar R package.

## Parameters

### MACS2 Peak Calling

Peaks are called using MACS2 with the following parameters:

```
-f BAM -B -q 0.01 --nomodel --extsize 100 --keep-dup all
```

High-confidence peaks with log-q-value < 10 are filtered and retained for downstream analysis.

### Peak Annotation

Peaks are annotated using the ChIPseeker R package with the following parameters:

```R
sameStrand = TRUE
tssRegion = c(-3000, 0)
genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", "Promoter", "Downstream", "Intergenic")
```

### Visualization

The distribution of peak sites on mRNA is visualized using the Guitar R package with default parameters.

## Dependencies

- MACS2 (v2.2.7.1)
- R (v4.1+)
- R packages:
  - ChIPseeker
  - TxDb.Mmusculus.UCSC.mm39.refGene (or appropriate genome annotation)
  - org.Mm.eg.db
  - Guitar
  - rtracklayer
  - GenomicFeatures

## Output Files

- `${sample_id}_peaks.narrowPeak`: MACS2 called peaks (filtered for high confidence)
- `${sample_id}_annotated_peaks.csv`: Annotated peaks with genomic features
- `${sample_id}_peak_distribution.pdf`: Visualization of peak distribution on mRNA

## Usage

The peak calling and annotation modules are automatically integrated into the main workflow. No additional parameters are required beyond the standard pipeline configuration.

### Command-line Integration

The peak calling module can be executed as part of the main pipeline or as a standalone step:

```bash
# As part of the main pipeline
./rna_mod.sh --input sample.bam --genome mm39 --output results/

# As a standalone module
./rna_mod.sh --module peak_calling --input sample.bam --genome mm39 --output results/
```

### Configuration

The default parameters can be customized in the configuration file `config/peak_calling.yaml`:

```yaml
peak_calling:
  macs2:
    q_value: 0.01
    extsize: 100
    high_confidence_threshold: 10  # -log10(q-value) threshold
  annotation:
    tss_region: [-3000, 0]
    priority: ["5UTR", "3UTR", "Exon", "Intron", "Promoter", "Downstream", "Intergenic"]
```

## Downstream Analysis

The peak calling results can be used for various downstream analyses:

### Motif Analysis

High-confidence peaks can be used for motif discovery using tools like HOMER or MEME:

```bash
# Example HOMER command
findMotifsGenome.pl ${sample_id}_peaks.narrowPeak mm39 motif_results/ -size 200 -mask
```

### Differential Methylation Analysis

Peak files can be used for differential methylation analysis between conditions:

```R
# Example R code for differential analysis using DiffBind
library(DiffBind)
sample_sheet <- read.csv("samples.csv")
dba <- dba(sampleSheet=sample_sheet)
dba <- dba.count(dba)
dba <- dba.contrast(dba, categories=DBA_CONDITION)
dba <- dba.analyze(dba)
diff_sites <- dba.report(dba)
```

### Integration with RNA-Seq Data

Peak annotations can be integrated with RNA-Seq data to correlate methylation with gene expression:

```R
# Example correlation analysis
library(GenomicFeatures)
library(rtracklayer)

# Load peak data and RNA-Seq data
peaks <- read.csv("${sample_id}_annotated_peaks.csv")
rna_seq <- read.csv("gene_expression.csv")

# Correlation analysis
merged_data <- merge(peaks, rna_seq, by="gene_id")
cor_result <- cor.test(merged_data$peak_score, merged_data$expression, method="spearman")
```

## Troubleshooting

### Common Issues

1. **No peaks detected**: 
   - Check input BAM file quality and alignment rate
   - Verify that the BAM file is properly indexed
   - Try relaxing the q-value threshold in the configuration

2. **Peak annotation errors**:
   - Ensure the correct genome annotation package is installed
   - Verify that chromosome names in BAM files match the annotation database

3. **Memory issues during peak calling**:
   - For large genomes or high-coverage data, increase available memory
   - Split the analysis by chromosomes and merge results

### Example Output Interpretation

The `${sample_id}_peaks.narrowPeak` file contains the following columns:

```
chromosome  start  end  name  score  strand  signalValue  pValue  qValue  peak
```

Example of annotated peaks in `${sample_id}_annotated_peaks.csv`:

```
peakID,chromosome,start,end,strand,annotation,geneID,geneName,distance,feature
Peak_1,chr1,12345678,12345778,+,Exon,ENSMUSG00000000001,Gene1,0,coding
Peak_2,chr5,78901234,78901334,-,5UTR,ENSMUSG00000000002,Gene2,0,non-coding
```