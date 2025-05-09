# RNA Methylation Analysis Pipeline

A Nextflow pipeline for analyzing RNA methylation (m6A-CT) data, including adapter trimming, genome alignment, insert size analysis, signal quantification, peak calling, and peak annotation.

## Peak Calling and Annotation

The pipeline includes modules for m6A-CT peak calling and annotation:

- **MACS2 Peak Calling**: Calls peaks using MACS2 with parameters `-f BAM -B -q 0.01 --nomodel --extsize 100 --keep-dup all`. High-confidence peaks with log-q-value < 10 are retained.

- **Peak Annotation**: Annotates peaks using ChIPseeker with parameters `sameStrand=TRUE, tssRegion=c(-3000, 0), genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", "Promoter", "Downstream", "Intergenic")`.

- **Peak Visualization**: Visualizes the distribution of peak sites on mRNA using the Guitar R package.

## Requirements

- Nextflow (>=21.04.0)
- Docker or Singularity
- Required tools (provided via containers):
  - Trim Galore
  - STAR (v2.7.2)
  - Picard
  - featureCounts
  - deepTools
  - R with DESeq2 package

## Input Data Requirements

- Paired-end RNA-seq reads (FASTQ format)
- Reference genome (mm9 or hg19)
- GTF annotation file

## Usage

1. Clone the repository:
```bash
git clone https://github.com/yourusername/RNA_mod.git
cd RNA_mod
```

2. Edit the configuration in `nextflow.config` to match your computing environment.

3. Prepare your input data and update the parameters in `nextflow.config`:
```groovy
params {
    reads = "path/to/reads/*_{1,2}.fastq.gz"
    genome = "mm9"  // or hg19
    gtf = "path/to/annotation.gtf"
    outdir = "results"
}
```

4. Run the pipeline:
```bash
nextflow run main.nf -profile docker
```

## Pipeline Steps

1. **Adapter Trimming** (Trim Galore)
   - Quality threshold: 30
   - Paired-end mode

2. **Genome Alignment** (STAR v2.7.2)
   - Parameters:
     - outFilterMatchNminOverLread: 0
     - outFilterScoreMinOverLread: 0
     - outFilterMatchNmin: 30
     - outFilterMismatchNmax: 3
     - outFilterMultimapNmax: 1

3. **Insert Size Analysis** (Picard)
   - CollectInsertSizeMetrics with default parameters

4. **m6A-CT Signal Quantification** (featureCounts)
   - Parameters:
     - Paired-end mode (-p)
     - Fragment counting (-P)
     - Both ends mapped (-B)
     - Minimum fragment length: 51 (-d 51)
     - Maximum fragment length: 200 (-D 200)
     - Unstranded (-s 0)
     - GTF format (-F GTF)
     - Feature type: exon (-t exon)

5. **Coverage Visualization** (deepTools)
   - bamCoverage with RPGC normalization
   - Optional INPUT library normalization using edgeR-computed size factors

## Output

The pipeline generates the following output directories:

- `results/trimmed/`: Trimmed FASTQ files
- `results/aligned/`: STAR alignment outputs (BAM files)
- `results/metrics/`: Insert size metrics and plots
- `results/counts/`: Feature counts output
- `results/bigwig/`: Normalized coverage tracks
- `results/qc/`: Quality control metrics summary (read counts, filtering rates, mapping rates)
- `results/`: Execution reports and logs

## Citation

If you use this pipeline, please cite the following tools:

- Trim Galore (https://github.com/FelixKrueger/TrimGalore)
- STAR: Dobin A, et al. (2013) STAR: ultrafast universal RNA-seq aligner
- Picard Tools (http://broadinstitute.github.io/picard)
- featureCounts: Liao Y, et al. (2014) featureCounts
- deepTools: Ramírez F, et al. (2016) deepTools2
- DESeq2: Love MI, et al. (2014) DESeq2

## License

This project is licensed under the MIT License - see the LICENSE file for details.