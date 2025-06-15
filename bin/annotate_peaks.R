#!/usr/bin/env Rscript

# This script is used by the PEAK_ANNOTATION process to annotate peaks and visualize their distribution

# Load required libraries
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(org.Mm.eg.db)
library(Guitar)
library(rtracklayer)
library(GenomicFeatures)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript annotate_peaks.R <peak_file> <gtf_file> <sample_id>")
}

peak_file <- args[1]
gtf_file <- args[2]
sample_id <- args[3]

# Create TxDb from GTF if provided, otherwise use mouse reference
if (file.exists(gtf_file)) {
  txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
  gff <- importGff(gtf_file)
  txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
} else {
  txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
}


# Read peak file
peaks <- readPeakFile(peak_file)

# Annotate peaks
peakAnno <- annotatePeak(peaks, 
                       tssRegion=c(-3000, 0),
                       TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       sameStrand=TRUE,
                       genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", "Promoter", "Downstream", "Intergenic"))

# Save annotation results
write.csv(as.data.frame(peakAnno), file=paste0(sample_id, "_annotated_peaks.csv"), row.names=FALSE)


###
pdf(str_c(sample_id,"_annotated_peaks.pdf"), width=12, height=8)

plot1 = GuitarPlot(stBedFiles = peak_file, txTxdb = txdb, 
                   stGroupName = sample_id,
                   pltTxType ="mrna",enableCI=F)

plot2 = plotAnnoPie(peakAnno, ndigit = 2)

dev.off()



# Generate summary statistics
summary_data <- as.data.frame(peakAnno)
cat("\nPeak Annotation Summary for", sample_id, "\n")
cat("Total peaks:", nrow(summary_data), "\n")
cat("Peaks by annotation:\n")
print(table(summary_data$annotation))