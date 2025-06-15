#!/usr/bin/env Rscript
library(RCAS)
library(data.table)
library(ggpubr)
library(argparse)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-b", "--bed_file", type = "character",
    help="BED file to be processed")
parser$add_argument("-g", "--gtf_file", type = "character",
    help="GTF file for annotation")
parser$add_argument("-v", "--genome_version", type="character",
    default = "mm39", 
    help="Genome version for the annotation (example: hg38, mm39 etc) [default: mm39]")

                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

sample_name = strsplit(args$bed_file, "/")[[1]][length(strsplit(args$bed_file, "\\.")[[1]])]
tmp_bed = paste0(sample_name,".bed")
reformated_bed =  importBed(filePath = args$bed_file, sampleN = 10000, keepStandardChr =F,
                          colnames = c('chrom', 'start', 'end', 'name', 'score' ,'strand'))
reformated_bed = sort(reformated_bed)
bed_df = as.data.frame(reformated_bed)
bed_df$strand = "+"
bed_df = bed_df[,c("seqnames", "start", "end", "name", "score", "strand")]
write.table(bed_df, file = tmp_bed, sep = '\t', quote = F, row.names = F, col.names = F)




runReport(queryFilePath = tmp_bed,
          genomeVersion = args$genome_version,
          gffFilePath = args$gtf_file,
          printProcessedTables = T,
          outputDir = paste0(sample_name, "_report"),
          ) 
