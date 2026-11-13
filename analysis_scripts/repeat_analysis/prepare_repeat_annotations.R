#!/usr/bin/env Rscript
# prepare_repeat_annotations.R
# Script to prepare repeat annotations from RepeatMasker output

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# Function to process RepeatMasker output
process_repeatmasker <- function(rm_file, output_file) {
  # Read RepeatMasker output
  # Standard RepeatMasker format has these columns:
  # SW_score, perc_div, perc_del, perc_ins, query_seq, query_start, query_end, 
  # query_left, strand, repeat_name, repeat_class, repeat_start, repeat_end, repeat_left, ID
  
  rm_data <- read.table(rm_file, header = FALSE, skip = 3, stringsAsFactors = FALSE)
  
  # Assign column names
  colnames(rm_data) <- c(
    "sw_score", "perc_div", "perc_del", "perc_ins", "chr", "start", "end",
    "query_left", "strand", "repeat_name", "repeat_class_family", 
    "repeat_start", "repeat_end", "repeat_left", "id"
  )
  
  # Clean chromosome names
  rm_data$chr <- gsub("chr", "", rm_data$chr)
  rm_data$chr <- paste0("chr", rm_data$chr)
  
  # Split repeat_class_family into class and family
  class_family <- strsplit(rm_data$repeat_class_family, "/")
  rm_data$repeat_class <- sapply(class_family, function(x) x[1])
  rm_data$repeat_family <- sapply(class_family, function(x) ifelse(length(x) > 1, x[2], x[1]))
  
  # Filter out low-quality annotations (optional)
  rm_data_filtered <- rm_data[rm_data$perc_div < 30, ]  # Less than 30% divergence
  
  # Create output format
  output_data <- data.frame(
    chr = rm_data_filtered$chr,
    start = rm_data_filtered$start,
    end = rm_data_filtered$end,
    repeat_name = rm_data_filtered$repeat_name,
    repeat_class = rm_data_filtered$repeat_class,
    repeat_family = rm_data_filtered$repeat_family,
    strand = rm_data_filtered$strand,
    perc_div = rm_data_filtered$perc_div,
    length = rm_data_filtered$end - rm_data_filtered$start + 1
  )
  
  # Write output
  write.table(output_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Processed", nrow(output_data), "repeat elements\n")
  cat("Output saved to:", output_file, "\n")
  
  return(output_data)
}

# Function to download and process UCSC repeat annotations
download_ucsc_repeats <- function(genome = "hg38", output_file = NULL) {
  if (is.null(output_file)) {
    output_file <- paste0(genome, "_repeats.txt")
  }
  
  # UCSC table browser URL for RepeatMasker track
  ucsc_url <- paste0(
    "http://hgdownload.cse.ucsc.edu/goldenPath/", genome, 
    "/database/rmsk.txt.gz"
  )
  
  # Download and process
  temp_file <- tempfile(fileext = ".txt.gz")
  download.file(ucsc_url, temp_file)
  
  # Read UCSC format
  ucsc_data <- read.table(gzfile(temp_file), header = FALSE, stringsAsFactors = FALSE)
  
  # UCSC rmsk format columns
  colnames(ucsc_data) <- c(
    "bin", "sw_score", "milliDiv", "milliDel", "milliIns", "chr", "start", "end",
    "genoLeft", "strand", "repeat_name", "repeat_class", "repeat_family",
    "repStart", "repEnd", "repLeft", "id"
  )
  
  # Process data
  output_data <- data.frame(
    chr = ucsc_data$chr,
    start = ucsc_data$start + 1,  # Convert to 1-based coordinates
    end = ucsc_data$end,
    repeat_name = ucsc_data$repeat_name,
    repeat_class = ucsc_data$repeat_class,
    repeat_family = ucsc_data$repeat_family,
    strand = ucsc_data$strand,
    perc_div = ucsc_data$milliDiv / 10,  # Convert to percentage
    length = ucsc_data$end - ucsc_data$start
  )
  
  # Filter standard chromosomes
  standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
  output_data <- output_data[output_data$chr %in% standard_chrs, ]
  
  # Write output
  write.table(output_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Downloaded and processed", nrow(output_data), "repeat elements for", genome, "\n")
  cat("Output saved to:", output_file, "\n")
  
  unlink(temp_file)
  return(output_data)
}

# Example usage
if (FALSE) {
  # Download repeat annotations for different genomes
  hg38_repeats <- download_ucsc_repeats("hg38", "data/annotations/hg38_repeats.txt")
  # Or process existing RepeatMasker files
  # process_repeatmasker("path/to/repeatmasker.out", "processed_repeats.txt")
}

cat("Repeat annotation preparation script loaded!\n")