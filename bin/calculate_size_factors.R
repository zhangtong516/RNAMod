#!/usr/bin/env Rscript

# Load required libraries
library(edgeR)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript calculate_size_factors.R <m6a_counts> <input_counts> <output_file>")
}

# Read count files
m6a_counts <- read.table(args[1], header=TRUE, row.names=1)
input_counts <- read.table(args[2], header=TRUE, row.names=1)

# Calculate size factors using edgeR
m6a_dge <- DGEList(counts = m6a_counts)
input_dge <- DGEList(counts = input_counts)

m6a_size_factor <- calcNormFactors(m6a_dge)
input_size_factor <- calcNormFactors(input_dge)

# Write size factors to output file
write.table(
    data.frame(
        sample = c(colnames(m6a_counts), colnames(input_counts)),
        size_factor = c(m6a_size_factor$samples$norm.factors, input_size_factor$samples$norm.factors)
    ),
    file = args[3],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)