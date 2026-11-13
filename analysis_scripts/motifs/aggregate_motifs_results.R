library(data.table)
library(stringr)
library(ggpubr)
library(xml2)
library(dplyr)

# Function to parse STREME XML output (adapted from parse_streme.R)
parse_streme_xml <- function(streme_file) {
  # Read XML file
  doc <- read_xml(streme_file)
  
  # Extract motifs
  motifs <- xml_find_all(doc, "//motif")
  
  motif_info <- list()
  
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
    
    # Extract motif attributes
    motif_id <- xml_attr(motif, "id")
    motif_alt <- xml_attr(motif, "alt")
    motif_width <- as.numeric(xml_attr(motif, "width"))
    motif_sites <- as.numeric(xml_attr(motif, "total_sites"))
    motif_llr <- as.numeric(xml_attr(motif, "test_log_pvalue"))
    motif_evalue <- as.numeric(xml_attr(motif, "test_evalue"))
    
    # Store motif information
    motif_info[[motif_id]] <- list(
      id = motif_id,
      alt = motif_alt,
      width = motif_width,
      sites = motif_sites,
      llr = motif_llr,
      evalue = motif_evalue
    )
  }
  
  return(motif_info)
}

# Function to create summary table from STREME data
create_streme_summary <- function(streme_data, sample_name) {
  if (length(streme_data) == 0) {
    return(data.frame())
  }
  
  summary_df <- do.call(rbind, lapply(names(streme_data), function(id) {
    info <- streme_data[[id]]
    data.frame(
      sample = sample_name,
      motif_id = id,
      alt_id = ifelse(is.null(info$alt) || info$alt == "", NA, info$alt),
      width = ifelse(is.null(info$width), NA, info$width),
      sites = ifelse(is.null(info$sites), NA, info$sites),
      log_pvalue = ifelse(is.null(info$llr), NA, info$llr),
      evalue = ifelse(is.null(info$evalue), NA, info$evalue),
      stringsAsFactors = FALSE
    )
  }))
  
  return(summary_df)
}

# Set the motif directory - update this path as needed
motif_dir <- "./"

# Find all STREME output directories
streme_dirs <- list.dirs(motif_dir, recursive = FALSE, full.names = TRUE)
streme_dirs <- streme_dirs[grepl("streme_", basename(streme_dirs))]

cat("Found", length(streme_dirs), "STREME output directories:\n")
print(basename(streme_dirs))

# Initialize list to store all motif data
all_motifs_list <- list()

# Process each STREME directory
for (streme_dir in streme_dirs) {
  # Extract sample name from directory name
  sample_name <- str_extract(basename(streme_dir), "(?<=streme_).+")
  
  # Look for streme.xml file
  streme_xml <- file.path(streme_dir, "streme.xml")
  
  if (file.exists(streme_xml)) {
    cat("Processing:", sample_name, "\n")
    
    # Parse STREME XML
    streme_data <- parse_streme_xml(streme_xml)
    
    # Create summary for this sample
    sample_summary <- create_streme_summary(streme_data, sample_name)
    
    if (nrow(sample_summary) > 0) {
      all_motifs_list[[sample_name]] <- sample_summary
    }
  } else {
    cat("Warning: streme.xml not found in", streme_dir, "\n")
  }
}

# Combine all motif summaries
if (length(all_motifs_list) > 0) {
  all_motifs <- rbindlist(all_motifs_list, fill = TRUE)
  
  # Sort by sample and E-value
  all_motifs <- all_motifs[order(sample, evalue)]
  
  # Print summary statistics
  cat("\n=== STREME Results Summary ===\n")
  cat("Total samples processed:", length(unique(all_motifs$sample)), "\n")
  cat("Total motifs found:", nrow(all_motifs), "\n")
  cat("Motifs per sample:\n")
  print(table(all_motifs$sample))
  
  # Save aggregated results
  output_file <- file.path(motif_dir, "aggregated_streme_results.csv")
  fwrite(all_motifs, output_file)
  cat("Aggregated results saved to:", output_file, "\n")
  
  # Create summary statistics by sample
  sample_stats <- all_motifs[, .(
    total_motifs = .N,
    avg_width = round(mean(width, na.rm = TRUE), 2),
    avg_sites = round(mean(sites, na.rm = TRUE), 2),
    min_evalue = min(evalue, na.rm = TRUE),
    max_evalue = max(evalue, na.rm = TRUE)
  ), by = sample]
  
  # Save sample statistics
  stats_file <- file.path(motif_dir, "streme_sample_statistics.csv")
  fwrite(sample_stats, stats_file)
  cat("Sample statistics saved to:", stats_file, "\n")
  
  # Create motif width distribution plot
  if (require(ggplot2, quietly = TRUE)) {
    width_plot <- ggplot(all_motifs, aes(x = factor(width), fill = sample)) +
      geom_bar(position = "dodge") +
      labs(title = "Motif Width Distribution by Sample",
           x = "Motif Width (bp)",
           y = "Number of Motifs") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plot_file <- file.path(motif_dir, "motif_width_distribution.pdf")
    ggsave(plot_file, width = 10, height = 6)
    cat("Width distribution plot saved to:", plot_file, "\n")
    
    # Create E-value distribution plot
    evalue_plot <- ggplot(all_motifs[!is.na(evalue)], aes(x = log10(evalue), fill = sample)) +
      geom_histogram(alpha = 0.7, bins = 30) +
      labs(title = "Motif E-value Distribution by Sample",
           x = "Log10(E-value)",
           y = "Number of Motifs") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    evalue_plot_file <- file.path(motif_dir, "motif_evalue_distribution.pdf")
    ggsave(evalue_plot_file, width = 10, height = 6)
    cat("E-value distribution plot saved to:", evalue_plot_file, "\n")
  }
  
  # Print top motifs by E-value
  cat("\n=== Top 10 Most Significant Motifs ===\n")
  top_motifs <- all_motifs[order(evalue)][1:min(10, nrow(all_motifs))]
  print(top_motifs[, .(sample, motif_id, width, sites, evalue)])
  
} else {
  cat("No STREME results found to aggregate.\n")
  cat("Please check that:\n")
  cat("1. STREME output directories exist in:", motif_dir, "\n")
  cat("2. Each directory contains a 'streme.xml' file\n")
  cat("3. Directory names follow the pattern 'streme_<sample_name>'\n")
}
