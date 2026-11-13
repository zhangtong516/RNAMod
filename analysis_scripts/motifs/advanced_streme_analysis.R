library(data.table)
library(stringr)
library(ggplot2)
library(xml2)
library(dplyr)
library(pheatmap)

# Set default motif directory
motif_dir <- "./"

# Source the main aggregation script functions
source("aggregate_motifs_results.R")

# Function to extract position weight matrices from STREME XML
extract_streme_pwms <- function(streme_file) {
  doc <- read_xml(streme_file)
  motifs <- xml_find_all(doc, "//motif")
  
  pwms <- list()
  
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
    motif_id <- xml_attr(motif, "id")
    
    # Extract position-specific scoring matrix
    pos_nodes <- xml_find_all(motif, ".//pos")
    
    if (length(pos_nodes) > 0) {
      # Initialize matrix
      pwm_matrix <- matrix(0, nrow = 4, ncol = length(pos_nodes))
      rownames(pwm_matrix) <- c("A", "C", "G", "U")
      
      for (j in seq_along(pos_nodes)) {
        pos <- pos_nodes[[j]]
        pwm_matrix["A", j] <- as.numeric(xml_attr(pos, "A"))
        pwm_matrix["C", j] <- as.numeric(xml_attr(pos, "C"))
        pwm_matrix["G", j] <- as.numeric(xml_attr(pos, "G"))
        pwm_matrix["U", j] <- as.numeric(xml_attr(pos, "U"))
      }
      
      pwms[[motif_id]] <- pwm_matrix
    }
  }
  
  return(pwms)
}

# Function to compare motifs across samples
compare_motifs_across_samples <- function(motif_dir = "./") {
  # Get all STREME directories
  streme_dirs <- list.dirs(motif_dir, recursive = FALSE, full.names = TRUE)
  streme_dirs <- streme_dirs[grepl("streme_", basename(streme_dirs))]
  
  # Extract all PWMs
  all_pwms <- list()
  sample_motif_map <- list()
  
  for (streme_dir in streme_dirs) {
    sample_name <- str_extract(basename(streme_dir), "(?<=streme_).+")
    streme_xml <- file.path(streme_dir, "streme.xml")
    
    if (file.exists(streme_xml)) {
      pwms <- extract_streme_pwms(streme_xml)
      
      for (motif_id in names(pwms)) {
        full_id <- paste(sample_name, motif_id, sep = "_")
        all_pwms[[full_id]] <- pwms[[motif_id]]
        sample_motif_map[[full_id]] <- sample_name
      }
    }
  }
  
  return(list(pwms = all_pwms, sample_map = sample_motif_map))
}

# Function to create motif consensus sequences
create_consensus_sequences <- function(pwms) {
  consensus_seqs <- sapply(pwms, function(pwm) {
    consensus <- apply(pwm, 2, function(col) {
      names(col)[which.max(col)]
    })
    paste(consensus, collapse = "")
  })
  
  return(consensus_seqs)
}

# Function to cluster similar motifs
cluster_similar_motifs <- function(motif_dir = "./", 
                                   similarity_threshold = 0.7) {
  
  # Get motif data
  motif_data <- compare_motifs_across_samples(motif_dir)
  pwms <- motif_data$pwms
  sample_map <- motif_data$sample_map
  
  if (length(pwms) < 2) {
    cat("Need at least 2 motifs for clustering\n")
    return(NULL)
  }
  
  # Create consensus sequences
  consensus_seqs <- create_consensus_sequences(pwms)
  
  # Calculate similarity matrix based on consensus sequences
  n_motifs <- length(consensus_seqs)
  motif_names <- names(consensus_seqs)
  
  sim_matrix <- matrix(0, nrow = n_motifs, ncol = n_motifs)
  rownames(sim_matrix) <- colnames(sim_matrix) <- motif_names
  
  for (i in 1:n_motifs) {
    for (j in i:n_motifs) {
      if (i == j) {
        sim_matrix[i, j] <- 1.0
      } else {
        seq1 <- consensus_seqs[i]
        seq2 <- consensus_seqs[j]
        
        # Simple similarity based on common subsequences
        if (nchar(seq1) == nchar(seq2)) {
          matches <- sum(strsplit(seq1, "")[[1]] == strsplit(seq2, "")[[1]])
          sim <- matches / nchar(seq1)
        } else {
          # For different lengths, use longest common subsequence approach
          sim <- 0  # Simplified for now
        }
        
        sim_matrix[i, j] <- sim_matrix[j, i] <- sim
      }
    }
  }
  
  # Create clustering plot
  if (n_motifs > 1) {
    # Add sample information as annotation
    sample_annotation <- data.frame(
      Sample = unlist(sample_map[motif_names]),
      row.names = motif_names
    )
    
    # Create heatmap
    pheatmap(sim_matrix,
             annotation_row = sample_annotation,
             annotation_col = sample_annotation,
             main = "Motif Similarity Matrix",
             filename = file.path(motif_dir, "motif_similarity_heatmap.pdf"),
             width = 10, height = 8)
    
    cat("Motif similarity heatmap saved to:", 
        file.path(motif_dir, "motif_similarity_heatmap.pdf"), "\n")
  }
  
  return(list(
    similarity_matrix = sim_matrix,
    consensus_sequences = consensus_seqs,
    sample_mapping = sample_map
  ))
}

# Function to generate comprehensive motif report
generate_motif_report <- function(motif_dir = "./") {
  
  cat("=== Generating Comprehensive STREME Motif Report ===\n\n")
  
  # 1. Basic aggregation
  cat("1. Aggregating basic motif information...\n")
  # This will run the main aggregation script
  
  # 2. Advanced clustering analysis
  cat("2. Performing motif similarity analysis...\n")
  clustering_results <- cluster_similar_motifs(motif_dir)
  
  # 3. Create summary report
  report_file <- file.path(motif_dir, "comprehensive_motif_report.txt")
  
  sink(report_file)
  cat("COMPREHENSIVE STREME MOTIF ANALYSIS REPORT\n")
  cat("==========================================\n\n")
  cat("Generated on:", Sys.time(), "\n\n")
  
  # Read aggregated results if available
  agg_file <- file.path(motif_dir, "aggregated_streme_results.csv")
  if (file.exists(agg_file)) {
    motif_data <- fread(agg_file)
    
    cat("SUMMARY STATISTICS:\n")
    cat("Total samples:", length(unique(motif_data$sample)), "\n")
    cat("Total motifs:", nrow(motif_data), "\n")
    cat("Average motifs per sample:", round(nrow(motif_data) / length(unique(motif_data$sample)), 2), "\n\n")
    
    cat("MOTIF WIDTH DISTRIBUTION:\n")
    print(table(motif_data$width))
    cat("\n")
    
    cat("TOP 5 MOST SIGNIFICANT MOTIFS:\n")
    top_motifs <- motif_data[order(evalue)][1:5]
    print(top_motifs[, .(sample, motif_id, width, sites, evalue)])
    cat("\n")
  }
  
  if (!is.null(clustering_results)) {
    cat("MOTIF SIMILARITY ANALYSIS:\n")
    cat("Number of motifs compared:", nrow(clustering_results$similarity_matrix), "\n")
    cat("Average similarity:", round(mean(clustering_results$similarity_matrix[upper.tri(clustering_results$similarity_matrix)]), 3), "\n\n")
    
    cat("CONSENSUS SEQUENCES:\n")
    for (i in 1:min(10, length(clustering_results$consensus_sequences))) {
      motif_name <- names(clustering_results$consensus_sequences)[i]
      consensus <- clustering_results$consensus_sequences[i]
      sample <- clustering_results$sample_mapping[[motif_name]]
      cat(sprintf("%-20s %-10s %s\n", motif_name, sample, consensus))
    }
  }
  
  sink()
  
  cat("Comprehensive report saved to:", report_file, "\n")
  cat("\n=== Analysis Complete ===\n")
}

# Run the comprehensive analysis if script is executed directly
if (!interactive()) {
  generate_motif_report()
}