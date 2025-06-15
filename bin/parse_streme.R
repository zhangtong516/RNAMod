# R Script to Parse STREME Output
# Required libraries
library(seqLogo)
library(Biostrings)
library(ggplot2)
library(ggseqlogo)
library(xml2)
library(jsonlite)

# Function to parse STREME XML output
parse_streme_xml <- function(streme_file) {
  # Read XML file
  doc <- read_xml(streme_file)
  
  # Extract motifs
  motifs <- xml_find_all(doc, "//motif")
  
  pwms <- list()
  motif_info <- list()
  
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
    
    # Extract motif attributes
    motif_id <- xml_attr(motif, "id")
    motif_alt <- xml_attr(motif, "alt")
    motif_width <- as.numeric(xml_attr(motif, "width"))
    motif_sites <- as.numeric(xml_attr(motif, "sites"))
    motif_llr <- as.numeric(xml_attr(motif, "llr"))
    motif_evalue <- as.numeric(xml_attr(motif, "e"))
    
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
  }
  
  return(list(pwms = pwms, info = motif_info))
}


# Function to plot sequence logo with STREME info
plot_streme_logo <- function(pwm, motif_info, title_prefix = "STREME Motif") {
  info_text <- ""
  if (!is.null(motif_info)) {
    info_parts <- c()
    if (!is.na(motif_info$sites)) info_parts <- c(info_parts, paste("Sites:", motif_info$sites))
    if (!is.na(motif_info$evalue)) info_parts <- c(info_parts, paste("E-value:", sprintf("%.2e", motif_info$evalue)))
    if (!is.na(motif_info$width)) info_parts <- c(info_parts, paste("Width:", motif_info$width))
    
    if (length(info_parts) > 0) {
      info_text <- paste(" (", paste(info_parts, collapse = ", "), ")", sep = "")
    }
  }
  
    main_title <- paste(title_prefix, motif_info$id, info_text, sep = " ")
    p <- ggseqlogo(pwm_df, method = 'bits') +
        ggtitle(main_title) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    return(p)
}

# Function to process STREME file and create plots
process_streme_and_plot <- function(streme_file, output_prefix, use_ggplot = TRUE) {
  # Create output directory

  # Parse STREME output
  cat("Parsing STREME output file:", streme_file, "\n")
  streme_data <- parse_streme_xml(streme_file)
  
  if (length(streme_data$pwms) == 0) {
    stop("No motifs found in STREME output file")
  }
  
  cat("Found", length(streme_data$pwms), "motifs\n")
  
  # Plot each motif
  for (motif_id in names(streme_data$pwms)) {
    pwm <- streme_data$pwms[[motif_id]]
    motif_info <- streme_data$info[[motif_id]]
    
    cat("Plotting motif:", motif_id, "\n")
    
    # Create filename
    filename <- paste0(output_prefix, "_plot_motifs.pdf")
    
    # Open graphics device
    pdf(filename, width = 800, height = 400)
    
    plot_streme_logo(pwm, motif_info, "STREME")
    
    dev.off()
    cat("Saved plot to:", filename, "\n")
  }
  
  return(streme_data)
}

# Function to create summary table of STREME motifs
create_streme_summary <- function(streme_data) {
  if (length(streme_data$info) == 0) {
    return(data.frame())
  }
  
  summary_df <- do.call(rbind, lapply(names(streme_data$info), function(id) {
    info <- streme_data$info[[id]]
    data.frame(
      Motif_ID = id,
      Alt_ID = ifelse(is.null(info$alt) || info$alt == "", NA, info$alt),
      Width = ifelse(is.null(info$width), NA, info$width),
      Sites = ifelse(is.null(info$sites), NA, info$sites),
      E_value = ifelse(is.null(info$evalue), NA, info$evalue),
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort by E-value
  if ("E_value" %in% colnames(summary_df) && !all(is.na(summary_df$E_value))) {
    summary_df <- summary_df[order(summary_df$E_value, na.last = TRUE), ]
  }
  
  return(summary_df)
}

# Function to export PWMs to various formats
export_streme_pwms <- function(streme_data, output_prefix) {
  
  # Export as CSV files
  for (motif_id in names(streme_data$pwms)) {
    pwm <- streme_data$pwms[[motif_id]]
    filename <- paste0(output_predix, "_", motif_id, "_pwm.csv")
    write.csv(pwm, filename, row.names = TRUE)
  }
  
  # Export as RData
  streme_pwms <- streme_data$pwms
  save(streme_pwms, file =paste0(output_predix, "_pwm.RData"))
  
  cat("PWMs exported to file:", paste0(output_predix, "_pwm.RData"), "\n")
}

# Example usage:
# streme_data <- process_streme_and_plot("streme_output.xml")
# summary_table <- create_streme_summary(streme_data)
# print(summary_table)
# export_streme_pwms(streme_data)

# Function to compare motifs (basic similarity)
compare_streme_motifs <- function(streme_data, method = "correlation") {
  motif_ids <- names(streme_data$pwms)
  n_motifs <- length(motif_ids)
  
  if (n_motifs < 2) {
    cat("Need at least 2 motifs for comparison\n")
    return(NULL)
  }
  
  # Create similarity matrix
  sim_matrix <- matrix(0, nrow = n_motifs, ncol = n_motifs)
  rownames(sim_matrix) <- colnames(sim_matrix) <- motif_ids
  
  for (i in 1:n_motifs) {
    for (j in i:n_motifs) {
      if (i == j) {
        sim_matrix[i, j] <- 1.0
      } else {
        pwm1 <- streme_data$pwms[[motif_ids[i]]]
        pwm2 <- streme_data$pwms[[motif_ids[j]]]
        
        # Simple correlation-based similarity (for same-width motifs)
        if (ncol(pwm1) == ncol(pwm2)) {
          if (method == "correlation") {
            sim <- cor(as.vector(pwm1), as.vector(pwm2))
          } else {
            # Euclidean distance
            sim <- 1 / (1 + sqrt(sum((pwm1 - pwm2)^2)))
          }
          sim_matrix[i, j] <- sim_matrix[j, i] <- sim
        }
      }
    }
  }
  
  return(sim_matrix)
}


###
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript compare_streme_motifs.R <streme_output_file> <output_prefix>\n")
  q()
}

streme_output_file = args[1] 
output_prefix = args[2] 
streme_data <- process_streme_and_plot(streme_output_file, output_prefix) 
summary_table <- create_streme_summary(streme_data)

write.table(summary_table, file = paste0(output_prefix, "_streme_summary_table.txt"), sep = "\t", quote = F, row.names = F)
similarity <- compare_streme_motifs(streme_data)
print(similarity)

write.table(similarity, file = paste0(output_prefix, "_streme_motif_similarity_matrix.txt"), sep = "\t", quote = F, row.names = T)
