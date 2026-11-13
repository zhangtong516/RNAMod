#!/usr/bin/env Rscript
# RNA_Mod_in_repeat.R
# Analysis of m6A modifications on repetitive elements across species

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(EnrichedHeatmap)
library(genomation)

# Function to load repeat annotations
load_repeat_annotations <- function(repeat_file) {
  # Load RepeatMasker output or similar repeat annotation
  # Expected format: chr, start, end, repeat_name, repeat_class, repeat_family
  repeats <- read.table(repeat_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Convert to GRanges
  repeat_gr <- GRanges(
    seqnames = repeats$chr,
    ranges = IRanges(start = repeats$start, end = repeats$end),
    repeat_name = repeats$repeat_name,
    repeat_class = repeats$repeat_class,
    repeat_family = repeats$repeat_family
  )
  
  return(repeat_gr)
}

# Function to load m6A peaks
load_m6a_peaks <- function(peak_file) {
  peaks <- read.table(peak_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(peaks) <- c("chr", "start", "end", "name", "score")
  
  peak_gr <- GRanges(
    seqnames = peaks$chr,
    ranges = IRanges(start = peaks$start, end = peaks$end),
    name = peaks$name,
    score = peaks$score
  )
  
  return(peak_gr)
}

# Function to load bigWig signal files
load_signal_data <- function(bigwig_file) {
  # Load bigWig file containing m6A signal
  signal <- import(bigwig_file, format = "BigWig")
  return(signal)
}

# Function to classify repeat families
classify_repeat_families <- function(repeat_gr) {
  # Create major repeat family classifications
  repeat_gr$major_family <- case_when(
    grepl("LINE", repeat_gr$repeat_class, ignore.case = TRUE) ~ "LINE",
    grepl("SINE", repeat_gr$repeat_class, ignore.case = TRUE) ~ "SINE",
    grepl("LTR|ERV", repeat_gr$repeat_class, ignore.case = TRUE) ~ "LTR/ERV",
    grepl("DNA", repeat_gr$repeat_class, ignore.case = TRUE) ~ "DNA_transposon",
    grepl("Satellite", repeat_gr$repeat_class, ignore.case = TRUE) ~ "Satellite",
    grepl("Simple_repeat", repeat_gr$repeat_class, ignore.case = TRUE) ~ "Simple_repeat",
    grepl("Low_complexity", repeat_gr$repeat_class, ignore.case = TRUE) ~ "Low_complexity",
    TRUE ~ "Other"
  )
  
  return(repeat_gr)
}

# Function to calculate overlap between peaks and repeats
calculate_repeat_overlap <- function(peaks_gr, repeats_gr) {
  # Find overlaps between peaks and repeats
  overlaps <- findOverlaps(peaks_gr, repeats_gr)
  
  # Create summary data frame
  overlap_data <- data.frame(
    peak_idx = queryHits(overlaps),
    repeat_idx = subjectHits(overlaps),
    peak_chr = as.character(seqnames(peaks_gr[queryHits(overlaps)])),
    peak_start = start(peaks_gr[queryHits(overlaps)]),
    peak_end = end(peaks_gr[queryHits(overlaps)]),
    peak_score = peaks_gr$score[queryHits(overlaps)],
    repeat_class = repeats_gr$repeat_class[subjectHits(overlaps)],
    repeat_family = repeats_gr$repeat_family[subjectHits(overlaps)],
    major_family = repeats_gr$major_family[subjectHits(overlaps)]
  )
  
  return(overlap_data)
}

# deepTools-based signal matrix creation with batch processing
create_signal_matrix_deeptools_batch <- function(signal_file, repeats_window, 
                                                upstream = 2000, downstream = 2000, 
                                                bin_size = 50, cores = 8, 
                                                batch_size = 10000000) {
  
  # Validate inputs
  if (!file.exists(signal_file)) {
    stop("Signal file does not exist: ", signal_file)
  }
  
  if (length(repeats_window) == 0) {
    stop("No repeat regions provided")
  }
  
  # Helper function for deepTools processing
  create_signal_matrix_deeptools <- function(signal_file, repeats_window, 
                                           upstream = 2000, downstream = 2000, 
                                           bin_size = 50, cores = 8) {
    
    # Validate batch input
    if (length(repeats_window) == 0) {
      stop("Empty batch of repeats provided")
    }
    
    # Create temporary BED file for repeat regions
    temp_bed <- tempfile(fileext = ".bed")
    
    # Convert GRanges to data frame for BED file with unique identifiers
    unique_ids <- paste0("region_", seq_along(repeats_window))
    bed_data <- data.frame(
      chr = as.character(seqnames(repeats_window)),
      start = pmax(0, start(repeats_window) - 1),  # Convert to 0-based
      end = end(repeats_window),
      name = unique_ids,  # Use unique sequential IDs
      score = 0,
      strand = "+"
    )
    
    # Store the mapping between sequential IDs and original indices
    id_mapping <- data.frame(
      region_id = unique_ids,
      original_index = seq_along(repeats_window),
      stringsAsFactors = FALSE
    )
    
    # Validate coordinates
    invalid_coords <- bed_data$start >= bed_data$end
    if (any(invalid_coords)) {
      cat("Warning: Removing", sum(invalid_coords), "regions with invalid coordinates\n")
      bed_data <- bed_data[!invalid_coords, ]
      repeats_window <- repeats_window[!invalid_coords]
    }
    
    if (nrow(bed_data) == 0) {
      stop("No valid regions after coordinate validation")
    }
    
    write.table(bed_data, temp_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Verify BED file was created and has content
    if (!file.exists(temp_bed) || file.size(temp_bed) == 0) {
      stop("Failed to create valid BED file")
    }
    
    # Create temporary output file
    temp_matrix <- tempfile(fileext = ".gz")
    
    # Build deepTools computeMatrix command with additional options
    cmd <- sprintf(
      "computeMatrix reference-point --referencePoint center -b %d -a %d --binSize %d -p %d -S %s -R %s -o %s --outFileNameMatrix %s.tab",
      upstream, downstream, bin_size, cores, signal_file, temp_bed, temp_matrix, temp_matrix
    )
    
    # Execute deepTools command
    cat("Running deepTools computeMatrix...\n")
    cat("Command:", cmd, "\n")
    system_result <- system(cmd, intern = FALSE)
    
    if (system_result != 0) {
      # Clean up temp files before stopping
      unlink(c(temp_bed, temp_matrix, paste0(temp_matrix, ".tab")))
      stop("deepTools computeMatrix failed with exit code: ", system_result)
    }
    
    # Read the matrix output
    matrix_file <- paste0(temp_matrix, ".tab")
    if (file.exists(matrix_file) && file.size(matrix_file) > 0) {
      if (!requireNamespace("data.table", quietly = TRUE)) {
        signal_matrix <- read.table(matrix_file, skip = 1, header = FALSE, sep = "\t")
      } else {
        library(data.table)
        signal_matrix <- fread(matrix_file, skip = 1, header = FALSE)
      }
      
      # Check if matrix has data
      if (nrow(signal_matrix) == 0) {
        unlink(c(temp_bed, temp_matrix, matrix_file))
        stop("deepTools produced empty matrix - no overlaps found")
      }
      
      # Extract region names from the matrix (4th column in deepTools output)
      returned_regions <- signal_matrix[[4]]
      
      # Remove the first 6 columns (BED info) to get only signal values
      signal_data <- as.matrix(signal_matrix[, 7:ncol(signal_matrix)])
      
      # Create a full matrix with zeros for missing regions
      n_bins <- ncol(signal_data)
      full_matrix <- matrix(0, nrow = length(repeats_window), ncol = n_bins)
      
      # Map returned regions back to their original positions
      for (i in seq_along(returned_regions)) {
        region_name <- returned_regions[i]
        original_idx <- id_mapping$original_index[id_mapping$region_id == region_name]
        if (length(original_idx) > 0) {
          full_matrix[original_idx, ] <- signal_data[i, ]
        }
      }
      
      # Clean up temporary files
      unlink(c(temp_bed, temp_matrix, matrix_file))
      
      return(full_matrix)
    } else {
      unlink(c(temp_bed, temp_matrix, matrix_file))
      stop("deepTools output file not found or empty")
    }
  }
  
  n_repeats <- length(repeats_window)
  n_batches <- ceiling(n_repeats / batch_size)
  
  cat(sprintf("Processing %d repeats in %d batches...\n", n_repeats, n_batches))
  
  all_matrices <- list()
  successful_batches <- 0
  
  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_repeats)
    
    cat(sprintf("Processing batch %d/%d (rows %d-%d)...\n", i, n_batches, start_idx, end_idx))
    
    batch_repeats <- repeats_window[start_idx:end_idx]
    
    tryCatch({
      batch_matrix <- create_signal_matrix_deeptools(signal_file, batch_repeats, 
                                                    upstream, downstream, bin_size, cores)
      all_matrices[[length(all_matrices) + 1]] <- batch_matrix
      successful_batches <- successful_batches + 1
    }, error = function(e) {
      cat("Warning: Batch", i, "failed:", e$message, "\n")
      # Create empty matrix for failed batch to maintain row order
      n_bins <- (upstream + downstream) / bin_size
      empty_matrix <- matrix(0, nrow = length(batch_repeats), ncol = n_bins)
      all_matrices[[length(all_matrices) + 1]] <- empty_matrix
    })
  }
  
  if (successful_batches == 0) {
    stop("All batches failed - check signal file and repeat coordinates")
  }
  
  # Combine all matrices
  final_matrix <- do.call(rbind, all_matrices)
  final_matrix = data.frame(final_matrix)
  # Set rownames for the final matrix using the original repeat names
  # Create unique row identifiers using genomic coordinates
  unique_ids <- paste0(seqnames(repeats_window), ":", start(repeats_window), "-", end(repeats_window))
  rownames(final_matrix) <- unique_ids
  
  # Create metadata from the original repeats_window GRanges object
  metadata <- data.frame(
    repeat_name = repeats_window$repeat_name,  # Keep original repeat_name as a column
    repeat_class = mcols(repeats_window)$repeat_class,
    repeat_family = mcols(repeats_window)$repeat_family,
    major_family = mcols(repeats_window)$major_family,
    repeat_length = width(repeats_window)
  )
  
  # Set rownames for metadata to match the matrix using unique identifiers
  rownames(metadata) <- unique_ids
  

  cat(sprintf("Successfully processed %d/%d batches\n", successful_batches, n_batches))
  
  # Return list with both matrix and metadata
  return(list(matrix = final_matrix, metadata = metadata))
}

# Function to plot signal heatmap by repeat family
plot_repeat_heatmap <- function(signal_data, species_name, output_prefix) {
  signal_matrix <- signal_data$matrix
  metadata <- signal_data$metadata
  
  # Remove rows with all zeros
  non_zero_rows <- rowSums(signal_matrix) > 0
  signal_matrix <- signal_matrix[non_zero_rows, ]
  metadata <- metadata[non_zero_rows, ]
  
  # Order by repeat family
  family_order <- order(metadata$major_family)
  signal_matrix <- signal_matrix[family_order, ]
  metadata <- metadata[family_order, ]
  
  # Create color annotation for repeat families
  family_colors <- RColorBrewer::brewer.pal(length(unique(metadata$major_family)), "Set3")
  names(family_colors) <- unique(metadata$major_family)
  
  # Create annotation
  row_annotation <- rowAnnotation(
    Family = metadata$major_family,
    col = list(Family = family_colors),
    annotation_name_side = "bottom"
  )
  
  # Create heatmap
  ht <- Heatmap(
    signal_matrix,
    name = "m6A Signal",
    col = colorRamp2(c(0, quantile(signal_matrix, 0.95)), c("white", "red")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    left_annotation = row_annotation,
    column_title = paste("m6A Signal on Repeat Elements -", species_name),
    heatmap_legend_param = list(title = "Signal Intensity")
  )
  
  # Save heatmap
  pdf(paste0(output_prefix, "_", species_name, "_repeat_heatmap.pdf"), width = 12, height = 8)
  draw(ht)
  dev.off()
  
  return(ht)
}

# Function to plot average signal profiles by repeat family
plot_signal_profiles <- function(signal_data, species_name, output_prefix, upstream = 2000, downstream = 2000) {
  signal_matrix <- signal_data$matrix
  metadata <- signal_data$metadata
  
  # Calculate average signal for each repeat family
  family_profiles <- list()
  
  for (family in unique(metadata$major_family)) {
    family_idx <- metadata$major_family == family
    if (sum(family_idx) > 0) {
      family_matrix <- signal_matrix[family_idx, , drop = FALSE]
      avg_profile <- colMeans(family_matrix, na.rm = TRUE)
      family_profiles[[family]] <- avg_profile
    }
  }
  
  # Create data frame for plotting
  bin_size <- (upstream + downstream) / ncol(signal_matrix)
  positions <- seq(-upstream, downstream - bin_size, by = bin_size)
  
  profile_df <- data.frame(
    position = rep(positions, length(family_profiles)),
    signal = unlist(family_profiles),
    family = rep(names(family_profiles), each = length(positions))
  )
  
  # Create plot
  p <- ggplot(profile_df, aes(x = position, y = signal, color = family)) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("m6A Signal Profiles on Repeat Elements -", species_name),
      x = "Distance from Repeat Center (bp)",
      y = "Average m6A Signal",
      color = "Repeat Family"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_brewer(type = "qual", palette = "Set2")
  
  # Save plot
  ggsave(paste0(output_prefix, "_", species_name, "_signal_profiles.pdf"), 
         plot = p, width = 10, height = 6)
  
  return(p)
}

# Function to compare repeat enrichment across species
compare_repeat_enrichment <- function(species_data, output_prefix) {
  # Calculate enrichment for each species and repeat family
  enrichment_data <- list()
  
  for (species in names(species_data)) {
    overlap_data <- species_data[[species]]$overlaps
    
    # Calculate total peaks and peaks in each repeat family
    total_peaks <- nrow(overlap_data)
    family_counts <- table(overlap_data$major_family)
    
    enrichment_data[[species]] <- data.frame(
      species = species,
      family = names(family_counts),
      count = as.numeric(family_counts),
      proportion = as.numeric(family_counts) / total_peaks
    )
  }
  
  # Combine data
  enrichment_df <- do.call(rbind, enrichment_data)
  
  # Create comparison plot
  p <- ggplot(enrichment_df, aes(x = family, y = proportion, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "m6A Peak Enrichment in Repeat Families Across Species",
      x = "Repeat Family",
      y = "Proportion of m6A Peaks",
      fill = "Species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_fill_brewer(type = "qual", palette = "Set1")
  
  ggsave(paste0(output_prefix, "_species_comparison.pdf"), 
         plot = p, width = 12, height = 6)
  
  return(p)
}

# Main analysis function
analyze_m6a_repeats <- function(config_file) {
  # Load configuration
  config <- read.table(config_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Expected columns: species, peak_file, signal_file, repeat_file
  
  species_data <- list()
  
  for (i in seq_len(nrow(config))) {
    species <- config$species[i]
    peak_file <- config$peak_file[i]
    signal_file <- config$signal_file[i]
    repeat_file <- config$repeat_file[i]
    
    cat("Processing", species, "...\n")
    
    # Load data
    peaks_gr <- load_m6a_peaks(peak_file)
    repeats_gr <- load_repeat_annotations(repeat_file)
    signal_gr <- load_signal_data(signal_file)
    
    # Classify repeat families
    repeats_gr <- classify_repeat_families(repeats_gr)
    
    # Calculate overlaps
    overlaps <- calculate_repeat_overlap(peaks_gr, repeats_gr)
    
    # Create signal matrix using deepTools batch processing
    signal_data <- create_signal_matrix_deeptools_batch(
      signal_file = signal_file,
      repeats_window = repeats_gr,
      upstream = 2000,
      downstream = 2000,
      bin_size = 50,
      cores = 8,
      batch_size = 10000
    )
    
    # Generate plots
    plot_repeat_heatmap(signal_data, species, "m6a_repeat_analysis")
    plot_signal_profiles(signal_data, species, "m6a_repeat_analysis")
    
    # Store data
    species_data[[species]] <- list(
      peaks = peaks_gr,
      repeats = repeats_gr,
      signal = signal_gr,
      overlaps = overlaps,
      signal_matrix = signal_data
    )

    saveRDS(species_data[[species]], paste0("results/rds/", species, "_m6A_repeat_data.rds"))
  }
  
  # Cross-species comparison
  compare_repeat_enrichment(species_data, "m6a_repeat_analysis")
  
  return(species_data)
}

# Example usage and configuration

# Create example configuration file
config_example <- data.frame(
  species = c(
    "H9_G",
    "H9_I",
    "LG_I",
    "OU_G",
    "OU_I",
    "PM_G",
    "PM_I"
  ),
  peak_file = c(
    "data/peaks/H9_G__m6A_peaks.narrowPeak",
    "data/peaks/H9_I__m6A_peaks.narrowPeak",
    "data/peaks/LG_I__m6A_peaks.narrowPeak",
    "data/peaks/OU_G__m6A_peaks.narrowPeak",
    "data/peaks/OU_I__m6A_peaks.narrowPeak",
    "data/peaks/PM_G__m6A_peaks.narrowPeak",
    "data/peaks/PM_I__m6A_peaks.narrowPeak"
  ),
  signal_file = c(
    'data/bigwig/H9_G__m6A_FE.bw',
    'data/bigwig/H9_I__m6A_FE.bw',
    'data/bigwig/LG_I__m6A_FE.bw',
    'data/bigwig/OU_G__m6A_FE.bw',
    'data/bigwig/OU_I__m6A_FE.bw',
    'data/bigwig/PM_G__m6A_FE.bw',
    'data/bigwig/PM_I__m6A_FE.bw'
    ),
  repeat_file = c(
    "data/repeats/hg38_repeat_clean.tsv",
    "data/repeats/hg38_repeat_clean.tsv",
    "data/repeats/nomLeu3_repeat_clean.tsv",
    "data/repeats/ponAbe3_repeat_clean.tsv",
    "data/repeats/ponAbe3_repeat_clean.tsv",
    "data/repeats/nasLar1_repeat_clean.tsv",
    "data/repeats/nasLar1_repeat_clean.tsv"
  )
)
  
write.table(config_example, "m6a_repeat_config.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
  
# Run analysis
results <- analyze_m6a_repeats("m6a_repeat_config.txt")

saveRDS(results, "species_data_m6A.Rds")
