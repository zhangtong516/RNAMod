library(RCAS)
library(ggpubr)

generateKmers <- function(k, letters = c("A", "C", "G", "T")) {
  kmer <- c()
  for(i in 1:k){
    kmer <- unlist(lapply(letters, function(x){paste(kmer, x, sep="")}))
  }
  return(kmer)
}

countPattern <- function(seqs, patterns, maxMismatch = 0, nCores = 1) {
  
  cl <- parallel::makeCluster(nCores)
  parallel::clusterExport(cl = cl, varlist = c('seqs', 'patterns', 'maxMismatch'), 
                          envir = environment())
  M <- do.call(rbind, pbapply::pblapply(cl = cl, X = patterns, 
                                        FUN = function(x) {
                                          Biostrings::vcountPattern(x, seqs, 
                                                                    max.mismatch = maxMismatch)
                                        }))
  parallel::stopCluster(cl)
  rownames(M) <- patterns
  return(t(M))
}

extractMatches <- function(seqs, patterns, minMismatch = 0, maxMismatch = 0) {
  matches <- lapply(patterns, function(x) {
    # Get matches - this returns different S4 classes on different platforms
    m_raw <- Biostrings::vmatchPattern(pattern = x, subject = seqs, max.mismatch = maxMismatch)
    
    # Convert to consistent IRangesList format regardless of platform
    if(class(m_raw)[1] == "ByPos_MIndex") {
      # Linux: convert ByPos_MIndex to IRangesList
      m_list <- as(m_raw, "IRangesList")
    } else {
      # Windows: already IRangesList or similar
      m_list <- m_raw
    }
    
    # Filter out empty matches
    m_list <- m_list[lengths(m_list) > 0]
    
    if(length(m_list) == 0) {
      return(NULL)
    }
    
    # Convert to flat IRanges with sequence names
    m_flat <- unlist(m_list)
    
    if(length(m_flat) == 0) {
      return(NULL)
    }
    
    # Get sequence names for each match
    seq_names <- names(m_flat)
    
    # Remove out of bounds matches
    valid_indices <- start(m_flat) >= 1
    
    # Check end bounds
    if(length(seq_names) > 0) {
      seq_widths <- Biostrings::width(seqs[unique(seq_names)])
      names(seq_widths) <- names(seqs[unique(seq_names)])
      
      # Match widths to each match position
      match_widths <- seq_widths[seq_names]
      valid_indices <- valid_indices & (end(m_flat) <= match_widths)
    }
    
    # Filter matches
    m_flat <- m_flat[valid_indices]
    seq_names <- seq_names[valid_indices]
    
    if(length(m_flat) == 0) {
      return(NULL)
    }
    
    # Group back by sequence
    m_grouped <- IRanges::IRangesList(split(m_flat, seq_names)) 
    
    # Get common sequences
    common <- intersect(names(m_grouped), names(seqs))
    if(length(common) == 0) {
      return(NULL)
    }
    
    # Extract sequences
    extracted <- try({
      # Extract the sequences at the match positions
      extracted_seqs <- Biostrings::extractAt(seqs[common], m_grouped[common])
      
      # Convert to character vector
      char_seqs <- as.character(unlist(extracted_seqs))
      
      return(char_seqs)
    }, silent = TRUE)
    
    if(inherits(extracted, "try-error")) {
      return(NULL)
    }
    
    return(extracted)
  })
  names(matches) <- patterns
  return(matches)
}

findDifferentialMotifs <- function(querySeqs, 
                                   controlSeqs, 
                                   motifWidth = 6,
                                   motifN = 5, 
                                   nCores = 1, 
                                   maxMismatch = 1) {
  
  kmers <- generateKmers(k = motifWidth)
  
  selected <- seq(1:length(querySeqs))
  if(length(querySeqs) > 1000) {
    selected <- sample(1:length(querySeqs), 1000)
  }
  
  query <- countPattern(querySeqs[selected], kmers, maxMismatch, nCores)
  ctrl <- countPattern(controlSeqs[selected], kmers, maxMismatch, nCores)
  
  # only keep patterns that are more frequent in the query
  queryHits <- apply(query, 2, function(x) sum(x > 0))
  controlHits <- apply(ctrl, 2, function(x) sum(x > 0))
  
  candidates <- names(which(log2((queryHits + 1) / (controlHits+1)) > 0))
  
  if(length(candidates) == 0) {
    warning("Couldn't find any motifs that occur more often in the query sequences
            compared to the background. Returning NULL. ")
    return(NULL)
  }
  
  query <- query[,candidates]
  ctrl <- ctrl[,candidates]
  
  # train a random forest model 
  df <- as.data.frame(rbind(query[,candidates], ctrl[,candidates]))
  df$label <- as.factor(c(rep("query", nrow(query)), 
                          rep("ctrl", nrow(ctrl))))
  
  fit <- ranger::ranger(label ~ ., df, importance = 'impurity')
  var.imp <- sort(ranger::importance(fit), decreasing = T)
  #get top variables 
  max <- ifelse(motifN < length(candidates), motifN, length(candidates))
  top <- names(var.imp[1:max])
  
  results <- list("counts_query" = query[,top, drop = F], 
                  "counts_ctrl" = ctrl[,top,drop = F], 
                  "matches_query" = extractMatches(seqs = querySeqs, 
                                                   patterns = top, 
                                                   maxMismatch = maxMismatch), 
                  "matches_ctrl" = extractMatches(seqs = controlSeqs, 
                                                  patterns = top, 
                                                  maxMismatch = maxMismatch))
  return(results)
}

getPWM <- function(sequences, letters = c('A', 'C', 'G', 'T')) {
  chars <- strsplit(sequences, '')
  if(length(unique(lengths(chars))) > 1) {
    stop("getPWM: All sequences must be of the same width.")
  }
  m <- do.call(rbind, chars) # matrix of chars 
  # for each position, count occurrences 
  counts <- apply(m, 2, function(y) {
    sapply(letters, function(l) {
      sum(y == l)
    })
  })
  # convert counts to fractions
  pwm <- apply(counts, 2, function(x) x / sum(x))
  return(pwm)
}

##### actual run 

run_motif = function(query_fasta, control_fasta,
                     motifWidth = 6, 
                     motifN = 5,
                     maxMismatch = 1, 
                     nCores = 2) {
  querySeqs = Biostrings::readDNAStringSet(query_fasta)
  controlSeqs = Biostrings::readDNAStringSet(control_fasta)
  
  motifResults <- findDifferentialMotifs(querySeqs = querySeqs, 
                                         controlSeqs = controlSeqs, 
                                         motifWidth = motifWidth, 
                                         motifN = motifN,
                                         maxMismatch = maxMismatch, 
                                         nCores = nCores)
  
  if(is.null(motifResults)) {
    return(NULL)
  }
  
  df <- data.frame('patterns' = names(motifResults$matches_query),
                   'queryHits' = colSums(motifResults$counts_query), # total number of motifs in sequences
                   'controlHits' = colSums(motifResults$counts_ctrl), 
                   'querySeqs' = colSums(motifResults$counts_query > 0), #number of seqs with motifs 
                   'controlSeqs' = colSums(motifResults$counts_ctrl > 0),
                   'queryFraction' = round(colSums(motifResults$counts_query > 0) / 
                                             nrow(motifResults$counts_query), 2),
                   'controlFraction' = round(colSums(motifResults$counts_ctrl > 0) / 
                                               nrow(motifResults$counts_ctrl), 2))
  
  
  df <- cbind(df, do.call(rbind, apply(df, 1, function(r) {
    querySeqCount <- as.numeric(r[4])
    ctrlSeqCount <- as.numeric(r[5])
    t <- fisher.test(matrix(data = c(querySeqCount, ctrlSeqCount, 
                                     nrow(motifResults$counts_query) - querySeqCount, 
                                     nrow(motifResults$counts_ctrl) - ctrlSeqCount), 
                            nrow = 2), alternative = 'greater')
    return(data.frame('oddsRatio' = round(t$estimate[[1]], 2), 
                      'pvalue' = t$p.value))
  })))
  return(list("result"= motifResults,
              "summary"= df))
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript annotate_peaks.R <query_fasta> <control_fasta> <sample_id>")
}

query_fasta <- args[1]
control_fasta <- args[2]
sample_id <- args[3]

motifs = run_motif(query_fasta, control_fasta)

if(is.null(motifs)) {
  cat("No differential motifs found. Exiting.\n")
  quit(status = 0)
}

motif_summary = motifs[['summary']]
data.table::fwrite(motif_summary, paste0(sample_id, "_motif_summary.csv"))

# Create plots with error handling
create_logo_plot <- function(matches, title) {
  if(is.null(matches) || length(matches) == 0) {
    return(ggplot() + ggtitle(paste(title, "(No matches found)")) + theme_void())
  }
  tryCatch({
    ggseqlogo::ggseqlogo(matches) + theme_bw() + ggtitle(title)
  }, error = function(e) {
    ggplot() + ggtitle(paste(title, "(Error creating logo)")) + theme_void()
  })
}

# Create plots for up to 5 motifs
plot_list <- list()
for(i in 1:min(5, length(motifs$result$matches_query))) {
  plot_list[[i]] <- create_logo_plot(
    motifs$result$matches_query[[i]], 
    names(motifs$result$matches_query)[i]
  )
}

# Fill remaining slots if less than 5 motifs
while(length(plot_list) < 5) {
  plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
}

# Create summary table
p6 <- ggtexttable(motif_summary[,c("patterns", "oddsRatio", "pvalue")])

all_plots <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                      plot_list[[4]], plot_list[[5]], p6, 
                      ncol=2, nrow=3)

pdf(paste0(sample_id, "_motifs_logoPlot.pdf"), height=6, width=8)
print(all_plots)
dev.off()

cat("Analysis completed successfully!\n")