library(argparse)
library(RCAS)
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-b", "--bed_file", type = "character",
    help="BED file to be processed")
parser$add_argument("-g", "--genome_fasta", type = "character",
    help="The reference genome fasta file.")
parser$add_argument("-o", "--output_prefid", type="character", 
    help="The Prefix of output files")
parser$add_argument("-m", "--motif_width", type="integer", 
    default = 6,
    help="The width of motifs to be discovered. Default = 6")

# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

### helper functions below:
createControlRegions <- function (bed_file, genome_fasta, n_peaks=10000) {
  ## make sure genome_fasta is indexed with .fai file exist.
  genome_fasta_index = paste(genome_fasta, ".fai", sep="")
  if(!file.exists(genome_fasta_index)) {
    message("Indexing genome fasta file...")
    cmd= paste0("samtools faidx ", genome_fasta)
    system(cmd)
  }
  bedtools_cmd = paste0("bedtools random -l 200 -n ",
                        n_peaks, " -g ", genome_fasta_index, 
                        " > cur_control.bed")
  system(bedtools_cmd)
  controlRegions <- RCAS::importBed(filePath = "cur_control.bed", sampleN = 10000, keepStandardChr =F) 
                          #colnames = c('chrom','start', 'end', 'name','width','strand'))
  controlRegions <- sort(controlRegions)
  return(controlRegions)

}


extractSequences <- function (queryRegions, genome) {
  sequences <- Biostrings::getSeq(genome, queryRegions)
  names(sequences) <- paste('seq', 1:length(sequences), sep='_')
  return (sequences)
}

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
    m <- Biostrings::vmatchPattern(pattern = x, subject = seqs, max.mismatch = maxMismatch)
    m <- unlist(m[lengths(m) > 0])
    if(length(m) == 0) {
      return(NULL)
    }
    # remove out of bounds matches
    m <- m[start(m) >= 1]
    m <- m[end(m) <= IRanges::width(seqs[names(m)])]
    # m <- split(m, names(m))
    common <- intersect(names(m), names(seqs))
    if(length(common) == 0) {
      return(NULL)
    }
    # Fix: Convert the extracted sequences to character strings properly
    extracted <- lapply(common, function(x) {
      tmp = Biostrings::extractAt(seqs[x], m[x])
      return(tmp)
    })
    # Convert the extracted DNAStringSet to character strings
    result <- sapply(extracted, function(x) paste(as.character(unlist(x)), collapse=""))
    paste(result, collapse="")
  })
  names(matches) <- patterns
  return(matches)
}

findDifferentialMotifs <- function(querySeqs, 
                                   controlSeqs, 
                                   motifWidth = 6,
                                   motifN = 1, 
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
  
  candidates <- names(which(log2((queryHits + 1) / (controlHits+1)) > 1 ))
  
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


discoverFeatureSpecificMotifs <- function(queryRegions, txdbFeatures, ...) {
  results <- lapply(names(txdbFeatures), function(f) {
    message("Looking for motifs in feature:",f)
    featureCoords <- txdbFeatures[[f]]
    #find query regions that overlap the target features
    q <- queryRegions[unique(queryHits(findOverlaps(queryRegions, featureCoords)))]
    
    if(length(q) > 0) {
      motifResults <- runMotifDiscovery(queryRegions = q, ...)
      return(motifResults)
    }
  })
  names(results) <- names(txdbFeatures)
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

### main functions below:

runMotifDiscovery <- function (bed_file, resizeN = 0, motifWidth = 6,
                        sampleN = 0, genome_fasta, maxMismatch = 1,
                        motifN = 5, nCores = 1) {

  if(sampleN > 0 && length(queryRegions) > sampleN) {
    message("Randomly sampling query regions for motif analysis. 
            Downsampling to ",sampleN," regions")
    queryRegions <- sample(queryRegions, sampleN)
  } 
  
  if(resizeN > 0) {
    resizeIntervals <- width(queryRegions) < resizeN
    if(sum(resizeIntervals) > 0){
      message("Found ",sum(resizeIntervals)," query regions shorter than ", 
              resizeN, " bps. Resizing those regions to ",resizeN,' bps')
      queryRegions[resizeIntervals] <- GenomicRanges::resize(
        x = queryRegions[resizeIntervals], 
        width = 15, 
        fix = 'center')
    }
  }

  genome <- Biostrings::readDNAStringSet(genome_fasta)
  message('extracting sequences from fasta..')
  queryRegions <- RCAS::importBed(filePath = bed_file, sampleN = 10000, keepStandardChr =F,
                          colnames = c('chrom','start', 'end', 'name','score','strand'))
  querySeqs <- extractSequences(queryRegions, genome)

  message('extracting background sequences from fasta..')
  controlRegions <- createControlRegions(bed_file, genome_fasta) 
  controlSeqs <- extractSequences(controlRegions, genome)
  remove(genome) 

  message('running motif discovery ... ')
  motifResults <- findDifferentialMotifs(querySeqs = querySeqs, 
                                         controlSeqs = controlSeqs, 
                                         motifWidth = motifWidth, 
                                         motifN = motifN,
                                         maxMismatch = maxMismatch, 
                                         nCores = nCores)
  return(motifResults)
}


getMotifSummaryTable <- function(motifResults){
  if(!is.null(motifResults)) {
    df <- data.frame('patterns' = names(motifResults$matches_query),
               'queryHits' = colSums(motifResults$counts_query), # total number of motifs in sequences
               'controlHits' = colSums(motifResults$counts_ctrl), 
               'querySeqs' = colSums(motifResults$counts_query > 0), #number of seqs with motifs 
               'controlSeqs' = colSums(motifResults$counts_ctrl > 0),
               'queryFraction' = round(colSums(motifResults$counts_query > 0) / 
                                         nrow(motifResults$counts_query), 2),
               'controlFraction' = round(colSums(motifResults$counts_ctrl > 0) / 
                                           nrow(motifResults$counts_ctrl), 2))
    # for each motif pattern, calculate p-value and odds ratio of the motif
    # occurrence in query sequences compared to the control sequences. 
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
    return(df)
  } else {
    return(data.frame('patterns' = 'NONE',
                      'queryHits' = 0,
                      'controlHits' = 0,
                      'querySeqs' = 0,
                      'controlSeqs' = 0,
                      'queryFraction' = 0,
                      'controlFraction' = 0, 
                      'oddsRatio' = 1,
                      'pvalue' = 1))
  }
}

##
runMotifDiscovery(queryRegions = args$bed_file , genome_fasta = args$genome_fasta,
                  motidN = args$motif_width)