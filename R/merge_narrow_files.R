library(dplyr)

union_peaks <- function(peakdf){
  
  peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = peakdf[, 1], 
                                                              start = peakdf[, 2], 
                                                              end = peakdf[, 3], 
                                                              strand = ".", keep.extra.columns = TRUE))
  
  peaks <- sort(peaks)
  
  while (!(GenomicRanges::isDisjoint(peaks))) {
    
    # Fast variable access
    chr_names <- as.character(GenomicRanges::seqnames(peaks))
    starts <- GenomicRanges::start(peaks)
    ends <- GenomicRanges::end(peaks)
    scores <- GenomicRanges::mcols(peaks)$score
    
    # See if consecutive peaks are overlapping
    overlap_next <- intersect(
      which(chr_names[1:(length(chr_names) - 1)] == chr_names[2:(length(chr_names))]), 
      which(ends[1:(length(ends) - 1)] >= starts[2:(length(starts))] )
    )
    
    # Compare consectuive peaks
    overlap_previous <- overlap_next + 1
    
    # overlap_comparison <- scores[keep_peaks[overlap_previous]] > scores[keep_peaks[overlap_next]]
    # discard <- keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
    
    starts[overlap_next][starts[overlap_next] >= starts[overlap_previous]] <- starts[overlap_previous][starts[overlap_previous] <= starts[overlap_next]]
    starts[overlap_previous][starts[overlap_previous] >= starts[overlap_next]] <- starts[overlap_next][starts[overlap_next] <= starts[overlap_previous]]
    
    ends[overlap_next][ends[overlap_next] <= ends[overlap_previous]] <- ends[overlap_previous][ends[overlap_previous] >= ends[overlap_next]]
    ends[overlap_previous][ends[overlap_previous] <= ends[overlap_next]] <- ends[overlap_next][ends[overlap_next] >= ends[overlap_previous]]
    
    
    new_peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chr_names, 
                                                                    start = starts,
                                                                    end = ends,
                                                                    strand = ".", keep.extra.columns = TRUE))
    peaks <- new_peaks[-which(duplicated(new_peaks))]
  }
  peaks <- GenomeInfoDb::sortSeqlevels(peaks)
  peaks <- sort(peaks)
  
  return(peaks)
}

merge_peaks <- function(peakdf, fdr_threshold  = 0.05){
  
  
  if(max(peakdf$qvalue) <= 1){

    peakdf <- peakdf %>% filter(qvalue < fdr_threshold) %>% 
                         mutate(qvalue = if_else(qvalue != 0, true = -log10(qvalue), false = 400))
  
  }else{
  
    peakdf <- peakdf %>% filter(qvalue > -log10(fdr_threshold))
  
  }

  peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = peakdf[, 1], 
                                                              start = peakdf[, 2], 
                                                              end = peakdf[, 3], 
                                                              strand = ".",
                                                              score = peakdf$qvalue), keep.extra.columns = TRUE)
  
  peaks <- sort(peaks)
  keep_peaks <- 1:length(peaks)
  
  while (!(GenomicRanges::isDisjoint(peaks[keep_peaks]))) {
    
    # Fast variable access
    chr_names <- as.character(GenomicRanges::seqnames(peaks[keep_peaks]))
    starts <- GenomicRanges::start(peaks[keep_peaks])
    ends <- GenomicRanges::end(peaks[keep_peaks])
    scores <- GenomicRanges::mcols(peaks)$score
    
    # See if consecutive peaks are overlapping
    overlap_next <- intersect(
      which(chr_names[1:(length(keep_peaks) - 1)] == chr_names[2:(length(keep_peaks))]), 
      which(ends[1:(length(keep_peaks) - 1)] >= starts[2:(length(keep_peaks))] )
    )
    
    # Compare consectuive peaks
    overlap_previous <- overlap_next + 1
    overlap_comparison <- scores[keep_peaks[overlap_previous]] > scores[keep_peaks[overlap_next]]
    discard <- keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
    keep_peaks <- keep_peaks[!keep_peaks %in% discard]
  }
  peaks <- GenomeInfoDb::sortSeqlevels(peaks)
  peaks <- sort(peaks)
  
  # Export the final result by making a data frame; getting the top (or as many) n peaks
  # based on the score and then resort based on genomic position.
  fP <- data.frame(peaks[keep_peaks], rank = 1:length(keep_peaks))
  nout <- min(Inf, dim(fP)[1])
  odf <- head(fP[order(fP$score, decreasing = TRUE),], nout)
  
  return(odf)
}

narrow_files <- list(snakemake@input[[1]], snakemake@input[[2]])
output_file <- snakemake@output[[1]]


if(!file.exists(dirname(output_file))){
  
  dir.create(out_dir, rec = T)
}

nf_list <- sapply(narrow_files, function(file_name){
    
    df <- readr::read_delim(file_name, col_names =F)
    colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "pvalue", "qvalue", "peak")
    df <- df %>% mutate(name = paste0(chr, start, "_", end))
    df
  }, simplify = F)


peakdf <- do.call("rbind", nf_list)
merged_peaks <- merge_peaks(peakdf) %>% dplyr::arrange(start)

rtracklayer::export(merged_peaks, output_file)
