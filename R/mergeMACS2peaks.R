library(dplyr)

merge_peaks <- function(peakdf, fdr_threshold  = 0.05){


	peakdf <- peakdf[peakdf$qvalue > -1*log10(fdr_threshold), ]
	peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = peakdf[, 1], 
																start = peakdf[, 2] + 1, 
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

create_gtf_entries <- function(peakdf, feat_type = "gene"){
	
	g_df <- data.frame(chr = as.character(peakdf[, 1]), 
					  start = peakdf[, 2], 
					  end = peakdf[, 3], 
					  strand = ".",
					  type = feat_type,
					  gene_id = paste0("peak_", peakdf[, 2], "_", peakdf[, 3]),
					  transcript_id = NA,
					  exon_id = NA)

	if(feat_type == "transcript"){

		g_df <- g_df %>% mutate(transcript_id = paste0("tx_peak_", peakdf[, 2], "_", peakdf[, 3]))

	}

	if(feat_type == "exon"){

		g_df <- g_df %>% mutate(transcript_id = paste0("tx_peak_", peakdf[, 2], "_", peakdf[, 3]),
							    exon_id = paste0("ex_peak_", peakdf[, 2], "_", peakdf[, 3]))

	}

	g_df
}

# data_dir <- "/vol/projects/agabel/NSP9_new/results/callpeaks_stranded/only_sars_cov2_woFilter"
data_dir <- "~/LRIB/eClip/Analysis/20042023/callpeaks_stranded/only_sars_cov2_woFilter"

out_dir <- "/vol/projects/agabel/NSP9_new/results/united_peaks/callpeaks_stranded/only_sars_cov2_woFilter/"

if(!file.exists(out_dir)){

	dir.create(out_dir, rec = T)
}

sel_dirs_list <- list(ctrl = c("cRIP_fourth_NSP9_8h_CTRL2-IP-SM", "cRIP_fourth_NSP9_12h_CTRL2-IP-SM"),
                			ko6 = c("cRIP_fourth_NSP9_8h_KO6_plus_eV-IP-SM", "cRIP_fourth_NSP9_12h_KO6_plus_eV-IP-SM"))#,
# 				 rescue = c("cRIP_fourth_NSP9_8h_KO6_plus_SND1-IP-SM", "cRIP_fourth_NSP9_12h_KO6_plus_SND1-IP-SM"))

# sel_dirs_list <- list("eCLIP_cnbp-IP-SM")

# sel_dirs_list <- list(snd1 = c("eCLIP_mixSG_snd1_huh7-IP-SM", "eCLIP_mixYW_snd1_huh7-IP-SM"),
					  # cnbp = "eCLIP_cnbp-IP-SM")

sel_dirs_list <- list("8h" = c("cRIP_fourth_NSP9_8h_CTRL2-IP-SM", "cRIP_fourth_NSP9_8h_KO6_plus_eV-IP-SM"),
                      "12h" = c("cRIP_fourth_NSP9_12h_CTRL2-IP-SM", "cRIP_fourth_NSP9_12h_KO6_plus_eV-IP-SM"))#

fdr_threshold <- 0.05

condition_peaks <- lapply(sel_dirs_list, function(sel_dirs){

	narrow_files <- file.path(data_dir, sel_dirs, "fwd_peaks.narrowPeak")

	nf_list <- sapply(narrow_files, function(file_name){

		rep <- gsub(x = gsub(x = file_name, pattern = "(.*)/(.*)_(.*?)\\-(.*)/(.*)", repl = "\\2_\\3"), pattern = "_(8h|fourth)", repl = "")
		df <- readr::read_delim(file_name, col_names =F)
		colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")
		df <- df %>% mutate(name = paste0(rep, start, "_", end))
		df
	}, simplify = F)


	peakdf <- do.call("rbind", nf_list)
	merge_peaks(peakdf)

})

united_peaks <- union_peaks(peakdf = do.call("rbind", condition_peaks))

peakdf <- as.data.frame(united_peaks)


united_peaks_df <- do.call("rbind", sapply(c("gene", "transcript", "exon"), FUN = create_gtf_entries, peakdf = peakdf, simplify = FALSE))


united_peaks_gtf <- GenomicRanges::sort(GenomicRanges::makeGRangesFromDataFrame(united_peaks_df, keep.extra.columns = TRUE))

rtracklayer::export(united_peaks_gtf, file.path(out_dir, "united_peaks_cRIP_8h_12h_fwd.gtf"))
rtracklayer::export(united_peaks, file.path(out_dir, "united_peaks_cRIP_8h_12h_fwd.bed"))


rtracklayer::export(condition_peaks$snd1, file.path(out_dir, "united_peaks_snd1_SG_YW_fwd.bed"))
rtracklayer::export(condition_peaks$cnbp, file.path(out_dir, "united_peaks_cnbp_fwd.bed"))



data_dir <- "~/LRIB/eClip/Analysis/20042023/confirmed_peaks/only_sars_cov2_woFilter"

out_dir <- "/vol/projects/agabel/NSP9_new/results/united_peaks/callpeaks_stranded/only_sars_cov2_woFilter/"

if(!file.exists(out_dir)){
  
  dir.create(out_dir, rec = T)
}


sel_dirs_list <- list("8h" = c("cRIP_fourth_NSP9_8h_CTRL2-IP-SM", "cRIP_fourth_NSP9_8h_KO6_plus_eV-IP-SM"),
                      "12h" = c("cRIP_fourth_NSP9_12h_CTRL2-IP-SM", "cRIP_fourth_NSP9_12h_KO6_plus_eV-IP-SM"))

fdr_threshold <- 0.05

condition_peaks <- lapply(sel_dirs_list, function(sel_dirs){
  
  narrow_files <- file.path(data_dir, sel_dirs, paste0(gsub(x = sel_dirs, pattern = "\\-", replacement = "_"),"_fwd_signif.narrowPeak"))
  
  nf_list <- sapply(narrow_files, function(file_name){
    
    rep <- gsub(x = gsub(x = file_name, pattern = "(.*)/(.*)_(.*?)\\-(.*)/(.*)", repl = "\\2_\\3"), pattern = "_(8|12h|fourth)", repl = "")
    df <- readr::read_delim(file_name, col_names =F)
    colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "pvalue", "qvalue", "peak")
    df <- df %>% mutate(name = paste0(rep, start, "_", end))
    df
  }, simplify = F)
  
  
  peakdf <- do.call("rbind", nf_list)
  merged_peks <- merge_peaks(peakdf)
  merged_peks
})

united_peaks <- union_peaks(peakdf = do.call("rbind", condition_peaks))
