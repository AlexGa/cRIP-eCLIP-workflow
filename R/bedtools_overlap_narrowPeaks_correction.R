library(dplyr)
library(rtracklayer)


merge_peaks <- function(peakdf, fdr_threshold  = 0.05){


	peakdf <- peakdf[peakdf$qvalue > -1*log10(fdr_threshold), ]
	peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = peakdf[, 1], 
																start = peakdf[, 2] + 1, 
																end = peakdf[, 3], 
																strand = peakdf$strand,
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
    
    	# Compare consectutive peaks
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
  	odf <- fP[order(fP$start, decreasing = FALSE),]

  	return(odf)
}

#vol/projects/agabel/
bed_dir_in <- "~/LRIB/hzi_projects/NSP9_new/results/united_peaks/callpeaks_stranded/only_sars_cov2_woFilter/bedtools_overlap_ctrl_ko_25102022"

bed_dir_out <- "~/LRIB/hzi_projects/NSP9_new/results/united_peaks/callpeaks_stranded/only_sars_cov2_woFilter/overlap_ko_ctrl_25102022/merged_peaks"

if(!file.exists(bed_dir_out)){

	dir.create(bed_dir_out)
}

bed_files <- dir(bed_dir_in)

# to be removed from the analysis
bed_dir_in <- "~/LRIB/eClip/Analysis/20042023/peak_intersect_raw/only_sars_cov2_woFilter/"
bed_file_name <- "cRIP_fourth_NSP9_12h-IP-SM_fwd_peaks.bed"
for(bed_file_name in bed_files){

	peak_orient = "rev"
	if(grepl(x = bed_file_name, pattern = "fwd"))
		peak_orient = "fwd"

	bed_file <- readr::read_delim(file.path(bed_dir_in,bed_file_name), col_names = NULL)
	colnames(bed_file) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")
	bed_file <- bed_file %>% mutate(name = paste0(peak_orient,"_", start, "_", end))
	new_bed_file <- merge_peaks(bed_file)

	if(! file.exists(file.path(bed_dir_out, bed_file_name))){

		write.table(new_bed_file, file.path(bed_dir_out, bed_file_name), col.names = F, row.names = F, quote = F, sep ="\t")

	}
}

