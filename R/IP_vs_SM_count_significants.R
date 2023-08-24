library(dplyr)
library(reshape2)
library(ggplot2)

perform_fisher_test <- function(gr1_peaks, gr2_peaks, gr1_name = "fwd", gr2_name = "rev", alternative_hyp = "two.sided", rm_gr1_zero = FALSE, rm_gr2_zero = FALSE, count_threshold = 10, select_peaks = NULL,
	libSize_df.1 = NULL, libSize_df.2 = NULL, filter_padj = FALSE){

	count_df <- dplyr::inner_join(
									data.frame(reshape2::melt(gr1_peaks, 
															  value.name = "count", 
															  variable.name = "sample")),
									data.frame(reshape2::melt(gr2_peaks, 
															  value.name = "count", 
															  variable.name = "sample")), 
									by = c("PeakId" = "PeakId", "sample" = "sample"), 
									suffix = paste0(".", 1:2)
								  )

	if(!is.null(select_peaks)){

		sample_names <- names(select_peaks)
		count_list <- list()

		for(i in seq_along(sample_names)){

			s_name <- names(select_peaks)[i]
			count_list[[i]] <- count_df %>% filter(sample %in% s_name, PeakId %in% select_peaks[[s_name]])
		}

		count_df <- do.call("rbind", count_list)

	}

	if(!is.null(libSize_df.1) & !is.null(libSize_df.2)){

		libSize_df <- dplyr::inner_join(libSize_df.1, libSize_df.2, 
									by = c("sample" = "sample"), 
									suffix = paste0(".", 1:2)
								  )

		count_df <- count_df %>% inner_join(libSize_df, by = c("sample" = "sample")) %>% 
					group_by(sample) %>%
					mutate(nonPeak.1 = libSize.1 - count.1, 
					   	   nonPeak.2 = libSize.2 - count.2)

	}else{

		count_df <- count_df %>% group_by(sample) %>% 
				mutate(nonPeak.1 = sum(count.1) - count.1, 
					   nonPeak.2 = sum(count.2) - count.2,
					   libSize.1 = sum(count.1), 
					   libSize.2 = sum(count.2))
	}
	

	count_df <- count_df %>% ungroup() %>% rowwise() %>% 
			 	mutate(p.value = fisher.test(matrix(data = c(count.1, count.2, 
			 										  		 nonPeak.1, nonPeak.2), 
			 									 	ncol = 2),
			 								 alternative = alternative_hyp)$p.value,
			 		   norm_count.1 = count.1 * (libSize.2/libSize.1),
			 		   norm_count.2 = count.2,
			 		   oddsRatio = fisher.test(matrix(data = c(count.1, count.2, 
			 										  		   nonPeak.1, nonPeak.2), 
			 									 	  ncol = 2))$estimate,
			 			`log2FC(gr1/gr2)` = log2(norm_count.1/norm_count.2)) %>% 
			 	ungroup() %>% relocate(PeakId, sample, matches("count"), matches("nonPeak"),
			 						   matches("libSize"), matches("norm"), matches("log"), oddsRatio)

	if(alternative_hyp  == "greater"){

		count_df <- count_df %>% filter(oddsRatio > 1)

	}else if(alternative_hyp  == "less"){

		count_df <- count_df %>% filter(oddsRatio < 1)
	}

	if(rm_gr1_zero){

		count_df <- count_df %>% filter(count.1 > count_threshold)
	}

	if(rm_gr2_zero){

		count_df <- count_df %>% filter(count.2 > count_threshold)
	}

	count_df <- cbind(count_df, p.adj = p.adjust(count_df$p.value, method = "BY")) %>% 
				as_tibble()

	if(filter_padj){

		count_df <- count_df %>% filter(p.adj < 0.05)

	}

	colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.1")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.1")], pattern = "\\.1", replacement = paste0(".", gr1_name))
	colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.2")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.2")], pattern = "\\.2", replacement = paste0(".", gr2_name))

	colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr1")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr1")], pattern = "gr1", replacement = paste0("gr.", gr1_name))
	colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr2")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr2")], pattern = "gr2", replacement = paste0("gr.", gr2_name))

	list_peak_stats <- split(count_df, count_df$sample)

	return(list_peak_stats)
}

perform_phyper_fwd_rev_test <- function(gr1_peaks, gr2_peaks, gr1_name = "fwd", gr2_name = "rev", rm_gr1_zero = FALSE, rm_gr2_zero = FALSE, count_threshold = 10, select_peaks = NULL, libSize_df.1 = NULL, libSize_df.2 = NULL, norm_libSizes = TRUE){

	count_df <- dplyr::inner_join(
									data.frame(reshape2::melt(gr1_peaks, 
															  value.name = "count", 
															  variable.name = "sample")),
									data.frame(reshape2::melt(gr2_peaks, 
															  value.name = "count", 
															  variable.name = "sample")), 
									by = c("PeakId" = "PeakId", "sample" = "sample"), 
									suffix = paste0(".", 1:2)
								  )

	if(!is.null(select_peaks)){

		sample_names <- names(select_peaks)
		count_list <- list()

		for(i in seq_along(sample_names)){

			s_name <- names(select_peaks)[i]
			count_list[[i]] <- count_df %>% filter(sample %in% s_name, PeakId %in% select_peaks[[s_name]])
		}

		count_df <- do.call("rbind", count_list)

	}


	if(!is.null(libSize_df.1) & !is.null(libSize_df.2)){

		libSize_df <- dplyr::inner_join(libSize_df.1, libSize_df.2, 
									by = c("sample" = "sample"), 
									suffix = paste0(".", 1:2)
								  )

		count_df <- count_df %>% inner_join(libSize_df, by = c("sample" = "sample"))

	}else{

		count_df <- count_df %>% group_by(sample) %>% 
				mutate(nonPeak.1 = sum(count.1) - count.1, 
					   nonPeak.2 = sum(count.2) - count.2,
					   libSize.1 = sum(count.1), 
					   libSize.2 = sum(count.2))
	}

	count_df <- count_df %>% ungroup() %>% rowwise() %>% 
			 	mutate(Peak_prob = (count.1 + count.2)/(libSize.1 + libSize.2),
			 		   lambda = if_else(norm_libSizes, count.2 * (libSize.1/libSize.2), count.2),
			 		   CPM.1 = if_else(norm_libSizes, count.1/libSize.1 * 1e6, count.1/(libSize.1 + libSize.2)* 1e6),
			 		   CPM.2 = if_else(norm_libSizes, count.2/libSize.2 * 1e6, count.2/(libSize.1 + libSize.2)* 1e6),
			 		   `log2FC(gr1/gr2)` = log2(CPM.1/CPM.2),
			 		   fraction.1 = CPM.1/(CPM.1+CPM.2),
			 		   fraction.2 = CPM.2/(CPM.1+CPM.2),
			 		   # p.value = phyper(q = count.1, 
			 					# 		m = libSize.1,
			 					# 		n = libSize.2, 
			 					# 		k = count.1 + count.2,
			 		   # 					lower.tail = F),
			 		   p.value_poisson = ppois(q = count.1,
			 		   				   		   lambda = lambda,
			 		   				   		   lower.tail = FALSE)) %>%
			 		   # p.value_binomial = pbinom(q = count.2, 
			 		   # 					 size = libSize.2, 
			 		   # 					 prob = count.1/libSize.1,
			 		   # 				   	 lower.tail = FALSE)) %>% 
			 	ungroup() 

	if(rm_gr1_zero){

		count_df <- count_df %>% filter(count.1 > count_threshold)
	}

	if(rm_gr2_zero){

		count_df <- count_df %>% filter(count.2 > count_threshold)
	}

	count_df <- cbind(count_df, p.adj = p.adjust(count_df$p.value_poisson, method = "BY")) %>% 
				as_tibble()

	# count_df <- cbind(count_df, p.adj_binomial = p.adjust(count_df$p.value_binomial, method = "BY")) %>% 
	# 			as_tibble()

	colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.1")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.1")], pattern = "\\.1", replacement = paste0(".", gr1_name))
	colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.2")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "\\.2")], pattern = "\\.2", replacement = paste0(".", gr2_name))

	colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr1")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr1")], pattern = "gr1", replacement = paste0("gr.", gr1_name))
	colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr2")] <- gsub(x = colnames(count_df)[grepl(x = colnames(count_df), pattern = "gr2")], pattern = "gr2", replacement = paste0("gr.", gr2_name))

	list_peak_stats <- split(count_df, count_df$sample)

	return(list_peak_stats)
}

super_res_dir <- "/vol/projects/agabel/NSP9_new/results/single_peaks/callpeaks_stranded/only_sars_cov2_woFilter/ko_vs_ctrl_const_libSizes_20102022"

count_dir <- "/vol/projects/agabel/NSP9_new/results/single_peaks/callpeaks_stranded/only_sars_cov2_woFilter/count_matrices_ko_ctrl"

lib_Size_df <- readr::read_delim("/vol/projects/agabel/NSP9_new/results/dedup/_sars_cov2_woFilter_mapped_fwd_rev.tsv")
colnames(lib_Size_df)[1] <- "sample"


const_libSizes <- TRUE

ip_sm_files <- dir(count_dir, pattern = "IP_SM")

res_dir <- file.path(super_res_dir, "IP_vs_SM")

if(!file.exists(res_dir)){
	dir.create(res_dir, rec = T)
}

ip_sm_list <- list()
ip_sm_list$fwd <- list()
ip_sm_list$rev <- list()

for(i in seq_along(ip_sm_files)){

	count_file <- ip_sm_files[i]

	peaks_df <- readr::read_delim(file.path(count_dir, count_file))

	colnames(peaks_df)[1] <- "PeakId"

	peaks_IP <- peaks_df %>% select(PeakId, starts_with("IP"))
	colnames(peaks_IP) <- gsub(x = colnames(peaks_IP), pattern = "IP-", repl = "")

	peaks_SM <- peaks_df %>% select(PeakId, starts_with("SM"))
	colnames(peaks_SM) <- gsub(x = colnames(peaks_SM), pattern = "SM-", repl = "")

	strand <- unique(gsub(x = count_file, pattern = "(.*)_(fwd|rev)(.*)", replacement = "\\2"))

	
	if(const_libSizes){

		lib_sizes_ip <- lib_Size_df %>% filter(paste0(sample, "_", strand) %in% colnames(peaks_df)[-1],
											grepl(x = sample, pattern = "^IP")) %>% 
					 select(sample, matches(strand)) %>% 
					 mutate(sample = paste0(gsub(x = sample, pattern = "IP-", repl = ""), "_",strand)) %>%
					 rename(libSize = `strand`)

		lib_sizes_sm <- lib_Size_df %>% filter(paste0(sample, "_", strand) %in% colnames(peaks_df)[-1],
											grepl(x = sample, pattern = "^SM")) %>% 
					 select(sample, matches(strand))%>% 
					 mutate(sample = paste0(gsub(x = sample, pattern = "SM-", repl = ""), "_",strand)) %>%
					 rename(libSize = `strand`)

		list_ip_sm_peaks <- perform_fisher_test(gr1_peaks = peaks_IP, 
												gr2_peaks = peaks_SM, 
												gr1_name = "IP", gr2_name = "SM", 
												alternative_hyp = "greater",
												rm_gr1_zero = TRUE,
												libSize_df.1 = lib_sizes_ip,
												libSize_df.2 = lib_sizes_sm,
												filter_padj = TRUE)
	}else{
		list_ip_sm_peaks <- perform_fisher_test(gr1_peaks = peaks_IP, 
											gr2_peaks = peaks_SM, 
											gr1_name = "IP", gr2_name = "SM", 
											alternative_hyp = "greater",
											rm_gr1_zero = TRUE,
											filter_padj = TRUE)
	}
	
	ip_sm_list[[strand]] <- append(ip_sm_list[[strand]], list_ip_sm_peaks)

	names(list_ip_sm_peaks) <- gsub(x = names(list_ip_sm_peaks), pattern = "_(fwd|rev)", replacement = "")

	names(list_ip_sm_peaks) <- gsub(x = names(list_ip_sm_peaks), pattern = "_fourth", replacement = "")

	list_ip_sm_peaks <- list_ip_sm_peaks[order(names(list_ip_sm_peaks))]
	outfile <- paste0("Fisher_Test_", gsub(x = count_file, pattern = ".tsv", replacement = ".xlsx"))

	writexl::write_xlsx(x = list_ip_sm_peaks, 
                    path = file.path(res_dir, outfile))


	gr_list <- lapply(list_ip_sm_peaks, function(mat){

		peak_coord <- do.call("rbind", lapply(strsplit(mat$PeakId, split = "_"), function(elem){c(elem[2:4])}))

		peak_names <- apply(peak_coord, 1, function(elem) paste0("Peak_", paste0(elem[-1], collapse = "_")))
		
		peak_out_df <- data.frame(chromosome = peak_coord[,1], 
								  start = peak_coord[,2],
								  end = peak_coord[,3],
								  name = peak_names,
								  score = ifelse(mat$p.adj > 0, -log10(mat$p.adj), 400),
								  strand = ".",
								  signalValue = ifelse(mat$p.value > 0, -log10(mat$p.value), 400),
								  pValue = ifelse(mat$p.value > 0, -log10(mat$p.value), 400),
								  qValue = ifelse(mat$p.adj > 0, -log10(mat$p.adj), 400),
								  peak = seq_along(peak_coord[,1]))
		# peak_out_gr <- GenomicRanges::makeGRangesFromDataFrame(peak_out_df, keep.extra.columns = T)		
		peak_out_df

		})

	for(i in seq_along(gr_list)){

		sample_name <- names(gr_list)[i]
		write.table(x = gr_list[[i]], file = file.path(res_dir, paste0(sample_name, "_", strand, "_p_adjusted.narrowPeak")), sep ="\t", dec = ".", row.names = F, col.names = F, quote = F)
	}

	gr_list <- lapply(list_ip_sm_peaks, function(mat){

		peak_coord <- do.call("rbind", lapply(strsplit(mat$PeakId, split = "_"), function(elem){c(elem[2:4])}))

		peak_names <- apply(peak_coord, 1, function(elem) paste0("Peak_", paste0(elem[-1], collapse = "_")))
		
		peak_out_df <- data.frame(chromosome = peak_coord[,1], 
								  start = peak_coord[,2],
								  end = peak_coord[,3],
								  name = peak_names,
								  score = mat$`log2FC(gr.IP/gr.SM)`,
								  strand = ".",
								  signalvalue = mat$`log2FC(gr.IP/gr.SM)`,
								  pValue = ifelse(mat$p.value > 0, -log10(mat$p.value), 400),
								  qValue = ifelse(mat$p.adj > 0, -log10(mat$p.adj), 400),
								  peak = seq_along(peak_coord[,1]))

		# peak_out_gr <- GenomicRanges::makeGRangesFromDataFrame(peak_out_df, keep.extra.columns = T)		
		peak_out_df

		})

	for(i in seq_along(gr_list)){

		sample_name <- names(gr_list)[i]
		write.table(x = gr_list[[i]], file = file.path(res_dir, paste0(sample_name,"_", strand, "_log2FC.narrowPeak")), sep ="\t", dec = ".", row.names = F, col.names = F, quote = F)
	}

	gr_list <- lapply(list_ip_sm_peaks, function(mat){

		peak_coord <- do.call("rbind", lapply(strsplit(mat$PeakId, split = "_"), function(elem){c(elem[2:4])}))

		peak_names <- apply(peak_coord, 1, function(elem) paste0("Peak_", paste0(elem[-1], collapse = "_")))
		
		peak_out_df <- data.frame(chromosome = peak_coord[,1], 
								  start = peak_coord[,2],
								  end = peak_coord[,3],
								  name = peak_names,
								  score = mat$oddsRatio,
								  strand = ".",
								  signalvalue = mat$oddsRatio,
								  pValue = ifelse(mat$p.value > 0, -log10(mat$p.value), 400),
								  qValue = ifelse(mat$p.adj > 0, -log10(mat$p.adj), 400),
								  peak = seq_along(peak_coord[,1]))
		# peak_out_gr <- GenomicRanges::makeGRangesFromDataFrame(peak_out_df, keep.extra.columns = T)		
		peak_out_df

		})

	for(i in seq_along(gr_list)){

		sample_name <- names(gr_list)[i]
		write.table(x = gr_list[[i]], file = file.path(res_dir, paste0(sample_name, "_", strand, "_oddsRatio.narrowPeak")), sep ="\t", dec = ".", row.names = F, col.names = F, quote = F)
	}
}
