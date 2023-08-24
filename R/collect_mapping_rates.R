organisms <- dir(".")[-2]
list_all <- list()

i <- 1
for(org in organisms){
	
	print(org)
	exps <- dir(org)

		for(exp_dir in exps){

			if(file.exists(file.path(org, exp_dir, "Log.final.out"))){
				stat_df <- read.delim(file.path(org, exp_dir, "Log.final.out"), sep ="|", header = F)

			stat_df[,2] <- as.numeric(gsub(stat_df[,2], pattern = "\t|%", repl=""))

			x <- stat_df[c(5, 8, 23, 25, 28, 30, 32),2]

			list_all[[i]] <- data.frame(organism = org, sample = exp_dir, input = x[1], 
									unique = x[2], multi_mapped = x[3], 
									unmapped_multi = x[4], unmapped_mismatch = x[5], 
									too_short = x[6], unmapped_other = x[7])
			i <<- i + 1
			}
			
		}
}

df_all <- do.call("rbind", list_all)


write.table(df_all, "../../mapping_rates_separate_newPipe.tsv", col.names = T, row.names = F, quote = F)



mapping_rates <- readr::read_table("~/LRIB/hzi_projects/NSP9_new/mapping_rates_all_old_and_new_060822.tsv")

# mapping_rates <- mapping_rates %>% filter(!sample %in% c("IP-eCLIP_second_NSP9_KO6_RESCUE", "SM-eCLIP_first_NSP9_KO6"))

input_reads <- mapping_rates %>% filter(organism == "repbase") %>% select(organism, sample, input)
mapped_reads <- mapping_rates %>% mutate(mapped = unique + multi_mapped,
                                         unmapped = unmapped_multi + unmapped_mismatch + too_short + unmapped_other) %>% select(organism, sample, input, mapped, unmapped)

mapped_reads_df <- rbind(mapped_reads, mapped_reads %>% filter(organism == "homo_sapiens_sars_cov2") %>% 
                 mutate(input = input,
                        mapped = unmapped,
                        unmapped= 0,
                        organism = "unknown")) %>%
                  mutate(treatment = tolower(gsub(x = gsub(x = sample, pattern = "(IP|SM)-(.*?)_(.*)", replacement = "\\3"), pattern = "(first|second)_", replacement = "")),
                         exp = factor(gsub(x = gsub(x = sample, pattern = "(IP|SM)-(.*?)_(.*)", replacement = "\\2"), pattern = "eRIP", replacement = "dRIP"), levels = c("dRIP", "eCLIP")),
                         type = factor(gsub(x = sample, pattern = "(IP|SM)-(.*?)_(.*)", replacement = "\\1"), levels = c("IP",  "SM")),
                         treatment = if_else(treatment == "nsp9_ko6", true = "nsp9_ko6_ev", false = treatment),
                         replicate = factor(case_when(grepl(x = sample, pattern = "first") ~ "rep1",
                                               grepl(x = sample, pattern = "eRIP") ~ "rep2",
                                               grepl(x = sample, pattern = "second") ~ "rep2",
                                               TRUE ~ "rep3",
                         ), levels = paste0("rep", 1:3)))

mapped_reads_df$sample <- factor(mapped_reads_df$sample, levels = unique(mapped_reads_df %>% arrange(treatment, type, exp) %>% pull(sample)))

plot_mapping_rate <- function(reads_df, y = "mapped", method = "abs", ylab = "#Mapped Reads", fill_type = "organism"){
  
    plt <- reads_df %>% ggplot2::ggplot(ggplot2::aes_string(x = "replicate", y = y, fill = fill_type))
    
    if(method == "abs"){
     
      plt <- plt + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge2()) +
        ggplot2::facet_wrap(~treatment + type)+
        ggsci::scale_fill_aaas() +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                       legend.position = "bottom") +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, title = NULL)) + 
        ggplot2::xlab("Replicate") + 
        ggplot2::ylab(ylab)
    }
    if(method == "rel"){
      
      plt <- plt + ggplot2::geom_bar(stat = "identity") +
        ggplot2::facet_wrap(~treatment + type)+
        ggsci::scale_fill_aaas() +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                       legend.position = "bottom") +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, title = NULL)) + 
        ggplot2::xlab("Replicate") + 
        ggplot2::ylab(ylab)
    }
   
    plt
}

plot_dedup_rate <- function(reads_df, y = "dedup_rate", ylab = "duplication rate", fill_type = "type"){
  
  plt <- reads_df %>% ggplot2::ggplot(ggplot2::aes_string(x = "replicate", y = y, fill = fill_type)) +  
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge2()) +
      ggplot2::facet_wrap(~treatment + type) +
      ggsci::scale_fill_aaas() +
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                     legend.position = "bottom") +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, title = NULL)) + 
      ggplot2::xlab("Replicate") + 
      ggplot2::ylab(ylab)
  
  plt
}


plt_dir <- "~/LRIB/eClip/Analysis/Mapping_rates"
dir.create(plt_dir)

plot_mapping_rate(reads_df = mapped_reads_df %>% filter(exp %in% "dRIP")) + ggplot2::ggtitle("dRIP - Samples")
ggplot2::ggsave(filename = file.path(plt_dir, "dRIP_mapping_absolute.pdf"), device = "pdf", height =10, width = 10)
plot_mapping_rate(mapped_reads_df %>% filter(exp %in% "eCLIP")) + ggplot2::ggtitle("eCLIP - Samples")
ggplot2::ggsave(filename = file.path(plt_dir, "eCLIP_mapping_absolute.pdf"), device = "pdf", height =10, width = 10)


plot_mapping_rate(mapped_reads_df %>% group_by(sample) %>% mutate(rel_mapped = mapped / sum(mapped)) %>% filter(exp %in% "dRIP"), 
                  y = "rel_mapped", method = "rel", ylab = "Relative # mapped reads") + ggplot2::ggtitle("dRIP - Samples")
ggplot2::ggsave(filename = file.path(plt_dir, "dRIP_mapping_relative.pdf"), device = "pdf", height =10, width = 10)
plot_mapping_rate(mapped_reads_df %>% group_by(sample) %>% mutate(rel_mapped = mapped / sum(mapped)) %>% filter(exp %in% "eCLIP"), 
                  y = "rel_mapped", method = "rel", ylab = "Relative # mapped reads") + ggplot2::ggtitle("eCLIP - Samples")
ggplot2::ggsave(filename = file.path(plt_dir, "eCLIP_mapping_relative.pdf"), device = "pdf", height =10, width = 10)



dedup_old_new <- read.table("~/LRIB/hzi_projects/NSP9_new/results/stats/dedup/homo_sapiens_sars_cov2/dedup_rate_old_and_new.csv", sep = ";", header =F)
colnames(dedup_old_new) <- c("sample", "dedup_rate")

dedup_old_new <- dedup_old_new %>% mutate(treatment = tolower(gsub(x = gsub(x = sample, pattern = "(IP|SM)-(.*?)_(.*)", replacement = "\\3"), pattern = "(first|second)_", replacement = "")),
         exp = factor(gsub(x = gsub(x = sample, pattern = "(IP|SM)-(.*?)_(.*)", replacement = "\\2"), pattern = "eRIP", replacement = "dRIP"), levels = c("dRIP", "eCLIP")),
         type = factor(gsub(x = sample, pattern = "(IP|SM)-(.*?)_(.*)", replacement = "\\1"), levels = c("IP",  "SM")),
         treatment = if_else(treatment == "nsp9_ko6", true = "nsp9_ko6_ev", false = treatment),
         replicate = case_when(grepl(x = sample, pattern = "first") ~ "rep1",
                               grepl(x = sample, pattern = "eRIP") ~ "rep2",
                               grepl(x = sample, pattern = "second") ~ "rep2",
                               TRUE ~ "rep3",
         )) %>% unique()

plot_dedup_rate(dedup_old_new %>% filter(exp %in% "dRIP"), 
                  y = "dedup_rate", ylab = "duplication rate", fill_type = "replicate") + ggplot2::ggtitle("dRIP - Samples")
ggplot2::ggsave(filename = file.path(plt_dir, "dRIP_duplication_rate.pdf"), device = "pdf", height =10, width = 10)
plot_dedup_rate(dedup_old_new %>% filter(exp %in% "eCLIP"), 
                y = "dedup_rate", ylab = "duplication rate", fill_type = "type") + ggplot2::ggtitle("eCLIP - Samples")
ggplot2::ggsave(filename = file.path(plt_dir, "eCLIP_duplication_rate.pdf"), device = "pdf", height =10, width = 10)
