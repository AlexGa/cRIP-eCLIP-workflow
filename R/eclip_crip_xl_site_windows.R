library(dplyr)
library(ggplot2)

add_count_vals <- function(bed_df, max_val = 29903, value = 0, chr = "MN908947.3"){
  
  all_pos = data.frame(chr = "MN908947.3", start = (1:max_val) - 1, end = 1:max_val, cov = value)
  missing_start_pos = all_pos[-bed_df$end, ]

  all_sites = rbind(bed_df, missing_start_pos) %>% arrange(start)
  return(all_sites)
}


conv2bed <- function(tsv_df, chr = "MN908947.3"){
  
  new_df <- tsv_df[,c(2, 2, 3)] 
  new_df[,2] <- new_df[,2] + 1
  colnames(new_df) <- c("start", "end", "cov")
  new_df <- new_df %>% mutate(chr = chr) %>% dplyr::relocate(chr, start, end, cov)
  new_df
}

window_sizes <- function(bed_df){
  
  pos <- diff(c(0, add_count_vals(bed_df)$cov))
  win_df <- data.frame(start = which(pos == 1), end = which(pos == -1)) %>% mutate(width = end - start)
  return(win_df)
}


####################################################################################
####################################################################################
##
##
## eCLIP VS cRIP based signifcant XL sites
##
##
####################################################################################
####################################################################################

file_eclip_1 = "/Users/alg22/LRIB/eClip/Analysis/20042023/XL_sites_from_peaks/only_sars_cov2_woFilter/eCLIP_N_rep1-IP-SM/eCLIP_N_rep1_IP_SM_fwd_signif.tsv"
file_eclip_2 = "/Users/alg22/LRIB/eClip/Analysis/20042023/XL_sites_from_peaks/only_sars_cov2_woFilter/eCLIP_N_rep2-IP-SM/eCLIP_N_rep2_IP_SM_fwd_signif.tsv"
file_crip_1 = "/Users/alg22/LRIB/eClip/Analysis/20042023/XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_N_rep1-IP-SM/cRIP_N_rep1_IP_SM_fwd_signif.tsv"
file_crip_2 = "/Users/alg22/LRIB/eClip/Analysis/20042023/XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_N_rep2-IP-SM/cRIP_N_rep2_IP_SM_fwd_signif.tsv"

xl_eclip_1 = conv2bed(readr::read_delim(file_eclip_1, delim ="\t", col_names = T))
xl_eclip_2 = conv2bed(readr::read_delim(file_eclip_2, delim ="\t", col_names = T))
xl_crip_1 = conv2bed(readr::read_delim(file_crip_1, delim ="\t", col_names = T))
xl_crip_2 = conv2bed(readr::read_delim(file_crip_2, delim ="\t", col_names = T))

##### Überlapp eCLIP und cRIP
#
# 
#
#############################

xl_eclip <- xl_eclip_1 %>% dplyr::inner_join(xl_eclip_2, by = c("chr", "start", "end"), suffix = c(".eclip1", ".eclip2") )
xl_crip <- xl_crip_1 %>% dplyr::inner_join(xl_crip_2, by = c("chr", "start", "end"), suffix = c(".crip1", ".crip2") )

xl_crip_eclip <- xl_eclip %>% dplyr::left_join(xl_crip, by = c("chr", "start", "end")) %>% mutate(cov.crip1 = if_else(is.na(cov.crip1), false = cov.crip1, true = 0),
                                                                                                  cov.crip2 = if_else(is.na(cov.crip2), false = cov.crip2, true = 0))

xl_crip_eclip_df <- xl_crip_eclip %>% dplyr::select(matches("cov"))

xl_crip_eclip_df %>% ggplot2::ggplot(ggplot2::aes(cov.crip1, cov.eclip1)) +
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(col = "red", se = F, method = "lm") + 
  ggplot2::theme_bw() + 
  ggplot2::theme()

xl_crip_eclip_df %>% ggplot2::ggplot(ggplot2::aes(cov.eclip1, cov.crip2)) +
  ggplot2::geom_point()

xl_crip_eclip_df %>% ggplot2::ggplot(ggplot2::aes(cov.eclip2, cov.crip1)) +
  ggplot2::geom_point()

xl_crip_eclip_df %>% ggplot2::ggplot(ggplot2::aes(cov.eclip2, cov.crip2)) +
  ggplot2::geom_point()


xl_crip_eclip_df

qqplot(xl_eclip_1$cov, xl_eclip_2$cov)
qqplot(xl_crip_1$cov, xl_crip_2$cov)

colnames(xl_eclip_1) <- colnames(xl_eclip_2) <- colnames(xl_crip_1) <- colnames(xl_crip_2) <- c("chr", "start", "end", "cov")
xl_eclip_1$cov <- xl_eclip_2$cov <- xl_crip_1$cov <- xl_crip_2$cov <- 1

xl_eclip_1_wins <- window_sizes(xl_eclip_1)
xl_eclip_2_wins <- window_sizes(xl_eclip_2)
xl_crip_1_wins <- window_sizes(xl_crip_1)
xl_crip_2_wins <- window_sizes(xl_crip_2)

win_sizes_eclip_crip <- data.frame("experiment" = c(rep("eclip1", times = length(xl_eclip_1_wins$width)), 
                                                    rep("eclip2", times = length(xl_eclip_2_wins$width)), 
                                                    rep("crip1", times = length(xl_crip_1_wins$width)), 
                                                    rep("crip2", times = length(xl_crip_2_wins$width))),
                                   width = c(xl_eclip_1_wins$width, xl_eclip_2_wins$width, xl_crip_1_wins$width, xl_crip_2_wins$width))

win_size_stat <- win_sizes_eclip_crip  %>% filter(width > 0) %>% group_by(experiment) %>% 
  summarise(medianWidth = median(width),
            N = n(),
            minCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[2],
            maxCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[3]) %>% 
  mutate(experiment = factor(experiment, levels = c("eclip1", "eclip2", "crip1", "crip2")))

parwise_wilcox <- win_sizes_eclip_crip %>% rstatix::pairwise_wilcox_test(width~experiment) %>% rstatix::add_xy_position()

win_size_stat %>% ggplot2::ggplot(ggplot2::aes(x = experiment, y = medianWidth)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = minCI, ymax = maxCI), width = 0.5, lwd = 1.1)+
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::ylab("size neighboring XL cluster in nt") + 
  ggplot2::coord_cartesian(ylim = c(1, 5)) + 
  ggplot2::xlab("") + 
  ggplot2::geom_text(ggplot2::aes(x = experiment, y = maxCI + 0.1, label = paste0("N = ", N))) 

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/Median_XL_cluster_signif_replicates.pdf", device = "pdf", dpi = 300, height = 5)

win_sizes_eclip_crip %>% filter(width > 0) %>% ggplot2::ggplot(ggplot2::aes(x = width, col = experiment)) + 
  ggplot2::stat_ecdf(lwd = 1.1) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::xlab("size neighboring XL cluster in nt") + 
  ggplot2::ylab("") + 
  ggsci::scale_color_jama() 

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/ECDF_XLsite_cluster_sizes_replicates_signif.pdf", device = "pdf", dpi = 300, height = 5)


win_sizes_eclip_crip <- data.frame("experiment" = c(rep("eCLIP", times = length(xl_eclip_1_wins$width)), 
                                                    rep("eCLIP", times = length(xl_eclip_2_wins$width)), 
                                                    rep("cRIP", times = length(xl_crip_1_wins$width)), 
                                                    rep("cRIP", times = length(xl_crip_2_wins$width))),
                                   width = c(xl_eclip_1_wins$width, xl_eclip_2_wins$width, xl_crip_1_wins$width, xl_crip_2_wins$width))

win_sizes_eclip_crip %>% filter(width > 0) %>% ggplot2::ggplot(ggplot2::aes(x = width, col = experiment)) + 
  ggplot2::stat_ecdf(lwd = 1.1) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::xlab("size neighboring XL cluster in nt") + 
  ggplot2::ylab("") + 
  ggsci::scale_color_jama() 

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/ECDF_XLsite_cluster_sizes_merge_signif.pdf", device = "pdf", dpi = 300, height = 5)


win_size_stat <- win_sizes_eclip_crip  %>% filter(width > 0) %>% group_by(experiment) %>% 
  summarise(medianWidth = median(width),
            N = n(),
            minCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[2],
            maxCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[3]) %>% 
  mutate(experiment = factor(experiment, levels = c("eCLIP","cRIP")))

p_val <- wilcox.test(win_sizes_eclip_crip  %>% filter(width > 0, experiment == "eCLIP") %>% pull(width), win_sizes_eclip_crip  %>% filter(width > 0, experiment == "cRIP") %>% pull(width), alternative = "less")$p.value

win_size_stat %>% ggplot2::ggplot(ggplot2::aes(x = experiment, y = medianWidth)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = minCI, ymax = maxCI), width = 0.5, lwd = 1.1)+
  ggplot2::geom_point(size = 5) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::ylab("size neighboring XL cluster in nt") + 
  ggplot2::xlab("") + 
  ggplot2::coord_cartesian(ylim = c(0, 3)) + 
  ggplot2::scale_y_continuous(expand = c(0,0)) + 
  ggsignif::geom_signif(comparisons = list(c("eCLIP", "cRIP")), y_position = 2.15, vjust = .1, annotations = paste0("p = ", format(p_val, digits = 2, scientific = T)), lwd = 1.1, textsize = 5)

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/Median_XL_cluster_signif_merge.pdf", device = "pdf", dpi = 300, height = 5)

####################################################################################
####################################################################################
##
##
## eCLIP VS cRIP based all XL sites
##
##
####################################################################################
####################################################################################

file_eclip_1 = "/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/bedgraph_data/IP-eCLIP_N_rep1/IP-eCLIP_N_rep1_fwd.bedgraph"
file_eclip_2 = "/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/bedgraph_data/IP-eCLIP_N_rep2/IP-eCLIP_N_rep2_fwd.bedgraph"
file_crip_1 = "/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/bedgraph_data/IP-cRIP_N_rep1/IP-cRIP_N_rep1_fwd.bedgraph"
file_crip_2 = "/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/bedgraph_data/IP-cRIP_N_rep2/IP-cRIP_N_rep2_fwd.bedgraph"

xl_eclip_1 = readr::read_delim(file_eclip_1, delim ="\t", col_names = F)
xl_eclip_2 = readr::read_delim(file_eclip_2, delim ="\t", col_names = F)
xl_crip_1 = readr::read_delim(file_crip_1, delim ="\t", col_names = F)
xl_crip_2 = readr::read_delim(file_crip_2, delim ="\t", col_names = F)

colnames(xl_eclip_1) <- colnames(xl_eclip_2) <- colnames(xl_crip_1) <- colnames(xl_crip_2) <- c("chr", "start", "end", "cov")
xl_eclip_1$cov <- xl_eclip_2$cov <- xl_crip_1$cov <- xl_crip_2$cov <- 1

xl_eclip_1_wins <- window_sizes(xl_eclip_1)
xl_eclip_2_wins <- window_sizes(xl_eclip_2)
xl_crip_1_wins <- window_sizes(xl_crip_1)
xl_crip_2_wins <- window_sizes(xl_crip_2)

win_sizes_eclip_crip <- data.frame("experiment" = c(rep("eclip1", times = length(xl_eclip_1_wins$width)), 
                                                    rep("eclip2", times = length(xl_eclip_2_wins$width)), 
                                                    rep("crip1", times = length(xl_crip_1_wins$width)), 
                                                    rep("crip2", times = length(xl_crip_2_wins$width))),
                                   width = c(xl_eclip_1_wins$width, xl_eclip_2_wins$width, xl_crip_1_wins$width, xl_crip_2_wins$width))

win_size_stat <- win_sizes_eclip_crip  %>% filter(width > 0) %>% group_by(experiment) %>% 
  summarise(medianWidth = median(width),
            N = n(),
            minCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[2],
            maxCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[3]) %>% 
  mutate(experiment = factor(experiment, levels = c("eclip1", "eclip2", "crip1", "crip2")))

parwise_wilcox <- win_sizes_eclip_crip %>% rstatix::pairwise_wilcox_test(width~experiment) %>% rstatix::add_xy_position()

win_size_stat %>% ggplot2::ggplot(ggplot2::aes(x = experiment, y = medianWidth)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = minCI, ymax = maxCI), width = 0.5, lwd = 1.1)+
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::ylab("size neighboring XL cluster in nt") + 
  ggplot2::xlab("") + 
  ggplot2::coord_cartesian(ylim = c(0, 20))  + 
  ggplot2::geom_text(ggplot2::aes(x = experiment, y = maxCI + 1, label = paste0("N = ", N))) 

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/Median_XL_cluster_raw_replicates.pdf", device = "pdf", dpi = 300, height = 5)

win_sizes_eclip_crip %>% filter(width > 0) %>% ggplot2::ggplot(ggplot2::aes(x = width, col = experiment)) + 
  ggplot2::stat_ecdf(lwd = 1.1) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::xlab("size neighboring XL cluster in nt") + 
  ggplot2::ylab("") + 
  ggsci::scale_color_jama() 

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/ECDF_XLsite_cluster_sizes_replicates_raw.pdf", device = "pdf", dpi = 300, height = 5)

win_sizes_eclip_crip <- data.frame("experiment" = c(rep("eCLIP", times = length(xl_eclip_1_wins$width)), 
                                                    rep("eCLIP", times = length(xl_eclip_2_wins$width)), 
                                                    rep("cRIP", times = length(xl_crip_1_wins$width)), 
                                                    rep("cRIP", times = length(xl_crip_2_wins$width))),
                                   width = c(xl_eclip_1_wins$width, xl_eclip_2_wins$width, xl_crip_1_wins$width, xl_crip_2_wins$width))

win_size_stat <- win_sizes_eclip_crip  %>% filter(width > 0) %>% group_by(experiment) %>% 
  summarise(medianWidth = median(width),
            N = n(),
            minCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[2],
            maxCI = DescTools::MedianCI(width, sides = "two.sided", method = "exact")[3]) %>% 
  mutate(experiment = factor(experiment, levels = c("eCLIP","cRIP")))

p_val <- wilcox.test(win_sizes_eclip_crip  %>% filter(width > 0, experiment == "eCLIP") %>% pull(width), win_sizes_eclip_crip  %>% filter(width > 0, experiment == "cRIP") %>% pull(width), alternative = "less")$p.value

win_size_stat %>% ggplot2::ggplot(ggplot2::aes(x = experiment, y = medianWidth)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = minCI, ymax = maxCI), width = 0.5, lwd = 1.1) +
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::ylab("size neighboring XL cluster in nt") + 
  ggplot2::xlab("") + 
  ggplot2::coord_cartesian(ylim = c(3, 15)) + 
  ggplot2::scale_y_continuous(expand = c(0,0)) + 
  ggsignif::geom_signif(comparisons = list(c("eCLIP", "cRIP")), y_position = 14, vjust = .1, annotations = paste0("p = ", format(p_val, digits = 2, scientific = T)), lwd = 1.1, textsize = 5)

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/Median_XL_cluster_raw_merge.pdf", device = "pdf", dpi = 300, height = 5)

win_sizes_eclip_crip %>% filter(width > 0) %>% ggplot2::ggplot(ggplot2::aes(x = width, col = experiment)) + 
  ggplot2::stat_ecdf(lwd = 1.1) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::xlab("size neighboring XL cluster in nt") + 
  ggplot2::ylab("") + 
  ggsci::scale_color_jama() 

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/ECDF_XLsite_cluster_sizes_merge_raw.pdf", device = "pdf", dpi = 300, height = 5)



colnames(xl_eclip_1) <- colnames(xl_eclip_2) <- colnames(xl_crip_1) <- colnames(xl_crip_2) <- c("chr", "start", "end", "cov")
xl_eclip_1$cov <- xl_eclip_2$cov <- xl_crip_1$cov <- xl_crip_2$cov <- 1

window_sizes <- function(bed_df){
  
  pos <- diff(c(0, add_count_vals(bed_df)$cov))
  win_df <- data.frame(start = which(pos == 1), end = which(pos == -1)) %>% mutate(width = end - start)
  return(win_df)
}

xl_eclip_1_wins <- window_sizes(xl_eclip_1)
xl_eclip_2_wins <- window_sizes(xl_eclip_2)
xl_crip_1_wins <- window_sizes(xl_crip_1)
xl_crip_2_wins <- window_sizes(xl_crip_2)

win_sizes_eclip_crip <- data.frame("experiment" = c(rep("eclip1", times = length(xl_eclip_1_wins$width)), 
      rep("eclip2", times = length(xl_eclip_2_wins$width)), 
      rep("crip1", times = length(xl_crip_1_wins$width)), 
      rep("crip2", times = length(xl_crip_2_wins$width))),
      width = c(xl_eclip_1_wins$width, xl_eclip_2_wins$width, xl_crip_1_wins$width, xl_crip_2_wins$width))

win_size_stat <- win_sizes_eclip_crip  %>% filter(width > 0) %>% group_by(experiment) %>% 
                                                summarise(meanWidth = mean(width),
                                                          N = n(),
                                                          sdWidth = sd(width)/sqrt(N)) %>% 
  mutate(experiment = factor(experiment, levels = c("eclip1", "eclip2", "crip1", "crip2")))

win_size_stat %>% ggplot2::ggplot(ggplot2::aes(x = experiment, y = meanWidth)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = meanWidth - sdWidth, ymax = meanWidth + sdWidth), width = 0.5, lwd = 1.1)+
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::ylab("size XL cluster in nt") + 
  ggplot2::xlab("") + 
  ggplot2::coord_cartesian(ylim = c(10, 50))

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/Mean_SE_XL_cluster_raw_single.pdf", device = "pdf", dpi = 300, height = 5)

win_sizes_eclip_crip <- data.frame("experiment" = c(rep("eCLIP", times = length(xl_eclip_1_wins$width)), 
                                                    rep("eCLIP", times = length(xl_eclip_2_wins$width)), 
                                                    rep("cRIP", times = length(xl_crip_1_wins$width)), 
                                                    rep("cRIP", times = length(xl_crip_2_wins$width))),
                                   width = c(xl_eclip_1_wins$width, xl_eclip_2_wins$width, xl_crip_1_wins$width, xl_crip_2_wins$width))

win_size_stat <- win_sizes_eclip_crip  %>% filter(width > 0) %>% group_by(experiment) %>% 
  summarise(meanWidth = mean(width),
            N = n(),
            sdWidth = sd(width)/sqrt(N)) %>% 
  mutate(experiment = factor(experiment, levels = c("eCLIP","cRIP")))

win_size_stat %>% ggplot2::ggplot(ggplot2::aes(x = experiment, y = meanWidth)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = meanWidth - sdWidth, ymax = meanWidth + sdWidth), width = 0.5, lwd = 1.1)+
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                 title = ggplot2::element_text(size = 14, color = "black", face = "bold")) +
  ggplot2::ylab("size XL cluster in nt") + 
  ggplot2::xlab("") + 
  ggplot2::coord_cartesian(ylim = c(15, 40)) + 
  ggplot2::scale_y_continuous(expand = c(0,0))

ggplot2::ggsave("/Users/alg22/LRIB/eClip/Analysis/cRIP_eCLIP_analysis/xl_comparisons/Mean_SE_XL_cluster_raw_merge.pdf", device = "pdf", dpi = 300, height = 5)
