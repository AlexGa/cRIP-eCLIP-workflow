library(dplyr)
library(ggplot2)
library(gggenes)

setwd("~/LRIB/eClip/Analysis/DBA_XL_sites020322/")

xl_sites_fwd <- readr::read_delim("XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h_CTRL2-IP-SM/cRIP_fifth_rep1_NSP9_12h_CTRL2_IP_SM_fwd_raw.tsv")[,-1] %>% filter(oddsRatio > 1, count.ip > 50)
xl_sites_rev <- readr::read_delim("XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h_CTRL2-IP-SM/cRIP_fifth_rep1_NSP9_12h_CTRL2_IP_SM_rev_raw.tsv")[,-1]%>% filter(oddsRatio > 1, count.ip > 50)


xl_sites_fwd <- readr::read_delim("DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fourth_NSP9_12h/cRIP_fourth_NSP9_12h_fwd_raw.tsv")[,-1] %>% 
                dplyr::rename(count.ip = count.ctrl,
                              `log2FC(ip/sm)` = `log2FC(ctrl/ko)`)
xl_sites_rev <- readr::read_delim("DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fourth_NSP9_12h/cRIP_fourth_NSP9_12h_rev_raw.tsv")[,-1] %>% 
                dplyr::rename(count.ip = count.ctrl,
                              `log2FC(ip/sm)` = `log2FC(ctrl/ko)`)

xl_sites_fwd <- xl_sites_fwd %>% mutate(p.adj = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj),
                                        `log2FC(ip/sm)` = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))

xl_sites_rev <- xl_sites_rev %>% mutate(p.adj = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj),
                                        `log2FC(ip/sm)` = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))

anno_scov <- rtracklayer::import("~/LRIB/eClip/references/SARS_CoV2/Sars_cov_2.ASM985889v3.101.gtf")
anno_scov <- anno_scov %>% as.data.frame() %>% filter(type %in% "gene") %>% dplyr::select(seqnames, start, end, strand, gene_name) %>% filter(end != 13483)
anno_scov <- anno_scov %>% dplyr::add_row(seqnames = "MN908947.3", start = c(0,29675), end = c(265, 29903), strand = "+", gene_name = c("5'", "3'")) %>% arrange(start)

anno_scov <- anno_scov %>% dplyr::add_row(seqnames = "MN908947.3", start = anno_scov[c(-1, -2, -12, -13),]$end+1, end = anno_scov[c(-1,-2, -3, -13),]$start-1, strand = "+", gene_name = "intergenic") %>% filter(start < end) %>% arrange(start)

anno_scov$gene_name <- factor(anno_scov$gene_name,
                              levels = c("5'", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", 
                                         "ORF7b", "ORF8", "N", "ORF10", "3'", "intergenic"))

anno_scov_gr <- GenomicRanges::makeGRangesFromDataFrame(anno_scov, keep.extra.columns = T)

xl_sites_fwd_gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(start = xl_sites_fwd$xl_site,
                                                                   end = xl_sites_fwd$xl_site,
                                                                   chr = "MN908947.3"))

xl_sites_rev_gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(start = xl_sites_rev$xl_site,
                                                                      end = xl_sites_rev$xl_site,
                                                                      chr = "MN908947.3"))

ov_fwd <- GenomicRanges::findOverlaps(query = anno_scov_gr, subject = xl_sites_fwd_gr)
ov_rev <- GenomicRanges::findOverlaps(query = anno_scov_gr, subject = xl_sites_rev_gr)

xl_sites_fwd <- xl_sites_fwd %>% mutate(gene_name = "intergenic")
xl_sites_rev <- xl_sites_rev %>% mutate(gene_name = "intergenic")

xl_sites_fwd[ov_fwd@to,]$gene_name <- anno_scov$gene_name[ov_fwd@from]
xl_sites_rev[ov_rev@to,]$gene_name <- anno_scov$gene_name[ov_rev@from]

xl_sites_fwd$gene_name <- factor(xl_sites_fwd$gene_name, levels = levels(anno_scov_gr$gene_name))
xl_sites_rev$gene_name <- factor(xl_sites_rev$gene_name, levels = levels(anno_scov_gr$gene_name))

xl_sites <- xl_sites_fwd %>% cbind(strand = "fwd") %>% 
  dplyr::bind_rows(xl_sites_rev %>% cbind(strand = "fwd"))

seg_colors <- c( "#F6313E", "#00441B", "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#FFA300", "#C390D4" ,"#FF5A00", "#0081C9" ,"#001588", "#490C65",  "#8F1336", "#4B4B4B", "#BA7FD0")
plt_fwd_rev <- ggplot(data = xl_sites) +
  ggplot2::geom_point(data = xl_sites_fwd %>% filter(p.adj < 0.05), aes(x = xl_site, y = -log10(p.adj), 
                                                                        col = factor(gene_name, levels = levels(anno_scov$gene_name))), size = 1, show.legend = F) +
  ggplot2::geom_point(data = xl_sites_rev %>% filter(p.adj < 0.05), aes(x = xl_site, y = log10(p.adj), 
                                                                        col = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                      size = 1, show.legend = F) + 
  ggplot2::geom_point(data = xl_sites_fwd %>% filter(p.adj >= 0.05), aes(x = xl_site, y = -log10(p.adj)), col = "gray", 
                      size = 1, show.legend = F) +
  ggplot2::geom_point(data = xl_sites_rev %>% filter(p.adj >= 0.05), aes(x = xl_site, y = log10(p.adj)), col = "gray", 
                      size = 1, show.legend = F) +
  ggplot2::geom_point(data = xl_sites_fwd %>% filter(gene_name == "N"), aes(x = xl_site, y = -log10(p.adj)), col = "black", pch = 1,
                      size = 1, show.legend = F) +
  ggplot2::geom_hline(yintercept = c(-log10(0.05), log10(0.05)), lty = 2) +
  geom_gene_arrow(data = anno_scov , 
                   aes(xmin = start, xmax = end, fill = gene_name, y = -120), 
                   arrow_body_height = ggplot2::unit(20, "points"), 
                   arrowhead_height = ggplot2::unit(20, "points")) +  
  # geom_gene_arrow(data = anno_scov %>% filter(!grepl(x = gene_name, pattern = "^(3|5)|intergenic")),
  #                 aes(xmin = start, xmax = end, fill = gene_name, y = -120),
  #                 arrow_body_height = ggplot2::unit(20, "points"),
  #                 arrowhead_height = ggplot2::unit(20, "points"), alpha = 0.5) +
  ggplot2::geom_rect(data = anno_scov %>% filter(grepl(x = gene_name, pattern = "^(3|5)")),
                     aes(xmin = start, xmax = end, fill = factor(gene_name, levels = levels(anno_scov$gene_name))), ymin = -128, ymax = -112, 
                     alpha = 1, col = "black", lwd = 0.2) + 
  geom_gene_label(data = anno_scov, aes(xmin = start, xmax = end, label = gene_name, y = -120), 
                  align = "left",min.size = "30pt") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.title = ggplot2::element_blank()) + 
  ggplot2::xlim(c(0, 29903)) +
  ggplot2::ylim(c(-120, 100)) +
  ggplot2::scale_fill_manual(breaks = levels(anno_scov$gene_name), values = alpha(seg_colors, alpha = 0.5)) +
  ggplot2::scale_color_manual(breaks = levels(anno_scov$gene_name), values = alpha(seg_colors, alpha = 1))


ggplot2::ggsave(plot = plt_fwd_rev, filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fourth_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)
  
plt_fwd_rev + 
  ggrepel::geom_label_repel(data = xl_sites_fwd %>% filter(p.adj < 10^-20, gene_name %in% c("N", "ORF7a", "intergenic")),
                      aes(x = xl_site, y = -log10(p.adj), label = xl_site), min.segment.length = 0.1)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fourth_ctrl_ko_with_label.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

adjust_xl_sites_ctrl_ko <- function(xl_sites_fwd_path, xl_sites_rev_path, anno_scov){
  
  xl_sites_fwd <- readr::read_delim(xl_sites_fwd_path)[,-1] %>% dplyr::rename(count.ip = count.ctrl, `log2FC(ip/sm)` = `log2FC(ctrl/ko)`)
  xl_sites_rev <- readr::read_delim(xl_sites_rev_path)[,-1] %>% dplyr::rename(count.ip = count.ctrl, `log2FC(ip/sm)` = `log2FC(ctrl/ko)`)
  
  xl_sites_fwd <- xl_sites_fwd %>% mutate(p.adj = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj),
                                          `log2FC(ip/sm)` = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))
  
  xl_sites_rev <- xl_sites_rev %>% mutate(p.adj = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj),
                                          `log2FC(ip/sm)` = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))
  
  anno_scov_gr <- GenomicRanges::makeGRangesFromDataFrame(anno_scov, keep.extra.columns = T)
  
  xl_sites_fwd_gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(start = xl_sites_fwd$xl_site,
                                                                        end = xl_sites_fwd$xl_site,
                                                                        chr = "MN908947.3"))
  
  xl_sites_rev_gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(start = xl_sites_rev$xl_site,
                                                                        end = xl_sites_rev$xl_site,
                                                                        chr = "MN908947.3"))
  
  ov_fwd <- GenomicRanges::findOverlaps(query = anno_scov_gr, subject = xl_sites_fwd_gr)
  ov_rev <- GenomicRanges::findOverlaps(query = anno_scov_gr, subject = xl_sites_rev_gr)
  
  xl_sites_fwd <- xl_sites_fwd %>% mutate(gene_name = "intergenic")
  xl_sites_rev <- xl_sites_rev %>% mutate(gene_name = "intergenic")
  
  xl_sites_fwd[ov_fwd@to,]$gene_name <- anno_scov$gene_name[ov_fwd@from]
  xl_sites_rev[ov_rev@to,]$gene_name <- anno_scov$gene_name[ov_rev@from]
  
  xl_sites_fwd$gene_name <- factor(xl_sites_fwd$gene_name, levels = levels(anno_scov_gr$gene_name))
  xl_sites_rev$gene_name <- factor(xl_sites_rev$gene_name, levels = levels(anno_scov_gr$gene_name))
  
  xl_sites <- xl_sites_fwd %>% cbind(strand = "fwd") %>% 
    dplyr::bind_rows(xl_sites_rev %>% cbind(strand = "fwd"))
  
  return(list(xl_sites, xl_sites_fwd, xl_sites_rev))
}

plot_manhattan <- function(xl_sites, xl_sites_fwd, xl_sites_rev, anno_scov){
  
  ggplot(data = xl_sites) +
    ggplot2::geom_point(data = xl_sites_fwd %>% filter(p.adj < 0.05), aes(x = xl_site, y = -log10(p.adj), 
                                                                          col = factor(gene_name, levels = levels(anno_scov$gene_name))), size = 1, show.legend = F) +
    ggplot2::geom_point(data = xl_sites_rev %>% filter(p.adj < 0.05), aes(x = xl_site, y = log10(p.adj), 
                                                                          col = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                        size = 1, show.legend = F) + 
    ggplot2::geom_point(data = xl_sites_fwd %>% filter(p.adj >= 0.05), aes(x = xl_site, y = -log10(p.adj)), col = "gray", 
                        size = 1, show.legend = F) +
    ggplot2::geom_point(data = xl_sites_rev %>% filter(p.adj >= 0.05), aes(x = xl_site, y = log10(p.adj)), col = "gray", 
                        size = 1, show.legend = F) +
    ggplot2::geom_hline(yintercept = c(-log10(0.05), log10(0.05)), lty = 2) +
    geom_gene_arrow(data = anno_scov , 
                    aes(xmin = start, xmax = end, fill = gene_name, y = -120), 
                    arrow_body_height = ggplot2::unit(20, "points"), 
                    arrowhead_height = ggplot2::unit(20, "points")) +  
    # geom_gene_arrow(data = anno_scov %>% filter(!grepl(x = gene_name, pattern = "^(3|5)|intergenic")),
    #                 aes(xmin = start, xmax = end, fill = gene_name, y = -120),
    #                 arrow_body_height = ggplot2::unit(20, "points"),
    #                 arrowhead_height = ggplot2::unit(20, "points"), alpha = 0.5) +
    ggplot2::geom_rect(data = anno_scov %>% filter(grepl(x = gene_name, pattern = "^(3|5)")),
                       aes(xmin = start, xmax = end, fill = factor(gene_name, levels = levels(anno_scov$gene_name))), ymin = -128, ymax = -112, 
                       alpha = 1, col = "black", lwd = 0.2) + 
    geom_gene_label(data = anno_scov, aes(xmin = start, xmax = end, label = gene_name, y = -120), 
                    align = "left",min.size = "30pt") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(face = "bold"),
                   axis.text = ggplot2::element_text(size = 12, color = "black")) + 
    ggplot2::xlim(c(0, 29903)) +
    ggplot2::ylim(c(-120, 100)) +
    ggplot2::xlab("genomic location") + 
    ggplot2::ylab(expression(paste(-log[10],'(adjusted p value)'))) + 
    ggplot2::scale_fill_manual(breaks = levels(anno_scov$gene_name), values = alpha(seg_colors, alpha = 0.5)) +
    ggplot2::scale_color_manual(breaks = levels(anno_scov$gene_name), values = alpha(seg_colors, alpha = 1))
}

list_xl_sites <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_fwd_raw.tsv",
                                         xl_sites_rev_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_rev_raw.tsv",
                                         anno_scov = anno_scov)

plot_manhattan(list_xl_sites[[1]], list_xl_sites[[2]], list_xl_sites[[3]], anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fourth_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio > 1), list_xl_sites[[2]]  %>% filter(oddsRatio > 1), list_xl_sites[[3]]  %>% filter(oddsRatio > 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fourth_ctrl_ko_OddsRatio_greater_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio < 1), list_xl_sites[[2]]%>% filter(`log2FC(ip/sm)` < 1), list_xl_sites[[3]] %>% filter(oddsRatio < 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fourth_ctrl_ko_OddsRatio_less_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)


plot_manhattan(xl_sites, xl_sites_fwd, xl_sites_rev, anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fourth_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(xl_sites %>% filter(oddsRatio > 1), xl_sites_fwd  %>% filter(oddsRatio > 1), xl_sites_rev  %>% filter(oddsRatio > 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fourth_ctrl_ko_OddsRatio_greater_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(xl_sites %>% filter(oddsRatio < 1), xl_sites_fwd %>% filter(oddsRatio < 1), xl_sites_rev %>% filter(`oddsRatio` < 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fourth_ctrl_ko_OddsRatio_less_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)



list_xl_sites <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_fwd_raw.tsv",
                        xl_sites_rev_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_rev_raw.tsv",
                        anno_scov = anno_scov)

plot_manhattan(list_xl_sites[[1]], list_xl_sites[[2]], list_xl_sites[[3]], anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fifth_rep1_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio > 1), list_xl_sites[[2]]  %>% filter(oddsRatio > 1), list_xl_sites[[3]]  %>% filter(oddsRatio > 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fifth_rep1_ctrl_ko_OddsRatio_greater_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio < 1), list_xl_sites[[2]]%>% filter(`log2FC(ip/sm)` < 1), list_xl_sites[[3]] %>% filter(oddsRatio < 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fifth_rep1_ctrl_ko_OddsRatio_less_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)


list_xl_sites <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep2_NSP9_12h/cRIP_fifth_rep2_NSP9_12h_fwd_raw.tsv",
                                         xl_sites_rev_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep2_NSP9_12h/cRIP_fifth_rep2_NSP9_12h_rev_raw.tsv",
                                         anno_scov = anno_scov)

plot_manhattan(list_xl_sites[[1]], list_xl_sites[[2]], list_xl_sites[[3]], anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fifth_rep2_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio > 1), list_xl_sites[[2]]  %>% filter(oddsRatio > 1), list_xl_sites[[3]]  %>% filter(oddsRatio > 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fifth_rep2_ctrl_ko_OddsRatio_greater_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio < 1), list_xl_sites[[2]]%>% filter(`log2FC(ip/sm)` < 1), list_xl_sites[[3]] %>% filter(oddsRatio < 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_12h_NSP9_fifth_rep2_ctrl_ko_OddsRatio_less_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)





list_xl_sites <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_8h/cRIP_fifth_rep1_NSP9_8h_fwd_raw.tsv",
                                         xl_sites_rev_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_8h/cRIP_fifth_rep1_NSP9_8h_rev_raw.tsv",
                                         anno_scov = anno_scov)

plot_manhattan(list_xl_sites[[1]], list_xl_sites[[2]], list_xl_sites[[3]], anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fifth_rep1_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio > 1), list_xl_sites[[2]]  %>% filter(oddsRatio > 1), list_xl_sites[[3]]  %>% filter(oddsRatio > 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fifth_rep1_ctrl_ko_OddsRatio_greater_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio < 1), list_xl_sites[[2]]%>% filter(`log2FC(ip/sm)` < 1), list_xl_sites[[3]] %>% filter(oddsRatio < 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fifth_rep1_ctrl_ko_OddsRatio_less_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)


list_xl_sites <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep2_NSP9_8h/cRIP_fifth_rep2_NSP9_8h_fwd_raw.tsv",
                                         xl_sites_rev_path = "DEG_XL_sites/only_sars_cov2_woFilter/cRIP_fifth_rep2_NSP9_8h/cRIP_fifth_rep2_NSP9_8h_rev_raw.tsv",
                                         anno_scov = anno_scov)

plot_manhattan(list_xl_sites[[1]], list_xl_sites[[2]], list_xl_sites[[3]], anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fifth_rep2_ctrl_ko.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio > 1), list_xl_sites[[2]]  %>% filter(oddsRatio > 1), list_xl_sites[[3]]  %>% filter(oddsRatio > 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fifth_rep2_ctrl_ko_OddsRatio_greater_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)

plot_manhattan(list_xl_sites[[1]] %>% filter(oddsRatio < 1), list_xl_sites[[2]]%>% filter(`log2FC(ip/sm)` < 1), list_xl_sites[[3]] %>% filter(oddsRatio < 1), anno_scov)
ggplot2::ggsave(filename = "~/LRIB/eClip/Analysis/DBA_XL_sites020322/plots/XL_sites_8h_NSP9_fifth_rep2_ctrl_ko_OddsRatio_less_1.pdf", device = "pdf",
                dpi = 300, width = 10, height = 5)



names(seg_colors) <- c("5'", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3'", "intergenic")
seg_colors <- c( "#F6313E", "#00441B", "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#FFA300", "#C390D4" ,"#FF5A00", "#0081C9" ,"#001588", "#490C65", "#BA7FD0", "#8F1336", "#4B4B4B")



anno_scov
