library(dplyr)
library(ggplot2)
library(gggenes)

setwd("~/LRIB/eClip/Analysis/20042023/")

segment_colors_pub <- c("#9ab7db", "#7fa18d", "#a2cf9f", "#7fd7cc", "#ffe0bc", "#cbecf4", "#ffd17f", "#c390d4", "#ffac7e", "#7ec0e4", "#7f89c3", "#490d65", "#f693c8","#4B4B4B")

plot_manhattan <- function(xl_sites, anno_scov, x_lim = c(0, 29903), y_lim = c(-120, 100), 
                           x_lab = "genomic location", y_lab = expression(paste(-log[10],'(adjusted p value)')),
                           border_colors = c("#9ab7db", "#7fa18d", "#a2cf9f", "#7fd7cc", "#ffe0bc", "#cbecf4", "#ffd17f", 
                                          "#c390d4", "#ffac7e", "#7ec0e4", "#7f89c3", "#490d65", "#f693c8","#4B4B4B"),
                           seg_colors = c("#1f3674", "#00441B", "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#FFA300", 
                                            "#C390D4" ,"#FF5A00", "#0081C9" ,"#001588", "#490C65",  "#8F1336", "#c9007d", "#4B4B4B"),
                           highlight_ends = TRUE,
                           label_ends = TRUE,
                           hightlight_ntop = 1){
  
  manhattan_plt <- ggplot(data = xl_sites) +
    ggplot2::geom_point(data = xl_sites %>% filter(p.adj < 0.05, strand %in% "fwd"), 
                        aes(x = xl_site, y = -log10(p.adj), 
                            col = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                        size = 1, show.legend = F) +
    ggplot2::geom_point(data = xl_sites %>% filter(p.adj < 0.05, strand %in% "rev"), 
                        aes(x = xl_site, y = log10(p.adj), 
                            col = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                        size = 1, show.legend = F) + 
    ggplot2::geom_point(data = xl_sites %>% filter(p.adj >= 0.05, strand %in% "fwd"), 
                        aes(x = xl_site, y = -log10(p.adj)), col = "gray", 
                        size = 1, show.legend = F) +
    ggplot2::geom_point(data = xl_sites %>% filter(p.adj >= 0.05, strand %in% "rev"), 
                        aes(x = xl_site, y = log10(p.adj)), col = "gray", 
                        size = 1, show.legend = F) +
    ggplot2::geom_hline(yintercept = c(-log10(0.05), log10(0.05)), lty = 2) +
    geom_gene_arrow(data = anno_scov , 
                    aes(xmin = start, xmax = end, fill = gene_name, y = -120), 
                    arrow_body_height = ggplot2::unit(20, "points"), 
                    arrowhead_height = ggplot2::unit(20, "points")) +  
    ggplot2::geom_rect(data = anno_scov %>% filter(grepl(x = gene_name, pattern = "^(3|5)")),
                       aes(xmin = start, xmax = end, fill = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                       ymin = -124, ymax = -116, 
                       alpha = 1, col = "black", lwd = 0.2) + 
    geom_gene_label(data = anno_scov, aes(xmin = start, xmax = end, label = gene_name, y = -120), 
                    align = "left",min.size = "30pt") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(face = "bold"),
                   axis.text = ggplot2::element_text(size = 12, color = "black")) + 
    ggplot2::xlim(x_lim) +
    ggplot2::ylim(y_lim) +
    ggplot2::xlab(x_lab) + 
    ggplot2::ylab(y_lab) + 
    ggplot2::scale_fill_manual(breaks = levels(anno_scov$gene_name), values = seg_colors) + # #alpha(seg_colors, alpha = 0.5)) +
    ggplot2::scale_color_manual(breaks = levels(anno_scov$gene_name), values = border_colors) #alpha(seg_colors, alpha = 1))
  
  
  if(highlight_ends){
    
    highlight_df <- xl_sites %>% filter(p.adj < 0.05, grepl(x = gene_name, pattern = "(3|5)\\'")) %>%
      group_by(strand, gene_name) %>% dplyr::slice_max(order_by = -log10(p.adj), n = hightlight_ntop) %>% ungroup()
    
    manhattan_plt <- manhattan_plt + 
      ggplot2::geom_point(data = highlight_df %>% filter(strand == "fwd"), 
                          aes(x = xl_site, y = -log10(p.adj), 
                              col = factor(gene_name, levels = levels(anno_scov$gene_name)), 
                              fill = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                          size = 4, show.legend = F, pch = 21) +
      ggplot2::geom_point(data = highlight_df %>% filter(strand == "rev"), 
                          aes(x = xl_site, y = log10(p.adj), 
                              col = factor(gene_name, levels = levels(anno_scov$gene_name)), 
                              fill = factor(gene_name, levels = levels(anno_scov$gene_name))), 
                          size = 4, show.legend = F, pch = 21) +
    ggplot2::geom_point(data = highlight_df %>% filter(strand == "fwd"), 
                        aes(x = xl_site, y = -log10(p.adj)), 
                        col = "black", 
                        size = 3, show.legend = F, pch = 1) +
      ggplot2::geom_point(data = highlight_df %>% filter(strand == "rev"), 
                          aes(x = xl_site, y = log10(p.adj)), 
                          col = "black", 
                          size = 3, show.legend = F, pch = 1) 
    
    if(label_ends){
      manhattan_plt <- manhattan_plt + 
        ggrepel::geom_label_repel(data = highlight_df %>% filter(strand == "fwd"), 
                            aes(x = xl_site, y = -log10(p.adj), label = xl_site), 
                            size = 3, show.legend = F, nudge_x = 100,
                            min.segment.length = 0, seed = 42, box.padding = 1.5) +
        ggrepel::geom_label_repel(data = highlight_df %>% filter(strand == "rev"), 
                                  aes(x = xl_site, y = log10(p.adj), label = xl_site), 
                                  size = 3, show.legend = F, nudge_x = 100, nudge_y = -20,
                                  min.segment.length = 0, seed = 42, box.padding = 1.5)
      
    }
                    
  }
  manhattan_plt
}


adjust_xl_sites_ctrl_ko <- function(xl_sites_fwd_path, xl_sites_rev_path, anno_scov, log_column_name = "log2FC(ctrl/ko)"){
  
  xl_sites_fwd <- readr::read_delim(xl_sites_fwd_path)
  xl_sites_rev <- readr::read_delim(xl_sites_rev_path)
  
  anno_scov_gr <- GenomicRanges::makeGRangesFromDataFrame(anno_scov, keep.extra.columns = T)
  
  if(length(xl_sites_fwd) != 0){
    
    # xl_sites_fwd <- xl_sites_fwd[,-1] %>% 
    #   dplyr::rename(count.ip = count.ctrl, `log2FC(ip/sm)` = `log2FC(ctrl/ko)`)
    
    # xl_sites_fwd <- xl_sites_fwd %>% dplyr::rename(count.ip = count.ctrl, `log2FC(ip/sm)` = !!sym(log_column_name))
    
    xl_sites_fwd <- xl_sites_fwd %>% mutate(p.adj = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))#,
                                            # `log2FC(ip/sm)` = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))
    
    xl_sites_fwd_gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(start = xl_sites_fwd$xl_site,
                                                                          end = xl_sites_fwd$xl_site,
                                                                          chr = "MN908947.3"))
    
    ov_fwd <- GenomicRanges::findOverlaps(query = anno_scov_gr, subject = xl_sites_fwd_gr)
    
    xl_sites_fwd <- xl_sites_fwd %>% mutate(gene_name = "intergenic")
    
    xl_sites_fwd[ov_fwd@to,]$gene_name <- anno_scov$gene_name[ov_fwd@from]
    
    xl_sites_fwd$gene_name <- factor(xl_sites_fwd$gene_name, levels = levels(anno_scov_gr$gene_name))
  }
 
  if(length(xl_sites_rev) != 0){
    
    # xl_sites_rev <- xl_sites_rev[,-1] %>% 
    #   dplyr::rename(count.ip = count.ctrl, `log2FC(ip/sm)` = `log2FC(ctrl/ko)`)
    
    # xl_sites_rev <- xl_sites_rev %>% dplyr::rename(count.ip = count.ctrl, `log2FC(ip/sm)` = !!sym(log_column_name))
    
    xl_sites_rev <- xl_sites_rev %>% mutate(p.adj = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))
                                            # `log2FC(ip/sm)` = if_else(p.adj < 1e-100, true = 1e-100, false = p.adj))
    
    xl_sites_rev_gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(start = xl_sites_rev$xl_site,
                                                                          end = xl_sites_rev$xl_site,
                                                                          chr = "MN908947.3"))
    
    ov_rev <- GenomicRanges::findOverlaps(query = anno_scov_gr, subject = xl_sites_rev_gr)
    
    xl_sites_rev <- xl_sites_rev %>% mutate(gene_name = "intergenic")
    
    xl_sites_rev[ov_rev@to,]$gene_name <- anno_scov$gene_name[ov_rev@from]
    
    xl_sites_rev$gene_name <- factor(xl_sites_rev$gene_name, levels = levels(anno_scov_gr$gene_name))
  }
  
  xl_sites <- c()
  if(length(xl_sites_fwd) == 0){
    
    xl_sites <- xl_sites_rev %>% cbind(strand = "rev")
    
  }else if(length(xl_sites_rev) == 0){
   
     xl_sites <- xl_sites_fwd %>% cbind(strand = "fwd")
    
  }else{
    
    xl_sites <- xl_sites_fwd %>% cbind(strand = "fwd") %>% 
      dplyr::bind_rows(xl_sites_rev %>% cbind(strand = "rev"))
  }
    
  
  return(list(xl_sites, xl_sites_fwd, xl_sites_rev))
}


anno_scov <- rtracklayer::import("~/LRIB/eClip/references/SARS_CoV2/Sars_cov_2.ASM985889v3.101.gtf")
anno_scov <- anno_scov %>% as.data.frame() %>%
             filter(type %in% "gene") %>% 
             dplyr::select(seqnames, start, end, strand, gene_name) %>%
             filter(end != 13483)
anno_scov <- anno_scov %>% 
             dplyr::add_row(seqnames = "MN908947.3", start = c(0,29675), end = c(265, 29903), strand = "+", gene_name = c("5'", "3'")) %>% 
             arrange(start)
anno_scov <- anno_scov %>% 
             dplyr::add_row(seqnames = "MN908947.3", start = anno_scov[c(-1, -2, -12, -13),]$end+1, 
                            end = anno_scov[c(-1,-2, -3, -13),]$start-1, strand = "+", gene_name = "intergenic") %>% 
             filter(start < end) %>% 
             arrange(start)
anno_scov$gene_name <- factor(anno_scov$gene_name,
                              levels = c("5'", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", 
                                         "ORF7b", "ORF8", "N", "ORF10", "3'", "intergenic"))


xl_sites_fwd_path <- "DEG_XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_fwd_signif.tsv"
xl_sites_rev_path <- "DEG_XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_rev_signif.tsv"
xl_site_list <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = xl_sites_fwd_path, xl_sites_rev_path = xl_sites_rev_path, anno_scov = anno_scov)
plot_manhattan(xl_sites = xl_site_list[[1]], anno_scov = anno_scov, hightlight_ntop = 1)

setwd("~/LRIB/eClip/Analysis/25052023/")
sub_dir <- "DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/"
comps <- dir(sub_dir)
plot_dir <- "DEG_XL_sites_from_pub_peaks_Manhattan_plots"

if(!file.exists(file.path(plot_dir, "pdf"))){
  dir.create(file.path(plot_dir, "pdf"), recursive = T)  
  dir.create(file.path(plot_dir, "svg"), recursive = T)  
}


for(i in seq_along(comps)){
  
  xl_sites_fwd_path <- file.path(sub_dir, comps[i], paste0(comps[i], "_fwd_signif.tsv"))
  xl_sites_rev_path <- file.path(sub_dir, comps[i], paste0(comps[i], "_rev_signif.tsv"))
  
  xl_site_list <- adjust_xl_sites_ctrl_ko(xl_sites_fwd_path = xl_sites_fwd_path, 
                                          xl_sites_rev_path = xl_sites_rev_path, 
                                          anno_scov = anno_scov,
                                          log_column_name = "log2FC(ko/ctrl)")
  
  for(ntop in 1:5){
    plot_manhattan(xl_sites = xl_site_list[[1]], anno_scov = anno_scov, hightlight_ntop = ntop, seg_colors = segment_colors_pub)
    ggplot2::ggsave(filename = file.path(plot_dir, "pdf", paste0(comps[i],"_TOP_",ntop,"_entries", ".pdf")), device = "pdf",
                    dpi = 300, width = 10, height = 5)
    ggplot2::ggsave(filename = file.path(plot_dir, "svg", paste0(comps[i],"_TOP_",ntop,"_entries", ".svg")), device = "svg",
                    dpi = 300, width = 10, height = 5)
  }
}


