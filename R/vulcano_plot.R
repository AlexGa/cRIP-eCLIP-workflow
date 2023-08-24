vulcano_plot <- function(df_tab, lfc_string, ylim = c(-10, 100), xlim = c(-2.5, 2.5), main_title = NULL, two_sided = FALSE, lfc_threshold = 1, p_threshold = 0.05, sc_limit = c(0, 30000)){
  
  x <- c('#0F3341', '#1563AA', '#0B99E6', '#3DCDFD', '#F7F7F7', '#EB9457', '#D1551F', '#B02F1B', '#8D1616')
  brick_palette <- structure(colorRampPalette(x)(3000), class = "palette", name = "ocean_brick")
  # BuenColors::jdb_palette("dark_blue")
  constant_sc_gradientn <- ggplot2::scale_colour_gradientn(colours = brick_palette, limits = sc_limit)
  
  df_tab <- df_tab %>% dplyr::rename(log2FC = all_of(lfc_string)) %>% 
    mutate(start = as.numeric(gsub(x = PeakId, pattern = ".*_(.*)_(.*)", replacement = "\\1")),
           p.adj = if_else(p.adj < 10^-100, true = 10^-100, false = p.adj),
           log2FC = if_else(is.infinite(log2FC), true = 10, false = log2FC)) 
  
  signif_peaks <- df_tab %>% filter(p.adj < 0.5, (start < 50 | start > 29800)) %>% mutate(label = gsub(x = PeakId, pattern = "peak.*_(.*)_(.*)", replacement = "[\\1 - \\2]"))
  
  v_plt <- ggplot2::ggplot(df_tab, ggplot2::aes(x = log2FC, y = -log10(p.adj), color = start)) + 
    ggplot2::geom_point(size = 5, alpha = 0.3) + 
    constant_sc_gradientn + 
    ggplot2::geom_point(data = df_tab %>% filter(p.adj < p_threshold),
                        aes(color = start), size = 5, alpha = 1) +
    ggplot2::geom_vline(xintercept = lfc_threshold, lty = 2, alpha = 0.5) + 
    ggplot2::geom_hline(yintercept = -log10(p_threshold), lty = 2, alpha = 0.5) +
    ggplot2::xlab(bquote(bold(log[2]~"fold change"))) + 
    ggplot2::ylab(bquote(bold(-log[10]~"("~adj.~p~value~")"))) + 
    ggplot2::labs(col = "genomic position") + 
    ggplot2::theme_classic() + 
    ggplot2::coord_cartesian(ylim = all_of(ylim), xlim = all_of(xlim)) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                   axis.title = ggplot2::element_text(size = 15, face = "bold", color = "black"),
                   legend.position=c(.15,.9),
                   legend.background = element_rect(fill='transparent'),
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 5)),
                   plot.title = ggplot2::element_text(size = 18, color = "black", face = "bold", hjust = 0.5)) + 
    ggplot2::geom_text(data = signif_peaks, ggplot2::aes(label = label), col = "black", hjust=-0.2, vjust = 0.5)
  # ggrepel::geom_text_repel(data = signif_peaks, 
  #                          ggplot2::aes(label = label), seed = 2,
  #                          size = 4, col = "black", max.overlaps = 50)
  
  if(!is.null(main_title)){
    
    v_plt <- v_plt + ggplot2::ggtitle(main_title)
  }
  
  if(two_sided){
    v_plt <- v_plt + ggplot2::geom_vline(xintercept = -lfc_threshold, lty = 2, alpha = 0.5) + 
      ggplot2::geom_point(data = df_tab %>% filter(p.adj < p_threshold, log2FC < -lfc_threshold), 
                          ggplot2::aes(color = start), size = 5, alpha = 1) + 
      ggplot2::geom_point(data = df_tab %>% filter(p.adj < p_threshold),
                          aes(color = start), size = 5, alpha = 1)
  }
  return(v_plt)
  
}
