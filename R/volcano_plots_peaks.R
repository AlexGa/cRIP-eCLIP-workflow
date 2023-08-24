library(dplyr)
library(reshape2)
library(ggplot2)
library(BuenColors)
# devtools::install_github("caleblareau/BuenColors")

ocean_brick = c('#0F3341', '#1563AA', '#0B99E6', '#3DCDFD', '#F7F7F7', '#EB9457', '#D1551F', '#B02F1B', '#8D1616')

print.palette <- function(x, ...) {
  n <- length(x)
  old <- par(mar = c(0.5, 0.5, 0.5, 0.5))
  on.exit(par(old))
  
  image(1:n, 1, as.matrix(1:n), col = x,
        ylab = "", xaxt = "n", yaxt = "n")
  
  rect(0, 0.9, n + 1, 1.1, col = rgb(1, 1, 1, 0.8), border = NA)
  text((n + 1) / 2, 1, labels = attr(x, "name"), cex = 1, family = "serif")
}

jdb_palette <- function(name, n, type = c("discrete", "continuous")) {
  type <- match.arg(type)
  
  pal <- jdb_palettes[[name]]
  
  if (is.null(pal)) stop("Palette not found.")
  if (missing(n)) n <- length(pal)
  
  if (type == "discrete" && n > length(pal)) {
    stop("Number of requested colors greater than what palette can offer")
  }
  
  out <- switch(type,
                continuous = colorRampPalette(pal)(1000),
                discrete = pal[1:n]
  )
  structure(out, class = "palette", name = name)
}

vulcano_plot <- function(df_tab, lfc_string, ylim = c(-10, 100), xlim = c(-2.5, 2.5), main_title = NULL, 
                         two_sided = FALSE, lfc_threshold = 1, p_threshold = 0.05, sc_limit = c(0, 30000)){
  
  x <- c('#0F3341', '#1563AA', '#0B99E6', '#3DCDFD', '#F7F7F7', '#EB9457', '#D1551F', '#B02F1B', '#8D1616')
  brick_palette <- structure(colorRampPalette(x)(3000), class = "palette", name = "ocean_brick")
  constant_sc_gradientn <- ggplot2::scale_colour_gradientn(colours = brick_palette, limits = sc_limit)
  
  df_tab <- df_tab %>% dplyr::rename(log2FC = all_of(lfc_string)) %>% 
    mutate(start = as.numeric(gsub(x = peak_name, pattern = ".*_(.*)_(.*)", replacement = "\\1")),
           p.adj = if_else(p.adj < 10^-100, true = 10^-100, false = p.adj),
           log2FC = if_else(is.infinite(log2FC), true = 10, false = log2FC)) 
  
  signif_peaks <- df_tab %>% filter(p.adj < 0.5, (start < 250 | start > 29800)) %>% mutate(label = gsub(x = peak_name, pattern = "peak.*_(.*)_(.*)", replacement = "[\\1 - \\2]"))
  
  v_plt <- ggplot2::ggplot(df_tab, ggplot2::aes(x = log2FC, y = -log10(p.adj), color = start)) + 
    ggplot2::geom_point(size = 5, alpha = 0.3) + 
    constant_sc_gradientn + 
    ggplot2::geom_point(data = df_tab %>% filter(p.adj < p_threshold),
                        aes(color = start), size = 5, alpha = 1) +
    ggplot2::geom_vline(xintercept = lfc_threshold, lty = 2, alpha = 0.5) + 
    ggplot2::geom_hline(yintercept = -log10(p_threshold), lty = 2, alpha = 0.5) +
    ggplot2::xlab(bquote(bold(log[2]~"fold change"))) + 
    ggplot2::ylab(bquote(-Log[10]~"("~adjusted~italic(P)~value~")")) + 
    ggplot2::labs(col = "genomic position") + 
    ggplot2::theme_classic() + 
    ggplot2::coord_cartesian(ylim = all_of(ylim), xlim = all_of(xlim)) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14, color = "black"),
                   axis.title = ggplot2::element_text(size = 15, face = "bold", color = "black"),
                   legend.position=c(.15,.9),
                   legend.background = element_rect(fill='transparent'),
                   legend.text = ggplot2::element_text(size = 13),
                   legend.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 5)),
                   plot.title = ggplot2::element_text(size = 18, color = "black", face = "bold", hjust = 0.5)) 
  
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
  
  v_plt <- v_plt + 
    # ggplot2::geom_text(data = signif_peaks, ggplot2::aes(label = label), col = "black", hjust=-0.2, vjust = 0.5)
    ggplot2::geom_point(data = signif_peaks, ggplot2::aes(x = log2FC, y = -log10(p.adj), color = start), 
                        col = "red", size = 7, pch = 1, alpha = 1) +
    ggrepel::geom_text_repel(data = signif_peaks, alpha = 1,
                             ggplot2::aes(label = label), seed = 2, 
                             size = 4, col = "black", #max.overlaps = 50, 
                             fontface = 'bold', 
                             box.padding = unit(1, "lines"),
                             point.padding = unit(0.5, "lines"),
                             segment.color = 'black')
  return(v_plt)
  
}


setwd("~/LRIB/eClip/Analysis/25052023/")

plot_dir <- "~/LRIB/Labmeeting/Data_meeting_15032023/plots"

if(!file.exists(plot_dir)){
  
  dir.create(plot_dir)  
}
new_dir <- "XL_sites_strand_from_peaks_corrected"

if(!file.exists(new_dir)){
  
  dir.create(new_dir)
}

sub_dirs <- dir("XL_sites_strand_from_peaks//")

for(i in seq_along(sub_dirs)){
  
  if(grepl(pattern = "plots|NSP12", x = sub_dirs[i])){
    next
  }
  for(strand in c("fwd", "rev")){
    
    file <- dir(file.path("XL_sites_strand_from_peaks", sub_dirs[i]), pattern = paste0(strand, "_XL_sites_raw.tsv"))
    if(length(file) == 0) {next;}
  # df_tab <- readr::read_delim(file.path("XL_sites_all", sub_dirs[i], paste0(sub_dirs[i],"_",strand,"_IP_SM_raw.tsv")))[,-1]
  df_tab <- readr::read_delim(file.path("XL_sites_strand_from_peaks", sub_dirs[i], file))[,-1]
  # df_tab <- df_tab %>% mutate(`log2FC(ctrl/ko)` = log2((count.ctrl/(count.ctrl + noPeak_ctrl))/(count.ko/(count.ko + noPeak_ko))))
  if(strand == "fwd"){

    df_tab <- df_tab %>% mutate(`log2FC(fwd/rev)` = log2((count_fwd/(libSize.fwd + libSize.rev))/(count_rev/(libSize.fwd + libSize.rev))))
  }
  if(strand == "rev"){

    df_tab <- df_tab %>% mutate(`log2FC(fwd/rev)` = log2((count_fwd/(libSize.fwd + libSize.rev))/(count_rev/(libSize.fwd + libSize.rev))))
  }
  # df_tab <- df_tab %>% mutate(`log2FC(ctrl/ko)` = log2((count.ctrl/(count.ctrl + no_xl_ctrl))/(count.ko/(count.ko + no_xl_ko))))
  # df_tab <- df_tab %>% mutate(`log2FC(ctrl/ko)` = log2((count.ctrl/(count.ctrl + no_xl_ctrl))/(count.ko/(count.ko + no_xl_ko))))
  # df_tab <- df_tab %>% mutate(`log2FC(ip/sm)` = log2((count.ip/(count.ip + no_xl_ip))/(count.sm/(count.sm + no_xl_sm))))
  
  if(!file.exists(file.path(new_dir, sub_dirs[i]))){
    
    dir.create(file.path(new_dir, sub_dirs[i]))
  }
  
  write.table(x = df_tab, file = file.path(new_dir, sub_dirs[i], paste0(sub_dirs[i],"_",strand,"_XL_sites_raw.tsv")), 
              row.names = F, quote = F, dec = ".", sep ="\t")
  
  write.table(x = df_tab %>% filter(`p.adj(BY)` < 0.05), 
              file = file.path(new_dir, sub_dirs[i], paste0(sub_dirs[i],"_",strand,"_XL_sites_signif.tsv")), 
              row.names = F, quote = F, dec = ".", sep ="\t")
  }
}

dba_dir <- "DBA_pub"
sub_dirs <- dir(file.path(dba_dir, "only_sars_cov2_woFilter"))
plot_dir <- paste0(dba_dir,"_volcano")

if(!file.exists(file.path(plot_dir, "pdf"))){
  dir.create(file.path(plot_dir, "pdf"), recursive = T)
  dir.create(file.path(plot_dir, "svg"), recursive = T)
}
for(i in seq_along(sub_dirs)){
  
  df_tab <- readr::read_delim(file.path(dba_dir, "only_sars_cov2_woFilter", sub_dirs[i], paste0(sub_dirs[i],"_rev_raw.tsv")))
  main_title <- gsub(x = sub_dirs[i], pattern = "_", repl = " ") 
  vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", 
               two_sided = T, xlim = c(-2.5, 2.5)) +
    ggplot2::xlab(bquote(Log[2]~"fold change"~"(SND1 KO/CTRL)"))
  ggsave(filename = file.path(plot_dir, "svg", paste0(sub_dirs[i], "_rev.svg")), device = "svg", width = 6, height = 7, dpi = 300)
  ggsave(filename = file.path(plot_dir, "pdf", paste0(sub_dirs[i], "_rev.pdf")), device = "pdf", width = 6, height = 7, dpi = 300)
  
  df_tab <- readr::read_delim(file.path(dba_dir, "only_sars_cov2_woFilter", sub_dirs[i], paste0(sub_dirs[i],"_fwd_raw.tsv")))
  main_title <- gsub(x = sub_dirs[i], pattern = "_", repl = " ") 
  vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", 
               two_sided = T, xlim = c(-2.5, 2.5)) +
    ggplot2::xlab(bquote(Log[2]~"fold change"~"(SND1 KO/CTRL)"))
  ggsave(filename = file.path(plot_dir, "svg", paste0(sub_dirs[i], "_fwd.svg")), device = "svg", width = 6, height = 7, dpi = 300)
  ggsave(filename = file.path(plot_dir, "pdf", paste0(sub_dirs[i], "_fwd.pdf")), device = "pdf", width = 6, height = 7, dpi = 300)
}

peak_tab <- readr::read_delim(file = "confirmed_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_12h_CTRL2-IP-SM/cRIP_fourth_NSP9_12h_CTRL2_IP_SM_fwd_raw.tsv")

peak_tab <- peak_tab %>% rowwise() %>% mutate(pval_new = fisher.test(matrix(data = c(count.IP, count.SM, 
                                                          noPeak_IP, noPeak_SM), 
                                                 ncol = 2),
                                          alternative = "greater")$p.value)
peak_tab <- cbind(peak_tab, "p.adj_new" = p.adjust(peak_tab$p.value, method = "BY"))


main_title <- "cRIP rep0 NSP9 12h (fwd)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_fourth_NSP9_12h_fwd_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)


df_tab <- read_delim("DBA/cRIP_fourth_NSP9_12h/cRIP_fourth_NSP9_12h_rev_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep0 NSP9 12h (rev)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_fourth_NSP9_12h_rev_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)

df_tab <- read_delim("DBA/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_fwd_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep0 NSP9 8h (fwd)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_fourth_NSP9_8h_fwd_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)


df_tab <- read_delim("DBA/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_rev_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep0 NSP9 8h (rev)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_fourth_NSP9_8h_rev_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)



df_tab <- read_delim("DBA/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_fwd_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep1 NSP9 12h (fwd)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep1_NSP9_12h_fwd_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)

df_tab <- read_delim("DBA/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_rev_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep1 NSP9 12h (rev)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep1_NSP9_12h_rev_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)


df_tab <- read_delim("DBA/cRIP_fifth_rep1_NSP9_8h/cRIP_fifth_rep1_NSP9_8h_fwd_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep1 NSP9 8h (fwd)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep1_NSP9_8h_fwd_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)

df_tab <- read_delim("DBA/cRIP_fifth_rep1_NSP9_8h/cRIP_fifth_rep1_NSP9_8h_rev_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep1 NSP9 8h (rev)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep1_NSP9_8h_rev_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)


df_tab <- read_delim("DBA/cRIP_fifth_rep2_NSP9_12h/cRIP_fifth_rep2_NSP9_12h_fwd_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep2 NSP9 12h (fwd)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep2_NSP9_12h_fwd_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)

df_tab <- read_delim("DBA/cRIP_fifth_rep2_NSP9_12h/cRIP_fifth_rep2_NSP9_12h_rev_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep2 NSP9 12h (rev)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep2_NSP9_12h_rev_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)


df_tab <- read_delim("DBA/cRIP_fifth_rep2_NSP9_8h/cRIP_fifth_rep2_NSP9_8h_fwd_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep2 NSP9 8h (fwd)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep2_NSP9_8h_fwd_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)

df_tab <- read_delim("DBA/cRIP_fifth_rep2_NSP9_8h/cRIP_fifth_rep2_NSP9_8h_rev_raw.tsv")
df_tab <- df_tab %>% mutate(`log2FC(ko/ctrl)` = log2((count.ko/(count.ko + noPeak_ko))/(count.ctrl/(count.ctrl + noPeak_ctrl))))

main_title <- "cRIP rep2 NSP9 8h (rev)"
vulcano_plot(df_tab = df_tab, main_title = main_title, lfc_string = "log2FC(ko/ctrl)", two_sided = T, xlim = c(-5, 5)) +
  ggplot2::xlab(bquote(bold(log[2]~"fold change"~"[ko/ctrl]")))
ggsave(filename = file.path(plot_dir, "cRIP_rep2_NSP9_8h_rev_raw.pdf"), device = "pdf", width = 6, height = 7, dpi = 300)

