library(dplyr)
library(GenomicRanges)
norm_df <- function(bw_cont, N){
  
  extend_df <- function(bw_ranges){
    
    bw_df <- bw_ranges %>% as.data.frame()
    j <- 0
    ext_coord_list <- apply(bw_df, 1, function(row_elem){
      
      j <<- j + 1
      data.frame(id = j, start = seq(from = as.numeric(row_elem[2]),to = as.numeric(row_elem[3]), by = 1), 
                 score = as.numeric(row_elem[6]))
    })
    
    ext_df <- do.call(rbind, ext_coord_list)
    ext_df
  }
  
  tmp_df <- extend_df(bw_cont)
  tmp_df$score <- as.numeric(tmp_df$score)/N * 1e6
  
  return(tmp_df)
}

norm_fw_rev <- function(bw_cont_fwd, bw_cont_rev, N_counts_fwd_rev){
  
  extend_df <- function(bw_ranges){
    
    bw_df <- bw_ranges %>% as.data.frame()
    j <- 0
    ext_coord_list <- apply(bw_df, 1, function(row_elem){
      
      j <<- j + 1
      data.frame(id = j, start = seq(from = as.numeric(row_elem[2]),to = as.numeric(row_elem[3]), by = 1), 
                 score = as.numeric(row_elem[6]))
    })
    
    ext_df <- do.call(rbind, ext_coord_list)
    ext_df
  }
  
  df_fwd <- extend_df(bw_cont_fwd)
  df_rev <- extend_df(bw_cont_rev)
  
  df_fwd$score <- as.numeric(df_fwd$score)/N_counts_fwd_rev * 1e6
  df_rev$score <- as.numeric(df_rev$score)/N_counts_fwd_rev * 1e6
  
  return(list(fwd = df_fwd, rev = df_rev))
}

df_to_bw <- function(cov_df, seq_length = 29903, chr_name = "MN908947.3"){
  
  cov_bw <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chromosome = chr_name,
                                                               start = cov_df$start,
                                                               end = cov_df$start,
                                                               strand = "*",
                                                               score = cov_df$score), keep.extra.columns = T)
  seqlengths(cov_bw) <- seq_length
  cov_bw
}


bw_dir <- "~/LRIB/eClip/Analysis/25102022/only_sars_cov2_woFilter_bw/"

libsizes <- readr::read_delim("~/LRIB/eClip/Analysis/25102022/_sars_cov2_woFilter_mapped_fwd_rev.tsv")
colnames(libsizes)[1] <- "id"

libsizes %>% filter(grepl(x = id, pattern = "^IP-(cRIP_fourth|eCLIP_(mix|N|cnbp)|SND1_A)")) %>% arrange(id) %>% rowwise() %>% mutate(N = fwd + rev)

libsizes <- libsizes %>% tidyr::pivot_longer(cols = c(fwd, rev), names_to = "strand", values_to = "libsize") %>% mutate(sample = paste0(id, "_", strand))

bigwig_files <- sort(dir(bw_dir, include.dirs = F, pattern = "*.bw$")) 
bigwig_files <- bigwig_files %>% as_tibble() %>% dplyr::rename(file_name = value) %>% filter(grepl(x = file_name, pattern = "^IP-(cRIP_fourth|eCLIP_(mix|N|cnbp)|SND1_A).*(fwd|rev).bw")) %>% as.data.frame()

bigwig_files <- bigwig_files %>% mutate(sample = gsub(x = file_name, pattern = ".bw", repl = ""))

bw_files <- bigwig_files %>% inner_join(libsizes, by = c("sample" = "sample")) %>% relocate(id, sample, file_name, libsize)

# Generate RPM tracks over both strands
bw_out_dir <- "~/LRIB/eClip/Analysis/25102022/RPM_norm"
if(!file.exists(bw_out_dir)){
  dir.create(bw_out_dir, recursive = T)
}

lapply(split(bw_files, bw_files$id), function(bw_strand_pair){
  
  cond <- unique(bw_strand_pair$id)
  
  fwd_file <- bw_strand_pair %>% filter(grepl(x = sample, pattern = "fwd")) %>% pull(file_name)
  rev_file <-  bw_strand_pair %>% filter(grepl(x = sample, pattern = "rev")) %>% pull(file_name)
  
  bw_fwd <- rtracklayer::import(file.path(bw_dir, fwd_file))
  bw_rev<- rtracklayer::import(file.path(bw_dir, rev_file))
  
  N_fwd_rev <- sum(bw_strand_pair %>% pull(libsize))
  rpm_fwd_rev <- norm_fw_rev(bw_cont_fwd = bw_fwd, bw_cont_rev = bw_rev, N_counts_fwd_rev = N_fwd_rev)

  fwd_rev_scaled_bw_list <- lapply(rpm_fwd_rev, function(df){
      df_to_bw(cov_df = df)
  })

  for(i in seq_along(fwd_rev_scaled_bw_list)){
      
      strand <- names(fwd_rev_scaled_bw_list)[i]
      bw_out_file_name <- paste0(cond, "_", strand,"_rpm_norm_fwd_rev.bw")
      
      rtracklayer::export.bw(object = df_to_bw(cov_df = fwd_rev_scaled_bw_list[[i]] %>% as.data.frame()), file.path(bw_out_dir, bw_out_file_name))
  }
  
  
})


bw_out_dir <- "~/LRIB/eClip/Analysis/25102022/RPM_norm_per_strand"
if(!file.exists(bw_out_dir)){
  dir.create(bw_out_dir)
}

lapply(split(bw_files, bw_files$sample), function(bw_strand_pair){

  bw_strand <- rtracklayer::import(file.path(bw_dir, bw_strand_pair$file_name))
  
  norm_df_strand <- norm_df(bw_cont = bw_strand, N = bw_strand_pair$libsize)
  
  file_name <- paste0(bw_strand_pair$sample,"_rpm_norm.bw")
  rtracklayer::export.bw(object = df_to_bw(cov_df = norm_df_strand %>% as.data.frame()), file.path(bw_out_dir, file_name))
  
})




bw_fwd_ctrl_12 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_12h_CTRL2_fwd.bw"))
bw_rev_ctrl_12 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_12h_CTRL2_rev.bw"))

bw_fwd_ko_12 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_12h_KO6_plus_eV_fwd.bw"))
bw_rev_ko_12 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_12h_KO6_plus_eV_rev.bw"))

bw_fwd_ctrl_8 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_8h_CTRL2_fwd.bw"))
bw_rev_ctrl_8 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_8h_CTRL2_rev.bw"))

bw_fwd_ko_8 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_8h_KO6_plus_eV_fwd.bw"))
bw_rev_ko_8 <- rtracklayer::import(file.path(bw_dir, "IP-cRIP_fourth_NSP9_8h_KO6_plus_eV_rev.bw"))

bw_fwd_snd1_sg <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixSG_snd1_huh7_fwd.bw"))
bw_rev_snd1_sg <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixSG_snd1_huh7_rev.bw"))

bw_fwd_snd1_yw <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_fwd.bw"))
bw_rev_snd1_yw <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_rev.bw"))

bw_fwd_N_ctrl <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_fwd.bw"))
bw_rev_N_ctrl <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_rev.bw"))

bw_fwd_N_ko <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_fwd.bw"))
bw_rev_N_ko <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_rev.bw"))

bw_fwd_cnbp <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_fwd.bw"))
bw_rev_cnbp <- rtracklayer::import(file.path(bw_dir, "IP-eCLIP_mixYW_snd1_huh7_rev.bw"))


bw_list <- list("12h_ctrl" = list(fwd = bw_fwd_ctrl_12, rev = bw_rev_ctrl_12, N = 628620, fwdN = 602748, revN = 25872),
                "12h_ko" = list(fwd = bw_fwd_ko_12, rev = bw_rev_ko_12, N = 244764, fwdN = 200356, revN = 44408),
                "8h_ctrl" = list(fwd = bw_fwd_ctrl_8, rev = bw_rev_ctrl_8, N = 52640, fwdN = 46576, revN = 6064),
                "8h_ko" = list(fwd = bw_fwd_ko_8, rev = bw_rev_ko_8, N = 14352, fwdN = 12592, revN = 1760),
                "snd1_sg" = list(fwd = bw_fwd_snd1_sg, rev = bw_rev_snd1_sg, N = 64292, fwdN = 23608, revN = 40684),
                "snd1_yw" = list(fwd = bw_fwd_snd1_yw, rev = bw_rev_snd1_yw, N = 13176, fwdN = 1732, revN = 11444))

rpm_list <- lapply(bw_list, function(lst_elem){
  norm_fw_rev(bw_cont_fwd = lst_elem$fwd, bw_cont_rev = lst_elem$rev, N_counts_fwd_rev = lst_elem$N)
})

rpm_gr_lst <- lapply(rpm_list, function(lst_elem){
  lapply(lst_elem, function(df){
    df_to_bw(cov_df = df)
  })
})

bw_out_dir <- "~/LRIB/eClip/Analysis/21102022/RPM_norm"
if(!file.exists(bw_out_dir)){
  dir.create(bw_out_dir)
}
for(i in seq_along(rpm_gr_lst)){
  cond <- names(rpm_gr_lst)[i]
  
  for(j in seq_along(rpm_gr_lst[[i]])){
    
    strand <- names(rpm_gr_lst[[i]])[j]
    file_name <- paste0(cond, "_", strand,"_rpm_norm_fwd_rev.bw")
    
    rtracklayer::export.bw(object = df_to_bw(cov_df = rpm_gr_lst[[i]][[j]] %>% as.data.frame()), file.path(bw_out_dir,file_name))
  }
}


rpm_list <- lapply(bw_list, function(lst_elem){
    fwd <- norm_df(bw_cont = lst_elem$fwd, N = lst_elem$fwdN)
    rev <- norm_df(bw_cont = lst_elem$rev, N = lst_elem$revN)
    list(fwd = fwd, rev = rev)
})

bw_out_dir <- "~/LRIB/eClip/Analysis/21102022/RPM_norm_per_strand"
if(!file.exists(bw_out_dir)){
  dir.create(bw_out_dir)
}
for(i in seq_along(rpm_list)){
  cond <- names(rpm_list)[i]
  
  for(j in seq_along(rpm_list[[i]])){
    
    strand <- names(rpm_list[[i]])[j]
    file_name <- paste0(cond, "_", strand,"_rpm_norm.bw")
    
    rtracklayer::export.bw(object = df_to_bw(cov_df = rpm_list[[i]][[j]] %>% as.data.frame()), file.path(bw_out_dir,file_name))
  }
}

