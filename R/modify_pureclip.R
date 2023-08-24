library(dplyr)
library(readr)
# pureclip_file <- "~/LRIB/hzi_projects/NSP9_new/results/pureclip/homo_sapiens_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_CTRL2-IP-SM/cRIP_fifth_rep12_NSP9_12h_CTRL2_IP_SM_XL_sites_fwd.bed"
# chromosomes <- c("MN908947.3")

pureclip_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

pureclip_df <- readr::read_delim(pureclip_file, col_names = c("chr", "start", "end", "name", "score", "strand", "anno"))

if(any(grepl(x = names(snakemake@params), pattern = "region"))){
  
  chromosomes <- snakemake@params$region
  pureclip_df <- pureclip_df %>% filter(chr %in% chromosomes)
}

pureclip_df_sub <- pureclip_df %>% dplyr::mutate(end = if_else(start == -1, true = end + 1, false = end),
                                              start = if_else(start == -1, true = 0, false = start)) %>% 
                                       dplyr::group_by(chr, start, end) %>% 
                                       dplyr::filter(score == max(score)) %>% ungroup() %>% 
                                       dplyr::arrange(chr, start)

readr::write_delim(x = pureclip_df_sub, file = output_file, delim = "\t", col_names = FALSE, quote = "none")