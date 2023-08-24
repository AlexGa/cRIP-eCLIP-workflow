library(dplyr)

units <- read.delim("units.tsv")
units <- units %>% filter(sample %in% c("SM", "IP"))
i_list <- split(units[,1:2], units$sample)

i_list[[1]] %>% inner_join(i_list[[2]], by = c("unit" = "unit")) %>% 
				rename(ip=sample.x,sm=sample.y) %>% 
				relocate(unit) %>% 
				write.table("peak_comparisons.tsv", sep = "\t", row.names = F, 
							 col.names = T, quote = F)