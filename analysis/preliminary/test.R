library(data.table)
library(tidyverse)

map<-fread("../naming/MAIN_MAPPING_bysample.tsv",sep="\t",header=T,data.table = F)

# Mapping stats against the catalog
stats.metaG<-fread("mapping_stats.tsv",sep="\t",header=T,data.table = F)
if(all(stats.metaG$sample %in% map$old_name)==F) cat("Not all samples match the map file")
stats.metaG$metaG_metaT<-"metaG"

stats.metaT<-fread("mapping_stats_T.tsv",sep="\t",header=T,data.table = F)
if(all(stats.metaT$sample %in% map$old_name)==F) cat("Not all samples match the map file")
stats.metaT$metaG_metaT<-"metaT"

# Rearrange
stats.all<-bind_rows(stats.metaG,stats.metaT) %>%
  left_join(map,by=c("sample"="old_name")) %>%
  dplyr::select(sample,PANGAEA_ID,metaG_metaT,total_inserts:db_average_entry_length)

# Get only the 180 + 187
stats.all.red<-stats.all[grep("0.22-[3,1.6]",stats.all$sample),]

# Save data
fwrite(stats.all,file = "mapping_stats_all.tsv",sep="\t",quote = F)
fwrite(stats.all.red,file = "mapping_stats_prok_size_fraction.tsv",sep="\t",quote = F)