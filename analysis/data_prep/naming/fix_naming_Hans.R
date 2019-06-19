library(data.table)
library(tidyverse)

# Load the correct mapping
mapping.correct<-fread("MAIN_MAPPING_bysample.tsv",sep="\t",header=T,data.table = F)

# Get only the metaG and remove the "_G"
mapping.correct<-mapping.correct %>%
  filter(grepl("_G",old_name)) %>%
  mutate(old_name=gsub("_G","",old_name))

# Load the mapping used by Hans
mapping.wrong<-fread("Table_S1.tsv",sep="\t",header=T,data.table = F)
mapping.wrong<-mapping.wrong %>%
  select(old_name=Sample_name,PANGAEA_ID_incorrect=`PANGAEA sample id`)

# Match both
correct.incorrect.map<-left_join(mapping.correct,mapping.wrong,by="old_name") %>%
  mutate(equal=PANGAEA_ID==PANGAEA_ID_incorrect)

# Save it
fwrite(correct.incorrect.map,file="correcct_incorrect_PANGAEA_ID_map.tsv",quote = F,sep = "\t")