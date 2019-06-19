library(data.table)
library(tidyverse)

mgs.genes<-fread("MGs.tab",sep="\t",header=F)

# Compute median of MGs in metaG and save
mgs.metaG<-fread("MG_metaG.tsv",sep="\t",header=T)

mgs.gath<-mgs.genes %>%
  select(OMRGC_ID=V2,KO=V3) %>%
  left_join(mgs.metaG,by="OMRGC_ID") %>%
  gather(key = "sample",value="value",-OMRGC_ID,-KO)

tmp<-mgs.gath %>%
  group_by(KO,sample) %>%
  summarise(value=sum(value))

tmp<-tmp %>%
  group_by(sample) %>%
  summarise(mg_median=median(value))
fwrite(tmp,"MG_metaG_median.tsv",sep="\t")

# Compute median of MGs in metaT and save
mgs.metaG<-fread("MG_metaT.tsv",sep="\t",header=T)

mgs.gath<-mgs.genes %>%
  select(OMRGC_ID=V2,KO=V3) %>%
  left_join(mgs.metaG,by="OMRGC_ID") %>%
  gather(key = "sample",value="value",-OMRGC_ID,-KO)

tmp<-mgs.gath %>%
  group_by(KO,sample) %>%
  summarise(value=sum(value))

tmp<-tmp %>%
  group_by(sample) %>%
  summarise(mg_median=median(value))
fwrite(tmp,"MG_metaT_median.tsv",sep="\t")