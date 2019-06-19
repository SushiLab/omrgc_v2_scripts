library(data.table)
library(tidyverse)
library(patchwork)

# MetaG
mgs.metaG<-fread("MG_metaG_median.tsv",sep="\t",header=T)
nifh.metaG<-fread("nifH_metaG.tsv",sep="\t",header=T)

nifh.metaG.gath<-nifh.metaG %>%
  gather(key="sample",value="nifh_ab",-OMRGC_ID) %>%
  left_join(mgs.metaG,by="sample") %>%
  mutate(nifh_ab_norm=nifh_ab/mg_median)

# MetaT
mgs.metaT<-fread("MG_metaT_median.tsv",sep="\t",header=T)
nifh.metaT<-fread("nifH_metaT.tsv",sep="\t",header=T)

nifh.metaT.gath<-nifh.metaT %>%
  gather(key="sample",value="nifh_tra",-OMRGC_ID) %>%
  left_join(mgs.metaT,by="sample") %>%
  mutate(nifh_tra_norm=nifh_tra/mg_median)

# Merge metaG and metaT
env.tab<-fread("../data/processed/KO_env.mat.match.txt.gz",sep="\t",header=T)
env.tab$metaG_sample<-sapply(strsplit(env.tab$Barcode,"-"),"[[",1)
env.tab$metaT_sample<-sapply(strsplit(env.tab$Barcode,"-"),"[[",2)

nifh.metaG.gath <- nifh.metaG.gath %>%
  filter(sample %in% env.tab$metaG_sample) %>%
  select(OMRGC_ID,sample,nifh_ab_norm) %>%
  left_join(env.tab,by=c("sample"="metaG_sample")) %>%
  select(OMRGC_ID,nifh_ab_norm,Barcode)
nifh.metaT.gath <- nifh.metaT.gath %>%
  filter(sample %in% env.tab$metaT_sample) %>%
  select(OMRGC_ID,sample,nifh_tra_norm) %>%
  left_join(env.tab,by=c("sample"="metaT_sample")) %>%
  select(OMRGC_ID,nifh_tra_norm,Barcode)

final.df<-nifh.metaG.gath %>%
  left_join(nifh.metaT.gath,by=c("Barcode","OMRGC_ID")) %>%
  left_join(env.tab,by="Barcode") %>%
  mutate(expression=log2(nifh_tra_norm+1e-10)-log2(nifh_ab_norm+1e-10))
final.df$expression[final.df$nifh_tra_norm==0 & final.df$nifh_ab_norm==0]<-NA

g1<-ggplot(data=final.df,aes(x=OMRGC_ID,y=nifh_ab_norm,fill=OMRGC_ID)) +
  geom_boxplot() +
  scale_y_sqrt()
g2<-ggplot(data=final.df,aes(x=OMRGC_ID,y=nifh_tra_norm,fill=OMRGC_ID)) +
  geom_boxplot() +
  scale_y_sqrt()
g1 / g2


ggplot(data=final.df,aes(x=log2(nifh_ab_norm+1e-07),y=log2(nifh_tra_norm+1e-07),col=OMRGC_ID)) +
  geom_point() +
  facet_wrap(~OMRGC_ID) +
  geom_smooth(method="loess") +
  geom_abline()
  
ggplot(data=final.df,aes(x=OMRGC_ID,y=expression,fill=OMRGC_ID)) +
  geom_boxplot() +
  #geom_violin(draw_quantiles = 0.5,scale = "width") +
  coord_flip()
