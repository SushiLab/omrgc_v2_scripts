library(data.table)
library(ggplot2)
library(tidyverse)

# Load metaG
metaG.norm<-fread("../data/processed/NOG_metaG.norm.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm)<-metaG.norm$V1
metaG.norm<-metaG.norm[,-1]
env.mat.metaG<-fread("../data/processed/NOG_env.mat.metaG.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.metaG)<-env.mat.metaG$V1
env.mat.metaG<-env.mat.metaG[,-1]
env.mat.metaG$epi<-env.mat.metaG$Layer
levels(env.mat.metaG$epi)<-c("EPI","MES","EPI","EPI")

domain<-fread("../data/processed/miTags/to_release/mitags_tab_domain.tsv",sep="\t",header=T,data.table = F)
domain<-domain %>%
  mutate(Eukaryota.prop=Eukaryota/(Archaea+Bacteria+Eukaryota+unclassified))

env.mat.metaG<-env.mat.metaG %>%
  left_join(domain,by=c("Barcode"="V1"))

env.mat.metaG$gsize<-apply(metaG.norm,1,sum)

p<-ggplot(data=env.mat.metaG,aes(x=ChlorophyllA,y=gsize,shape=epi)) +
  geom_point() +
  #geom_smooth(method="lm",se=F,col="gray") +
  theme_bw() +
  xlab("Clorophyll a concentration") +
  ylab("Genome size (number of OGs)") +
  facet_grid(.~paste("Upper filter:",upper.size.fraction))

to.share<-env.mat.metaG %>%
  dplyr::select(Barcode,Sample_name_original=Sample_name,Sample_name_goodlayer,lower.size.fraction,upper.size.fraction,Genome_size=gsize,Latitude,Longitude,Depth.nominal,Layer,Temperature,Salinity,ChlorophyllA,Biogeographical.province,Eukaryota.prop)

ggsave("../results/Genome_size.pdf",p,width = 9,height=4)
fwrite(to.share,file = "../results/Genome_size.tsv",sep="\t",quote = F)
