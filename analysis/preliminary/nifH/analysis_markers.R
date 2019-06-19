library(data.table)
library(tidyverse)
library(patchwork)
library(vegan)
setwd('nifH')
source("../data/lib/eveness_functions.R")

find_hull<-function(dat,x.name,y.name,lev.name.1,lev.name.2){
  hull<-NULL
  for (i in levels(as.factor(dat[,lev.name.1]))){
    for (j in levels(as.factor(dat[,lev.name.2]))){
      X<-na.exclude(dat[,x.name][dat[,lev.name.1]==i & dat[,lev.name.2]==j])
      Y<-na.exclude(dat[,y.name][dat[,lev.name.1]==i & dat[,lev.name.2]==j])
      if(length(X)>=1 & length(Y)>=1) pos<-chull(X,Y)
      if(length(X)>=1 & length(Y)>=1) hull<-rbind(hull,data.frame(x=X[pos],y=Y[pos],lev1=rep(i,length(pos)),lev2=rep(j,length(pos))))
    }
  }
  colnames(hull)<-c(x.name,y.name,lev.name.1,lev.name.2)
  hull
}

# MetaG
mgs.metaG<-fread("MG_metaG_median.tsv",sep="\t",header=T)
markers.metaG<-fread("marker_KOs_metaG.tsv",sep="\t",header=T)

markers.metaG.gath<-markers.metaG %>%
  gather(key="sample",value="markers_ab",-OMRGC_ID) %>%
  left_join(mgs.metaG,by="sample") %>%
  mutate(markers_ab_norm=markers_ab/mg_median)

# MetaT
mgs.metaT<-fread("MG_metaT_median.tsv",sep="\t",header=T)
markers.metaT<-fread("marker_KOs_metaT.tsv",sep="\t",header=T)

markers.metaT.gath<-markers.metaT %>%
  gather(key="sample",value="markers_tra",-OMRGC_ID) %>%
  left_join(mgs.metaT,by="sample") %>%
  mutate(markers_tra_norm=markers_tra/mg_median)

# Merge metaG and metaT
env.tab<-fread("../data/processed/KO_env.mat.match.txt.gz",sep="\t",header=T)
env.tab$metaG_sample<-sapply(strsplit(env.tab$Barcode,"-"),"[[",1)
env.tab$metaT_sample<-sapply(strsplit(env.tab$Barcode,"-"),"[[",2)
env.tab$epi<-as.factor(env.tab$Layer)
levels(env.tab$epi)<-c("EPI","MES","EPI","EPI")



markers.metaG.gath <- markers.metaG.gath %>%
  filter(sample %in% env.tab$metaG_sample) %>%
  select(OMRGC_ID,sample,markers_ab_norm) %>%
  left_join(env.tab,by=c("sample"="metaG_sample")) %>%
  select(OMRGC_ID,markers_ab_norm,Barcode)
markers.metaT.gath <- markers.metaT.gath %>%
  filter(sample %in% env.tab$metaT_sample) %>%
  select(OMRGC_ID,sample,markers_tra_norm) %>%
  left_join(env.tab,by=c("sample"="metaT_sample")) %>%
  select(OMRGC_ID,markers_tra_norm,Barcode)

final.df<-markers.metaG.gath %>%
  left_join(markers.metaT.gath,by=c("Barcode","OMRGC_ID")) %>%
  left_join(env.tab,by="Barcode") %>%
  mutate(expression=log2(markers_tra_norm+1e-10)-log2(markers_ab_norm+1e-10),markers_tra_norm_log2=log2(markers_tra_norm+1e-10),markers_ab_norm_log2=log2(markers_ab_norm+1e-10))
final.df$expression[final.df$markers_tra_norm==0 & final.df$markers_ab_norm==0]<-NA

marker.info<-fread("marker_KOs_catalog.tab",sep="\t",header=F)
marker.info<-marker.info %>%
  select(OMRGC_ID=V2,KO=V3)
final.df<-final.df %>%
  left_join(marker.info,by="OMRGC_ID")

# Plots nifH
# ggplot() +
#   geom_line(data=final.df %>% filter(KO=="K02588"),aes(x=fct_reorder(Sample_name_goodlayer,markers_ab_norm,.fun = median),y=markers_ab_norm,col="gray",group=1)) +
#   geom_line(data=final.df %>% filter(KO=="K02588"),aes(x=fct_reorder(Sample_name_goodlayer,markers_ab_norm,.fun = median),y=markers_tra_norm,col="red",group=1)) +
#   theme_bw() +
#   facet_wrap(~OMRGC_ID) +
#   scale_y_log10()

ggplot(data=final.df %>% filter(KO=="K02588" & !is.na(expression) & markers_ab_norm>0 & markers_tra_norm>0),aes(x=NO2NO3,y=expression),col="black",alpha=0.5) +
  geom_point() +
  theme_bw() +
  geom_smooth(method="lm",formula=y~sqrt(x),se=F) +
  geom_hline(yintercept = 0,linetype=2) +
  facet_wrap(~OMRGC_ID)


ggplot() +
  geom_point(data=final.df %>% filter(KO=="K02588"),aes(x=abs(Latitude),y=markers_ab_norm),col="black",alpha=0.5) +
  geom_point(data=final.df %>% filter(KO=="K02588"),aes(x=abs(Latitude),y=markers_tra_norm),col="red",alpha=0.5) +
  #geom_point() +
  #geom_line() +
  theme_bw() +
  facet_wrap(~OMRGC_ID) +
  scale_y_log10() +
  geom_vline(xintercept = 60,linetype=2)

ggplot() +
  geom_point(data=final.df %>% filter(KO=="K02588"),aes(x=abs(Latitude),y=markers_ab_norm),col="black",alpha=0.5) +
  geom_point(data=final.df %>% filter(KO=="K02588"),aes(x=abs(Latitude),y=markers_tra_norm),col="red",alpha=0.5) +
  #geom_point() +
  #geom_line() +
  theme_bw() +
  facet_grid(epi~OMRGC_ID) +
  scale_y_log10() +
  geom_vline(xintercept = 60,linetype=2)

ggplot(data=final.df %>% filter(KO=="K02588" & !is.na(expression)),aes(x=fct_reorder(OMRGC_ID,expression,.fun = median),y=expression)) +
  geom_boxplot(outlier.size = 0.1) +
  #geom_violin(draw_quantiles = 0.5,scale = "width") +
  coord_flip() +
  geom_hline(yintercept = 0,linetype=2)

########################################################################################################################



rich.df<-final.df %>%
  group_by(Barcode,KO) %>%
  summarise(richness.tra=sum(markers_tra_norm>0),richness.ab=sum(markers_ab_norm>0),
            shannon.ab=diversity(markers_ab_norm),shannon.tra=diversity(markers_tra_norm),
            eveness.ab=PielouE(markers_ab_norm),eveness.tra=PielouE(markers_tra_norm)) %>%
  gather(key="richness_type",value="richness",-KO,-Barcode,-shannon.ab,-shannon.tra,-eveness.ab,-eveness.tra) %>%
  gather(key="shannon_type",value="shannon",-KO,-Barcode,-richness_type,-richness,-eveness.ab,-eveness.tra) %>%
  gather(key="eveness_type",value="eveness",-KO,-Barcode,-richness_type,-richness,-shannon_type,-shannon) %>%
  left_join(env.tab,by="Barcode")

ggplot(data=rich.df,aes(x=fct_reorder(KO,richness,.fun=median),y=richness,fill=richness_type)) +
  geom_boxplot() +
  scale_y_sqrt() +
  coord_flip()
ggplot(data=rich.df,aes(x=fct_reorder(KO,richness,.fun=median),y=shannon,fill=shannon_type)) +
  geom_boxplot() +
  scale_y_sqrt() +
  coord_flip()


ggplot(data=rich.df %>% spread(key="eveness_type",value="eveness") %>% arrange(epi),aes(x=eveness.ab,y=eveness.tra,col=epi,fill=epi)) +
  geom_point(alpha=0.5) +
  #geom_density_2d() +
  geom_polygon(data=find_hull(rich.df %>% spread(key="eveness_type",value="eveness") %>% arrange(epi) %>% as.data.frame(),"eveness.ab","eveness.tra","epi","KO"),alpha=0.5) +
  facet_wrap(~KO) +
  geom_abline() +
  coord_fixed()

ggplot(data=rich.df %>% spread(key="shannon_type",value="shannon") %>% arrange(epi),aes(x=shannon.ab,y=shannon.tra,col=epi,fill=epi)) +
  geom_point(alpha=0.5) +
  #geom_density_2d() +
  geom_polygon(data=find_hull(rich.df %>% spread(key="shannon_type",value="shannon") %>% arrange(epi) %>% as.data.frame(),"shannon.ab","shannon.tra","epi","KO"),alpha=0.5) +
  facet_wrap(~KO) +
  geom_abline() +
  coord_fixed()

ggplot(data=rich.df %>% spread(key="richness_type",value="richness") %>% arrange(Biogeographical.province),aes(x=richness.ab,y=richness.tra,col=Biogeographical.province,fill=Biogeographical.province)) +
  geom_point(alpha=0.5) +
  geom_polygon(data=find_hull(rich.df %>% spread(key="richness_type",value="richness") %>% arrange(Biogeographical.province) %>% as.data.frame(),"richness.ab","richness.tra","Biogeographical.province","KO"),alpha=0.5) +
  facet_wrap(~KO) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline() +
  coord_fixed() +
  theme(legend.position = "none")
########################################################
p<-ggplot() +
  geom_point(data=final.df %>% filter(KO=="K01602"),aes(x=abs(Latitude),y=markers_ab_norm),col="black",alpha=0.5) +
  geom_point(data=final.df %>% filter(KO=="K01602"),aes(x=abs(Latitude),y=markers_tra_norm),col="red",alpha=0.5) +
  #geom_point() +
  #geom_line() +
  theme_bw() +
  facet_wrap(~OMRGC_ID) +
  scale_y_log10() +
  geom_vline(xintercept = 60,linetype=2)
#ggsave("/Users/guillemsalazar/Desktop/prova.pdf",p,width = 100,height=100,limitsize = F)
p
