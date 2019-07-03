library(data.table)
library(tidyverse)
library(patchwork)
library(vegan)

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
mgs.metaG<-fread("../../analysis/preliminary/nifH/MG_metaG_median.tsv",sep="\t",header=T)
markers.metaG<-fread("../../analysis/preliminary/nifH/marker_KOs_metaG.tsv.gz",sep="\t",header=T)

markers.metaG.gath<-markers.metaG %>%
  gather(key="sample",value="markers_ab",-OMRGC_ID) %>%
  left_join(mgs.metaG,by="sample") %>%
  mutate(markers_ab_norm=markers_ab/mg_median)

# MetaT
mgs.metaT<-fread("../../analysis/preliminary/nifH/MG_metaT_median.tsv",sep="\t",header=T)
markers.metaT<-fread("../../analysis/preliminary/nifH/marker_KOs_metaT.tsv.gz",sep="\t",header=T)

markers.metaT.gath<-markers.metaT %>%
  gather(key="sample",value="markers_tra",-OMRGC_ID) %>%
  left_join(mgs.metaT,by="sample") %>%
  mutate(markers_tra_norm=markers_tra/mg_median)

# Merge metaG and metaT
env.tab<-fread("../../data/processed/KO_env.mat.match.txt.gz",sep="\t",header=T)
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

marker.info<-fread("../../analysis/preliminary/nifH/marker_KOs_catalog.tab.gz",sep="\t",header=F)
marker.info<-marker.info %>%
  select(OMRGC_ID=V2,KO=V3,V6:V13) %>%
  mutate(tax=gsub(";;","",paste(V6,V7,V8,V9,V10,V11,V12,V12,V13,sep=";")))
final.df<-final.df %>%
  left_join(marker.info,by="OMRGC_ID")

# Plots nifH
naming<-fread("../../analysis/preliminary/nifH/nifH_vs_Zehr_DB_HBDs_nr_processed.m8",sep="\t",header=T,data.table = F)
naming$V1<-sapply(strsplit(as.character(naming$V1),"_"),"[[",7)
naming$V2.y<-sapply(strsplit(as.character(naming$V2.y)," "),"[[",1)
naming.clean<-naming %>% filter(V15>=90 & V3>=90)


toplot<-final.df %>%
  filter(KO %in% c("K02588")) %>%
  select(OMRGC_ID,markers_ab_norm_log2,markers_tra_norm_log2,expression,Latitude,epi,polar) %>%
  gather(key="measure",value="value",-OMRGC_ID,-Latitude,-epi,-polar)

toplot<-final.df %>%
  filter(KO %in% c("K02588")) %>%
  select(OMRGC_ID,markers_ab_norm,markers_tra_norm,expression,Latitude,epi,polar,Depth.nominal) %>%
  gather(key="measure",value="value",-OMRGC_ID,-Latitude,-epi,-polar,-Depth.nominal)
order.df<-toplot %>% filter(measure=="markers_ab_norm") %>% filter(value>0) %>% group_by(OMRGC_ID) %>% summarise(median=median(value,na.rm=T)) %>% arrange(median)
toplot$OMRGC_ID<-factor(toplot$OMRGC_ID,levels = c(order.df$OMRGC_ID),ordered = T)

toplot<-toplot %>% left_join(naming.clean,by=c("OMRGC_ID"="V1"))
toplot$V2.y[is.na(toplot$V2.y)]<-"no annotation"
toplot$goodname<-paste(toplot$OMRGC_ID,"\n(",toplot$V2.y,")",sep="")
toplot$goodname[toplot$goodname=="OM-RGC.v2.009268969\n(ACJ53724.1|uncultured_cyanobacterium)"]<-"OM-RGC.v2.009268969\n(ACJ53724.1|UCYN-A)"

g1<-ggplot(data=toplot %>% filter(measure!="expression" & value>0),aes(x=goodname,y=log2(value),col=measure,fill=measure)) +
  geom_boxplot(alpha=0.3,position = position_dodge(width = 0.9),outlier.size = 0) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2)) +
  #geom_violin(scale = "width",draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.position = "top",axis.text.x = element_blank()) +
  xlab(NULL) +
  ylab("Gene / Transcript\nabundance [log2]") +
  scale_fill_manual(name="",labels=c("Gene abundance","Transcript abundance"),values=c("#EAAD52","#35CCFF")) +
  scale_color_manual(name="",labels=c("Gene abundance","Transcript abundance"),values=c("#EAAD52","#35CCFF"))

g2<-ggplot(data=toplot %>% filter(measure!="expression" & value>0),aes(x=goodname,y=Latitude,col=measure,fill=measure)) +
  geom_hline(yintercept = c(-60,60),linetype=2) +
  geom_boxplot(alpha=0.3,position = position_dodge(width = 0.9),outlier.size = 0) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2)) +
  #geom_violin(scale = "width",draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.position = "none",axis.text.x = element_blank()) +
  xlab(NULL) +
  ylab("Absolute latitude") +
  scale_fill_manual(values=c("#EAAD52","#35CCFF")) +
  scale_color_manual(values=c("#EAAD52","#35CCFF"))

g3<-ggplot(data=toplot %>% filter(measure!="expression" & value>0),aes(x=goodname,y=Depth.nominal,col=measure,fill=measure)) +
  geom_hline(yintercept = 200,linetype=2) +
  geom_boxplot(alpha=0.3,position = position_dodge(width = 0.9),outlier.size = 0) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2)) +
  #geom_violin(scale = "width",draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.position = "none",axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab(NULL) +
  ylab("Depth [m]") +
  scale_fill_manual(values=c("#EAAD52","#35CCFF")) +
  scale_color_manual(values=c("#EAAD52","#35CCFF")) +
  scale_y_reverse()

g<-g1 / g2 / g3
ggsave("../../results/paper_figures/review/FigX_nifH_genelevel.pdf",g,width=8,height=12)

