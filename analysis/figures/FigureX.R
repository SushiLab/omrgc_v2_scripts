library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
library(geosphere)
source("../data/lib/sushipal.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])

###############
## LOAD DATA ##
###############

metaG.norm.match.log2<-fread("zcat < ../data/processed/NOG_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("zcat < ../data/processed/NOG_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("zcat < ../data/processed/NOG_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("zcat < ../data/processed/NOG_metaT.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("zcat < ../data/processed/NOG_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]

env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match$epi<-env.mat.match$Layer
levels(env.mat.match$epi)<-c("EPI","MES","EPI","EPI")
env.mat.match$date<-as.POSIXct(substr(as.character(env.mat.match$Event.date),1,10))
env.mat.match$doy<-yday(env.mat.match$date)
env.mat.match$daylength<-as.numeric(daylength(lat = env.mat.match$Latitude,doy=env.mat.match$doy))
env.mat.match$hour<-as.numeric(substr(env.mat.match$Event.date,12,13))
env.mat.match<-env.mat.match[match(rownames(metaG.norm.match),env.mat.match$Barcode),]

########################

metaG.norm.match.log2.gath<-metaG.norm.match.log2 %>% 
  rownames_to_column(var="sample") %>%
  gather(key="OG",value="log_abundance",-sample)
metaG.norm.match.gath<-metaG.norm.match %>% 
  rownames_to_column(var="sample") %>%
  gather(key="OG",value="abundance",-sample)
metaT.norm.match.log2.gath<-metaT.norm.match.log2 %>% 
  rownames_to_column(var="sample") %>%
  gather(key="OG",value="log_transcription",-sample)
metaT.norm.match.gath<-metaT.norm.match %>% 
  rownames_to_column(var="sample") %>%
  gather(key="OG",value="transcription",-sample)
ratio.mat.gath<-ratio.mat %>% 
  rownames_to_column(var="sample") %>%
  gather(key="OG",value="expression",-sample)
toplot<-left_join(metaG.norm.match.log2.gath,ratio.mat.gath,by = c("sample", "OG")) %>%
  left_join(metaG.norm.match.gath,by = c("sample", "OG")) %>%
  left_join(metaT.norm.match.log2.gath,by = c("sample", "OG")) %>%
  left_join(metaT.norm.match.gath,by = c("sample", "OG")) %>%
  left_join(env.mat.match,by=c("sample"="Barcode")) %>%
  filter(abundance>0 | transcription>0) %>%
  filter(OG!="unknown")
occ.df<-data.frame(OG=colnames(metaG.norm.match),occ=apply(metaG.norm.match,2,function(x){length(x[x>0])}))
occ.df<-occ.df %>% filter(OG!="unknown")
toplot<-toplot %>% left_join(occ.df,by="OG")
MGs<-fread("../data/lib/40_MGs.tsv",sep="\t",header=T)
toplot<-toplot %>%
  left_join(MGs,by=c("OG"="COG")) %>%
  mutate(expression.type=cut(expression,c(-Inf,-2,2,Inf))) %>%
  mutate(expression.type=fct_recode(expression.type,"Underexpressed"="(-Inf,-2]","Overexpressed"="(2, Inf]","Not differentially expressed"="(-2,2]"))
toplot<-toplot %>%  mutate(core=cut(occ,c(-1,128,129),labels = c("non-core","core")))

#############################################
nog.info<-fread("../data/lib/NOG.annotations.tsv")
nog.cat<-fread("../data/lib/NOG.categories.tsv",header=F)

nog.info<-nog.info %>% select(OG=V2,OG.cat=V5,OG.description=V6) %>%
  left_join(nog.cat,by=c("OG.cat"="V1")) %>%
  rename(OG.cat.description=V2)

tmp<-toplot %>%
  select(OG,log_abundance,expression,log_transcription) %>%
  group_by(OG) %>%
  summarise_all(mean) %>% head()
  gather(key = "metric",value="value") %>%
  left_join(occ.df,by="OG") %>%
  mutate(core=cut(occ,c(-1,120,129),labels = c("non-core","core"))) %>%
  left_join(nog.info,by="OG")
  

ggplot(data=tmp,aes(x=OG.cat.description,y=log_abundance,fill=core)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  coord_flip()

ggplot(data=tmp,aes(x=OG.cat.description,y=log_transcription,fill=core)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  coord_flip()

ggplot(data=tmp,aes(x=OG.cat.description,y=expression,fill=core)) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_violin(draw_quantiles = 0.5,scale = "width",alpha=0.6) +
  coord_flip() +
  theme_bw()

#############################################
#limits<-range(c(range(toplot$log_transcription),range(toplot$log_abundance)))
#for (i in unique(toplot$sample)){
#  tit<-unique(toplot$Sample_name_goodlayer[toplot$sample==i])
#  
#  main.plot<-ggplot(toplot %>% filter(sample==i),aes(x=log_abundance,y=log_transcription,col=cut(expression,c(-10,-2,2,10)))) +
#    geom_point(alpha=0.2) +
#    geom_density_2d(col="gray40") +
#    #geom_point(data=toplot %>% filter(sample==i & !is.na(type)),aes(x=log_abundance,y=log_transcription,col=type)) +
#    geom_smooth(method="loess",se=F,span=0.5,col="black") +
#    #facet_grid(core~.) +
#    theme_bw() +
#    geom_abline(linetype=2) +
#    #scale_color_manual(values=c("forestgreen","blue4")) +
#    theme(legend.position = "bottom",legend.title = element_blank()) +
#    xlim(limits) +
#    ylim(limits) +
#    coord_fixed()
#  ggsave(paste("../results/abundance_vs_transcription/",tit,".pdf",sep=""),main.plot)
#} 
  
#############################################

tmp<-toplot %>%
  mutate(expression.type=cut(expression,c(-Inf,-2,2,Inf))) %>%
  group_by(sample,expression.type) %>%
  summarise(nOGs=n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(propOGs=nOGs/sum(nOGs)) %>%
  left_join(env.mat.match,by=c("sample"="Barcode"))

ggplot(tmp,aes(x=epi,y=propOGs,fill=polar)) +
  geom_boxplot() +
  facet_wrap(~expression.type)


ggplot(tmp%>% filter(Layer %in% c("SRF","DCM")),aes(x=Temperature,y=propOGs,col=expression.type,shape=Layer)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm",se=F)

ggplot(tmp%>% filter(Layer %in% c("SRF","DCM")),aes(x=abs(Latitude),y=propOGs,col=expression.type,shape=Layer)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess",se=F)

ggplot(tmp%>% filter(Layer %in% c("SRF","DCM")),aes(x=ChlorophyllA,y=propOGs,col=expression.type)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm",se=F,formula=y~log10(x))

lm.epi<-lm(tmp$propOGs[tmp$Layer %in% c("SRF","DCM") & tmp$expression.type=="(2, Inf]"]~log10(tmp$ChlorophyllA[tmp$Layer %in% c("SRF","DCM")  & tmp$expression.type=="(2, Inf]"]))
summary(lm.epi)

#########################
## Logistic regression ##
#########################


toplot<-toplot %>%
  mutate(expression.type=cut(expression,c(-Inf,-2,2,Inf))) %>%
  mutate(expression.type=fct_recode(expression.type,"Underexpressed"="(-Inf,-2]","Overexpressed"="(2, Inf]","Not differentially expressed"="(-2,2]"))
ggplot(data=toplot %>% filter(sample=="TARA_A100000164-TARA_A100000165"),aes(x=log_abundance,fill=expression.type.2levels)) +
  geom_density(alpha=0.5)


res<-NULL
for (i in unique(toplot$sample)){
  tmp<-toplot %>%
    filter(sample==i)
  mod<-glm(prob~abundance,family = "binomial",data=tmp)
  res<-rbind(res,summary(mod)$coefficients[c(2,8)])
}
colnames(res)<-c("slope","Prob")
res<-as.data.frame(res)
res$sample<-unique(toplot$sample)
res<-res %>%
  left_join(env.mat.match,by=c("sample"="Barcode"))

ggplot(data=res,aes(x=Latitude,y=slope,col=cut(Prob,c(-Inf,0.05,Inf)))) +
  geom_point() +
  facet_wrap(~epi)



ggplot(data=tmp,aes(x=abundance,y=prob)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

ggplot(data=tmp,aes(x=log_abundance,fill=expression.type.2levels)) +
  geom_density(alpha=0.5)
