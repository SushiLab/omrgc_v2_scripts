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
# Correlation analysis #
########################
corr.T<-function(x){
  require(ppcor)
  cat("Correlation nº: ",x,"\n")
  tmp<-data.frame(metaT=as.numeric(metaT.norm.match.log2[x,]),metaG=as.numeric(metaG.norm.match.log2[x,]),Exp=as.numeric(ratio.mat[x,]))
  tmp.res.cor<-cor(tmp,method = "spearman")
  tmp.res.pcor<-pcor(tmp,method = "spearman")
  res<-c(tmp.res.cor[1,2],tmp.res.cor[1,3],tmp.res.pcor$estimate[1,2],tmp.res.pcor$estimate[1,3])
  names(res)<-c("TGcorr","TEcorr","TGpcorr","TEpcorr")
  res
}
corr.res<-sapply(1:nrow(metaG.norm.match.log2),corr.T)
corr.res<-cbind(t(corr.res),env.mat.match)
corr.res$Layer<-factor(corr.res$Layer,levels=c("SRF","DCM","MES","other"))

ggplot(data=corr.res,aes(x=Temperature,y=TEpcorr/TGpcorr,col=epi,label=gsub("TARA_","",Sample_name_goodlayer))) +
  geom_point(size=4,alpha=0.5,aes(shape=polar)) +
  geom_smooth(data=corr.res %>% filter(epi %in% c("EPI")),aes(x=Temperature,y=TEpcorr/TGpcorr,col=epi),method="lm",alpha=0.6,se=F) +
  facet_grid(.~epi) +
  theme_bw() +
  scale_color_manual(values = sushi.palette(alpha=0.7)[c(2,1,14)]) +
  ylab("Expression / Abundance") +
  geom_text_repel(size=4,color="black")
ggsave(file="../results/Fig4A_Varpart_corr_EPI.pdf",width=15,height=10)

summary(lm(corr.res$TEpcorr[corr.res$epi=="EPI"]/corr.res$TGpcorr[corr.res$epi=="EPI"]~corr.res$Temperature[corr.res$epi=="EPI"]))
summary(lm(corr.res$TEpcorr[corr.res$epi=="MES"]/corr.res$TGpcorr[corr.res$epi=="MES"]~corr.res$Temperature[corr.res$epi=="MES"]))
     

g0<-ggplot(data=corr.res,aes(x=daylength,y=TEpcorr/TGpcorr,label=sapply(strsplit(as.character(Sample),"_"),"[[",2))) +
  geom_point(size=4,alpha=0.5,aes(shape=polar)) +
  geom_smooth(data=corr.res %>% filter(epi %in% c("EPI")),aes(x=daylength,y=TEpcorr/TGpcorr),col="black",method="lm",alpha=0.6,se=F) +
  facet_grid(.~epi) +
  theme_bw() +
  scale_color_manual(values = sushi.palette(alpha=0.7)[c(2,1,14)]) +
  ylab("Expression / Abundance") #+
#ggsave(file="../results/SFig7A_Varpart_corr_EPI_daylength.pdf",width=10,height=5)

summary(lm(corr.res$TEpcorr[corr.res$epi=="EPI"]/corr.res$TGpcorr[corr.res$epi=="EPI"]~corr.res$daylength[corr.res$epi=="EPI"]))
summary(lm(corr.res$TEpcorr[corr.res$epi=="MES"]/corr.res$TGpcorr[corr.res$epi=="MES"]~corr.res$daylength[corr.res$epi=="MES"]))

## Test for differences in polar/non-polar after correcting for temperature after considering daylength
mod1<-lm(corr.res$TEpcorr[corr.res$epi=="EPI"]/corr.res$TGpcorr[corr.res$epi=="EPI"]~corr.res$Temperature[corr.res$epi=="EPI"])
summary(mod1)
mod2<-lm(corr.res$TEpcorr[corr.res$epi=="EPI"]/corr.res$TGpcorr[corr.res$epi=="EPI"]~corr.res$daylength[corr.res$epi=="EPI"])
summary(mod2)

#mod3<-lm(residuals(mod1)~corr.res$polar[corr.res$epi=="EPI"])
#summary(mod3)
#mod4<-lm(residuals(mod2)~corr.res$polar[corr.res$epi=="EPI"])
#summary(mod4)

t.test(residuals(mod1)~corr.res$polar[corr.res$epi=="EPI"])
t.test(residuals(mod2)~corr.res$polar[corr.res$epi=="EPI"])


toplot<-corr.res[corr.res$epi=="EPI",]
toplot$residuals.temp<-residuals(mod1)
toplot$residuals.daylength<-residuals(mod2)

g1<-ggplot(data=toplot,aes(x=polar,y=residuals.temp)) +
  geom_violin(draw_quantiles = 0.5,scale = "width",fill="gray") +
  theme_bw() +
  xlab("") +
  ylab("Residuals") +
  labs(title="Relative contribution expression ~ Temperature")
g2<-ggplot(data=toplot,aes(x=polar,y=residuals.daylength)) +
  geom_violin(draw_quantiles = 0.5,scale = "width",fill="gray") +
  theme_bw() +
  xlab("") +
  ylab("Residuals") +
  labs(title="Relative contribution expression ~ Day length")

g0 / (g1 | g2)
ggsave(file="../results/SFig7A_Varpart_temp_vs_daylength.pdf",width=10,height=10)
