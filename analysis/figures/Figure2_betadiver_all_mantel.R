library(tidyverse)
library(data.table)
library(patchwork)
source("../data/lib/sushipal.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])
library(vegan)

# Load environmental data
env.mat<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat)<-env.mat$V1
env.mat<-env.mat[,-1]
env.mat$epi<-env.mat$Layer
levels(env.mat$epi)<-c("EPI","MES","EPI","EPI")
levels(env.mat$Layer)<-c("DCM","MES","DCM","SRF")
env.mat$Absolute.Latitude<-abs(env.mat$Latitude)

# Load marine distances
mar.dist.mat<-fread("zcat < ../data/processed/mardist.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(mar.dist.mat)<-mar.dist.mat$V1
mar.dist.mat<-mar.dist.mat[,-1]
#mar.dist.mat<-mar.dist.mat[match(env.mat$Sample_name,rownames(mar.dist.mat)),match(env.mat$Sample_name,rownames(mar.dist.mat))]

# metaG
metaG.norm<-fread("zcat < ../data/processed/NOG_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm)<-metaG.norm$V1
metaG.norm<-metaG.norm[,-1]
BCdist.metaG<-dist(metaG.norm)

# metaT
metaT.norm<-fread("zcat < ../data/processed/NOG_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm)<-metaT.norm$V1
metaT.norm<-metaT.norm[,-1]
BCdist.metaT<-dist(metaT.norm)

# Expression
exp.norm<-fread("zcat < ../data/processed/NOG_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(exp.norm)<-exp.norm$V1
exp.norm<-exp.norm[,-1]
BCdist.exp<-dist(exp.norm)

# Reorder env mat
pos<-match(rownames(metaG.norm),env.mat$Barcode)
env.mat<-env.mat[pos,]

# Reorder marine distances
pos<-match(env.mat$Barcode,rownames(mar.dist.mat))
mar.dist.mat<-mar.dist.mat[pos,pos]

# miTags
otu.table<-fread("zcat < ../data/processed/miTags/OTU.tab.TARA180.noeuks.log2.txt.gz",sep="\t",header=T,data.table = F)
rownames(otu.table)<-otu.table$V1
otu.table<-otu.table[,-1]
otu.table<-otu.table[match(gsub("-",".",env.mat$Sample_name),rownames(otu.table)),]
BCdist.tax<-dist(otu.table)

env.mat %>%
  filter(Layer %in% c("SRF","DCM")) %>%
  mutate(lat.fct=cut(abs(Latitude),c(0,20,40,60,80))) %>%
  group_by(lat.fct) %>%
  summarise(median(Temperature,na.rm=T))

###########################
# Latitude vs Temperature #
###########################
ggplot(data=env.mat %>% mutate(Layer=fct_relevel(Layer,"SRF","DCM","MES")) %>% filter(Layer!="MES"),aes(x=Latitude,y=Temperature,col=Layer)) +
  geom_point(alpha=0.5,size=1) +
  #geom_line() +
  geom_smooth(se=F,span=0.4,size=0.5) +
  xlab("Latitude (deg)") +
  ylab("Temperature (ºC)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values=c("#0cb14c","#b82327","#136fba"))
ggsave("../results/FigSX_Lat_Temp.pdf",width=unit(5,"cm"),height=unit(2.5,"cm"))

ggplot(data=env.mat %>% mutate(Layer=fct_relevel(Layer,"SRF","DCM","MES")) %>% filter(Layer!="MES"),aes(x=cut(abs(Latitude),c(0,20,40,60,80)),y=Temperature,fill=Layer)) +
  geom_violin(scale = "width",draw_quantiles = 0.5) +
  #geom_point(alpha=0.5,size=1) +
  #geom_line() +
  #geom_smooth(se=F,span=0.4,size=0.5) +
  xlab("Absolute latitude (deg)") +
  ylab("Temperature (ºC)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#0cb14c","#b82327","#136fba"))
ggsave("../results/FigSX_Lat_Temp_fct.pdf",width=unit(5,"cm"),height=unit(2.5,"cm"))

##################################
## Correlation between datasets ##
##################################

toplot.df<-data.frame(metaG.dist=c(BCdist.metaG),metaT.dist=c(BCdist.metaT),tax.dist=c(BCdist.tax))

test<-mantel(BCdist.tax,BCdist.metaG)
pa<-ggplot(data=toplot.df,aes(x=tax.dist,y=metaG.dist)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  xlab("Euclidean distance (Taxonomic composition)") +
  ylab("Euclidean distance (Metagenomic composition)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

test<-mantel(BCdist.tax,BCdist.metaT)
pb<-ggplot(data=toplot.df,aes(x=tax.dist,y=metaT.dist)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  xlab("Euclidean distance (Taxonomic composition)") +
  ylab("Euclidean distance (Metatranscriptomic composition)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

test<-mantel(BCdist.metaG,BCdist.metaT)
pc<-ggplot(data=toplot.df,aes(x=metaG.dist,y=metaT.dist)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  xlab("Euclidean distance (Metagenomic composition)") +
  ylab("Euclidean distance (Metatranscriptomic composition)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

p<-pa | pb | pc
ggsave("../results/FigSX2_distances.pdf",width = unit(12,"cm"),height = unit(5,"cm"))


############################
# CORE-NONCORE exploration #
############################
# Load metaG 
metaG<-fread("zcat < ../data/processed/NOG_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG)<-metaG$V1
metaG<-metaG[,-1]
if (all(rownames(metaG)==rownames(metaG.norm))==F) stop("metaG and metaG.norm samples are not in the same order")
if (all(colnames(metaG)==colnames(metaG.norm))==F) stop("metaG and metaG.norm OGs are not in the same order")
occ<-apply(metaG,2,function(x){length(x[x>0])})
core<-base::factor(cut(occ,c(0,128,129),include.lowest = T))
levels(core)<-c("non core","core")

og.annot<-fread("../data/lib/NOG.annotations.tsv",sep="\t")
og.annot$V2<-gsub("ENOG41","",og.annot$V2)


scg<-fread("../data/lib/40_MGs.tsv") %>%
  mutate(SCG="YES") %>%
  rename(OG=COG)
tmp.df<-data.frame(OG=colnames(metaG.norm),mean.G=apply(metaG.norm,2,mean),variance.G=apply(metaG.norm,2,var),mean.T=apply(metaT.norm,2,mean),variance.T=apply(metaT.norm,2,var),mean.E=apply(exp.norm,2,mean),variance.E=apply(exp.norm,2,var)) %>%
  left_join(scg,by="OG") %>%
  left_join(data.frame(OG=colnames(metaG),occ,core),by="OG") %>%
  left_join(og.annot,by=c("OG"="V2"))

ggplot(data=tmp.df,aes(x=mean.G,y=variance.G,col=core)) +
  geom_point(alpha=0.5) +
  facet_wrap(~SCG) +
  theme_bw()

ggplot(data=tmp.df %>% arrange(core),aes(x=mean.G,y=variance.G,col=cut(occ,c(0,120,129),include.lowest = T))) +
  geom_point(alpha=0.2)

ggplot(data=tmp.df %>% arrange(core),aes(x=mean.G,y=variance.E,col=cut(occ,c(0,120,129),include.lowest = T))) +
  geom_point(alpha=0.2)


######################################################
## Correlation between datasets (core and non-core) ##
######################################################

# Compute distances
BCdist.metaG.core<-dist(metaG.norm[,core=="core"])
BCdist.metaG.noncore<-dist(metaG.norm[,core=="non core"])
BCdist.metaT.core<-dist(metaT.norm[,core=="core"])
BCdist.metaT.noncore<-dist(metaT.norm[,core=="non core"])

toplot.df<-data.frame(metaG.dist.core=c(BCdist.metaG.core),
                      metaG.dist.noncore=c(BCdist.metaG.noncore),
                      metaT.dist.core=c(BCdist.metaT.core),
                      metaT.dist.noncore=c(BCdist.metaT.noncore),
                      tax.dist=c(BCdist.tax))

# GC - GNC
test<-mantel(BCdist.metaG.core,BCdist.metaG.noncore)
p1<-ggplot(data=toplot.df,aes(x=metaG.dist.core,y=metaG.dist.noncore)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic core)") +
  ylab("Euclidean distance\n(Metagenomic non core)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# GC - TC
test<-mantel(BCdist.metaG.core,BCdist.metaT.core)
p2<-ggplot(data=toplot.df,aes(x=metaG.dist.core,y=metaT.dist.core)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic core)") +
  ylab("Euclidean distance\n(Metatranscriptomic core)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# GC - TNC
test<-mantel(BCdist.metaG.core,BCdist.metaT.noncore)
p3<-ggplot(data=toplot.df,aes(x=metaG.dist.core,y=metaT.dist.noncore)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic core)") +
  ylab("Euclidean distance\n(Metatranscriptomic non core)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# GC - tax
test<-mantel(BCdist.metaG.core,BCdist.tax)
p4<-ggplot(data=toplot.df,aes(x=metaG.dist.core,y=tax.dist)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic core)") +
  ylab("Euclidean distance\n(Taxonomic)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# GNC - TC
test<-mantel(BCdist.metaG.noncore,BCdist.metaT.core)
p5<-ggplot(data=toplot.df,aes(x=metaG.dist.noncore,y=metaT.dist.core)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic non core)") +
  ylab("Euclidean distance\n(Metatranscriptomic core)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# GNC - TNC
test<-mantel(BCdist.metaG.noncore,BCdist.metaT.noncore)
p6<-ggplot(data=toplot.df,aes(x=metaG.dist.noncore,y=metaT.dist.noncore)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic non core)") +
  ylab("Euclidean distance\n(Metatranscriptomic non core)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# GC - tax
test<-mantel(BCdist.metaG.noncore,BCdist.tax)
p7<-ggplot(data=toplot.df,aes(x=metaG.dist.noncore,y=tax.dist)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metagenomic non core)") +
  ylab("Euclidean distance\n(Taxonomic)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# TC - TNC
test<-mantel(BCdist.metaT.core,BCdist.metaT.noncore)
p8<-ggplot(data=toplot.df,aes(x=metaT.dist.core,y=metaT.dist.noncore)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metatranscriptomic core)") +
  ylab("Euclidean distance\n(Metatranscriptomic non core)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# TC - tax
test<-mantel(BCdist.metaT.core,BCdist.tax)
p9<-ggplot(data=toplot.df,aes(x=metaT.dist.core,y=tax.dist)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metatranscriptomic core)") +
  ylab("Euclidean distance\n(Taxonomic)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

# TNC - tax
test<-mantel(BCdist.metaT.noncore,BCdist.tax)
p10<-ggplot(data=toplot.df,aes(x=metaT.dist.noncore,y=tax.dist)) +
  geom_point(alpha=0.5,col="gray") +
  theme_bw() +
  xlab("Euclidean distance\n(Metatranscriptomic non core)") +
  ylab("Euclidean distance\n(Taxonomic)") +
  annotate(geom = 'text', label=paste("Mantel r: ",round(test$statistic,2),"\n","P-value < ",round(test$signif,3)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5)

p<-wrap_plots(list(plot_spacer(),p1,p2,p3,p4,
                   plot_spacer(),plot_spacer(),p5,p6,p7,
                   plot_spacer(),plot_spacer(),plot_spacer(),p8,p9,
                   plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),p10,
                   plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),ncol=5)
ggsave("../results/FigSX_core_noncore.pdf",p,width = unit(12,"cm"),height = unit(12,"cm"))


######################
## Analysis for EPI ##
######################
library(ggcorrplot)
# Select variables
vars<-c("HCO3","CO3","Carbon.total","Alkalinity.total","NO2","PO4","NO2NO3","Si","Temperature","Salinity","Density","Oxygen","NO3","ChlorophyllA","Fluorescence","PAR.TO","Iron.5m","Ammonium.5m","Gradient.Surface.temp(SST)","Okubo.Weiss","Lyapunov","Residence.time","Depth.Mixed.Layer","Brunt.Väisälä","Depth.Max.O2","Depth.Min.O2","Nitracline")

# Compute the mantel tests
multimantel<-function(distance,env.df,geo.dist){
  BCdist<-distance
  statistic<-NULL
  pval<-NULL
  n.obs<-NULL
  for (i in 1:ncol(env.df)){
    na.pos<-which(is.na(env.df[,i]))
    if (length(na.pos)>0) tmp<-mantel.partial(as.dist(as.matrix(BCdist)[-c(na.pos),-c(na.pos)]),dist(env.df[-c(na.pos),i]),as.dist(as.matrix(geo.dist)[-c(na.pos),-c(na.pos)]),method = "pearson",permutations = 1000) else tmp<-mantel.partial(BCdist,dist(env.df[,i]),geo.dist,method = "pearson",permutations = 1000)
    statistic<-c(statistic,tmp$statistic)
    pval<-c(pval,tmp$signif)
    n.obs<-c(n.obs,nrow(env.df)-length(na.pos))
  }
  data.frame(var=colnames(env.df),statistic,pval,p.corr=p.adjust(pval,method="bonferroni"),n.obs)
}

# miTags
res.tax<-multimantel(as.dist(as.matrix(BCdist.tax)[env.mat$epi=="EPI",env.mat$epi=="EPI"]),env.mat[env.mat$epi=="EPI",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="EPI",env.mat$epi=="EPI"]))
res.tax %>% arrange(abs(statistic))

# metaG
res.metaG<-multimantel(as.dist(as.matrix(BCdist.metaG)[env.mat$epi=="EPI",env.mat$epi=="EPI"]),env.mat[env.mat$epi=="EPI",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="EPI",env.mat$epi=="EPI"]))
res.metaG %>% arrange(abs(statistic))


# metaT
res.metaT<-multimantel(as.dist(as.matrix(BCdist.metaT)[env.mat$epi=="EPI",env.mat$epi=="EPI"]),env.mat[env.mat$epi=="EPI",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="EPI",env.mat$epi=="EPI"]))
res.metaT %>% arrange(abs(statistic))


# Expression
res.exp<-multimantel(as.dist(as.matrix(BCdist.exp)[env.mat$epi=="EPI",env.mat$epi=="EPI"]),env.mat[env.mat$epi=="EPI",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="EPI",env.mat$epi=="EPI"]))
res.exp %>% arrange(abs(statistic))

cor.mat<-data.frame(Taxonomy=res.tax$statistic,Function=res.metaG$statistic,Transcriptome=res.metaT$statistic,row.names = res.exp$var)
pcor.mat<-data.frame(Taxonomy=res.tax$pval,Function=res.metaG$p.corr,Expression=res.exp$pval,row.names = res.exp$var,Transcriptome=res.metaT$pval)
ordre<-order(apply(cor.mat[,1:3],1,mean),decreasing = T)
ggcorrplot(cor.mat[ordre,c(3,2,1)],p.mat=pcor.mat[ordre,c(3,2,1)],insig = "blank",sig.level = 0.05,method = "square",lab=T,lab_size = 2.5,colors=c("#2874b2","white","#ba2832"))
ggsave("../results/Fig2_Mantel_env_vars_EPI.pdf",width = unit(7,"cm"),height = unit(3.5,"cm"))

# Export it for Chris
fwrite(cor.mat,file="../results/cormat.bio_vs_env.tsv",quote = F,sep = "\t",row.names = T)
fwrite(pcor.mat,file="../results/pcormat.bio_vs_env.tsv",quote = F,sep = "\t",row.names = T)


# Construct the correlation plot
cor.mat<-cor(env.mat[env.mat$epi=="EPI",vars],use = "pairwise.complete.obs",method = "pearson")
cor.pmat<-cor_pmat(env.mat[env.mat$epi=="EPI",vars],use = "pairwise.complete.obs",method="pearson")
#hc<-hclust(as.dist(1-abs(cor.mat)))
ggcorrplot(cor.mat[ordre,ordre],p.mat = cor.pmat[ordre,ordre],insig = "blank",sig.level = 0.05,method = "square",lab=T,show.diag = F,lab_size = 0,colors=c("#2874b2","white","#ba2832"))
ggsave("../results/Fig2_env_vars_corr_EPI.pdf",width = unit(8,"cm"),height = unit(7,"cm"))

# Export it for Chris
fwrite(as.data.frame(cor.mat),file="../results/cormat.env_vs_env.tsv",quote = F,sep = "\t",row.names = T)
fwrite(as.data.frame(cor.pmat),file="../results/pcormat.env_vs_env.tsv",quote = F,sep = "\t",row.names = T)



######################
## Analysis for MES ##
######################

# Compute the mantel tests

# miTags
res.tax<-multimantel(as.dist(as.matrix(BCdist.tax)[env.mat$epi=="MES",env.mat$epi=="MES"]),env.mat[env.mat$epi=="MES",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="MES",env.mat$epi=="MES"]))
res.tax %>% arrange(abs(statistic))

# metaG
res.metaG<-multimantel(as.dist(as.matrix(BCdist.metaG)[env.mat$epi=="MES",env.mat$epi=="MES"]),env.mat[env.mat$epi=="MES",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="MES",env.mat$epi=="MES"]))
res.metaG %>% arrange(abs(statistic))


# metaT
res.metaT<-multimantel(as.dist(as.matrix(BCdist.metaT)[env.mat$epi=="MES",env.mat$epi=="MES"]),env.mat[env.mat$epi=="MES",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="MES",env.mat$epi=="MES"]))
res.metaT %>% arrange(abs(statistic))


# Expression
res.exp<-multimantel(as.dist(as.matrix(BCdist.exp)[env.mat$epi=="MES",env.mat$epi=="MES"]),env.mat[env.mat$epi=="MES",vars],as.dist(as.matrix(mar.dist.mat)[env.mat$epi=="MES",env.mat$epi=="MES"]))
res.exp %>% arrange(abs(statistic))

cor.mat<-data.frame(Taxonomy=res.tax$statistic,Function=res.metaG$statistic,Expression=res.exp$statistic,Transcriptome=res.metaT$statistic,row.names = res.exp$var)
pcor.mat<-data.frame(Taxonomy=res.tax$pval,Function=res.metaG$p.corr,Expression=res.exp$pval,row.names = res.exp$var,Transcriptome=res.metaT$pval)
cor.mat[is.na(cor.mat)]<-0
pcor.mat[is.na(pcor.mat)]<-0
vars.mes<-c("Temperature","Salinity","Density","Oxygen","NO3","Fluorescence","Gradient.Surface.temp(SST)","Okubo.Weiss","Lyapunov","Residence.time","Depth.Mixed.Layer","Brunt.Väisälä","Depth.Min.O2","Nitracline")
rownames(cor.mat)[which(rownames(cor.mat) %in% vars.mes==F)]

ggcorrplot(cor.mat[ordre,c(3,4,2,1)],p.mat=pcor.mat[ordre,c(3,4,2,1)],insig = "blank",sig.level = 0.05,method = "square",lab=T,lab_size = 2.5)
ggsave("../results/Fig2_Mantel_env_vars_MES.pdf",width = unit(7,"cm"),height = unit(3.5,"cm"))

