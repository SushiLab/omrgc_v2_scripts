# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure 2A: Beta diversity ==============================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(patchwork)
source("../lib/sushipal.R")
source("../lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])
library(vegan)

logeucdist<-function(mat){
  mat<-mat/rowSums(mat)
  pseudocount<-1*10^(floor(log10(min(mat[mat>0]))))
  cat("Default pseudocount: floor the minimum non-zero value to the same number of digits\nPseudocount used: ",pseudocount)
  dist(log10(mat+pseudocount))
}

# Load data ------------------------------------------------------------------------------

# Load environmental data
env.mat<-fread("../../data/OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_matched.tsv",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat)<-env.mat$`PANGAEA sample id`
env.mat<-env.mat[,-1]
env.mat$epi<-env.mat$Layer
levels(env.mat$epi)<-c("EPI","MES","EPI","EPI")
levels(env.mat$Layer)<-c("DCM","MES","DCM","SRF")
env.mat$Absolute.Latitude<-abs(env.mat$Latitude)

# Load marine distances
mar.dist.mat<-fread("../../data/marine_distances/mardist.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(mar.dist.mat)<-mar.dist.mat$V1
mar.dist.mat<-mar.dist.mat[,-1]
mar.dist.mat<-mar.dist.mat[match(rownames(env.mat),rownames(mar.dist.mat)),match(rownames(env.mat),rownames(mar.dist.mat))]

# metaG
metaG.norm<-fread("../../data/OM-RGC_v2_functional_profile_eggNOG/OG_metaG.norm.match.log2.tsv.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm)<-metaG.norm$V1
metaG.norm<-metaG.norm[,-1]
metaG.norm<-metaG.norm[match(rownames(env.mat),rownames(metaG.norm)),]
BCdist.metaG<-dist(metaG.norm)

# metaT
metaT.norm<-fread("../../data/OM-RGC_v2_functional_profile_eggNOG/OG_metaT.norm.match.log2.tsv.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm)<-metaT.norm$V1
metaT.norm<-metaT.norm[,-1]
metaT.norm<-metaT.norm[match(rownames(env.mat),rownames(metaT.norm)),]
BCdist.metaT<-dist(metaT.norm)

# Expression
exp.norm<-fread("../../data/OM-RGC_v2_functional_profile_eggNOG/OG_exp.norm.match.log2.tsv.gz",header=T,sep="\t",data.table = F)
rownames(exp.norm)<-exp.norm$V1
exp.norm<-exp.norm[,-1]
exp.norm<-exp.norm[match(rownames(env.mat),rownames(exp.norm)),]
BCdist.exp<-dist(exp.norm)

# miTags
otu.table<-fread("../../data/OM-RGC_v2_taxonomic_profiles/mitags_tab_otu_rr.tsv",sep="\t",header=T,data.table = F)
rownames(otu.table)<-otu.table$V1
otu.table<-otu.table[,-1]
otu.table<-otu.table[match(sapply(strsplit(rownames(env.mat),"-"),"[[",1),rownames(otu.table)),]
BCdist.tax<-logeucdist(otu.table)

######################
## Analysis for EPI ##
######################
library(ggcorrplot)
# Select variables
vars<-c("HCO3","CO3","Carbon.total","Alkalinity.total","NO2","PO4","NO2NO3","Si","Temperature","Salinity","Density","Oxygen","NO3","ChlorophyllA","Fluorescence","PAR.PC","Iron.5m","Ammonium.5m","Gradient.Surface.temp(SST)","Okubo.Weiss","Lyapunov","Residence.time","Depth.Mixed.Layer","Brunt.VÃ¤isÃ¤lÃ¤","Depth.Max.O2","Depth.Min.O2","Nitracline")

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
p<-ggcorrplot(cor.mat[ordre,c(3,2,1)],p.mat=pcor.mat[ordre,c(3,2,1)],insig = "blank",sig.level = 0.05,method = "square",lab=T,lab_size = 2.5,colors=c("#2874b2","white","#ba2832"))
ggsave("../../results/figures/Figure4AI.pdf",p,width = unit(7,"cm"),height = unit(3.5,"cm"))


# Construct the correlation plot
cor.mat<-cor(env.mat[env.mat$epi=="EPI",vars],use = "pairwise.complete.obs",method = "pearson")
cor.pmat<-cor_pmat(env.mat[env.mat$epi=="EPI",vars],use = "pairwise.complete.obs",method="pearson")
#hc<-hclust(as.dist(1-abs(cor.mat)))
p<-ggcorrplot(cor.mat[ordre,ordre],p.mat = cor.pmat[ordre,ordre],insig = "blank",sig.level = 0.05,method = "square",lab=T,show.diag = F,lab_size = 0,colors=c("#2874b2","white","#ba2832"))
ggsave("../../results/figures/Figure4AII.pdf",p,width = unit(8,"cm"),height = unit(7,"cm"))
