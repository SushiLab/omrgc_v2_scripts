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
ggsave("../../results/figures/SFigure5.pdf",width = unit(12,"cm"),p,height = unit(5,"cm"))

