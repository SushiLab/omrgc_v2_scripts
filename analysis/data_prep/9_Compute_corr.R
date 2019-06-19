library(vegan)
library(tidyverse)
library(data.table)

###############
## LOAD DATA ##
###############

metaG.norm.match.log2<-fread("zcat ../data/processed/NOGplusGF_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("zcat ../data/processed/NOGplusGF_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("zcat ../data/processed/NOGplusGF_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("zcat ../data/processed/NOGplusGF_metaT.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("zcat ../data/processed/NOGplusGF_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]
env.mat.match<-fread("zcat ../data/processed/NOGplusGF_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]

############################
# Compute the correlations #
############################
source("../data/lib/cor.cutoff.R")
occ<-function(x){length(x[x>0])/length(x)}

ratio.mat.red<-ratio.mat[,apply(metaG.norm.match,2,occ)>=0.1] # filter by occurrence

cor.mat<-cor.cutoff(ratio.mat.red,r.min = 0.6,blocksize = 1000) # compute correlation
fwrite(cor.mat,file = "../data/processed/cor.mat_0.1occ_0.6rcutoff.tsv",sep="\t",nThread = 1)
