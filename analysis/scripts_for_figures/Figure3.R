# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure 3: ecological differentiation ===================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(DESeq2) # BiocManager::install("DESeq2")
library(data.table)
library(EcolUtils) # devtools::install_github("GuillemSalazar/EcolUtils")
library(ggplot2)
library(tidyverse)

# Functions ------------------------------------------------------------------------------

logeucdist<-function(mat){
  mat<-mat/rowSums(mat)
  pseudocount<-1*10^(floor(log10(min(mat[mat>0]))))
  cat("Default pseudocount: floor the minimum non-zero value to the same number of digits\nPseudocount used: ",pseudocount)
  dist(log10(mat+pseudocount))
}

# Load data ------------------------------------------------------------------------------

# meta g
metaG.norm<-fread("../../data/OM-RGC_v2_functional_profile_eggNOG/OG_metaG.norm.tsv.gz",header=T,sep="\t",data.table = F)
env.mat.metaG<-fread("../../data/OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_metaG.tsv",header=T,sep="\t",data.table = F,stringsAsFactors = T)

# meta t
metaT.norm<-fread("../../data/OM-RGC_v2_functional_profile_eggNOG/OG_metat.norm.tsv.gz",header=T,sep="\t",data.table = F)
env.mat.metaT<-fread("../../data/OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_metaT.tsv",header=T,sep="\t",data.table = F,stringsAsFactors = T)

# ratios
ratio.norm<-fread("../../data/OM-RGC_v2_functional_profile_eggNOG/OG_exp.norm.match.log2.tsv.gz",header=T,sep="\t",data.table = F)
env.mat.ratio<-fread("../../data/OM-RGC_v2_auxiliary_data/OM-RGC_v2_auxiliary_data_matched.tsv",header=T,sep="\t",data.table = F,stringsAsFactors = T)

# taxo
miTags<-fread("../../data/OM-RGC_v2_taxonomic_profiles/mitags_tab_otu_rr.tsv",sep="\t",header=T,data.table = F)


# Transform data --------------------------------------------------------------------------

# meta g
rownames(metaG.norm)<-metaG.norm$V1
metaG.norm<-metaG.norm[,-1]
rownames(env.mat.metaG)<-env.mat.metaG$`PANGAEA sample id`
env.mat.metaG<-env.mat.metaG[,-1]
input.data.merged<-t(metaG.norm)
colnames(input.data.merged)<-paste(colnames(input.data.merged),"_G",sep="")
input.metadata.merged<-env.mat.metaG
to.integer<-function(taula){taula<-round(taula/(max(taula)/1e9))}
input.data.merged<-to.integer(input.data.merged)
rownames(input.metadata.merged)<-colnames(input.data.merged)
dds <- DESeqDataSetFromMatrix(countData = input.data.merged,colData = input.metadata.merged,design = ~ 1)
vsd <- vst(dds, blind=F)
res<-assay(vsd)
input.data.norm<-res[,1:nrow(metaG.norm)]
metaG.norm.log2<-t(input.data.norm)+log2((max(input.data.merged)/1e9))
rownames(metaG.norm.log2)<-gsub("_G","",rownames(metaG.norm))

# meta t
rownames(metaT.norm)<-metaT.norm$V1
metaT.norm<-metaT.norm[,-1]
rownames(env.mat.metaT)<-env.mat.metaT$`PANGAEA sample id`
env.mat.metaT<-env.mat.metaT[,-1]
input.data.merged<-t(metaT.norm)
colnames(input.data.merged)<-paste(colnames(input.data.merged),"_G",sep="")
input.metadata.merged<-env.mat.metaT
to.integer<-function(taula){taula<-round(taula/(max(taula)/1e9))}
input.data.merged<-to.integer(input.data.merged)
rownames(input.metadata.merged)<-colnames(input.data.merged)
dds <- DESeqDataSetFromMatrix(countData = input.data.merged,colData = input.metadata.merged,design = ~ 1)
vsd <- vst(dds, blind=F)
res<-assay(vsd)
input.data.norm<-res[,1:nrow(metaT.norm)]
metaT.norm.log2<-t(input.data.norm)+log2((max(input.data.merged)/1e9))
rownames(metaT.norm.log2)<-gsub("_G","",rownames(metaT.norm))

# ratios
rownames(ratio.norm)<-ratio.norm$V1
ratio.norm<-ratio.norm[,-1]
rownames(env.mat.ratio)<-env.mat.ratio$`PANGAEA sample id`
env.mat.ratio<-env.mat.ratio[,-1]


# Match tables
rownames(miTags)<-miTags$V1
miTags<-miTags[,-1]

metaG.norm.log2<-metaG.norm.log2[match(rownames(env.mat.metaG),rownames(metaG.norm.log2)),]
metaT.norm.log2<-metaT.norm.log2[match(rownames(env.mat.metaT),rownames(metaT.norm.log2)),]
ratio.norm<-ratio.norm[match(rownames(env.mat.ratio),rownames(ratio.norm)),]
miTags<-miTags[match(rownames(env.mat.metaG),rownames(miTags)),]

# Define window size for the sliding window analysis
w.size<-10

## SMWDA ##
# SRF
res.metaG<-smwda(dist(metaG.norm.log2[which(env.mat.metaG$Layer=="SRF"),]),env.mat.metaG$Latitude[which(env.mat.metaG$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.metaT<-smwda(dist(metaT.norm.log2[which(env.mat.metaT$Layer=="SRF"),]),env.mat.metaT$Latitude[which(env.mat.metaT$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.ratio<-smwda(dist(ratio.norm[which(env.mat.ratio$Layer=="SRF"),]),env.mat.ratio$Latitude[which(env.mat.ratio$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.tax.SRF<-smwda(logeucdist(miTags[which(env.mat.metaG$Layer=="SRF"),]),env.mat.metaG$Latitude[which(env.mat.metaG$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))

res.SRF<-rbind(res.metaG$windows,res.metaT$windows,res.ratio$windows,res.tax.SRF$windows)
res.SRF$data.source<-c(rep("metaG",nrow(res.metaG$windows)),rep("metaT",nrow(res.metaT$windows)),rep("Expression",nrow(res.ratio$windows)),rep("miTags",nrow(res.tax.SRF$windows)))

# DCM
res.metaG<-smwda(dist(metaG.norm.log2[which(env.mat.metaG$Layer=="DCM"),]),env.mat.metaG$Latitude[which(env.mat.metaG$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.metaT<-smwda(dist(metaT.norm.log2[which(env.mat.metaT$Layer=="DCM"),]),env.mat.metaT$Latitude[which(env.mat.metaT$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.ratio<-smwda(dist(ratio.norm[which(env.mat.ratio$Layer=="DCM"),]),env.mat.ratio$Latitude[which(env.mat.ratio$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.tax.DCM<-smwda(logeucdist(miTags[which(env.mat.metaG$Layer=="DCM"),]),env.mat.metaG$Latitude[which(env.mat.metaG$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))

res.DCM<-rbind(res.metaG$windows,res.metaT$windows,res.ratio$windows,res.tax.DCM$windows)
res.DCM$data.source<-c(rep("metaG",nrow(res.metaG$windows)),rep("metaT",nrow(res.metaT$windows)),rep("Expression",nrow(res.ratio$windows)),rep("miTags",nrow(res.tax.DCM$windows)))

# Merge
res<-rbind(res.SRF,res.DCM)
res$Layer<-c(rep("SRF",nrow(res.SRF)),rep("DCM",nrow(res.DCM)))
res$Layer<-factor(res$Layer,levels=c("SRF","DCM"))
res$data.source<-factor(res$data.source,levels=c("miTags","metaG","metaT","Expression"))

# Plot
res<-res %>% filter(data.source!="Expression")
p<-ggplot(data=res,aes(x=env.var.mean,y=stat.real.zscore,col=data.source)) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_point(aes(shape=sign),size=2) +
  geom_errorbarh(data=res %>% filter(sign=="sign"),aes(x=env.var.mean,y=stat.real.zscore,col=data.source,xmin=env.var.min,xmax=env.var.max),alpha=0.5,linetype=1,height=0.2) +
  #geom_line() +
  geom_smooth(span=0.4,se=F) +
  scale_shape_manual(values=c(3,19)) +
  theme_bw() +
  facet_grid(.~Layer) +
  xlab("Latitude") +
  ylab("Differentiation index") +
  scale_color_manual(values=c("gray","#ed9406","#33ccff"),labels=c("Taxonomic composition","Metagenomic composition","Metatranscriptomic composition")) +
  scale_x_continuous(breaks=c(-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80)) +
  coord_flip()
  

ggsave(file="../../results/figures/Figure3.pdf",p,width = 210,height = 140,units = "mm")


