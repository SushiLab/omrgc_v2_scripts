library(DESeq2)
library(data.table)
library(EcolUtils)
library(ggplot2)
library(tidyverse)

logeucdist<-function(mat){
  mat<-mat/rowSums(mat)
  pseudocount<-1*10^(floor(log10(min(mat[mat>0]))))
  cat("Default pseudocount: floor the minimum non-zero value to the same number of digits\nPseudocount used: ",pseudocount)
  dist(log10(mat+pseudocount))
}

# Load metaG
metaG.norm<-fread("../data/processed/NOG_metaG.norm.txt",header=T,sep="\t",data.table = F)
rownames(metaG.norm)<-metaG.norm$V1
metaG.norm<-metaG.norm[,-1]
env.mat.metaG<-fread("../data/processed/NOG_env.mat.metaG.txt",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.metaG)<-env.mat.metaG$V1
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

# Load metaT
metaT.norm<-fread("../data/processed/NOG_metaT.norm.txt",header=T,sep="\t",data.table = F)
rownames(metaT.norm)<-metaT.norm$V1
metaT.norm<-metaT.norm[,-1]
env.mat.metaT<-fread("../data/processed/NOG_env.mat.metaT.txt",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.metaT)<-env.mat.metaT$V1
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

# Load ratio
ratio.norm<-fread("../data/processed/NOG_ratio.mat.txt",header=T,sep="\t",data.table = F)
rownames(ratio.norm)<-ratio.norm$V1
ratio.norm<-ratio.norm[,-1]
env.mat.ratio<-fread("../data/processed/NOG_env.mat.match.txt",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.ratio)<-env.mat.ratio$V1
env.mat.ratio<-env.mat.ratio[,-1]
env.mat.ratio<-env.mat.ratio[match(rownames(ratio.norm),env.mat.ratio$Barcode),]

# Load tax
miTags<-fread("../data/processed/miTags/OTU.tab.rr.TARA180.noeuks.txt",sep="\t",header=T,data.table = F)
rownames(miTags)<-miTags$V1
miTags<-miTags[,-1]

# Load environmental data
env.mat<-fread("../data/processed/NOG_env.mat.metaG.txt",sep="\t",header=T,data.table = F,stringsAsFactors = T)
env.mat$Layer<-factor(env.mat$Layer,levels = c("SRF","DCM","MES","other"))

miTags<-miTags[match(gsub("-",".",env.mat$Sample_name),rownames(miTags)),]
wTINA<-as.matrix(logeucdist(miTags))


w.size<-10

## LATITUDE ##
# SRF
res.metaG<-smwda(dist(metaG.norm.log2[which(env.mat.metaG$Layer=="SRF"),]),env.mat.metaG$Latitude[which(env.mat.metaG$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.metaT<-smwda(dist(metaT.norm.log2[which(env.mat.metaT$Layer=="SRF"),]),env.mat.metaT$Latitude[which(env.mat.metaT$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.ratio<-smwda(dist(ratio.norm[which(env.mat.ratio$Layer=="SRF"),]),env.mat.ratio$Latitude[which(env.mat.ratio$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.tax.SRF<-smwda(wTINA[which(env.mat$Layer=="SRF"),which(env.mat$Layer=="SRF")],env.mat$Latitude[which(env.mat$Layer=="SRF")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))

res.SRF<-rbind(res.metaG$windows,res.metaT$windows,res.ratio$windows,res.tax.SRF$windows)
res.SRF$data.source<-c(rep("metaG",nrow(res.metaG$windows)),rep("metaT",nrow(res.metaT$windows)),rep("Expression",nrow(res.ratio$windows)),rep("miTags",nrow(res.tax.SRF$windows)))

# DCM
res.metaG<-smwda(dist(metaG.norm.log2[which(env.mat.metaG$Layer=="DCM"),]),env.mat.metaG$Latitude[which(env.mat.metaG$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.metaT<-smwda(dist(metaT.norm.log2[which(env.mat.metaT$Layer=="DCM"),]),env.mat.metaT$Latitude[which(env.mat.metaT$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.ratio<-smwda(dist(ratio.norm[which(env.mat.ratio$Layer=="DCM"),]),env.mat.ratio$Latitude[which(env.mat.ratio$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))
res.tax.DCM<-smwda(wTINA[which(env.mat$Layer=="DCM"),which(env.mat$Layer=="DCM")],env.mat$Latitude[which(env.mat$Layer=="DCM")],w.size = w.size,nrep = 10000,probs = c(0.005,0.995))

res.DCM<-rbind(res.metaG$windows,res.metaT$windows,res.ratio$windows,res.tax.DCM$windows)
res.DCM$data.source<-c(rep("metaG",nrow(res.metaG$windows)),rep("metaT",nrow(res.metaT$windows)),rep("Expression",nrow(res.ratio$windows)),rep("miTags",nrow(res.tax.DCM$windows)))

# Merge
res<-rbind(res.SRF,res.DCM)
res$Layer<-c(rep("SRF",nrow(res.SRF)),rep("DCM",nrow(res.DCM)))
res$Layer<-factor(res$Layer,levels=c("SRF","DCM"))
res$data.source<-factor(res$data.source,levels=c("miTags","metaG","metaT","Expression"))

# Plot
ggplot(data=res,aes(x=env.var.mean,y=stat.real.zscore,col=data.source)) +
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
  scale_color_manual(values=c("gray","#ed9406","#33ccff","#c6e9af"),labels=c("Taxonomic composition","Genomic composition","Transcriptomic composition","Expression")) +
  scale_x_continuous(breaks=c(-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80)) +
  coord_flip()
  

ggsave(file="../results/Fig2A.pdf",width = 210,height = 140,units = "mm")


