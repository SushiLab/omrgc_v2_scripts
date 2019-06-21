library(data.table)
library(dplyr)
library(tidyverse)

cat("Loading gene abundances (metaG)\n")
gene.ab<-fread("/nfs/cds/OceanMicrobiome_v2/OM-RGC_v2/profiles/tara.all.prok.180.gene.profile.screened.tara.adapters.on.OM-RGCv2.95.padded.solexaqa.allbest.l45.p95.base.mm.dist.among.unique.norm.gene",header=T,data.table = F)
gene.ab<-gene.ab %>% # remove the first row (-1)
  slice(2:n())
cat("Done\n\n")

rownames(gene.ab)<-gene.ab$gene
gene.ab<-gene.ab[,-c(1)]

# Load environmental data
env.mat<-fread("zcat /nfs/cds/scratch/guillems/Thermogenomics_release/data/processed/KO_env.mat.metaG.txt.gz",sep="\t", header=T)

all(gsub("_G","",colnames(gene.ab)) %in% env.mat$Sample_name)

env.mat.red<-env.mat[match(gsub("_G","",colnames(gene.ab)),env.mat$Sample_name),]

n.genes<-NULL
for (i in 1:nrow(env.mat.red)){
  cat(i,"\n")
  if (i==1) n.genes<-c(n.genes,length(which(gene.ab[,1:i]>0))) else n.genes<-c(n.genes,length(which(rowSums(gene.ab[,1:i])>0)))
}
names(n.genes)<-colnames(gene.ab)

fwrite(data.frame(sample=names(n.genes),n.genes),file="accum.genes.tsv",sep="\t",quote=F,row.names = F)

