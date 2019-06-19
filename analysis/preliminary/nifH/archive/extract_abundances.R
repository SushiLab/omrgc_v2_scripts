library(data.table)
library(tidyverse)

target.genes<-fread("nifH_catalog.tab",sep="\t",header=F,data.table = F)
mgs<-fread("MGs.tab",sep="\t",header=F,data.table = F)

# Load the metaG abundance ############################################################################################
gene.ab<-fread("/nfs/cds/share/OceanMicrobiome_v2/OM-RGC_v2_gene_profile_metaG.tsv.gz",sep="\t",header=T,data.table = F)

# Extract MGs abundances
if (all(mgs$V2 %in% gene.ab$OMRGC_ID)==F) stop("Not all MG genes have been found")
mgs.ab<-gene.ab %>% filter(OMRGC_ID %in% mgs$V2)
fwrite(mgs.ab,file = "MG_metaG.tsv",sep="\t")

# Extract nifH abundances
if (all(target.genes$V2 %in% gene.ab$OMRGC_ID)==F) stop("Not all genes have been found")
target.genes.ab<-gene.ab %>% filter(OMRGC_ID %in% target.genes$V2)
fwrite(target.genes.ab,file = "nifH_metaG.tsv",sep="\t")

# Load the metaT abundance ############################################################################################
gene.ab<-NULL
gene.ab<-fread("/nfs/cds/share/OceanMicrobiome_v2/OM-RGC_v2_gene_profile_metaT.tsv.gz",sep="\t",header=T,data.table = F)

# Extract MGs abundances
if (all(mgs$V2 %in% gene.ab$OMRGC_ID)==F) stop("Not all MG genes have been found")
mgs.ab<-gene.ab %>% filter(OMRGC_ID %in% mgs$V2)
fwrite(mgs.ab,file = "MG_metaT.tsv",sep="\t")

# Extract nifH abundances
if (all(target.genes$V2 %in% gene.ab$OMRGC_ID)==F) stop("Not all genes have been found")
target.genes.ab<-gene.ab %>% filter(OMRGC_ID %in% target.genes$V2)
fwrite(target.genes.ab,file = "nifH_metaT.tsv",sep="\t")

