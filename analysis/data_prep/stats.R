library(data.table)
library(dplyr)
library(tidyverse)

# Load the gene catalog
cat("Loading the gene catalog\n")
gene.cat<-fread("zcat OM-RGC_v2_release.tsv.gz",header=T,data.table = F)
cat("Done\n\n")

# Compute the number of functionally annotated genes
n.tot<-nrow(gene.cat)
n.ko<-length(which(gene.cat$KO!=""))
n.NOG<-length(which(gene.cat$NOG!=""))
n.GF<-length(which(gene.cat$GF!=""))

# Compute the number of taxonomically annotated genes
n.Domain<-gene.cat %>% group_by(Domain) %>% summarise(n=n()) %>% as.data.frame()
n.Phylum<-gene.cat %>% group_by(Phylum) %>% summarise(n=n()) %>% as.data.frame()
n.Class<-gene.cat %>% group_by(Class) %>% summarise(n=n()) %>% as.data.frame()

# Save results
fwrite(data.frame(annotation=c("gene","KO","NOG","GF"),n=c(n.tot,n.ko,n.NOG,n.GF)),file="./stats/func.annotation.stats",nThread=1,sep="\t",quote = F,row.names = F)
fwrite(n.Domain,file="./stats/tax.annotation.Domain.stats",nThread=1,sep="\t",quote = F,row.names = F)
fwrite(n.Phylum,file="./stats/tax.annotation.Phylum.stats",nThread=1,sep="\t",quote = F,row.names = F)
fwrite(n.Class,file="./stats/tax.annotation.Class.stats",nThread=1,sep="\t",quote = F,row.names = F)

