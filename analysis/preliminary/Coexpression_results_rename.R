library(data.table)
library(fastmatch)
cor.mat<-fread("../results/Coexpression_pairs.tsv",header=T,sep="\t",data.table = F)

gene.cat<-fread("zcat /nfs/cds/OceanMicrobiome_v2/OM-RGC_v2_data/OM-RGC_v2_release.tsv.gz",header=T,sep="\t",data.table = F)

pos<-fmatch(cor.mat$`GC/OG representative 1`,gene.cat$gene)
all(is.na(pos)==FALSE)
cor.mat$`GC/OG representative 1`<-gene.cat$OMRGC_ID[pos]

pos<-which(grepl("TARA",cor.mat$`GC/OG representative 2`))
pos2<-fmatch(cor.mat$`GC/OG representative 2`[pos],gene.cat$gene)
all(is.na(pos2)==FALSE)
cor.mat$`GC/OG representative 2`[pos]<-gene.cat$OMRGC_ID[pos2]
cor.mat<-cor.mat[,-4]
colnames(cor.mat)<-c("Representative 1","Representative 2","Co-expression Pearson r","Representative 2 type","Representative 2 annotation")
fwrite(cor.mat,file="../results/Coexpression_pairs_renamed.tsv",sep="\t",na = NA,quote = F)