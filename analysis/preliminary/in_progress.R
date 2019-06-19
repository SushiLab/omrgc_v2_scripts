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

#########################
# Mean betadiversity

tmp<-as.matrix(BCdist.tax)
diag(tmp)<-NA
tmp[lower.tri(tmp)]<-NA
  
tmp<-tmp %>% as.data.frame() %>%
  rownames_to_column(var="Sample1") %>%
  gather(key = "Sample2",value="value",-Sample1) %>%
  filter(!is.na(value))

tmp$epi1<-env.mat$epi[match(tmp$Sample1,gsub("-",".",env.mat$Sample_name))]
tmp$epi2<-env.mat$epi[match(tmp$Sample2,gsub("-",".",env.mat$Sample_name))]
tmp$polar1<-env.mat$polar[match(tmp$Sample1,gsub("-",".",env.mat$Sample_name))]
tmp$polar2<-env.mat$polar[match(tmp$Sample2,gsub("-",".",env.mat$Sample_name))]
tmp$epi<-interaction(tmp$epi1,tmp$epi2)
tmp$polar<-interaction(tmp$polar1,tmp$polar2)


ggplot(data=tmp %>% filter(epi=="EPI.EPI" & polar %in% c("polar.polar","non polar.non polar")),aes(x=polar,y=value)) +
  geom_violin(scale = "width",draw_quantiles = 0.5,fill="gray") +
  theme_bw()

t.test(tmp$value[tmp$epi=="EPI.EPI" & tmp$polar=="polar.polar"],tmp$value[tmp$epi=="EPI.EPI" & tmp$polar=="non polar.non polar"])

#########################
plot(BCdist.tax,BCdist.metaG)
mantel(BCdist.tax,BCdist.metaG)


plot(dist(env.mat$Latitude[env.mat.match$epi=="EPI"]),as.dist(as.matrix(BCdist.tax)[env.mat$epi=="EPI",env.mat$epi=="EPI"]))

