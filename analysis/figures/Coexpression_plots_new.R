library(propagate)
library(vegan)
library(tidyverse)
library(data.table)
library(patchwork)
#source("http://sunagawalab.ethz.ch/sushilab/resources/R/sushi.palette.R")
source("../data/lib/sushipal.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])

###############
## LOAD DATA ##
###############

metaG.norm.match.log2<-fread("zcat < ../data/processed/NOGplusGF_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("zcat < ../data/processed/NOGplusGF_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("zcat < ../data/processed/NOGplusGF_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("zcat < ../data/processed/NOGplusGF_metaT.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("zcat < ../data/processed/NOGplusGF_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]
env.mat.match<-fread("zcat < ../data/processed/NOGplusGF_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]

##############
## COG info ##
##############
cog.info<-read.table("../data/lib/cognames2003-2014.tab",sep="\t",quote="",comment.char="",header=T)
#cog.10MGs<-read.table("../data/lib/10MG_COG_ids",sep="\t",quote="",comment.char="",header=F)
#cog.10MGs<-cog.info[match(cog.10MGs$V1,cog.info$X..COG),]
nog.info<-fread("../data/lib/NOG.annotations.tsv",sep="\t",header=F,data.table = F)
nog.info<-nog.info[,c(2,5,6)]
nog.info<-nog.info %>%
  #filter(!grepl("COG",V2)) %>%
  mutate(V2=gsub("ENOG41","",V2))%>%
  dplyr::rename(OG=V2) %>%
  dplyr::rename(annotation=V6) %>%
  dplyr::rename(func=V5)
pos<-which((colnames(metaG.norm.match) %in% nog.info$OG)==F)
tmp2<-data.frame(OG=colnames(metaG.norm.match)[pos],func=rep("S",length(colnames(metaG.norm.match)[pos])),annotation=rep(NA,length(colnames(metaG.norm.match)[pos])))
annot<-bind_rows(tmp2,nog.info) %>%
  filter(OG %in% colnames(metaG.norm.match))
annot$known.status<-rep("Known known",nrow(annot))
annot$known.status[grepl("TARA",annot$OG)]<-"Unknown unknown (new GC)"
annot$known.status[which(!grepl("TARA",annot$OG) & is.na(annot$annotation))]<-"Known unknown"
annot$source<-"OG"
annot$source[grepl("TARA",annot$OG)]<-"GC"

# Overwrite annotation from COG
pos<-match(cog.info$X..COG,annot$OG)
na.pos<-which(is.na(pos))
annot$annotation[pos[-c(na.pos)]]<-as.character(cog.info$name[-c(na.pos)])

ggplot(data=annot,aes(x=source,fill=known.status)) +
  geom_bar()

###################
## Co-expression ##
###################
occ<-function(x){length(x[x>0])/length(x)}
metaG.norm.match.occ<-apply(metaG.norm.match,2,occ)
metaG.norm.match.log2.red<-metaG.norm.match.log2[,metaG.norm.match.occ>=0.1] # Change the occurrence threshold
annot.red<-annot[annot$OG %in% colnames(metaG.norm.match.log2.red),]

cor.mat<-fread("zcat < ../data/processed/cor.mat_0.1occ_0.6rcutoff.tsv.gz",sep="\t",data.table = F)

# Filter based on occurrence
cor.mat<-cor.mat %>%
  filter(x %in% colnames(metaG.norm.match.log2.red)) %>%
  filter(y %in% colnames(metaG.norm.match.log2.red))

# Get he best correlated for each OG/GC
tmp<-data.frame(x=cor.mat$y,y=cor.mat$x,r=cor.mat$r)
cor.mat<-bind_rows(cor.mat,tmp) %>%
  group_by(x) %>%
  filter(r==max(r)) %>%
  as.data.frame()
  
# Filter by r value
cor.mat<-cor.mat %>% filter(r>=0.86)

# Add annotation
cor.mat<-cor.mat %>%
  left_join(annot.red,by=c("x"="OG")) %>%
  left_join(annot.red,by=c("y"="OG"))

cor.mat$compar.known<-interaction(cor.mat$known.status.x,cor.mat$known.status.y,sep="-")
cor.mat$compar.source<-interaction(cor.mat$source.x,cor.mat$source.y,sep="-")

# Remove duplicated pairs
#cor.mat<-arrange(cor.mat,desc(r))
#cor.mat$pair<-apply(cor.mat,1,function(x){paste(sort(c(x[1],x[2])),collapse="-")})
#dupl<-duplicated(cor.mat$pair)
#cor.mat<-cor.mat[-c(which(dupl)),]

# Statistics
total.gc<-table(annot.red$source)[1]
total.og<-table(annot.red$source)[2]
gc.to.og<-length(which(cor.mat$source.x=="GC" & cor.mat$source.y=="OG"))
gc.to.og.kk<-length(which(cor.mat$source.x=="GC" & cor.mat$source.y=="OG" & cor.mat$known.status.y=="Known known"))
gc.to.og.ku<-length(which(cor.mat$source.x=="GC" & cor.mat$source.y=="OG" & cor.mat$known.status.y=="Known unknown"))
gc.to.gc<-length(which(cor.mat$source.x=="GC" & cor.mat$source.y=="GC"))
cat("Total OGs: ",total.og,"\nTotal GCs: ",total.gc,"\nGCs linked to OGs: ",gc.to.og," perc: ",100*gc.to.og/total.gc,"\nGCs linked to GCs: ",gc.to.gc," perc: ",100*gc.to.gc/total.gc,"\nGCs linked to known known OGs: ",gc.to.og.kk,"\nGCs linked to known unknown OGs: ",gc.to.og.ku)

# Save results
to.save<-cor.mat %>%
  filter(source.x=="GC" & source.y=="OG" | source.x=="GC" & source.y=="GC") %>%
  select(x,y,r,source.x,source.y,annotation.y) %>%
  arrange(desc(r))
colnames(to.save)<-c("GC/OG representative 1","GC/OG representative 2","Co-expressiom Pearson r","Representative 1 type","Representative 2 type","Representative 2 annotation")

fwrite(to.save,file="../results/Coexpression_pairs.tsv",sep="\t",na = NA)
