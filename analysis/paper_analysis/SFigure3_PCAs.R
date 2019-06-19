library(tidyverse)
library(data.table)
library(patchwork)
source("../data/lib/sushipal.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])
library(vegan)

find_hull<-function(dat,x.name,y.name,lev.name){
  hull<-NULL
  for (i in levels(as.factor(dat[,lev.name]))){
    X<-na.exclude(dat[,x.name][dat[,lev.name]==i])
    Y<-na.exclude(dat[,y.name][dat[,lev.name]==i])
    if(length(X)>=1 & length(Y)>=1) pos<-chull(X,Y)
    if(length(X)>=1 & length(Y)>=1) hull<-rbind(hull,data.frame(x=X[pos],y=Y[pos],lev=rep(i,length(pos))))
  }
  colnames(hull)<-c(x.name,y.name,lev.name)
  hull
}

# Load environmental data match
env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match$epi<-env.mat.match$Layer
levels(env.mat.match$epi)<-c("EPI","MES","EPI","EPI")
levels(env.mat.match$Layer)<-c("DCM","MES","DCM","SRF")
env.mat.match$Absolute.Latitude<-abs(env.mat.match$Latitude)

# Load environmental data metaG
env.mat.metaG<-fread("zcat < ../data/processed/NOG_env.mat.metaG.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.metaG)<-env.mat.metaG$V1
env.mat.metaG<-env.mat.metaG[,-1]
env.mat.metaG$epi<-env.mat.metaG$Layer
levels(env.mat.metaG$epi)<-c("EPI","MES","EPI","EPI")
levels(env.mat.metaG$Layer)<-c("DCM","MES","DCM","SRF")
env.mat.metaG$Absolute.Latitude<-abs(env.mat.metaG$Latitude)

# Load environmental data metaT
env.mat.metaT<-fread("zcat < ../data/processed/NOG_env.mat.metaT.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.metaT)<-env.mat.metaT$V1
env.mat.metaT<-env.mat.metaT[,-1]
env.mat.metaT$epi<-env.mat.metaT$Layer
levels(env.mat.metaT$epi)<-c("EPI","MES","EPI","EPI")
levels(env.mat.metaT$Layer)<-c("DCM","MES","DCM","SRF")
env.mat.metaT$Absolute.Latitude<-abs(env.mat.metaT$Latitude)

# metaG
metaG.norm<-fread("zcat < ../data/processed/NOG_metaG.norm.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm)<-metaG.norm$V1
metaG.norm<-log10(metaG.norm[,-1]+1E-12)
BCdist.metaG<-dist(metaG.norm)

# metaT
metaT.norm<-fread("zcat < ../data/processed/NOG_metaT.norm.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm)<-metaT.norm$V1
metaT.norm<-log10(metaT.norm[,-1]+1E-12)
metaT.norm<-metaT.norm[which(env.mat.metaT$Layer!="NA"),]
env.mat.metaT<-env.mat.metaT[which(env.mat.metaT$Layer!="NA"),]
BCdist.metaT<-dist(metaT.norm)


# Expression
exp.norm<-fread("zcat < ../data/processed/NOG_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(exp.norm)<-exp.norm$V1
exp.norm<-exp.norm[,-1]
BCdist.exp<-dist(exp.norm)

env.mat.match<-env.mat.match[match(rownames(exp.norm),env.mat.match$Barcode),]

# miTags
otu.table<-fread("zcat < ../data/processed/miTags/OTU.tab.TARA180.noeuks.log2.txt.gz",sep="\t",header=T,data.table = F)
rownames(otu.table)<-otu.table$V1
otu.table<-otu.table[,-1]
otu.table<-otu.table[match(gsub("-",".",env.mat.metaG$Sample_name),rownames(otu.table)),]
BCdist.tax<-dist(otu.table)

################################
## Test for depth differences ##
################################

env.mat.metaG$polar2<-as.character(env.mat.metaG$polar)
env.mat.metaG$polar2[env.mat.metaG$polar2=="polar" & env.mat.metaG$Latitude<0]<-"South.pole"
env.mat.metaG$polar2[env.mat.metaG$polar2=="polar" & env.mat.metaG$Latitude>0]<-"North.pole"
env.mat.metaG$polar2<-as.factor(env.mat.metaG$polar2)

env.mat.metaT$polar2<-as.character(env.mat.metaT$polar)
env.mat.metaT$polar2[env.mat.metaT$polar2=="polar" & env.mat.metaT$Latitude<0]<-"South.pole"
env.mat.metaT$polar2[env.mat.metaT$polar2=="polar" & env.mat.metaT$Latitude>0]<-"North.pole"
env.mat.metaT$polar2<-as.factor(env.mat.metaT$polar2)

env.mat.metaT$polar2<-as.character(env.mat.metaT$polar)
env.mat.metaT$polar2[env.mat.metaT$polar2=="polar" & env.mat.metaT$Latitude<0]<-"South.pole"
env.mat.metaT$polar2[env.mat.metaT$polar2=="polar" & env.mat.metaT$Latitude>0]<-"North.pole"
env.mat.metaT$polar2<-as.factor(env.mat.metaT$polar2)

env.mat.match$polar2<-as.character(env.mat.match$polar)
env.mat.match$polar2[env.mat.match$polar2=="polar" & env.mat.match$Latitude<0]<-"South.pole"
env.mat.match$polar2[env.mat.match$polar2=="polar" & env.mat.match$Latitude>0]<-"North.pole"
env.mat.match$polar2<-as.factor(env.mat.match$polar2)


# Taxonomy
adonis(BCdist.tax~env.mat.metaG$epi+env.mat.metaG$Layer,permutations = 10000)
adonis(BCdist.tax~env.mat.metaG$polar+env.mat.metaG$polar2,permutations = 10000)
adonis(BCdist.tax~env.mat.metaG$epi+env.mat.metaG$polar,permutations = 10000)

# metaG
adonis(BCdist.metaG~env.mat.metaG$epi+env.mat.metaG$Layer,permutations = 10000)
adonis(BCdist.metaG~env.mat.metaG$polar+env.mat.metaG$polar2,permutations = 10000)
adonis(BCdist.metaG~env.mat.metaG$epi+env.mat.metaG$polar,permutations = 10000)

# Expression
adonis(BCdist.exp~env.mat.match$epi+env.mat.match$Layer,permutations = 10000)
adonis(BCdist.exp~env.mat.match$polar+env.mat.match$polar2,permutations = 10000)
adonis(BCdist.exp~env.mat.match$epi+env.mat.match$polar,permutations = 10000)

# metaT
adonis(BCdist.metaT~env.mat.metaT$epi+env.mat.metaT$Layer,permutations = 10000)
adonis(BCdist.metaT~env.mat.metaT$polar+env.mat.metaT$polar2,permutations = 10000)
adonis(BCdist.metaT~env.mat.metaT$epi+env.mat.metaT$polar,permutations = 10000)

##########
## NMDs ##
##########

# Taxonomy
pca<-prcomp(otu.table,center = T,scale. = F)
R2<-summary(pca)$importance[2,1:2]
eigenvectors<-pca$x[,1:2]
tmp.df<-cbind(eigenvectors,env.mat.metaG)
tmp.df$Layer<-factor(tmp.df$Layer,levels=c("SRF","DCM","MES"))
tmp.df$inter<-interaction(tmp.df$Layer,tmp.df$polar)
g1<-ggplot(data=tmp.df,aes(x=PC1,y=PC2)) +
  geom_polygon(data=find_hull(tmp.df,"PC1","PC2","inter"),aes(x=PC1,y=PC2,group=inter),alpha=0.1,col="black") +
  geom_point(aes(col=Layer,shape=polar)) +
  theme_bw() +
  scale_color_manual(values=sushi.palette(alpha=0.7)[c(2,3,1)]) +
  xlab(paste("PC 1 (",round(100*R2[1],1)," % explained variance)",sep="")) +
  ylab(paste("PC 2 (",round(100*R2[2],1)," % explained variance)",sep="")) +
  labs(title="Taxonomic composition")

# metaG
pca<-prcomp(metaG.norm,center = T,scale. = F)
R2<-summary(pca)$importance[2,1:2]
eigenvectors<-pca$x[,1:2]
tmp.df<-cbind(eigenvectors,env.mat.metaG)
tmp.df$Layer<-factor(tmp.df$Layer,levels=c("SRF","DCM","MES"))
tmp.df$inter<-interaction(tmp.df$Layer,tmp.df$polar)
g2<-ggplot(data=tmp.df,aes(x=PC1,y=PC2)) +
  geom_polygon(data=find_hull(tmp.df,"PC1","PC2","inter"),aes(x=PC1,y=PC2,group=inter),alpha=0.1,col="black") +
  geom_point(aes(col=Layer,shape=polar)) +
  theme_bw() +
  scale_color_manual(values=sushi.palette(alpha=0.7)[c(2,3,1)]) +
  xlab(paste("PC 1 (",round(100*R2[1],1)," % explained variance)",sep="")) +
  ylab(paste("PC 2 (",round(100*R2[2],1)," % explained variance)",sep="")) +
  labs(title="Genomic composition")

# metaT
pca<-prcomp(metaT.norm,center = T,scale. = T)
R2<-summary(pca)$importance[2,1:2]
eigenvectors<-pca$x[,1:2]
tmp.df<-cbind(eigenvectors,env.mat.metaT)
tmp.df$Layer<-factor(tmp.df$Layer,levels=c("SRF","DCM","MES"))
tmp.df$inter<-interaction(tmp.df$Layer,tmp.df$polar)
g3<-ggplot(data=tmp.df %>% filter(!is.na(Layer)),aes(x=PC1,y=PC2)) +
  geom_polygon(data=find_hull(tmp.df,"PC1","PC2","inter"),aes(x=PC1,y=PC2,group=inter),alpha=0.1,col="black") +
  geom_point(aes(col=Layer,shape=polar)) +
  theme_bw() +
  scale_color_manual(values=sushi.palette(alpha=0.7)[c(2,3,1)]) +
  xlab(paste("PC 1 (",round(100*R2[1],1)," % explained variance)",sep="")) +
  ylab(paste("PC 2 (",round(100*R2[2],1)," % explained variance)",sep="")) +
  labs(title="Transcriptomic composition")

# Expression
pca<-prcomp(exp.norm,center = T,scale. = T)
R2<-summary(pca)$importance[2,1:2]
eigenvectors<-pca$x[,1:2]
tmp.df<-cbind(eigenvectors,env.mat.match)
tmp.df$Layer<-factor(tmp.df$Layer,levels=c("SRF","DCM","MES"))
tmp.df$inter<-interaction(tmp.df$Layer,tmp.df$polar)
g4<-ggplot(data=tmp.df,aes(x=PC1,y=PC2)) +
  geom_polygon(data=find_hull(tmp.df,"PC1","PC2","inter"),aes(x=PC1,y=PC2,group=inter),alpha=0.1,col="black") +
  geom_point(aes(col=Layer,shape=polar)) +
  theme_bw() +
  scale_color_manual(values=sushi.palette(alpha=0.7)[c(2,3,1)]) +
  xlab(paste("PC 1 (",round(100*R2[1],1)," % explained variance)",sep="")) +
  ylab(paste("PC 2 (",round(100*R2[2],1)," % explained variance)",sep="")) +
  labs(title="Expression")

(g1 | g2) / (g3 | g4)
ggsave("../results/SFig3_PCAs.pdf")


