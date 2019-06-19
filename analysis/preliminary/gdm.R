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

# miTags
otu.table<-fread("../data/processed/miTags/OTU.tab.TARA180.noeuks.txt",sep="\t",header=T,data.table = F)
rownames(otu.table)<-otu.table$V1
otu.table<-otu.table[,-1]
otu.table<-otu.table[match(gsub("-",".",env.mat.metaG$Sample_name),rownames(otu.table)),]
BCdist.tax<-vegdist(otu.table)

#########
## GDM ##
#########
library(gdm)
env.mat<-env.mat.match %>% dplyr::select(site=Barcode,Longitude,Latitude,Temperature,Layer,ChlorophyllA,Oxygen,Carbon.total,NO2NO3) %>% mutate(site=gsub("-",".",site),Carbon.total=sqrt(Carbon.total))
env.mat<-env.mat %>% filter(Layer %in% c("SRF","DCM"))
env.mat<-na.exclude(env.mat)
env.mat<-env.mat %>% dplyr::select(-Layer)

bio.dist<-as.data.frame(as.matrix(BCdist.metaG))
bio.dist<-bio.dist/max(bio.dist)
bio.dist<-bio.dist[match(gsub("\\.","-",env.mat$site),rownames(bio.dist)),match(gsub("\\.","-",env.mat$site),colnames(bio.dist))]
bio.dist<-cbind(site=rownames(bio.dist),bio.dist)

input.dat<-formatsitepair(bioData = bio.dist,predData = env.mat,bioFormat = 3,XColumn = "Longitude",YColumn = "Latitude",siteColumn = "site")
mod<-gdm(input.dat,geo = F)
summary(mod)
#gdm.varImp(input.dat,geo=T)

plot(mod, plot.layout=c(3,3))
