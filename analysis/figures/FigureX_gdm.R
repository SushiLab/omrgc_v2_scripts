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

# Load marine distances
mar.dist.mat<-fread("zcat < ../data/processed/mardist.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(mar.dist.mat)<-mar.dist.mat$V1
mar.dist.mat<-mar.dist.mat[,-1]
#mar.dist.mat<-mar.dist.mat[match(env.mat$Sample_name,rownames(mar.dist.mat)),match(env.mat$Sample_name,rownames(mar.dist.mat))]

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
rownames(otu.table)<-env.mat$Barcode
BCdist.tax<-dist(otu.table)

#########################################################
# Select variables
vars<-c("Temperature","HCO3","CO3","Carbon.total","Alkalinity.total","NO2","PO4","NO2NO3","Si","Salinity","Density","Oxygen","NO3","ChlorophyllA","Fluorescence","PAR.TO","Iron.5m","Ammonium.5m","Gradient.Surface.temp(SST)","Okubo.Weiss","Lyapunov","Residence.time","Depth.Mixed.Layer","Brunt.Väisälä","Depth.Max.O2","Depth.Min.O2","Nitracline")

# Convert Bray-Curtis into a matrix and add the rownames as a variable (this variable needs to exist)
target.data<-BCdist.tax
bio<-(as.matrix(target.data)/max(target.data)) %>%
  as.data.frame() %>%
  rownames_to_column(var="sample") %>%
  as.data.frame()

# Add the rownames as a variable in the environmental table (this variable needs to exist and be names the same way as above)
environt<-env.mat %>%
  filter(epi=="MES") %>%
  dplyr::select(sample=Barcode,vars,Latitude,Longitude) %>%
  as.data.frame()

# Remove NAs
num.nas<-apply(environt,2,function(x){length(which(is.na(x)))})
environt<-environt[,which(num.nas<10)]
environt<-na.exclude(environt)
bio<-bio[match(environt$sample,bio$sample),c(1,match(environt$sample,colnames(bio)))]
nrow(bio)

# Format the data for GDM
library(gdm)
input.dat<-formatsitepair(bioData = bio,predData = environt,bioFormat = 3,siteColumn = "sample",XColumn = "Longitude",YColumn = "Latitude")

# Run the model (without geography)
mod<-gdm(input.dat,geo = F)
summary(mod)

# Get the spline functions and importance of each variable
mod.ispline<-isplineExtract(mod)
mod.ispline.df<-cbind(gather(as.data.frame(mod.ispline$x),key="var",value="x"),dplyr::select(gather(as.data.frame(mod.ispline$y),key="var",value="function_of_x"),"function_of_x"))

# Get the variable's importance
# (percent change in deviance explained by the full model and the deviance explained by a model fit with that variable permuted)
#mod.imp<-gdm.varImp(input.dat,geo=F,fullModelOnly = T) # takes a while
#mod.imp.df<-data.frame(predictor=rownames(mod.imp[[2]]),importance=mod.imp[[2]][,1])

# Get the ecological (i.e. environmental) distance, the observed biological distance and the predicted biological distance (i.e. predicted through the model)
mod.dist.df<-data.frame(eco.dist=mod$ecological,bio.obs=mod$observed,bio.pred=mod$predicted)


# Plot model's fit and variable's importance
g1<-ggplot(data = mod.dist.df,aes(x=eco.dist,y=bio.obs)) +
  geom_point(alpha=0.7) +
  geom_line(data = mod.dist.df,aes(x=eco.dist,y=bio.pred),col="gray",size=1) +
  theme_bw() +
  labs(title = "Generalized Dissimilarity Model",subtitle = paste("Percentage of deviance explained: ",round(mod$explained,2),sep="")) +
  xlab("Ecological distance") +
  ylab("Observed Bray-Curtis distance")
g2<-ggplot(data = mod.dist.df,aes(x=bio.obs,y=bio.pred)) +
  geom_abline(linetype=2) +
  geom_point(alpha=0.7) +
  theme_bw() +
  labs(subtitle = paste("Pearson r: ",round(cor(mod.dist.df$bio.obs,mod.dist.df$bio.pred),2)," | ","Spearman r: ",round(cor(mod.dist.df$bio.obs,mod.dist.df$bio.pred,method = "spearman"),2),sep="")) +
  xlab("Observed Bray-Curtis distance") +
  ylab("Predicted Bray-Curtis distance") +
  xlim(range(mod.dist.df[,2:3])) +
  ylim(range(mod.dist.df[,2:3]))
#g3<-ggplot(data = mod.imp.df,aes(x=fct_reorder(predictor,importance,.desc = T),y=importance)) +
#  geom_bar(stat = "identity") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
#  xlab("") +
#  ylab("Variable importance\n(% change in deviance)")

#p<-(g1 | g2) / g3
#p

g1 | g2

# Plot the variable's spline transformations
q<-ggplot(data=mod.ispline.df,aes(x=x,y=function_of_x)) +
  geom_line() +
  facet_wrap(~var,scales = "free_x") +
  theme_bw()
q
