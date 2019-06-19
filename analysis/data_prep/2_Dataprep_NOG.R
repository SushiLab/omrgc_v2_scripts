library(patchwork)
library(tidyquant)
library(vegan)
library(DESeq2)
library(data.table)
library(tidyverse)
library(vsn)
################
## METAGENOME ##
################
# Load environmental data
env.mat.metaG<-fread("../data/raw/metaG_metadata.tsv",sep="\t",header=T,data.table = F)
#rownames(env.mat.metaG)<-env.mat.metaG$Sample_name_goodlayer # add Sample as rowname
env.mat.metaG$Depth.nominal<-as.numeric(as.character(env.mat.metaG$Depth.nominal)) # convert Depth to numeric
env.mat.metaG$Layer<-factor(env.mat.metaG$Layer)
levels(env.mat.metaG$Layer)<-c("DCM","MES","MIX","SRF","MIX")
env.mat.metaG$Layer<-factor(env.mat.metaG$Layer,levels=c("SRF","DCM","MES","MIX"),ordered=T)
env.mat.metaG$polar<-rep("non polar",nrow(env.mat.metaG)) # Code polar/nonpolar factor
env.mat.metaG$polar[abs(env.mat.metaG$Latitude)>50]<-"polar"
env.mat.metaG$polar<-as.factor(env.mat.metaG$polar)

# Load data
metaG<-fread("zcat < ../data/raw/OM-RGC_v2_release_profile_by_NOG_metaG.tsv.gz",sep="\t",header=T,data.table = F)
rownames(metaG)<-metaG$NOG # Move NOG variale to rownames
metaG<-metaG[,-1]
colnames(metaG)<-sub("_G","",colnames(metaG)) # Remove "_G" from sample names
#no.metadata<-which(colnames(metaG) %in% env.mat.metaG$Sample_name_goodlayer==F) # exclude 4 samples without environmental information
#metaG<-metaG[,-c(no.metadata)]
env.mat.metaG<-env.mat.metaG[match(colnames(metaG),env.mat.metaG$Sample_name),]

# Tests same order between environmental data matrix and metaG abundance matrix
all(colnames(metaG)==env.mat.metaG$Sample_name)
colnames(metaG)<-env.mat.metaG$Barcode

# Transpose and remove 0-abundant OGs
metaG<-t(metaG)
metaG<-metaG[,which(colSums(metaG)>0)]



#######################
## METATRANSCRIPTOME ##
#######################
# Load environmental data
env.mat.metaT<-fread("../data/raw/metaT_metadata.tsv",sep="\t",header=T,data.table = F)
#rownames(env.mat.metaT)<-env.mat.metaT$Sample_name_goodlayer # add Sample as rowname
env.mat.metaT$Depth.nominal<-as.numeric(as.character(env.mat.metaT$Depth.nominal)) # convert Depth to numeric
env.mat.metaT$Layer<-factor(env.mat.metaT$Layer)
levels(env.mat.metaT$Layer)<-c("DCM","MES","MIX","SRF","MIX")
env.mat.metaT$Layer<-factor(env.mat.metaT$Layer,levels=c("SRF","DCM","MES","MIX"),ordered=T)
env.mat.metaT$polar<-rep("non polar",nrow(env.mat.metaT)) # Code polar/nonpolar factor
env.mat.metaT$polar[abs(env.mat.metaT$Latitude)>50]<-"polar"
env.mat.metaT$polar<-as.factor(env.mat.metaT$polar)

# Load data
metaT<-fread("zcat < ../data/raw/OM-RGC_v2_release_profile_by_NOG_metaT.tsv.gz",sep="\t",header=T,data.table = F)
rownames(metaT)<-metaT$NOG # Move NOG variale to rownames
metaT<-metaT[,-1]
colnames(metaT)<-sub("_T","",colnames(metaT)) # Remove "_G" from sample names
#no.metadata<-which(colnames(metaT) %in% env.mat.metaT$Sample_name_goodlayer==F) # exclude 4 samples without environmental information
#metaT<-metaT[,-c(no.metadata)]
env.mat.metaT<-env.mat.metaT[match(colnames(metaT),env.mat.metaT$Sample_name),]

# Tests same order between environmental data matrix and metaT abundance matrix
all(colnames(metaT)==env.mat.metaT$Sample_name)
colnames(metaT)<-env.mat.metaT$Barcode

# Transpose and remove 0-abundant OGs
metaT<-t(metaT)
metaT<-metaT[,which(colSums(metaT)>0)]


##############
## COG info ##
##############
cog.info<-read.table("../data/lib/cognames2003-2014.tab",sep="\t",quote="",comment.char="",header=T)
cog.10MGs<-read.table("../data/lib/10MG_COG_ids",sep="\t",quote="",comment.char="",header=F)
cog.10MGs<-cog.info[match(cog.10MGs$V1,cog.info$X..COG),]

# Normalize metaG and metaT dividing by the median abundance of 10MGs and computing the relative abundance
metaG.norm<-metaG/apply(metaG[,match(cog.10MGs$X..COG,colnames(metaG))],1,median)
metaT.norm<-metaT/apply(metaT[,match(cog.10MGs$X..COG,colnames(metaT))],1,median)

# Save the sum of the MGs
write.table(apply(metaG[,match(cog.10MGs$X..COG,colnames(metaG))],1,median),file = "../data/processed/NOG_metaG.sumMGs.tsv",sep="\t",quote=F,col.names = F)
write.table(apply(metaT[,match(cog.10MGs$X..COG,colnames(metaT))],1,median),file = "../data/processed/NOG_metaT.sumMGs.tsv",sep="\t",quote=F,col.names = F)

#########################################
## MATCH metaG and meaT SHARED SAMPLES ##
#########################################
env.mat.match<-fread("../data/raw/match_metadata.tsv",sep="\t",header=T,data.table = F)
env.mat.match$Depth.nominal<-as.numeric(as.character(env.mat.match$Depth.nominal)) # convert Depth to numeric
env.mat.match$Layer<-factor(env.mat.match$Layer)
levels(env.mat.match$Layer)<-c("DCM","MES","MIX","SRF","MIX")
env.mat.match$Layer<-factor(env.mat.match$Layer,levels=c("SRF","DCM","MES","MIX"),ordered=T)
env.mat.match$polar<-rep("non polar",nrow(env.mat.match)) # Code polar/nonpolar factor
env.mat.match$polar[abs(env.mat.match$Latitude)>50]<-"polar"
env.mat.match$polar<-as.factor(env.mat.match$polar)

common.samples<-env.mat.match$Barcode

# Test that Latitude, Longitude and Depth are equal for the matching metaG and metaT samples
#for (i in 1:length(common.samples)){
#  res<-all(apply(rbind(env.mat.metaG[match(common.samples[i],env.mat.metaG$Sample_name_goodlayer),],env.mat.metaT[match(common.samples[i],env.mat.metaT$Sample_name_goodlayer),])[,c('Latitude','Longitude','Depth.nominal')],2,function(x){length(unique(x))==1}))
#  cat(common.samples[i]," -> ",res,"\n")
#}

metaG.red<-metaG[match(sapply(strsplit(env.mat.match$Barcode,"-"),"[[",1),rownames(metaG)),]
rownames(metaG.red)<-env.mat.match$Barcode
metaT.red<-metaT[match(sapply(strsplit(env.mat.match$Barcode,"-"),"[[",2),rownames(metaT)),]
rownames(metaT.red)<-env.mat.match$Barcode

metaG.norm.red<-metaG.norm[match(sapply(strsplit(env.mat.match$Barcode,"-"),"[[",1),rownames(metaG.norm)),]
rownames(metaG.norm.red)<-env.mat.match$Barcode
metaT.norm.red<-metaT.norm[match(sapply(strsplit(env.mat.match$Barcode,"-"),"[[",2),rownames(metaT.norm)),]
rownames(metaT.norm.red)<-env.mat.match$Barcode


# metaG and metaT
metaG.gath<-as.data.frame(metaG.red) %>%
  rownames_to_column(var="Sample") %>%
  gather(key = "KO",value = "ab",-Sample) %>%
  mutate(dataset="metaG")
metaT.gath<-as.data.frame(metaT.red) %>%
  rownames_to_column(var="Sample") %>%
  gather(key = "KO",value = "ab",-Sample) %>%
  mutate(dataset="metaT")
merged<-bind_rows(metaG.gath,metaT.gath) %>%
  spread(key = KO,value = ab,fill = 0) %>%
  filter(Sample %in% common.samples)
metaG.match<-merged %>%
  filter(dataset=="metaG") %>%
  select(-dataset)
metaT.match<-merged %>%
  filter(dataset=="metaT") %>%
  select(-dataset)
rownames(metaG.match)<-metaG.match$Sample
metaG.match<-metaG.match[,-1]
rownames(metaT.match)<-metaT.match$Sample
metaT.match<-metaT.match[,-1]

# normalized metaG and metaT
metaG.norm.gath<-as.data.frame(metaG.norm.red) %>%
  rownames_to_column(var="Sample") %>%
  gather(key = "KO",value = "ab",-Sample) %>%
  mutate(dataset="metaG")
metaT.norm.gath<-as.data.frame(metaT.norm.red) %>%
  rownames_to_column(var="Sample") %>%
  gather(key = "KO",value = "ab",-Sample) %>%
  mutate(dataset="metaT")
merged<-bind_rows(metaG.norm.gath,metaT.norm.gath) %>%
  spread(key = KO,value = ab,fill = 0) %>%
  filter(Sample %in% common.samples)
metaG.norm.match<-merged %>%
  filter(dataset=="metaG") %>%
  select(-dataset)
metaT.norm.match<-merged %>%
  filter(dataset=="metaT") %>%
  select(-dataset)
rownames(metaG.norm.match)<-metaG.norm.match$Sample
metaG.norm.match<-metaG.norm.match[,-1]
rownames(metaT.norm.match)<-metaT.norm.match$Sample
metaT.norm.match<-metaT.norm.match[,-1]

#########################
## MetaT / MetaG ratio ##
#########################
input.data<-t(metaG.norm.match)
colnames(input.data)<-paste(colnames(input.data),"_G",sep="")
input.dataT<-t(metaT.norm.match)
colnames(input.dataT)<-paste(colnames(input.data),"_T",sep="")
input.metadata<-env.mat.match
input.data.merged<-cbind(input.data,input.dataT)
to.integer<-function(taula){taula<-round(taula/(max(taula)/1e9))}
input.data.merged<-to.integer(input.data.merged)
input.metadata.merged<-rbind(input.metadata,input.metadata)
rownames(input.metadata.merged)<-colnames(input.data.merged)
dds <- DESeqDataSetFromMatrix(countData = input.data.merged,colData = input.metadata.merged,design = ~ 1)
vsd <- vst(dds, blind=F)
res<-assay(vsd)
res.untrans<-assay(dds)

toplot.df<-data.frame(ko=rownames(res),mean=apply(res,1,mean),variance=apply(res,1,var),transformation="Transformation: VST")
toplot.df$rank<-rank(toplot.df$mean)
toplot.df.untr<-data.frame(ko=rownames(res.untrans),mean=apply(log2(res.untrans+1),1,mean),variance=apply(log2(res.untrans+1),1,var),transformation="Transformation: none")
toplot.df.untr$rank<-rank(toplot.df.untr$mean)
limits.x<-range(c(range(toplot.df$mean),range(toplot.df.untr$mean)))
limits.y<-range(c(range(toplot.df$variance),range(toplot.df.untr$variance)))
toplot.merged<-bind_rows(toplot.df,toplot.df.untr)

p1<-ggplot(toplot.merged %>% filter(ko!="unknown") %>% arrange(mean),aes(x=mean,y=variance)) +
  geom_point(alpha=0.1,col="gray60") +
  geom_ma(ma_fun = DEMA,n = 500,linetype=1,col="black") +
  theme_bw() +
  ylim(limits.y) +
  facet_grid(.~transformation,scales = "free_x")

p2<-ggplot(toplot.merged %>% filter(ko!="unknown") %>% arrange(mean),aes(x=mean,y=variance)) +
  geom_hex() +
  geom_ma(ma_fun = DEMA,n = 500,linetype=1,col="black") +
  theme_bw() +
  ylim(limits.y) +
  scale_fill_gradient(low="lightblue",high="darkred",trans="sqrt",breaks = c(100,1000,5000,10000,20000), labels = c(100,1000,5000,10000,20000)) +
  facet_grid(.~transformation,scales = "free_x")
ggsave("../results/VST_evaluation_OG_dots.pdf",p1,width=6,height=4)
ggsave("../results/VST_evaluation_OG_hex.pdf",p2,width=6,height=4)

input.data.norm<-res[,1:nrow(metaG.norm.match)]
input.data.normT<-res[,(nrow(metaG.norm.match)+1):(nrow(metaG.norm.match)*2)]
metaG.norm.match.log2<-t(input.data.norm)+log2((max(cbind(input.data,input.dataT))/1e9))
rownames(metaG.norm.match.log2)<-gsub("_G","",rownames(metaG.norm.match.log2))
metaT.norm.match.log2<-t(input.data.normT)+log2((max(cbind(input.data,input.dataT))/1e9))
rownames(metaT.norm.match.log2)<-gsub("_T","",rownames(metaG.norm.match.log2))

# Renormalize with the 10 MGs
metaG.norm.match.log2<-metaG.norm.match.log2-apply(metaG.norm.match.log2[,which(colnames(metaG.norm.match.log2) %in% cog.10MGs$X..COG)],1,median)
metaT.norm.match.log2<-metaT.norm.match.log2-apply(metaG.norm.match.log2[,which(colnames(metaT.norm.match.log2) %in% cog.10MGs$X..COG)],1,median)
ratio.mat<-metaT.norm.match.log2-metaG.norm.match.log2


###############
## Save data ##
###############

write.table(env.mat.metaG,file="../data/processed/NOG_env.mat.metaG.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(env.mat.metaT,file="../data/processed/NOG_env.mat.metaT.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(env.mat.match,file="../data/processed/NOG_env.mat.match.txt",sep="\t",row.names=T,col.names=NA,quote=F)

write.table(metaG,file = "../data/processed/NOG_metaG.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.match,file = "../data/processed/NOG_metaG.match.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.norm,file = "../data/processed/NOG_metaG.norm.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.norm.match,file = "../data/processed/NOG_metaG.norm.match.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.norm.match.log2,file = "../data/processed/NOG_metaG.norm.match.log2.txt",sep="\t",row.names=T,col.names=NA,quote=F)

write.table(metaT,file = "../data/processed/NOG_metaT.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.match,file = "../data/processed/NOG_metaT.match.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.norm,file = "../data/processed/NOG_metaT.norm.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.norm.match,file = "../data/processed/NOG_metaT.norm.match.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.norm.match.log2,file = "../data/processed/NOG_metaT.norm.match.log2.txt",sep="\t",row.names=T,col.names=NA,quote=F)

write.table(ratio.mat,file = "../data/processed/NOG_ratio.mat.txt",sep="\t",row.names=T,col.names=NA,quote=F)

#####################################################
## Save data transformed to integer in BIOM format ##
## (input for sparCC)                              ##
#####################################################
#write.biom<-function(taula,file=NULL){
#  write.table(t(c("#OTU ID",colnames(t(taula)))),file=file,quote=F,sep="\t",col.names=F,row.names=F)
#  write.table(t(taula),file=file,quote=F,sep="\t",col.names=F,append = T)
#}

#occ.reduce<-function(taula,occ=0.2){
#  taula[,which(colSums(decostand(taula,"pa"))>occ*nrow(taula))]}

#metaG.norm.int<-occ.reduce(to.integer(metaG.norm))
#write.biom(metaG.norm.int,file="../data/processed/NOG_metaG.norm.int.biom")
#metaG.norm.match.int<-occ.reduce(to.integer(metaG.norm.match))
#write.biom(metaG.norm.match.int,file="../data/processed/NOG_metaG.norm.match.int.biom")
#metaT.norm.int<-occ.reduce(to.integer(metaT.norm))
#write.biom(metaT.norm.int,file="../data/processed/NOG_metaT.norm.int.biom")
#metaT.norm.match.int<-occ.reduce(to.integer(metaT.norm.match))
#write.biom(metaT.norm.match.int,file="../data/processed/NOG_metaT.norm.match.int.biom")

#ratio.mat.int<-occ.reduce(to.integer(ratio.mat+abs(min(ratio.mat))))
#write.biom(ratio.mat.int,file="../data/processed/NOG_ratio.mat.int.biom")

