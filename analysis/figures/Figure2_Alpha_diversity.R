# Load palette
library(patchwork)
library(ggplot2)
library(vegan)
library(data.table)
library(tidyverse)
source("../data/lib/sushipal.R")
palette(sushi.palette(alpha=0.7)[c(2,3,1,14)])

########################
## Taxonomic Richness ##
########################

# Load subsampled OTU table
otu.table<-fread("../data/processed/miTags/OTU.tab.rr.TARA180.noeuks.txt",sep="\t",header=T,data.table = F)
rownames(otu.table)<-otu.table$V1
otu.table<-otu.table[,-1]

# Load environmental data
env.mat<-fread("zcat < ../data/processed/NOG_env.mat.metaG.txt.gz",sep="\t",header=T,data.table = F,stringsAsFactors = T)
levels(env.mat$Layer)<-c("DCM","MES","DCM","SRF")
env.mat$Layer<-factor(env.mat$Layer,levels=c("SRF","DCM","MES"))

source("../data/lib/eveness_functions.R")
Richness<-specnumber(otu.table)
Shannon<-diversity(otu.table)
Eveness<-evenness(otu.table)[2,]
res<-data.frame(Richness,Shannon,Eveness,Smpl=names(Shannon)) %>% dplyr::select(Smpl,Richness,Shannon,Eveness) %>% gather(Metric,value,-Smpl)
res<-cbind(res,env.mat[match(res$Smpl,gsub("-",".",env.mat$Sample_name)),])
res$Metric<-factor(res$Metric,levels=c("Richness","Eveness","Shannon"))

p1<-ggplot(data=res %>% filter(Layer!="other" & Metric=="Richness") %>% mutate(facetting="metaG"),aes(x=Layer,y=value,fill=polar)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  theme_bw() +
  theme(legend.position="none",plot.margin = margin(0,0,0,0, "cm")) +
  scale_fill_manual(values=c("#FB8072","#80B1D3")) +
  ylab("Number of OTUs") +
  #facet_grid(.~facetting) +
  labs(title="Species richness")

test.data<-res[res$Metric=="Richness" & res$Layer %in% c("SRF","DCM","MES"),]
summary(aov(value~Layer+polar,data=test.data))


#######################
# Functional richness #
#######################
env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",sep="\t",header=T,data.table = F,stringsAsFactors = T)
env.mat.match$Layer<-factor(env.mat.match$Layer,levels = c("SRF","DCM","MES","other"))

# Common samples in both metaG and metaT
alpha.ko.metaG<-read.table("../data/processed/NOG_metaG.match.scaled.txt.rrmedian_alpha_diversity.tsv",header=T,sep="\t",skip=1)
alpha.ko.metaG<-alpha.ko.metaG %>% dplyr::select(Smpl,Richness,Shannon,Eveness) %>% gather(Metric,value,-Smpl)
alpha.ko.metaT<-read.table("../data/processed/NOG_metaT.match.scaled.txt.rrmedian_alpha_diversity.tsv",header=T,sep="\t",skip=1)
alpha.ko.metaT<-alpha.ko.metaT %>% dplyr::select(Smpl,Richness,Shannon,Eveness) %>% gather(Metric,value,-Smpl)

plot.df<-bind_rows(alpha.ko.metaG,alpha.ko.metaT)
plot.df$dataset<-c(rep("metaG",nrow(alpha.ko.metaG)),rep("metaT",nrow(alpha.ko.metaT)))
plot.df<-cbind(plot.df,env.mat.match[match(plot.df$Smpl,env.mat.match$Barcode),])
plot.df$Metric<-factor(plot.df$Metric,levels=c("Richness","Eveness","Shannon"))
plot.df$dataset<-as.factor(plot.df$dataset)
levels(plot.df$dataset)<-c("Genomic composition","Transcriptomic composition")

p2<-ggplot(data=plot.df %>% filter(Layer!="other" & Metric=="Richness"),aes(x=Layer,y=value/1000,fill=polar)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  theme_bw() +
  scale_fill_manual(values=c("#FB8072","#80B1D3")) +
  labs(title="Functional richness") +
  facet_grid(.~dataset,scales = "free_y") +
  ylab("Number of OGs (10^3)")


test.data<-plot.df[plot.df$Metric=="Richness" & plot.df$Layer %in% c("SRF","DCM","MES") & plot.df$dataset=="Genomic composition",]
summary(aov(value~Layer+polar,data=test.data))

test.data<-plot.df[plot.df$Metric=="Richness" & plot.df$Layer %in% c("SRF","DCM","MES") & plot.df$dataset=="Transcriptomic composition",]
summary(aov(value~Layer+polar,data=test.data))


tax<-res %>%
  dplyr::select("Sample_name","Metric","value") %>%
  filter(Metric=="Richness") %>%
  rename(tax.richness=value)
metaG<-plot.df %>%
  filter(dataset=="Genomic composition") %>%
  dplyr::select("Sample_name","Metric","value") %>%
  filter(Metric=="Richness") %>%
  rename(metaG.richness=value)
metaT<-plot.df %>%
  filter(dataset=="Transcriptomic composition") %>%
  dplyr::select("Sample_name","Metric","value") %>%
  filter(Metric=="Richness") %>%
  rename(metaT.richness=value)

all<-left_join(tax,metaG,by="Sample_name") %>% left_join(metaT,by="Sample_name") %>% left_join(env.mat.metaG,by="Sample_name")

p.a<-ggplot(data=all,aes(x=tax.richness,y=metaG.richness)) +
  geom_point(aes(col=polar,shape=epi)) +
  geom_smooth(method="lm",se=F,col="black") +
  theme_bw() +
  scale_fill_manual(values=c("#e59d8d","#81b0d2")) +
  xlab("Number of OTUs\n(Taxonomic composition)") +
  ylab("Number of OGs\n(Metagenomic composition)")

p.b<-ggplot(data=all,aes(x=tax.richness,y=metaT.richness)) +
  geom_point(aes(col=polar,shape=epi)) +
  geom_smooth(method="lm",se=F,col="black") +
  theme_bw() +
  scale_fill_manual(values=c("#e59d8d","#81b0d2")) +
  xlab("Number of OTUs\n(Taxonomic composition)") +
  ylab("Number of OGs\n(Metatranscriptomic composition)")

p.c<-ggplot(data=all,aes(x=metaG.richness,y=metaT.richness)) +
  geom_point(aes(col=polar,shape=epi)) +
  geom_smooth(method="lm",se=F,col="black") +
  theme_bw() +
  scale_fill_manual(values=c("#e59d8d","#81b0d2")) +
  xlab("Number of OGs\n(Metagenomic composition)") +
  ylab("Number of OGs\n(Metatranscriptomic composition)")

p<-p.a | p.b | p.c
ggsave(filename = "../results/FigSX2_Diversity_comparison.pdf",p,width=unit(12,"cm"),height=unit(3,"cm"))

cor.test(all$tax.richness,all$metaG.richness,use = "pairwise.complete.obs")
cor.test(all$tax.richness,all$metaT.richness,use = "pairwise.complete.obs")
cor.test(all$metaG.richness,all$metaT.richness,use = "pairwise.complete.obs")

#########
## VCR ##
#########

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


# Load data metaG
vir<-fread("zcat < ../data/raw/OM-RGC_v2_release_profile_by_viralcontig_metaG.tsv.gz",sep="\t",header=T,data.table = F)
rownames(vir)<-vir$viral_contig # Move NOG variale to rownames
vir<-vir[,-1]
colnames(vir)<-sub("_G","",colnames(vir)) # Remove "_G" from sample names


# Tests same order between environmental data matrix and metaG abundance matrix
all(colnames(vir)==env.mat.metaG$Sample_name)

# Compute the sum for all contigs (i.e. the VPR)
vir.sum<-colSums(vir)

res<-env.mat.metaG
res$vir.sum<-vir.sum

p3<-ggplot(data=res %>% filter(Layer %in% c("SRF","DCM","MES")),aes(x=Layer,y=vir.sum,fill=polar)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  scale_y_log10() +
  scale_fill_manual(values=c("#FB8072","#80B1D3")) +
  theme_bw() +
  theme(legend.position="none",plot.margin = margin(0,0,0,0, "cm")) +
  ylab("VCR") +
  labs(title="Virus-to-cell ratio (VCR)")


test.data<-res[res$Layer %in% c("SRF","DCM","MES"),]
summary(aov(vir.sum~Layer+polar,data=test.data))

res %>%
  group_by(Layer,polar) %>%
  summarise(mean(vir.sum))

p<-(p1 | p3) / p2
ggsave(filename = "../results/Fig2_Alpha_diversity.pdf",p,width=unit(5,"cm"),height=unit(5,"cm"))
