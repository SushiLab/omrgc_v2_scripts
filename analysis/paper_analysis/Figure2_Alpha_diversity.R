# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure 2: Aplha diversity ==============================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

## FIXME
## variable pi was not defined in correlatino plots
## I used SRF/DCM as EPI the rest as MES

# Libraries ------------------------------------------------------------------------------

library(patchwork)
library(ggplot2)
library(vegan)
library(data.table)
library(tidyverse)
# Load palette
source("lib/sushipal.R")
palette(sushi.palette(alpha=0.7)[c(2,3,1,14)])
# Source functions
source("lib/eveness_functions.R")

# Load data ------------------------------------------------------------------------------

# Load subsampled OTU table
otu.table<-fread("../data/processed/mitags.OTU.tab.rr.TARA180.noeuks.txt",sep="\t",header=T,data.table = F)
rownames(otu.table)<-otu.table$V1
otu.table<-otu.table[,-1]

# Load meta g and t profiles
# Common samples in both metaG and metaT
alpha.ko.metaG<-read.table("../data/processed/NOG_metaG.match.scaled.txt.rrmedian_alpha_diversity.tsv",header=T,sep="\t",skip=1)
alpha.ko.metaG<-alpha.ko.metaG %>% dplyr::select(Smpl,Richness,Shannon,Eveness) %>% gather(Metric,value,-Smpl)
alpha.ko.metaT<-read.table("../data/processed/NOG_metaT.match.scaled.txt.rrmedian_alpha_diversity.tsv",header=T,sep="\t",skip=1)
alpha.ko.metaT<-alpha.ko.metaT %>% dplyr::select(Smpl,Richness,Shannon,Eveness) %>% gather(Metric,value,-Smpl)

# Load environmental data
env.mat<-fread("zcat < ../data/processed/NOG_env.mat.metaG.txt.gz",sep="\t",header=T,data.table = F,stringsAsFactors = T)
levels(env.mat$Layer)<-c("DCM","MES","DCM","SRF")
env.mat$Layer<-factor(env.mat$Layer,levels=c("SRF","DCM","MES"))
env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",sep="\t",header=T,data.table = F,stringsAsFactors = T)
env.mat.match$Layer<-factor(env.mat.match$Layer,levels = c("SRF","DCM","MES","other"))

# Plot taxonomic richness ----------------------------------------------------------------

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

# Plot functional richness ---------------------------------------------------------------

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

p.viol <- p1 | p2

ggsave(filename = "../results/figures/Figure_diversity_violin.pdf",p.viol,width=unit(12,"cm"),height=unit(3,"cm"))

# Correlation plots ----------------------------------------------------------------------

tax<-res %>%
  dplyr::select("Sample_name","Metric","value") %>%
  filter(Metric=="Richness") %>%
  dplyr::rename(tax.richness=value)
metaG<-plot.df %>%
  filter(dataset=="Genomic composition") %>%
  dplyr::select("Sample_name","Metric","value") %>%
  filter(Metric=="Richness") %>%
  dplyr::rename(metaG.richness=value)
metaT<-plot.df %>%
  filter(dataset=="Transcriptomic composition") %>%
  dplyr::select("Sample_name","Metric","value") %>%
  filter(Metric=="Richness") %>%
  dplyr::rename(metaT.richness=value)

all<-left_join(tax,metaG,by="Sample_name") %>% left_join(metaT,by="Sample_name") %>% left_join(env.mat,by="Sample_name") %>%
  mutate(epi = ifelse(grepl('SRF|DCM', Layer), 'EPI', 'MES'))

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

p <- p.a | p.b | p.c

ggsave(filename = "../results/figures/Figure_diversity_comparison.pdf",p,width=unit(12,"cm"),height=unit(3,"cm"))


# Statistical tests ----------------------------------------------------------------------

sink('../results/tables/Table_stats_for_alpha_div.txt')
print('ANOVA')
print('mitags')
test.data<-res[res$Metric=="Richness" & res$Layer %in% c("SRF","DCM","MES"),]
summary(aov(value~Layer+polar,data=test.data))

print('metag')
test.data<-plot.df[plot.df$Metric=="Richness" & plot.df$Layer %in% c("SRF","DCM","MES") & plot.df$dataset=="Genomic composition",]
summary(aov(value~Layer+polar,data=test.data))

print('metat')
test.data<-plot.df[plot.df$Metric=="Richness" & plot.df$Layer %in% c("SRF","DCM","MES") & plot.df$dataset=="Transcriptomic composition",]
summary(aov(value~Layer+polar,data=test.data))

print('')
print('CORR')

cor.test(all$tax.richness,all$metaG.richness,use = "pairwise.complete.obs")

cor.test(all$tax.richness,all$metaT.richness,use = "pairwise.complete.obs")

cor.test(all$metaG.richness,all$metaT.richness,use = "pairwise.complete.obs")
sink()
