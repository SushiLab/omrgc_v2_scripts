# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure X: Virus to Cell Ratios (VCR) ===================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

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

# Load environmental data
env.mat.metaG<-fread("../data/raw/metaG_metadata.tsv",sep="\t",header=T,data.table = F)

# Prep data ------------------------------------------------------------------------------

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

# Plot data ------------------------------------------------------------------------------

p3<-ggplot(data=res %>% filter(Layer %in% c("SRF","DCM","MES")),aes(x=Layer,y=vir.sum,fill=polar)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  scale_y_log10() +
  scale_fill_manual(values=c("#FB8072","#80B1D3")) +
  theme_bw() +
  theme(legend.position="none",plot.margin = margin(0,0,0,0, "cm")) +
  ylab("VCR") +
  labs(title="Virus-to-cell ratio (VCR)")


ggsave(filename = "../results/figures/Figure_Virus_Cell_ratios.pdf",p3,width=unit(2.5,"cm"),height=unit(2.5,"cm"))

# Statistical tests ----------------------------------------------------------------------

sink('../results/tables/Table_stats_for_VCR.txt')
test.data<-res[res$Layer %in% c("SRF","DCM","MES"),]
summary(aov(vir.sum~Layer+polar,data=test.data))

res %>%
  group_by(Layer,polar) %>%
  summarise(mean(vir.sum))
sink()
