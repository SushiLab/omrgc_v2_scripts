library(tidyverse)
library(data.table)
library(patchwork)
source("../lib/sushipal.R")
source("../lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])
library(vegan)
library(EcolUtils)

# Define number of iterations:
iter<-1000

# Load environmental data
env.mat<-fread("../../data/processed/KO_env.mat.metaG.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat)<-env.mat$V1
env.mat<-env.mat[,-1]
env.mat$epi<-env.mat$Layer
levels(env.mat$epi)<-c("EPI","MES","EPI","EPI")
levels(env.mat$Layer)<-c("DCM","MES","DCM","SRF")
env.mat$Absolute.Latitude<-abs(env.mat$Latitude)

# Restrict to EPIPELAGIC
env.mat<-env.mat[env.mat$epi=="EPI",]

# ###############################################
# # Analysis for taxonomic composition (miTags) #
# #                genus level                  #
# ###############################################
# 
# # Load genus table
# otu.table<-fread("../../data/release/OM-RGC_v2_taxonomic_profiles/mitags_tab_genus.tsv.gz",sep="\t",header=T,data.table = F)
# otu.table<-otu.table %>% rename(Barcode=V1) %>% select(-unclassified)  %>% column_to_rownames(var = "Barcode")
# otu.table<-otu.table[match(env.mat$Barcode,rownames(otu.table)),]
# 
# # Mean relative abundance
# mean.ab<-otu.table %>%
#   decostand(method="total") %>%
#   rownames_to_column(var="Barcode") %>%
#   gather(key="genus",value="relab",-Barcode) %>%
#   group_by(genus) %>%
#   summarise(relab=mean(relab))
# occ<-otu.table %>%
#   rownames_to_column(var="Barcode") %>%
#   gather(key="genus",value="relab",-Barcode) %>%
#   filter(relab>0) %>%
#   group_by(genus) %>%
#   summarise(occ=n())
# 
# # Latitudinal niche value
# res.lat<-niche.val(otu.table,env.mat$Absolute.Latitude,n = iter)
# toplot.lat<-res.lat %>%
#   rownames_to_column(var="genus") %>%
#   left_join(mean.ab,by="genus") %>%
#   left_join(occ,by="genus") %>%
#   dplyr::arrange(desc(relab)) %>%
#   mutate(sign=fct_recode(sign,"POLAR"="HIGHER","NON-POLAR"="LOWER","N.S."="NON SIGNIFICANT")) %>%
#   filter(grepl(";Mitochondria",genus)==F)
# 
# # Proportion of families with association to polar /non-polar
# toplot.lat %>% filter(occ>10) %>% {table(.$sign)} %>% {prop.table(.)}
# 
# # Plots
# g<-ggplot() +
#   geom_errorbar(data=toplot.lat %>% head(n=30),aes(x=fct_reorder(genus,relab),ymax=uppCI,ymin=lowCI),col="gray40") +
#   geom_point(data=toplot.lat %>% head(n=30),aes(x=fct_reorder(genus,relab),y=mean.simulated),col="gray40") +
#   geom_point(data=toplot.lat %>% head(n=30),aes(x=fct_reorder(genus,relab),y=observed,fill=sign,size=100*relab),shape=21) +
#   coord_flip() +
#   theme_bw() +
#   theme(legend.position = "right",legend.direction = "vertical") +
#   scale_size(name="Mean relative abundance [%]",breaks=c(0.5,1,5,10),limits=c(0.5,20)) +
#   xlab(NULL) +
#   ylab("Abundance-weighted\nabsolute latitude") +
#   scale_fill_manual(name="Polar / Non-polar association",values=c("#1C6EAC","#F48B12","#BEBEBE")) +
#   ylim(0,90)
# 
# ggsave(filename = "../../results/paper_figures/review/FigSX_taxonomy_latitude_associations_genus.pdf",g,width=12,height=5)
# 
# ###############################################
# # Analysis for taxonomic composition (miTags) #
# #                family level                  #
# ###############################################
# 
# # Load family table
# otu.table<-fread("../../data/release/OM-RGC_v2_taxonomic_profiles/mitags_tab_family.tsv.gz",sep="\t",header=T,data.table = F)
# otu.table<-otu.table %>% rename(Barcode=V1) %>% select(-unclassified)  %>% column_to_rownames(var = "Barcode")
# otu.table<-otu.table[match(env.mat$Barcode,rownames(otu.table)),]
# 
# # Mean relative abundance
# mean.ab<-otu.table %>%
#   decostand(method="total") %>%
#   rownames_to_column(var="Barcode") %>%
#   gather(key="family",value="relab",-Barcode) %>%
#   group_by(family) %>%
#   summarise(relab=mean(relab))
# occ<-otu.table %>%
#      rownames_to_column(var="Barcode") %>%
#      gather(key="family",value="relab",-Barcode) %>%
#      filter(relab>0) %>%
#      group_by(family) %>%
#      summarise(occ=n())
# 
# # Latitudinal niche value
# res.lat<-niche.val(otu.table,env.mat$Absolute.Latitude,n = iter)
# toplot.lat<-res.lat %>%
#   rownames_to_column(var="family") %>%
#   left_join(mean.ab,by="family") %>%
#   left_join(occ,by="family") %>%
#   dplyr::arrange(desc(relab)) %>%
#   mutate(sign=fct_recode(sign,"POLAR"="HIGHER","NON-POLAR"="LOWER","N.S."="NON SIGNIFICANT")) %>%
#   filter(grepl(";Mitochondria",family)==F)
# 
# # Proportion of families with association to polar /non-polar
# toplot.lat %>% filter(occ>10) %>% {table(.$sign)} %>% {prop.table(.)}
# 
# # Plots
# g<-ggplot() +
#   geom_errorbar(data=toplot.lat %>% head(n=30),aes(x=fct_reorder(family,relab),ymax=uppCI,ymin=lowCI),col="gray40") +
#   geom_point(data=toplot.lat %>% head(n=30),aes(x=fct_reorder(family,relab),y=mean.simulated),col="gray40") +
#   geom_point(data=toplot.lat %>% head(n=30),aes(x=fct_reorder(family,relab),y=observed,fill=sign,size=100*relab),shape=21) +
#   coord_flip() +
#   theme_bw() +
#   theme(legend.position = "right",legend.direction = "vertical") +
#   scale_size(name="Mean relative abundance [%]",breaks=c(0.5,1,5,10),limits=c(0.5,20)) +
#   xlab(NULL) +
#   ylab("Abundance-weighted\nabsolute latitude") +
#   scale_fill_manual(name="Polar / Non-polar association",values=c("#1C6EAC","#F48B12","#BEBEBE")) +
#   ylim(0,90)
# 
# ggsave(filename = "../../results/paper_figures/review/FigSX_taxonomy_latitude_associations_family.pdf",g,width=12,height=5)

###############################################
# Analysis for taxonomic composition (miTags) #
#                otu level                  #
###############################################

# Load otu table
otu.table<-fread("../../data/release/OM-RGC_v2_taxonomic_profiles/mitags_tab_otu.tsv.gz",sep="\t",header=T,data.table = F)
otu.table<-otu.table %>% rename(Barcode=V1) %>% select(-unclassified)  %>% column_to_rownames(var = "Barcode")
otu.table<-otu.table[match(env.mat$Barcode,rownames(otu.table)),]

# Mean relative abundance
mean.ab<-otu.table %>%
  decostand(method="total") %>%
  rownames_to_column(var="Barcode") %>%
  gather(key="otu",value="relab",-Barcode) %>%
  group_by(otu) %>%
  summarise(relab=mean(relab))
occ<-otu.table %>%
  rownames_to_column(var="Barcode") %>%
  gather(key="otu",value="relab",-Barcode) %>%
  filter(relab>0) %>%
  group_by(otu) %>%
  summarise(occ=n())

# Latitudinal niche value
res.lat<-niche.val(otu.table,env.mat$Absolute.Latitude,n = iter)
toplot.lat<-res.lat %>%
  rownames_to_column(var="otu") %>%
  left_join(mean.ab,by="otu") %>%
  left_join(occ,by="otu") %>%
  dplyr::arrange(desc(relab)) %>%
  mutate(sign=fct_recode(sign,"POLAR"="HIGHER","NON-POLAR"="LOWER","N.S."="NON SIGNIFICANT")) %>%
  filter(grepl(";Mitochondria",otu)==F)

# Proportion of families with association to polar /non-polar
toplot.lat %>% filter(occ>10) %>% {table(.$sign)} %>% {prop.table(.)}

# Plots
g<-ggplot() +
  geom_errorbar(data=toplot.lat %>% head(n=60),aes(x=fct_reorder(otu,observed),ymax=uppCI,ymin=lowCI),col="gray40") +
  geom_point(data=toplot.lat %>% head(n=60),aes(x=fct_reorder(otu,observed),y=mean.simulated),col="gray40") +
  geom_point(data=toplot.lat %>% head(n=60),aes(x=fct_reorder(otu,observed),y=observed,fill=sign,size=100*relab),shape=21) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right",legend.direction = "vertical") +
  scale_size(name="Mean relative abundance [%]",breaks=c(0.05,0.1,0.5,1),limits=c(0,5),range = c(1,10)) +
  xlab(NULL) +
  ylab("Abundance-weighted\nabsolute latitude") +
  scale_fill_manual(name="Polar / Non-polar association",values=c("#1C6EAC","#F48B12","#BEBEBE")) +
  ylim(0,90)

ggsave(filename = "../../results/paper_figures/review/FigSX_taxonomy_latitude_associations_otu.pdf",g,width=15,height=10)

# ###########################################################
# # Analysis for functional composition (eggNOG abundances) #
# ###########################################################
# # Load OG abundances
# og.table<-fread("../../data/processed/NOG_metaG.norm.txt.gz",sep="\t",header=T,data.table = F)
# og.table<-og.table %>% rename(Barcode=V1) %>% column_to_rownames(var = "Barcode")
# og.table<-og.table[match(env.mat$Barcode,rownames(og.table)),]
# 
# # Get OG description
# og.desc<-fread("../lib/NOG.annotations.tsv",sep="\t",header=F,data.table = F)
# og.desc<-og.desc %>% select(OG=V2,description=V6) %>% mutate(OG=gsub("ENOG41","",OG))
# cog.desc<-fread("../lib/cognames2003-2014.tab",sep="\t",header=T,data.table = F)
# desc<-left_join(og.desc,cog.desc,by=c("OG"="# COG"))
# desc$name[is.na(desc$name)]<-desc$description[is.na(desc$name)]
# 
# # Mean per-cell abundance
# mean.ab<-og.table %>%
#   rownames_to_column(var="Barcode") %>%
#   gather(key="OG",value="ab",-Barcode) %>%
#   group_by(OG) %>%
#   summarise(ab=mean(ab))
# occ<-og.table %>%
#   rownames_to_column(var="Barcode") %>%
#   gather(key="OG",value="ab",-Barcode) %>%
#   filter(ab>0) %>%
#   group_by(OG) %>%
#   summarise(occ=n())
# 
# # Latitudinal niche value
# res.lat<-niche.val(og.table,env.mat$Absolute.Latitude,n = iter)
# toplot.lat<-res.lat %>%
#   rownames_to_column(var="OG") %>%
#   left_join(mean.ab,by="OG") %>%
#   left_join(occ,by="OG") %>%
#   dplyr::arrange(desc(ab)) %>%
#   mutate(sign=fct_recode(sign,"POLAR"="HIGHER","NON-POLAR"="LOWER","N.S."="NON SIGNIFICANT")) %>%
#   filter(OG!="unknown") %>%
#   left_join(desc,by="OG") %>%
#   mutate(OG=paste(OG,name))
# 
# 
# # Plots
# g<-ggplot() +
#   geom_errorbar(data=toplot.lat %>% filter(occ<142) %>% arrange(desc(observed)) %>% head(n=30),aes(x=fct_reorder(OG,ab),ymax=uppCI,ymin=lowCI),col="gray40") +
#   geom_point(data=toplot.lat %>% filter(occ<142)  %>% arrange(desc(observed)) %>% head(n=30),aes(x=fct_reorder(OG,ab),y=mean.simulated),col="gray40") +
#   geom_point(data=toplot.lat %>% filter(occ<142)  %>% arrange(desc(observed)) %>% head(n=30),aes(x=fct_reorder(OG,ab),y=observed,fill=sign,size=100*ab),shape=21) +
#   coord_flip() +
#   theme_bw() +
#   theme(legend.direction = "vertical") +
#   scale_size(name="Per-cell abundance") +
#   xlab(NULL) +
#   ylab("Abundance-weighted\nabsolute latitude") +
#   scale_fill_manual(values=c("#1C6EAC","#F48B12","#BEBEBE")) +
#   ylim(0,90)
# ggsave(filename = "../../results/paper_figures/review/FigSX_taxonomy_latitude_associations_OG.pdf",g,width=15,height=5)
