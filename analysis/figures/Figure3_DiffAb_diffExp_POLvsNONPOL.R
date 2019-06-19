#library(propagate)
library(vegan)
library(tidyverse)
library(ggrepel)
library(data.table)
library(patchwork)
#source("http://sunagawalab.ethz.ch/sushilab/resources/R/sushi.palette.R")
source("../data/lib/sushipal.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])

###############
## LOAD DATA ##
###############

metaG.norm.match.log2<-fread("../data/processed/KO_metaG.norm.match.log2.txt",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("../data/processed/KO_metaG.norm.match.txt",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("../data/processed/KO_metaT.norm.match.log2.txt",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("../data/processed/KO_metaT.norm.match.txt",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("../data/processed/KO_ratio.mat.txt",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]
env.mat.match<-fread("../data/processed/KO_env.mat.match.txt",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match$epi<-env.mat.match$Layer
levels(env.mat.match$epi)<-c("EPI","MES","EPI","EPI")
env.mat.match<-env.mat.match[match(rownames(metaG.norm.match),env.mat.match$Barcode),]

####################
## Load KEGG info ##
####################
# ko-reaction
ko.reaction.link<-read.table("../data/lib/kegg/ko_reaction.link",sep="\t",header=F)
ko.reaction.link$V1<-substr(ko.reaction.link$V1,4,100)
ko.reaction.link$V2<-substr(ko.reaction.link$V2,4,100)
ko.reaction.link<-droplevels(ko.reaction.link)


# ko-module
ko.module.link<-read.table("../data/lib/kegg/ko_module.link",sep="\t",header=F)
ko.module.link$V1<-substr(ko.module.link$V1,4,100)
ko.module.link$V2<-substr(ko.module.link$V2,4,100)
ko.module.link<-droplevels(ko.module.link)

# ko-pathway
ko.pathway.link<-read.table("../data/lib/kegg/ko_pathway.link",sep="\t",header=F)
ko.pathway.link$V1<-substr(ko.pathway.link$V1,6,100)
ko.pathway.link$V2<-substr(ko.pathway.link$V2,4,100)
ko.pathway.link<-ko.pathway.link[grep("map",ko.pathway.link$V1),]
ko.pathway.link<-droplevels(ko.pathway.link)

# module-pathway
module.pathway.link<-read.table("../data/lib/kegg/module_pathway.link",sep="\t",header=F)
module.pathway.link$V1<-substr(module.pathway.link$V1,6,100)
module.pathway.link$V2<-substr(module.pathway.link$V2,4,100)
module.pathway.link<-module.pathway.link[grep("map",module.pathway.link$V1),]
module.pathway.link<-droplevels(module.pathway.link)


# pathway info
pathway.info<-read.table("../data/lib/kegg/pathway.list",sep="\t",header=F,quote="")
pathway.info$V1<-substr(pathway.info$V1,6,100)

# module info
module.info<-read.table("../data/lib/kegg/module.list",sep="\t",header=F,quote="")
module.info$V1<-substr(module.info$V1,4,100)
module.info$V2<-as.character(module.info$V2)
dupl<-which(module.info$V2 %in% module.info$V2[which(duplicated(module.info$V2))])
module.info$V2[dupl]<-paste0(module.info$V2[dupl]," (",module.info$V1[dupl],")")
module.info$V2<-as.factor(module.info$V2)

# reaction info
reaction.info<-read.table("../data/lib/kegg/reaction.list",sep="\t",header=F,quote="")
reaction.info$V1<-substr(reaction.info$V1,4,100)
reaction.info$V2<-as.character(reaction.info$V2)
dupl<-which(reaction.info$V2 %in% reaction.info$V2[which(duplicated(reaction.info$V2))])
reaction.info$V2[dupl]<-paste0(reaction.info$V2[dupl]," (",reaction.info$V1[dupl],")")
reaction.info$V2<-as.factor(reaction.info$V2)


# ko info
ko.info<-read.table("../data/lib/kegg/ko.list",sep="\t",header=F,quote="")
ko.info$V1<-substr(ko.info$V1,4,100)

# add info
ko.pathway.link$ko.description<-ko.info$V2[match(ko.pathway.link$V2,ko.info$V1)]
ko.pathway.link$pathway.description<-pathway.info$V2[match(ko.pathway.link$V1,pathway.info$V1)]

ko.module.link$ko.description<-ko.info$V2[match(ko.module.link$V2,ko.info$V1)]
ko.module.link$module.description<-module.info$V2[match(ko.module.link$V1,module.info$V1)]

ko.reaction.link$ko.description<-ko.info$V2[match(ko.reaction.link$V2,ko.info$V1)]
ko.reaction.link$reaction.description<-reaction.info$V2[match(ko.reaction.link$V1,reaction.info$V1)]

# Build a ko to module link for the KOs univocally assigned to a single module
ko.module.tab<-ko.module.link %>%
  dplyr::select(V1,V2) %>%
  distinct() %>%
  table()
ko.module.link.uniquemod<-ko.module.link[match(names(which(colSums(ko.module.tab)==1)),ko.module.link$V2),]

# Build a ko to pathway link for the KOs univocally assigned to a single pathway
ko.pathway.tab<-ko.pathway.link %>%
  filter(!grepl("map011",V1)) %>% # remove global maps
  filter(!grepl("map012",V1)) %>% # remove overview maps
  dplyr::select(V1,V2) %>%
  distinct() %>%
  table()
ko.pathway.link.uniquemod<-ko.pathway.link[match(names(which(colSums(ko.pathway.tab)==1)),ko.pathway.link$V2),]

# Build a ko to reaction link for the KOs univocally assigned to a single reaction
ko.reaction.tab<-ko.reaction.link %>%
  dplyr::select(V1,V2) %>%
  distinct() %>%
  table()
ko.reaction.link.uniquemod<-ko.reaction.link[match(names(which(colSums(ko.reaction.tab)==1)),ko.reaction.link$V2),]

# Load selected pathways and modules
sel.path<-read.table("../data/lib/kegg/selected_pathways.txt",header=F,sep="\t")
sel.path<-fread("../data/lib/kegg/selected_pathways.txt",sep="\t",header=F,data.table = F)
sel.path<-sel.path %>%
  left_join(pathway.info,by=c("V1"="V1")) %>%
  dplyr::rename(map=V1,Category=V2.x,pathway.description=V2.y)

sel.mod<-read.table("../data/lib/kegg/selected_modules.txt",header=F,sep="\t")

# Load selected marker genes
sel.genes<-read.table("../data/lib/kegg/Selection_Genes.txt",header=T,sep="\t")
sel.genes<-sel.genes[which(sel.genes$KO %in% colnames(metaG.norm.match.log2)),]
sel.genes$Enzyme<-paste(paste(sel.genes$KO,sel.genes$Enzyme)," (",sel.genes$Module,")",sep="")
#sel.genes2<-read.table("../data/lib/kegg/List KO used in Shini paper Science 2015_modifiedbySilvia.txt",header=T,sep="\t")

##############
## Analysis ##
##############
# Want to reduce everything to EPI?
metaG.norm.match.log2<-metaG.norm.match.log2[env.mat.match$Layer %in% c("SRF","DCM"),]
metaT.norm.match.log2<-metaT.norm.match.log2[env.mat.match$Layer %in% c("SRF","DCM"),]
ratio.mat<-ratio.mat[env.mat.match$Layer %in% c("SRF","DCM"),]
env.mat.match<-env.mat.match[env.mat.match$Layer %in% c("SRF","DCM"),]


# Gather data
metaG.norm.match.core.gath<-metaG.norm.match.log2 %>%
  rownames_to_column(var="Barcode") %>%
  gather(key="OG",value="Abundance",-Barcode) %>%
  left_join(env.mat.match,by="Barcode")
metaT.norm.match.core.gath<-metaT.norm.match.log2 %>%
  rownames_to_column(var="Barcode") %>%
  gather(key="OG",value="Transcription",-Barcode) %>%
  left_join(env.mat.match,by="Barcode")
Exp.norm.match.core.gath<-ratio.mat %>%
  rownames_to_column(var="Barcode") %>%
  gather(key="OG",value="Expression",-Barcode) %>%
  left_join(env.mat.match,by="Barcode")

# Compute mean for polar and non-polar
metaG.norm.match.core.gath<-metaG.norm.match.core.gath %>%
  group_by(OG,polar) %>%
  mutate(mean_Abundance=mean(Abundance))
metaT.norm.match.core.gath<-metaT.norm.match.core.gath %>%
  group_by(OG,polar) %>%
  mutate(mean_Transcription=mean(Transcription))
Exp.norm.match.core.gath<-Exp.norm.match.core.gath %>%
  group_by(OG,polar) %>%
  mutate(mean_Expression=mean(Expression))

# Merge data
metaT.norm.match.core.gath<-metaT.norm.match.core.gath %>% dplyr::select(mean_Transcription,Transcription)
Exp.norm.match.core.gath<-Exp.norm.match.core.gath %>% dplyr::select(mean_Expression,Expression)
merged.core.gath<-bind_cols(metaG.norm.match.core.gath,metaT.norm.match.core.gath,Exp.norm.match.core.gath) %>%
  #left_join(ko.pathway.link,by=c("OG"="V2")) %>%
  filter(OG!="unknown") %>%
  as.data.frame()

# Compute mean(polar) - mean(non polar)
tmp.ab<-merged.core.gath %>%
  dplyr::select(OG,mean_Abundance,mean_Expression,mean_Transcription,polar) %>%
  distinct() %>%
  dplyr::select(OG,mean_Abundance,polar) %>%
  spread(key = "polar",value="mean_Abundance") %>%
  mutate(Ab.ratio=polar-`non polar`) %>%
  dplyr::select(OG,Ab.ratio)
tmp.exp<-merged.core.gath %>%
  dplyr::select(OG,mean_Abundance,mean_Expression,mean_Transcription,polar) %>%
  distinct() %>%
  dplyr::select(OG,mean_Expression,polar) %>%
  spread(key = "polar",value="mean_Expression") %>%
  mutate(Exp.ratio=polar-`non polar`) %>%
  dplyr::select(OG,Exp.ratio)
tmp.tra<-merged.core.gath %>%
  dplyr::select(OG,mean_Abundance,mean_Expression,mean_Transcription,polar) %>%
  distinct() %>%
  dplyr::select(OG,mean_Transcription,polar) %>%
  spread(key = "polar",value="mean_Transcription") %>%
  mutate(Tra.ratio=polar-`non polar`) %>%
  dplyr::select(OG,Tra.ratio)

# Compute differential abundance and expression test
fun<-function(x,env.var){
  res<-wilcox.test(x~env.var)
  tmp<-c(pval=res$p.value,mean.lev1=mean(x[env.var==levels(env.var)[1]]),mean.lev2=mean(x[env.var==levels(env.var)[2]]))
  names(tmp)<-c("pval",paste("mean",levels(env.var)[1],sep="."),paste("mean",levels(env.var)[2],sep="."))
  tmp
}

res.ab<-as.data.frame(t(apply(metaG.norm.match.log2,2,fun,env.var=env.mat.match$polar)))
res.ab$pval.fdr<-p.adjust(res.ab$pval,method = "bonferroni")
sign<-NULL
for (i in 1:nrow(res.ab)){
  if (is.na(res.ab[i,'pval.fdr'])) sign<-c(sign,"N.S.") else if (res.ab[i,'pval.fdr']<=0.05 & res.ab[i,'mean.polar']>res.ab[i,'mean.non polar']) sign<-c(sign,"polar") else if (res.ab[i,'pval.fdr']<=0.05 & res.ab[i,'mean.polar']<res.ab[i,'mean.non polar']) sign<-c(sign,"non polar") else if (res.ab[i,'pval.fdr']>0.05) sign<-c(sign,"N.S.")
}
res.ab$sign<-sign

res.exp<-as.data.frame(t(apply(ratio.mat,2,fun,env.var=env.mat.match$polar)))
res.exp$pval.fdr<-p.adjust(res.exp$pval,method = "bonferroni")
sign<-NULL
for (i in 1:nrow(res.exp)){
  if (is.na(res.exp[i,'pval.fdr'])) sign<-c(sign,"N.S.") else if (res.exp[i,'pval.fdr']<=0.05 & res.exp[i,'mean.polar']>res.exp[i,'mean.non polar']) sign<-c(sign,"polar") else if (res.exp[i,'pval.fdr']<=0.05 & res.exp[i,'mean.polar']<res.exp[i,'mean.non polar']) sign<-c(sign,"non polar") else if (res.exp[i,'pval.fdr']>0.05) sign<-c(sign,"N.S.")
}
res.exp$sign<-sign

res.tra<-as.data.frame(t(apply(metaT.norm.match.log2,2,fun,env.var=env.mat.match$polar)))
res.tra$pval.fdr<-p.adjust(res.tra$pval,method = "bonferroni")
sign<-NULL
for (i in 1:nrow(res.tra)){
  if (is.na(res.tra[i,'pval.fdr'])) sign<-c(sign,"N.S.") else if (res.tra[i,'pval.fdr']<=0.05 & res.tra[i,'mean.polar']>res.tra[i,'mean.non polar']) sign<-c(sign,"polar") else if (res.tra[i,'pval.fdr']<=0.05 & res.tra[i,'mean.polar']<res.tra[i,'mean.non polar']) sign<-c(sign,"non polar") else if (res.tra[i,'pval.fdr']>0.05) sign<-c(sign,"N.S.")
}
res.tra$sign<-sign

colnames(res.ab)<-paste("ab",colnames(res.ab),sep=".")
colnames(res.exp)<-paste("exp",colnames(res.exp),sep=".")
colnames(res.tra)<-paste("tra",colnames(res.tra),sep=".")

wilcox.df<-res.ab %>%
  rownames_to_column(var="OG") %>%
  bind_cols(res.exp) %>%
  dplyr::select(OG,ab.pval,ab.pval.fdr,ab.sign,exp.pval,exp.pval.fdr,exp.sign)

# Plot Selected Genes
toplot.df<-left_join(tmp.ab,tmp.exp,by="OG") %>%
  left_join(sel.genes,by=c("OG"="KO")) %>%
  left_join(wilcox.df,by="OG") %>%
  filter(!is.na(Enzyme))
toplot.df$toplot<-"no"
toplot.df$toplot[toplot.df$exp.sign!="N.S." | toplot.df$ab.sign!="N.S."]<-"yes"
toplot.df$toplot[abs(toplot.df$Ab.ratio)<1 & abs(toplot.df$Exp.ratio)<1]<-"no"

toplot.df<- toplot.df %>%
  filter(Metabolism %in% c("CARBON FIXATION","NITROGEN","PHOTOSYNTHESIS","SULFUR")) %>%
  droplevels()
toplot.df$Metabolism<-factor(toplot.df$Metabolism,levels=c("CARBON FIXATION","PHOTOSYNTHESIS","NITROGEN","SULFUR"))

g<-ggplot() +
  geom_hline(yintercept = 0,color="gray") +
  geom_vline(xintercept = 0,color="gray") +
  geom_point(data=toplot.df %>% filter(toplot=="no"),aes(x=Ab.ratio,y=Exp.ratio),size=4,col="gray") +
  geom_point(data=toplot.df %>% filter(toplot=="yes"),aes(x=Ab.ratio,y=Exp.ratio,col=exp.sign,shape=ab.sign),size=4) +
  geom_text_repel(data=toplot.df %>% filter(toplot=="yes"),aes(x=Ab.ratio,y=Exp.ratio,label=OG),size=2,force=10) +
  theme_bw() +
  facet_wrap(~Metabolism) +
  xlab("Difference in mean abundance (POLAR - NON POLAR)") +
  ylab("Difference in mean expression (POLAR - NON POLAR)") +
  scale_color_manual(name="Significance (expression)",values=c("darkgoldenrod1","#FB8072FF","#8DD3C7FF")) +
  scale_shape_manual(name="Significance (abundance)",values=c(19,15,17))
ggsave(filename = "../results/POL_NONPOL_diff.pdf",g,width =160,height=160,units = "mm")
