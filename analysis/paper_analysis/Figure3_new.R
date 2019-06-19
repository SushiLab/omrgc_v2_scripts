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

metaG.norm.match.log2<-fread("zcat < ../data/processed/KO_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("zcat < ../data/processed/KO_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("zcat < ../data/processed/KO_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("zcat < ../data/processed/KO_metaT.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("zcat < ../data/processed/KO_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]
env.mat.match<-fread("zcat < ../data/processed/KO_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match$epi<-env.mat.match$Layer
levels(env.mat.match$epi)<-c("EPI","MES","EPI","EPI")
env.mat.match<-env.mat.match[match(rownames(metaG.norm.match),env.mat.match$Barcode),]
env.mat.match$date<-as.POSIXct(sapply(strsplit(as.character(env.mat.match$Event.date),"T"),"[[",1),format="%Y-%m-%d")

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

# Compute mean for EPI and MES
metaG.norm.match.core.gath<-metaG.norm.match.core.gath %>%
  group_by(OG,epi) %>%
  mutate(mean_Abundance=mean(Abundance))
metaT.norm.match.core.gath<-metaT.norm.match.core.gath %>%
  group_by(OG,epi) %>%
  mutate(mean_Transcription=mean(Transcription))
Exp.norm.match.core.gath<-Exp.norm.match.core.gath %>%
  group_by(OG,epi) %>%
  mutate(mean_Expression=mean(Expression))

# Merge data
metaT.norm.match.core.gath<-metaT.norm.match.core.gath %>% dplyr::select(mean_Transcription,Transcription)
Exp.norm.match.core.gath<-Exp.norm.match.core.gath %>% dplyr::select(mean_Expression,Expression)
merged.core.gath<-bind_cols(metaG.norm.match.core.gath,metaT.norm.match.core.gath,Exp.norm.match.core.gath) %>%
  #left_join(ko.pathway.link,by=c("OG"="V2")) %>%
  filter(OG!="unknown") %>%
  as.data.frame()
##########################################################################################
toplot.df<-merged.core.gath %>%
  left_join(sel.genes,by=c("OG"="KO")) %>%
  filter(!is.na(Enzyme)) %>%
  filter(Metabolism %in% c("CARBON FIXATION","NITROGEN","PHOTOSYNTHESIS","SULFUR"))
tmp<-toplot.df %>%
  select(Sample_name,OG,Abundance,Expression,Transcription) %>%
  gather(key="data.type",value = "value",-OG,-Sample_name)
toplot.df<-tmp %>%
  left_join(sel.genes,by=c("OG"="KO")) %>%
  filter(!is.na(Enzyme)) %>%
  filter(Metabolism %in% c("CARBON FIXATION","NITROGEN","PHOTOSYNTHESIS","SULFUR")) %>%
  left_join(env.mat.match,by="Sample_name")

ggplot(toplot.df %>% filter(epi=="EPI"),aes(x=Enzyme,y=value,fill=polar)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  facet_grid(data.type~Metabolism,scales = "free",space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) +
  scale_fill_manual(values=c("#e59d8d","#81b0d2"))
ggsave(file="../results/SFig_extra_POLAR.pdf")

ggplot(toplot.df,aes(x=Enzyme,y=value,fill=epi)) +
  geom_violin(draw_quantiles = 0.5,scale = "width") +
  facet_grid(data.type~Metabolism,scales = "free",space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) +
  scale_fill_manual(values=c("#ab8b17","#136fba"))
ggsave(file="../results/SFig_extra_EPI.pdf")

###
ggplot(toplot.df %>% filter(OG %in% c("K02588","K02586")),aes(x=NO2NO3,y=value)) +
  geom_point(aes(col=epi),size=2,alpha=0.6) +
  geom_smooth(span=0.2,se=F,col="black") +
  facet_grid(data.type~str_wrap(Enzyme,20),scales = "free",space = "free_x") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) +
  scale_color_manual(values=c("#ab8b17","#136fba")) +
  scale_x_continuous(breaks=c(0.2,1,5,10,20,40),trans = "pseudo_log") +
  xlab("Nitrate and nitrite (uM/L)")
ggsave(file="../results/SFig_nifHD.pdf",scale = 0.6)

a<-toplot.df %>% filter(OG %in% c("K02588","K02586")) %>% filter(data.type=="Expression" & value>2.5 & epi=="EPI")
range(abs(a$Latitude))
median(abs(a$Latitude))
median(a$NO2NO3,na.rm = T)
median(env.mat.match$NO2NO3[env.mat.match$epi=="EPI"],na.rm = T)

##########################################################################################
# Compute mean(EPI) - mean(MES)
tmp.ab<-merged.core.gath %>%
  dplyr::select(OG,mean_Abundance,mean_Expression,mean_Transcription,epi) %>%
  distinct() %>%
  dplyr::select(OG,mean_Abundance,epi) %>%
  spread(key = "epi",value="mean_Abundance") %>%
  mutate(Ab.ratio=EPI-MES) %>%
  dplyr::select(OG,Ab.ratio)
tmp.exp<-merged.core.gath %>%
  dplyr::select(OG,mean_Abundance,mean_Expression,mean_Transcription,epi) %>%
  distinct() %>%
  dplyr::select(OG,mean_Expression,epi) %>%
  spread(key = "epi",value="mean_Expression") %>%
  mutate(Exp.ratio=EPI-MES) %>%
  dplyr::select(OG,Exp.ratio)

# Compute differential abundance and expression test
fun<-function(x,env.var){
  res<-wilcox.test(x~env.var)
  tmp<-c(pval=res$p.value,mean.lev1=mean(x[env.var==levels(env.var)[1]]),mean.lev2=mean(x[env.var==levels(env.var)[2]]))
  names(tmp)<-c("pval",paste("mean",levels(env.var)[1],sep="."),paste("mean",levels(env.var)[2],sep="."))
  tmp
}

#fun2<-function(x,env.var,block){
#  library(coin)
#  res<-wilcox_test(x~env.var | block)
#  tmp<-c(pval=pvalue(res),mean.lev1=mean(x[env.var==levels(env.var)[1]]),mean.lev2=mean(x[env.var==levels(env.var)[2]]))
#  names(tmp)<-c("pval",paste("mean",levels(env.var)[1],sep="."),paste("mean",levels(env.var)[2],sep="."))
#  tmp
#}

res.ab<-as.data.frame(t(apply(metaG.norm.match.log2,2,fun,env.var=env.mat.match$epi)))
#res.ab2<-as.data.frame(t(apply(metaG.norm.match.log2,2,fun2,env.var=env.mat.match$epi)))
res.ab$pval.fdr<-p.adjust(res.ab$pval,method = "bonferroni")
sign<-NULL
for (i in 1:nrow(res.ab)){
  if (is.na(res.ab[i,'pval.fdr'])) sign<-c(sign,"N.S.") else if (res.ab[i,'pval.fdr']<=0.05 & res.ab[i,'mean.EPI']>res.ab[i,'mean.MES']) sign<-c(sign,"EPI") else if (res.ab[i,'pval.fdr']<=0.05 & res.ab[i,'mean.EPI']<res.ab[i,'mean.MES']) sign<-c(sign,"MES") else if (res.ab[i,'pval.fdr']>0.05) sign<-c(sign,"N.S.")
}
res.ab$sign<-sign

res.exp<-as.data.frame(t(apply(ratio.mat,2,fun,env.var=env.mat.match$epi)))
res.exp$pval.fdr<-p.adjust(res.exp$pval,method = "bonferroni")
sign<-NULL
for (i in 1:nrow(res.exp)){
  if (is.na(res.exp[i,'pval.fdr'])) sign<-c(sign,"N.S.") else if (res.exp[i,'pval.fdr']<=0.05 & res.exp[i,'mean.EPI']>res.exp[i,'mean.MES']) sign<-c(sign,"EPI") else if (res.exp[i,'pval.fdr']<=0.05 & res.exp[i,'mean.EPI']<res.exp[i,'mean.MES']) sign<-c(sign,"MES") else if (res.exp[i,'pval.fdr']>0.05) sign<-c(sign,"N.S.")
}
res.exp$sign<-sign

colnames(res.ab)<-paste("ab",colnames(res.ab),sep=".")
colnames(res.exp)<-paste("exp",colnames(res.exp),sep=".")

wilcox.df<-res.ab %>%
  rownames_to_column(var="OG") %>%
  bind_cols(res.exp) %>%
  dplyr::select(OG,ab.pval,ab.pval.fdr,ab.sign,exp.pval,exp.pval.fdr,exp.sign)

# Compute occurrence and eveness
source("../data/lib/eveness_functions.R")
metaG.norm.match.occ<-apply(metaG.norm.match,2,function(x){length(which(x>0))/length(x)})
metaT.norm.match.occ<-apply(metaT.norm.match,2,function(x){length(which(x>0))/length(x)})
metaG.norm.match.even<-apply(metaG.norm.match,2,SimpE)
metaT.norm.match.even<-apply(metaT.norm.match,2,SimpE)

occ.df<-data.frame(OG=colnames(metaG.norm.match),metaG.occ=metaG.norm.match.occ,metaT.occ=metaT.norm.match.occ)

# Plot Selected Genes
toplot.df<-left_join(tmp.ab,tmp.exp,by="OG") %>%
  left_join(sel.genes,by=c("OG"="KO")) %>%
  left_join(wilcox.df,by="OG") %>%
  filter(!is.na(Enzyme)) %>%
  left_join(occ.df,by="OG") %>%
  mutate(sign.combo=interaction(ab.sign,exp.sign))
toplot.df$toplot<-"no"
toplot.df$toplot[toplot.df$exp.sign!="N.S." | toplot.df$ab.sign!="N.S."]<-"yes"
toplot.df$toplot[abs(toplot.df$Ab.ratio)<1 & abs(toplot.df$Exp.ratio)<1]<-"no"
toplot.df$ab.sign<-factor(toplot.df$ab.sign,levels = c("N.S.","EPI","MES"))
toplot.df$exp.sign<-factor(toplot.df$exp.sign,levels = c("N.S.","EPI","MES"))

toplot.df<- toplot.df %>%
  filter(Metabolism %in% c("CARBON FIXATION","NITROGEN","PHOTOSYNTHESIS","SULFUR")) %>%
  droplevels()
toplot.df$Metabolism<-factor(toplot.df$Metabolism,levels=c("CARBON FIXATION","PHOTOSYNTHESIS","NITROGEN","SULFUR"))


#g<-ggplot() +
#  geom_hline(yintercept = 0,color="gray") +
#  geom_vline(xintercept = 0,color="gray") +
#  geom_point(data=toplot.df %>% filter(toplot=="no"),aes(x=Ab.ratio,y=Exp.ratio),size=4,col="gray") +
#  geom_point(data=toplot.df %>% filter(toplot=="yes"),aes(x=Ab.ratio,y=Exp.ratio,col=exp.sign,shape=ab.sign),size=4) +
#  geom_text_repel(data=toplot.df %>% filter(toplot=="yes"),aes(x=Ab.ratio,y=Exp.ratio,label=OG),size=2,force=10) +
#  theme_bw() +
#  facet_wrap(~Metabolism) +
#  xlab("Difference in mean abundance (EPI - MES)") +
#  ylab("Difference in mean expression (EPI - MES)") +
#  scale_color_manual(name="Significance (expression)",values=c("orange",sushi.palette()[c(2,1)])) +
#  scale_shape_manual(name="Significance (abundance)",values=c(19,15,17))
#ggsave(filename = "../results/MES_EPI_diff.pdf",g,width =160,height=160,units = "mm")

lims<-range(c(toplot.df$Ab.ratio,toplot.df$Exp.ratio))
toplot.df<-toplot.df %>% separate(Enzyme,"nom",sep = ";",extra="drop") %>% mutate(nom=paste(nom,Module))
g1<-ggplot(data=toplot.df,aes(x=nom,y=Ab.ratio,col=ab.sign,xend=nom,yend=0)) +
  geom_hline(yintercept = 0) +
  geom_point(size=3) +
  geom_segment() +
  coord_flip() +
  facet_grid(Metabolism~.,scales = "free",space = "free") +
  theme_bw() +
  ylim(lims) +
  xlab(NULL) +
  ylab("Difference in mean gene abundance") +
  theme(legend.position = "none",strip.background = element_blank(),strip.text.y = element_blank()) +
  scale_color_manual(values=c("#BEBEBE","#AB8B17","#81B0D2"))
g2<-ggplot(data=toplot.df,aes(x=nom,y=Exp.ratio,col=exp.sign,xend=nom,yend=0)) +
  geom_hline(yintercept = 0) +
  geom_point(size=3) +
  geom_segment() +
  coord_flip() +
  facet_grid(Metabolism~.,scales = "free",space = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank(),strip.background = element_blank(),strip.text.y = element_blank()) +
  ylim(lims) +
  ylab("Difference in mean gene expression") +
  xlab(NULL) +
  theme(legend.position = "top") +
  scale_color_manual(values=c("#BEBEBE","#AB8B17","#81B0D2"))

tmp<-metaT.norm.match.log2[,which(colnames(metaT.norm.match.log2) %in% toplot.df$OG)] %>%
  rownames_to_column(var="Barcode") %>%
  gather(key = "OG",value="transcription",-Barcode) %>%
  left_join(env.mat.match,by="Barcode") %>%
  left_join(sel.genes,by = c("OG"="KO")) %>%
  mutate(Metabolism=fct_relevel(Metabolism,"CARBON FIXATION","PHOTOSYNTHESIS","NITROGEN","SULFUR"))

tmp<-tmp %>% separate(Enzyme,"nom",sep = ";",extra="drop") %>% mutate(nom=paste(nom,Module))

g3<-ggplot(data=tmp,aes(x=nom,y=transcription,fill=epi)) +
  geom_boxplot(outlier.size = 0.5) +
  #geom_violin(draw_quantiles = 0.5,scale = "width") +
  coord_flip() +
  facet_grid(Metabolism~.,scales = "free",space = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  ylab("Transcript abundance (log2-transformed)") +
  xlab(NULL) +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("#AB8B17","#81B0D2"))

  
g1 | g2 | g3
