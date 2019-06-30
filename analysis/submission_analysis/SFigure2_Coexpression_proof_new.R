library(caret)
library(tidyverse)
library(data.table)
library(vegan)
library(patchwork)
source("../data/lib/sushipal.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])
min.occ<-0.1
st<-0.02 # step for binning
occ<-function(x){length(x[x>0])/length(x)}
get.roc<-function(true.vec,r.vec,r.cutoff){
  pred.vec<-sapply(r.vec,function(x){x>=r.cutoff})
  tmp<-confusionMatrix(factor(pred.vec,levels=c("TRUE","FALSE")),factor(true.vec,levels=c("TRUE","FALSE")),positive = "TRUE")
  c(r.cutoff=r.cutoff,tmp$byClass['Sensitivity'],tmp$byClass['Specificity'])
}
###############
## LOAD DATA ##
###############

metaG.norm.match.log2<-fread("zcat < ../data/processed/NOG_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("zcat < ../data/processed/NOG_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("zcat < ../data/processed/NOG_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("zcat < ../data/processed/NOG_metaT.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("zcat < ../data/processed/NOG_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]

env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match<-env.mat.match[match(rownames(metaG.norm.match),env.mat.match$Barcode),]

#####################
## LOAD NOG-KO map ##
#####################
nog.ko.map<-fread("../data/lib/NOG.KO.map",header=T,sep="\t",data.table = F)

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

##################
## Module level ##
##################

# Add module to the map
nog.ko.map.module<-nog.ko.map %>%
  left_join(ko.module.link,by=c("KO"="V2")) %>%
    dplyr::select(NOG,KO,module=V1,ko.description,module.description) %>%
  filter(!is.na(module))

# Filter by abundance-based occurrence
ab.occ<-apply(metaG.norm.match,2,occ)
ab.mat<-metaG.norm.match.log2[,which(ab.occ>min.occ)]
exp.mat<-ratio.mat[,which(ab.occ>min.occ)]
tra.mat<-metaT.norm.match.log2[,which(ab.occ>min.occ)]

# Get only the NOGs associated to KOs
ab.mat<-ab.mat[,which(colnames(ab.mat) %in% nog.ko.map.module$NOG)]
exp.mat<-exp.mat[,which(colnames(exp.mat) %in% nog.ko.map.module$NOG)]
tra.mat<-tra.mat[,which(colnames(tra.mat) %in% nog.ko.map.module$NOG)]

# Compute the correlation
ab.mat.corr<-cor(ab.mat)
exp.mat.corr<-cor(exp.mat)
tra.mat.corr<-cor(tra.mat)
diag(ab.mat.corr)<-NA
diag(exp.mat.corr)<-NA
diag(tra.mat.corr)<-NA

# Gather correlation matrix
ab.mat.corr.gath<-ab.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)
exp.mat.corr.gath<-exp.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)
tra.mat.corr.gath<-tra.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)

# Get only the best correlated NOG for each NOG
ab.mat.corr.gath<-ab.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()
exp.mat.corr.gath<-exp.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()
tra.mat.corr.gath<-tra.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()

# Add a variable with common modules (if exist)
fun<-function(x){
  mod1<-nog.ko.map.module$module[which(nog.ko.map.module$NOG %in% x[1])]
  mod2<-nog.ko.map.module$module[which(nog.ko.map.module$NOG %in% x[2])]
  if (any(mod1 %in% mod2)) paste(mod1[which(mod1 %in% mod2)],collapse=",") else NA
}
ab.mat.corr.gath$common.modules<-apply(ab.mat.corr.gath,1,fun)
exp.mat.corr.gath$common.modules<-apply(exp.mat.corr.gath,1,fun)
tra.mat.corr.gath$common.modules<-apply(tra.mat.corr.gath,1,fun)

# Add a variable dichotomous variable (module.shared / not shared)
ab.mat.corr.gath$module.shared<-!is.na(ab.mat.corr.gath$common.modules)
exp.mat.corr.gath$module.shared<-!is.na(exp.mat.corr.gath$common.modules)
tra.mat.corr.gath$module.shared<-!is.na(tra.mat.corr.gath$common.modules)

# Compute Specificity 
roc.ab<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=ab.mat.corr.gath$module.shared,r.vec=ab.mat.corr.gath$r,r.cutoff=x)})),data.type="Abundance")
roc.exp<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=exp.mat.corr.gath$module.shared,r.vec=exp.mat.corr.gath$r,r.cutoff=x)})),data.type="Expression")
roc.tra<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=tra.mat.corr.gath$module.shared,r.vec=tra.mat.corr.gath$r,r.cutoff=x)})),data.type="Transcription")
toplot.module<-rbind(roc.ab,roc.exp,roc.tra)
toplot.module$kegg.level<-"Module"

##################
## Pathway level ##
##################

# Add pathway to the map
nog.ko.map.pathway<-nog.ko.map %>%
  left_join(ko.pathway.link,by=c("KO"="V2")) %>%
  dplyr::select(NOG,KO,pathway=V1,ko.description,pathway.description) %>%
  filter(!is.na(pathway))

# Filter by abundance-based occurrence
ab.occ<-apply(metaG.norm.match,2,occ)
ab.mat<-metaG.norm.match.log2[,which(ab.occ>min.occ)]
exp.mat<-ratio.mat[,which(ab.occ>min.occ)]
tra.mat<-metaT.norm.match.log2[,which(ab.occ>min.occ)]

# Get only the NOGs associated to KOs
ab.mat<-ab.mat[,which(colnames(ab.mat) %in% nog.ko.map.pathway$NOG)]
exp.mat<-exp.mat[,which(colnames(exp.mat) %in% nog.ko.map.pathway$NOG)]
tra.mat<-tra.mat[,which(colnames(tra.mat) %in% nog.ko.map.pathway$NOG)]

# Compute the correlation
ab.mat.corr<-cor(ab.mat)
exp.mat.corr<-cor(exp.mat)
tra.mat.corr<-cor(tra.mat)
diag(ab.mat.corr)<-NA
diag(exp.mat.corr)<-NA
diag(tra.mat.corr)<-NA

# Gather correlation matrix
ab.mat.corr.gath<-ab.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)
exp.mat.corr.gath<-exp.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)
tra.mat.corr.gath<-tra.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)

# Get only the best correlated NOG for each NOG
ab.mat.corr.gath<-ab.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()
exp.mat.corr.gath<-exp.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()
tra.mat.corr.gath<-tra.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()

# Add a variable with common pathways (if exist)
fun<-function(x){
  mod1<-nog.ko.map.pathway$pathway[which(nog.ko.map.pathway$NOG %in% x[1])]
  mod2<-nog.ko.map.pathway$pathway[which(nog.ko.map.pathway$NOG %in% x[2])]
  if (any(mod1 %in% mod2)) paste(mod1[which(mod1 %in% mod2)],collapse=",") else NA
}
ab.mat.corr.gath$common.pathways<-apply(ab.mat.corr.gath,1,fun)
exp.mat.corr.gath$common.pathways<-apply(exp.mat.corr.gath,1,fun)
tra.mat.corr.gath$common.pathways<-apply(tra.mat.corr.gath,1,fun)

# Add a variable dichotomous variable (pathway.shared / not shared)
ab.mat.corr.gath$pathway.shared<-!is.na(ab.mat.corr.gath$common.pathways)
exp.mat.corr.gath$pathway.shared<-!is.na(exp.mat.corr.gath$common.pathways)
tra.mat.corr.gath$pathway.shared<-!is.na(tra.mat.corr.gath$common.pathways)

# Compute Specificity 
roc.ab<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=ab.mat.corr.gath$pathway.shared,r.vec=ab.mat.corr.gath$r,r.cutoff=x)})),data.type="Abundance")
roc.exp<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=exp.mat.corr.gath$pathway.shared,r.vec=exp.mat.corr.gath$r,r.cutoff=x)})),data.type="Expression")
roc.tra<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=tra.mat.corr.gath$pathway.shared,r.vec=tra.mat.corr.gath$r,r.cutoff=x)})),data.type="Transcription")
toplot.pathway<-rbind(roc.ab,roc.exp,roc.tra)
toplot.pathway$kegg.level<-"Pathway"

##################
## reaction level ##
##################

# Add reaction to the map
nog.ko.map.reaction<-nog.ko.map %>%
  left_join(ko.reaction.link,by=c("KO"="V2")) %>%
  dplyr::select(NOG,KO,reaction=V1,ko.description,reaction.description) %>%
  filter(!is.na(reaction))

# Filter by abundance-based occurrence
ab.occ<-apply(metaG.norm.match,2,occ)
ab.mat<-metaG.norm.match.log2[,which(ab.occ>min.occ)]
exp.mat<-ratio.mat[,which(ab.occ>min.occ)]
tra.mat<-metaT.norm.match.log2[,which(ab.occ>min.occ)]

# Get only the NOGs associated to KOs
ab.mat<-ab.mat[,which(colnames(ab.mat) %in% nog.ko.map.reaction$NOG)]
exp.mat<-exp.mat[,which(colnames(exp.mat) %in% nog.ko.map.reaction$NOG)]
tra.mat<-tra.mat[,which(colnames(tra.mat) %in% nog.ko.map.reaction$NOG)]

# Compute the correlation
ab.mat.corr<-cor(ab.mat)
exp.mat.corr<-cor(exp.mat)
tra.mat.corr<-cor(tra.mat)
diag(ab.mat.corr)<-NA
diag(exp.mat.corr)<-NA
diag(tra.mat.corr)<-NA

# Gather correlation matrix
ab.mat.corr.gath<-ab.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)
exp.mat.corr.gath<-exp.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)
tra.mat.corr.gath<-tra.mat.corr %>%
  as.data.frame() %>%
  rownames_to_column(var = "NOG1") %>%
  gather(key = "NOG2",value="r",-NOG1) %>%
  filter(!is.na(r)) %>%
  filter(r>0)

# Get only the best correlated NOG for each NOG
ab.mat.corr.gath<-ab.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()
exp.mat.corr.gath<-exp.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()
tra.mat.corr.gath<-tra.mat.corr.gath %>%
  group_by(NOG1) %>%
  filter(r==max(r)) %>%
  as.data.frame()

# Add a variable with common reactions (if exist)
fun<-function(x){
  mod1<-nog.ko.map.reaction$reaction[which(nog.ko.map.reaction$NOG %in% x[1])]
  mod2<-nog.ko.map.reaction$reaction[which(nog.ko.map.reaction$NOG %in% x[2])]
  if (any(mod1 %in% mod2)) paste(mod1[which(mod1 %in% mod2)],collapse=",") else NA
}
ab.mat.corr.gath$common.reactions<-apply(ab.mat.corr.gath,1,fun)
exp.mat.corr.gath$common.reactions<-apply(exp.mat.corr.gath,1,fun)
tra.mat.corr.gath$common.reactions<-apply(tra.mat.corr.gath,1,fun)

# Add a variable dichotomous variable (reaction.shared / not shared)
ab.mat.corr.gath$reaction.shared<-!is.na(ab.mat.corr.gath$common.reactions)
exp.mat.corr.gath$reaction.shared<-!is.na(exp.mat.corr.gath$common.reactions)
tra.mat.corr.gath$reaction.shared<-!is.na(tra.mat.corr.gath$common.reactions)

# Compute Specificity 
roc.ab<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=ab.mat.corr.gath$reaction.shared,r.vec=ab.mat.corr.gath$r,r.cutoff=x)})),data.type="Abundance")
roc.exp<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=exp.mat.corr.gath$reaction.shared,r.vec=exp.mat.corr.gath$r,r.cutoff=x)})),data.type="Expression")
roc.tra<-data.frame(do.call(rbind,lapply(seq(0,1,by=0.01),function(x){get.roc(true.vec=tra.mat.corr.gath$reaction.shared,r.vec=tra.mat.corr.gath$r,r.cutoff=x)})),data.type="Transcription")
toplot.reaction<-rbind(roc.ab,roc.exp,roc.tra)
toplot.reaction$kegg.level<-"Reaction"

############################
## Merge results and plot ##
############################
toplot<-rbind(toplot.module,toplot.pathway,toplot.reaction) %>% mutate(fpr=1-Specificity)
toplot$kegg.level<-factor(toplot$kegg.level,levels=c("Pathway","Module","Reaction"))

toplot %>% filter(fpr<0.05 & data.type=="Expression") %>% arrange(kegg.level) 

g1<-ggplot(data=toplot %>% filter(r.cutoff<1),aes(x=r.cutoff,y=fpr,col=data.type,shape=kegg.level,linetype=kegg.level)) +
  geom_line() +
  #facet_grid(.~kegg.level) +
  theme_bw() +
  ylab("False positive rate") +
  scale_x_continuous(minor_breaks = seq(0,1,0.05),limits = c(0,1)) +
  scale_y_continuous(minor_breaks = seq(0,1,0.05)) +
  scale_color_manual(values=c("#ED9406","#BD68B1","#33CCFF")) +
  xlab("Pearson r cutoff (rmin)")

g2<-ggplot(data=toplot %>% filter(r.cutoff<1),aes(x=r.cutoff,y=Sensitivity,col=data.type,shape=kegg.level,linetype=kegg.level)) +
  #geom_point() +
  geom_line() +
  #facet_grid(.~kegg.level) +
  theme_bw() +
  xlab("Pearson r cutoff (rmin)") +
  scale_x_continuous(minor_breaks = seq(0,1,0.05),limits = c(0,1)) +
  scale_y_continuous(minor_breaks = seq(0,1,0.05)) +
  scale_color_manual(values=c("#ED9406","#BD68B1","#33CCFF")) +
  ylab("Sensitivity")

g3<-ggplot(data=toplot %>% filter(r.cutoff<1),aes(x=1-Specificity,y=Sensitivity,col=data.type,shape=kegg.level,linetype=kegg.level,text=r.cutoff)) +
  geom_abline(linetype=2) +
  #geom_point() +
  geom_line() +
  #facet_grid(.~kegg.level) +
  theme_bw() +
  xlab("False positive rate") +
  scale_x_continuous(minor_breaks = seq(0,1,0.05),limits = c(0,1)) +
  scale_y_continuous(minor_breaks = seq(0,1,0.05)) +
  scale_color_manual(values=c("#ED9406","#BD68B1","#33CCFF")) +
  ylab("Sensitivity")

g3 | g1 | g2

ggsave("../results/Coexpression_ROC.pdf",width=15,height=4)
