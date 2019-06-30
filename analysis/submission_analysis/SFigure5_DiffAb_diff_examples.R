#library(propagate)
library(vegan)
library(tidyverse)
library(ggrepel)
library(data.table)
library(patchwork)
#source("http://sunagawalab.ethz.ch/sushilab/resources/R/sushi.palette.R")
source("../data/lib/sushi.palette.R")
source("../data/lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])
library(geosphere)
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

####################



####################
ab.gath<-metaG.norm.match.log2 %>%
  rownames_to_column(var="Barcode") %>%
  gather(key = "KO",value="value",-Barcode) %>%
  mutate(data_type="Abundance")
exp.gath<-ratio.mat %>%
  rownames_to_column(var="Barcode") %>%
  gather(key = "KO",value="value",-Barcode) %>%
  mutate(data_type="Expression")
trans.gath<-metaT.norm.match.log2 %>%
  rownames_to_column(var="Barcode") %>%
  gather(key = "KO",value="value",-Barcode) %>%
  mutate(data_type="Transcription")
all.gath<-bind_rows(ab.gath,exp.gath,trans.gath) %>%
  left_join(env.mat.match,by="Barcode") %>%
  filter(KO %in% sel.genes$KO) %>%
  left_join(sel.genes,by="KO")
all.gath$date<-as.POSIXct(substr(as.character(all.gath$Event.date),1,10))
all.gath$doy<-yday(all.gath$date)
all.gath$daylength<-as.numeric(daylength(lat = all.gath$Latitude,doy=all.gath$doy))


# Rubisco
g1<-ggplot(all.gath %>% filter(KO %in% c("K01601","K01602")),aes(x=-Depth.nominal,y=value)) +
  geom_point(col="gray60") +
  geom_smooth(method="loess",se=F,col="black") +
  facet_grid(str_wrap(Enzyme,10)~data_type,scales = "free_x") +
  coord_flip() +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  ylab("log2-transformed value") +
  xlab("Depth (m)")

test1<-cor.test(metaT.norm.match.log2$K02703,metaT.norm.match.log2$K01601)
g2<-ggplot(cbind(metaT.norm.match.log2[match(env.mat.match$Barcode,rownames(metaT.norm.match.log2)),],env.mat.match),aes(x=K02703,y=K01601)) +
  geom_point(aes(shape=epi),col="gray60") +
  geom_smooth(method="lm",se=F,col="black") +
  theme_bw() +
  xlab(str_wrap(paste(sel.genes$Enzyme[grep("K02703",sel.genes$Enzyme)],"transcript abundance"),40)) +
  ylab(str_wrap(paste(sel.genes$Enzyme[grep("K01601",sel.genes$Enzyme)],"transcript abundance"),40)) +
  geom_label(aes(x=-5,y=7,label=paste("r =",round(test1$estimate,2),"P < 0.001")),hjust=0)

test2<-cor.test(metaT.norm.match.log2$K02703,metaT.norm.match.log2$K01602)
g3<-ggplot(cbind(metaT.norm.match.log2[match(env.mat.match$Barcode,rownames(metaT.norm.match.log2)),],env.mat.match),aes(x=K02703,y=K01602)) +
  geom_point(aes(shape=epi),col="gray60") +
  geom_smooth(method="lm",se=F,col="black") +
  theme_bw() +
  xlab(str_wrap(paste(sel.genes$Enzyme[grep("K02703",sel.genes$Enzyme)],"transcript abundance"),40)) +
  ylab(str_wrap(paste(sel.genes$Enzyme[grep("K01602",sel.genes$Enzyme)],"transcript abundance"),40)) +
  geom_label(aes(x=-5,y=6,label=paste("r =",round(test2$estimate,2),"P < 0.001")),hjust=0)

g1 | (g2 / g3)
ggsave("../results/SFig6_Rubisco.pdf")

####################
ggplot(all.gath %>% filter(Metabolism=="PHOTOSYNTHESIS" & epi=="EPI" & Layer!="MIX"),aes(x=Temperature,y=value)) +
  geom_point() +
  geom_smooth(method="loess",span=1,se=F) +
  facet_grid(data_type~str_wrap(Enzyme,30),scales = "free_y") +
  theme_bw() +
  theme(strip.text.x = element_text(angle=90))

ggplot(cbind(env.mat.match,metaT.norm.match.log2) %>% filter(Layer !="MIX"),aes(x=Depth.nominal,y=K01601,col=Layer)) + geom_point() + geom_smooth()

####################
metaG.norm.match.log2<-fread("../data/processed/NOG_metaG.norm.match.log2.txt",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("../data/processed/NOG_metaG.norm.match.txt",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("../data/processed/NOG_metaT.norm.match.log2.txt",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("../data/processed/NOG_metaT.norm.match.txt",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("../data/processed/NOG_ratio.mat.txt",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]
env.mat.match<-fread("../data/processed/NOG_env.mat.match.txt",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match$epi<-env.mat.match$Layer
levels(env.mat.match$epi)<-c("EPI","MES","EPI","EPI")

ab.gath.nog<-metaG.norm.match.log2 %>%
  rownames_to_column(var="Sample_name") %>%
  select(Sample_name,'0XS7T') %>%
  gather(key = "KO",value="value",-Sample_name) %>%
  mutate(data_type="Abundance")
exp.gath.nog<-ratio.mat %>%
  rownames_to_column(var="Sample_name") %>%
  select(Sample_name,'0XS7T') %>%
  gather(key = "KO",value="value",-Sample_name) %>%
  mutate(data_type="Expression")
trans.gath.nog<-metaT.norm.match.log2 %>%
  rownames_to_column(var="Sample_name") %>%
  select(Sample_name,'0XS7T') %>%
  gather(key = "KO",value="value",-Sample_name) %>%
  mutate(data_type="Transcription")
all.gath.nog<-bind_rows(ab.gath.nog,exp.gath.nog,trans.gath.nog) %>%
  filter(data_type=="Expression") %>%
  select(Sample_name,dmsp.synt=value) %>%
  unique()

tmp.df<-left_join(all.gath,all.gath.nog,by="Sample_name")


ggplot(tmp.df %>% filter(KO %in% c("K00956","K00957","K00390","K00381","K00380") & epi=="EPI" & data_type=="Expression"),aes(x=dmsp.synt,y=value)) +
  geom_point(col="gray60") +
  geom_smooth(method="lm",se=F,col="black") +
  facet_grid(str_wrap(Enzyme,10)~.,scales = "free_y") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0))

cor.test(tmp.df$dmsp.synt[tmp.df$KO=="K00956" & tmp.df$data_type=="Expression"],tmp.df$value[tmp.df$KO=="K00956" & tmp.df$data_type=="Expression"])
cor.test(tmp.df$dmsp.synt[tmp.df$KO=="K00957" & tmp.df$data_type=="Expression"],tmp.df$value[tmp.df$KO=="K00957" & tmp.df$data_type=="Expression"])
cor.test(tmp.df$dmsp.synt[tmp.df$KO=="K00390" & tmp.df$data_type=="Expression"],tmp.df$value[tmp.df$KO=="K00390" & tmp.df$data_type=="Expression"])
cor.test(tmp.df$dmsp.synt[tmp.df$KO=="K00381" & tmp.df$data_type=="Expression"],tmp.df$value[tmp.df$KO=="K00381" & tmp.df$data_type=="Expression"])
cor.test(tmp.df$dmsp.synt[tmp.df$KO=="K00380" & tmp.df$data_type=="Expression"],tmp.df$value[tmp.df$KO=="K00380" & tmp.df$data_type=="Expression"])

# Load OTUtab
otu.table<-fread("../data/processed/miTags/Merged_table_domain.txt",sep="\t",header=T,data.table=F)
rownames(otu.table)<-otu.table[,1]
otu.table<-otu.table[,-c(1)]
otu.table<-as.data.frame(t(otu.table[,grepl("TARA",colnames(otu.table))]))
otu.table$Sample<-as.character(sapply(strsplit(rownames(otu.table),"_"),function(x){paste(x[1:4],collapse="_")}))

otu.table<-otu.table %>%
  gather(key="Domain",value="Abundance",-Sample) %>%
  group_by(Sample,Domain) %>%
  summarise(Abundance=sum(Abundance)) %>%
  group_by(Sample) %>%
  mutate(RAbundance=Abundance/sum(Abundance)) %>%
  filter(Domain=="Eukaryota;") %>%
  select(Sample_name=Sample,Euk=RAbundance)

test.df<-all.gath %>%
  left_join(otu.table,by="Sample_name")

g1<-ggplot(test.df %>% filter(KO %in% c("K00956","K00957","K00390","K00381","K00380","K12339") & epi=="EPI" & data_type=="Expression"),aes(x=Euk,y=value)) +
  geom_point(col="gray60") +
  geom_smooth(method="lm",se=F,col="black") +
  facet_grid(str_wrap(Enzyme,10)~.,scales = "free_x") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0))

g2<-ggplot(test.df %>% filter(KO %in% c("K00956","K00957","K00390","K00381","K00380","K12339") & epi=="EPI" & data_type=="Expression"),aes(x=ChlorophyllA,y=value)) +
  geom_point(col="gray60") +
  geom_smooth(method="lm",se=F,col="black") +
  facet_grid(str_wrap(Enzyme,10)~.,scales = "free_x") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0))

g1 | g2

tmp.df<-test.df %>%
  filter(KO %in% c("K00956","K00957","K00390","K00381","K00380","K12339") & epi=="EPI" & data_type=="Expression") %>%
  select(Sample_name,Euk,ChlorophyllA,KO,value) %>%
  spread(key = KO,value = value)

cor.test(tmp.df$Euk,tmp.df[,4])
cor.test(tmp.df$Euk,tmp.df[,5])
cor.test(tmp.df$Euk,tmp.df[,6])
cor.test(tmp.df$Euk,tmp.df[,7])
cor.test(tmp.df$Euk,tmp.df[,8])
cor.test(tmp.df$Euk,tmp.df[,9])

cor.test(tmp.df$ChlorophyllA,tmp.df[,4])
cor.test(tmp.df$ChlorophyllA,tmp.df[,5])
cor.test(tmp.df$ChlorophyllA,tmp.df[,6])
cor.test(tmp.df$ChlorophyllA,tmp.df[,7])
cor.test(tmp.df$ChlorophyllA,tmp.df[,8])
cor.test(tmp.df$ChlorophyllA,tmp.df[,9])



# Nitrogen fixation
#vars<-c("NO2","PO4","NO2NO3","Temperature","Oxygen","NO3","Iron.5m","Ammonium.5m","Nitracline","ChlorophyllA","HCO3")
vars<-c("HCO3","CO3","Carbon.total","Alkalinity.total","NO2","PO4","NO2NO3","Si","Temperature","Salinity","Density","Oxygen","NO3","ChlorophyllA","Fluorescence","PAR.TO","Iron.5m","Ammonium.5m","Gradient.Surface.temp(SST)","Okubo.Weiss","Lyapunov","Residence.time","Depth.Mixed.Layer","Brunt.Väisälä","Depth.Max.O2","Depth.Min.O2","Nitracline")
tmp.df<-env.mat.match[,vars]
tmp.df<-decostand(tmp.df,method = "range",MARGIN = 2,na.rm = T)
tmp.df$K05376<-ratio.mat$K05376
tmp.df<-apply(tmp.df,2,function(x){cor(log10(x+0.00001),tmp.df$K05376,use = "complete.obs")})
sort(tmp.df)
tmp.df<-na.exclude(tmp.df)
library(MASS)
fit <- lm(K02586~.,data=tmp.df)
step <- stepAIC(fit, direction="backward",na.action="na.omit")
summary(step)
#step$anova # display results

plot(tmp.df$K02586,predict(step))
abline(0,1)


library(relaimpo)
calc.relimp(step,type=c("lmg","last","first","pratt"),rela=TRUE)

ggplot(all.gath %>% filter(KO %in% c("K02586","K02588") & epi=="EPI"),aes(x=NO2NO3,y=value,col=polar)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess",se=F) +
  facet_grid(data_type~str_wrap(Enzyme,10),scales = "free_y") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  xlab("Nitrate and nitrite concentration [uM / L]") +
  ylab("Expression") +
  scale_x_sqrt()

ggplot(all.gath %>% filter(KO %in% c("K02586","K02588") & epi=="EPI" & data_type=="Expression"),aes(x=Oxygen,y=NO2NO3,col=value)) +
  geom_point(size=2) +
  #geom_smooth(method="loess",se=F) +
  facet_grid(.~str_wrap(Enzyme,10)) +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  xlab("Nitrate and nitrite concentration [uM / L]") +
  ylab("Expression") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_color_continuous(low="blue",high = "red")



g1<-ggplot(cbind(env.mat.match,ratio.mat) %>% filter(epi=="EPI"),aes(x=NO2NO3,y=K02586,col=polar)) +
  geom_point() +
  geom_smooth(method="lm",formula=y~log10(x+0.0001),se=F) +
  scale_x_sqrt(breaks=c(0,0.1,1,2,5,10,20,30,40)) +
  theme_bw() +
  scale_color_manual(values=c("#FB8072FF","#8DD3C7FF")) +
  xlab("Nitrate and nitrite concentration [uM / L]") +
  ylab("nifD expression")
g2<-ggplot(cbind(env.mat.match,ratio.mat) %>% filter(epi=="EPI"),aes(x=NO2NO3,y=K02588,col=polar)) +
  geom_point() +
  geom_smooth(method="lm",formula=y~log10(x+0.0001),se=F) +
  scale_x_sqrt(breaks=c(0,0.1,1,2,5,10,20,30,40)) +
  theme_bw() +
  scale_color_manual(values=c("#FB8072FF","#8DD3C7FF")) +
  xlab("Nitrate and nitrite concentration [uM / L]") +
  ylab("nifH expression")
g1 | g2
ggsave("../results/nifDH_expression.pdf",width=9,height=4)

tmp<-cbind(env.mat.match[,c("NO2NO3","Layer","Oxygen","epi","polar")],ratio.mat[,c("K02586","K02588")]) %>% na.exclude()
m1<-lm(tmp$K02586~tmp$polar+log10(tmp$NO2NO3+0.0001))
summary(m1)
m2<-lm(tmp$K02588~tmp$polar+log10(tmp$NO2NO3+0.0001))
summary(m2)


# Maen 
ggplot(cbind(env.mat.match,metaG.norm.match.log2),aes(x=-Depth.nominal,y=K00376,size=Oxygen)) +
  geom_point() +
  geom_smooth(method="loess",se=F) +
  theme_bw() +
  scale_color_manual(values=c("#FB8072FF","#8DD3C7FF")) +
  coord_flip()

# Amonia monoxygenase
ggplot(cbind(env.mat.match,metaG.norm.match.log2),aes(x=-Depth.nominal,y=K10944)) +
  geom_point() +
  geom_smooth(method="loess",se=F) +
  theme_bw() +
  scale_color_manual(values=c("#FB8072FF","#8DD3C7FF")) +
  coord_flip()

cbind(env.mat.match,metaG.norm.match) %>%
  group_by(epi) %>%
  summarise(mean(K01601),mean(K01602),sd(K01601),sd(K01602))

cbind(env.mat.match,metaT.norm.match) %>%
  group_by(epi) %>%
  summarise(mean(K01601),mean(K01602))

cbind(env.mat.match,ratio.mat) %>%
  group_by(epi) %>%
  summarise(mean(K01601),mean(K01602))
