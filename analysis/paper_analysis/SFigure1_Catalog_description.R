library(data.table)
library(EcolUtils)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("../data/lib/varpart.sqr.euc_functions.R")
source("../data/lib/sushipal.R")
pal<-sushi.palette()[c(14,1:13,15)]

func.stats<-fread("../data/processed/gene.cat.stats/func.annotation.stats",header=T,sep="\t",data.table = F)

func.stats.ko<-func.stats[c(1,2),]
func.stats.ko[1,2]<-func.stats.ko[1,2]-func.stats.ko[2,2]
func.stats.ko[1,1]<-"unannotated"
func.stats.ko$type<-"KEGG"
func.stats.ko<-func.stats.ko %>% mutate(perc=100*n/sum(n))

func.stats.egg<-func.stats[c(1,3,4),]
func.stats.egg[1,2]<-func.stats.egg[1,2]-sum(func.stats.egg[c(2,3),2])
func.stats.egg[1,1]<-"unannotated"
func.stats.egg$type<-"eggNOG"
func.stats.egg<-func.stats.egg %>% mutate(perc=100*n/sum(n))

func.stats.merged<-rbind(func.stats.ko,func.stats.egg)
func.stats.merged<-func.stats.merged %>% mutate(annotation=paste(annotation," (",round(perc,2)," %)",sep=""))
func.stats.merged$annotation<-fct_relevel(func.stats.merged$annotation,"unannotated (76.41 %)","unannotated (17.26 %)")

p0<-ggplot(func.stats.merged,aes(x=type,y=n,fill=annotation)) +
  geom_bar(stat = "identity",position = position_stack()) +
  theme_bw() +
  ylim(0,50000000) +
  scale_fill_manual(values=pal[c(1,1,2,3,4)])
#ggsave("../results/gene.cat.stats.func.pdf",width=10,height=2)


tax.stats.dom<-fread("../data/processed/gene.cat.stats/tax.annotation.Domain.stats",header=T,sep="\t",data.table = F)
tax.stats.dom$type<-"Domain"
tax.stats.dom$Domain[which(tax.stats.dom$Domain=="")]<-"unannotated"
tax.stats.dom<-tax.stats.dom %>% mutate(perc=100*n/sum(n))
tax.stats.dom<-tax.stats.dom %>% mutate(Domain=paste(Domain," (",round(perc,2)," %)",sep=""))
tax.stats.dom$Domain<-fct_relevel(tax.stats.dom$Domain,"No annotation (27.48 %)","Bacteria (59.26 %)","Viruses (5.26 %)","LUCA (3.36 %)","Eukaryota (2.61 %)","Archaea (2.03 %)")

p1<-ggplot(tax.stats.dom,aes(x=type,y=n,fill=Domain)) +
  geom_bar(stat = "identity",position = position_stack()) +
  theme_bw() +
  ylim(0,50000000) +
  scale_fill_manual(values=pal)


tax.stats.phyl<-fread("../data/processed/gene.cat.stats/tax.annotation.Phylum.stats",header=T,sep="\t",data.table = F)
tax.stats.class<-fread("../data/processed/gene.cat.stats/tax.annotation.Class.stats",header=T,sep="\t",data.table = F)
colnames(tax.stats.class)[1]<-"Phylum"
tax.stats.phyl<-tax.stats.phyl[-c(which(tax.stats.phyl$Phylum=="Proteobacteria")),]
tax.stats.phyl<-rbind(tax.stats.phyl,tax.stats.class[grep("proteobacteria",tax.stats.class$Phylum),])
tax.stats.phyl$n[which(tax.stats.phyl$Phylum=="")]<-tax.stats.phyl$n[which(tax.stats.phyl$Phylum=="")]+(sum(tax.stats.class$n)-sum(tax.stats.phyl$n))
tax.stats.phyl$type<-"Phylum"
tax.stats.phyl$Phylum[which(tax.stats.phyl$Phylum=="")]<-"unannotated"



tax.stats.phyl$Phylum[tax.stats.phyl$n<200000]<-"Other"
tax.stats.phyl<-tax.stats.phyl %>%
  group_by(Phylum,type) %>%
  summarise(n=sum(n))
tax.stats.phyl$Phylum<-fct_reorder(tax.stats.phyl$Phylum,tax.stats.phyl$n)

tax.stats.phyl<-tax.stats.phyl %>% as.data.frame() %>% mutate(perc=100*n/sum(n))
tax.stats.phyl<-tax.stats.phyl %>% as.data.frame() %>% mutate(Phylum=paste(Phylum," (",round(perc,2)," %)",sep=""))
tax.stats.phyl$Phylum<-fct_reorder(tax.stats.phyl$Phylum,desc(tax.stats.phyl$n))
tax.stats.phyl$Phylum<-fct_relevel(tax.stats.phyl$Phylum,"unannotated (52.06 %)")

p2<-ggplot(tax.stats.phyl,aes(x=type,y=n,fill=Phylum)) +
  geom_bar(stat = "identity",position = position_stack()) +
  theme_bw() +
  ylim(0,50000000) +
  scale_fill_manual(values=pal)

## Gene accumulation ##
dades<-fread("../data/processed/gene.cat.stats/accum.genes.tsv",sep="\t",header=T,data.table = F)
dades$sample.num<-1:nrow(dades)
env.mat<-fread("zcat < ../data/processed/KO_env.mat.metaG.txt.gz",sep="\t", header=T)
dades<-cbind(dades,env.mat)
dades$polar[abs(dades$Latitude)<60]<-"non polar"
dades$polar[abs(dades$Latitude)>=60]<-"polar"

p3<-ggplot(dades,aes(x=sample.num,y=n.genes,col=polar)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = c(0.9,0.3)) +
  scale_color_manual(values=c("#FB8072FF","#80b1d3")) +
  xlab("Sample") +
  ylab("Number of genes") +
  ylim(0,50000000) +
  geom_vline(xintercept = 139.5,linetype=2) +
  geom_hline(yintercept = 47000000,linetype=2)

(p3 | p1 | p2 | p0) + plot_layout(ncol=4,widths = c(3,1,1,1))
#ggsave("../results/Fig1B.pdf",width=20,height=8)
