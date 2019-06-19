library(tidyverse)
library(data.table)
library(ggrepel)

assoc<-fread("../results/Coexpression_pairs_renamed.tsv",sep="\t",header=T,data.table = F)
OGcat<-fread("../data/lib/NOG.annotations.tsv",sep="\t",header=F,data.table = F)
OGcat$V2<-gsub("ENOG41","",OGcat$V2)
tmp<-fread("../data/lib/NOG.categories.tsv",sep="\t",header=F,data.table = F)
OGcat<-OGcat %>% left_join(tmp,by=c("V5"="V1"))
assoc<-assoc %>% left_join(OGcat,by=c("Representative 2"="V2.x"))

assoc.tosave<-assoc %>% select("Representative 1","Representative 2","Co-expression Pearson r","Representative 2 type","Representative 2 annotation",`Functional category`=V2.y)

assoc.og<-assoc %>% filter(`Representative 2 type`=="OG")
fwrite(assoc.tosave,file = "../results/Coexpression_pairs_renamed_withFunctionalCategory.tsv",sep="\t",quote = F,row.names = F)

exp.ma<-fread("zcat < ../data/processed/NOG_ratio.mat.txt.gz",sep="\t",header=T,data.table = F)
all.ogs<-data.frame(og=colnames(exp.ma)[-1])
all.ogs<-all.ogs %>% left_join(OGcat,by=c("og"="V2.x"))

all.tab<-table(all.ogs$V2.y) %>% as.data.frame()
assoc.tab<-table(assoc.og$V2.y) %>% as.data.frame()
colnames(assoc.tab)<-c("Var1","Freq.assoc")

merged<-left_join(all.tab,assoc.tab) %>% mutate_all(funs(replace(., is.na(.), 0)))

ggplot(data=merged,aes(x=Freq,y=Freq.assoc)) +
  geom_point(size=5,alpha=0.5) +
  geom_smooth(method="lm",se=F,col="black",alpha=0.5) +
  geom_label_repel(data=merged %>% filter(Freq.assoc>0),aes(x=Freq,y=Freq.assoc,label=str_wrap(Var1,30)),force=5,size=2.5,fontface = 'bold',point.padding = 0.5,box.padding = 0.9,fill="gray90",alpha=0.8) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("Total number of OGs") +
  ylab("Number of OGs involved in a GC-OG association")
ggsave(filename = "../results/GC_OG_assoc_prop.pdf",width = unit(8,"cm"),height = unit(8,"cm"))

cont.tab<-merged
rownames(cont.tab)<-cont.tab$Var1
cont.tab<-as.table(as.matrix(cont.tab[,-1]))

library(XNomial)
xmonte(merged$Freq.assoc,merged$Freq,detail = 3,histobins = T,statName = "Chisq")
