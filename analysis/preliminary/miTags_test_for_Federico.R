library(DESeq2)
library(data.table)
library(EcolUtils)
library(ggplot2)
library(tidyverse)
library(vegan)

# Load tax
miTags<-fread("../data/processed/miTags/OTU.tab.rr.TARA180.noeuks.txt",sep="\t",header=T,data.table = F)
rownames(miTags)<-miTags$V1
miTags<-miTags[,-1]

# Load environmental data
env.mat<-fread("zcat < ../data/processed/NOG_env.mat.metaG.txt.gz",sep="\t",header=T,data.table = F,stringsAsFactors = T)
env.mat$Layer<-factor(env.mat$Layer,levels = c("SRF","DCM","MES","other"))
miTags<-miTags[match(gsub("-",".",env.mat$Sample_name),rownames(miTags)),]

# Compare to de-novo clustering
#mat.16S<-read.table("/Users/guillemsalazar/polybox/ETH/PROJECTS/extras/metaphlan_mtags_mOTU_comparison/v3/data/otutab.all.final.txt",sep="\t",header=T,comment.char="",quote="",row.names=1)
#colnames(mat.16S)<-sub("_G","",colnames(mat.16S))
#to.keep<-rownames(mat.16S)
#mat.16S<-t(mat.16S)
#miTags.matched<-miTags[match(rownames(mat.16S),rownames(miTags)),]
#env.mat.matched<-env.mat[match(rownames(mat.16S),gsub("-",".",env.mat$Sample_name)),]
#miTags.matched<-miTags.matched/rowSums(miTags.matched)
#mat.16S<-mat.16S/rowSums(mat.16S)

#toplot.df<-data.frame(shan.mitags=diversity(miTags.matched),shan.denovo=diversity(mat.16S),env.mat.matched)

#ggplot(data=toplot.df %>% filter(!is.na(Layer)),aes(x=Latitude,y=shan.denovo,col=Layer)) +
#  geom_point() +
#  geom_smooth(method="loess") +
#  facet_grid(.~Layer)

# Remove chloroplasts
to.remove<-which(grepl(";Chloroplast",colnames(miTags)))
miTags<-miTags[,-to.remove]

# Normailze (relative abundance)
miTags<-miTags/rowSums(miTags)

# Split into Archaea (A), Bacteria phototrophs (P) and Bacteria heterotrophs (H)
arch<-which(grepl("Archaea;",colnames(miTags)))
miTags.A<-miTags[,arch]
miTags.B<-miTags[,-arch]
phot<-which(grepl(";Cyanobacteria",colnames(miTags.B)))
phot<-c(phot,which(grepl(";Chloroflexi",colnames(miTags.B))))
miTags.P<-miTags.B[,phot]
miTags.H<-miTags.B[,-phot]

# Compute the sum
sum.A<-rowSums(miTags.A)
sum.P<-rowSums(miTags.P)
sum.H<-rowSums(miTags.H)

# Plot
env.mat$Sample_name<-gsub("-",".",env.mat$Sample_name)
toplot.df<-data.frame(Sample_name=names(c(sum.A,sum.P,sum.H)),ra=c(sum.A,sum.P,sum.H),troph.g=rep(c("A","P","H"),each=nrow(miTags.A)))
toplot.df<-toplot.df %>% left_join(env.mat,by="Sample_name")
toplot.df$Latitude.bin<-cut(toplot.df$Latitude,breaks = c(-90,-70,-40,-20,0,20,40,60,80,90))

toplot.df<-toplot.df %>%
  select(Latitude.bin,ra,troph.g) %>%
  group_by(Latitude.bin,troph.g) %>%
  summarise(ra=sum(ra)) %>%
  mutate(freq = ra / sum(ra)) %>%
  as.data.frame()

ggplot(data=toplot.df,aes(x=Latitude.bin,y=freq*100,fill=troph.g)) +
  geom_bar(stat="identity") +
  coord_flip() +
  ylab("Relative abundance") +
  theme_bw()
ggsave(filename = "/Users/guillemsalazar/Desktop/Fig1B_for_Federico.pdf",width = 4,height = 6)
  

# Compute the richness and diversity
rich.A<-specnumber(miTags.A)
rich.P<-specnumber(miTags.P)
rich.H<-specnumber(miTags.H)

shan.A<-diversity(miTags.A)
shan.P<-diversity(miTags.P)
shan.H<-diversity(miTags.H)


toplot.df<-data.frame(Sample_name=names(c(shan.A,shan.P,shan.H)),D=c(shan.A,shan.P,shan.H),S=c(rich.A,rich.P,rich.H),troph.g=rep(c("A","P","H"),each=nrow(miTags.A)),env.mat)

g1<-ggplot(data=toplot.df %>% filter(Layer %in% c("SRF","DCM")),aes(x=Latitude,y=D,col=troph.g)) +
  geom_point() +
  geom_smooth(se=F) +
  theme_bw() +
  facet_grid(.~Layer) +
  labs(title = "A. Diversity") +
  ylab("Diversity (Shannon)")
  
g2<-ggplot(data=toplot.df %>% filter(Layer %in% c("SRF","DCM")),aes(x=Latitude,y=S,col=troph.g)) +
  geom_point() +
  geom_smooth(se=F) +
  theme_bw() +
  facet_grid(.~Layer) +
  scale_y_log10() +
  labs(title = "B. Richness") +
  ylab("Richness (observed OTUs)")

g3<-ggplot(data=toplot.df %>% filter(Layer %in% c("SRF","DCM")),aes(x=S,y=D,col=troph.g)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  theme_bw() +
  facet_grid(.~Layer) +
  scale_x_log10() +
  labs(title = "C. Richness vs. Diversity") +
  xlab("Richness (observed OTUs)") +
  ylab("Diversity (Shannon)")

library(patchwork)
g1 / g2 / g3 + plot_layout()
ggsave(filename = "/Users/guillemsalazar/Desktop/Richness_Diversity.pdf",width = 5,height = 9)
