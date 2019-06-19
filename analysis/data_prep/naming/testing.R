library(data.table)
library(tidyverse)
#mapping.correct.metaG<-mapping.correct[grepl("_G",mapping.correct$old_name),]
#mapping.correct.metaT<-mapping.correct[grepl("_T",mapping.correct$old_name),]

dat<-fread("/Users/guillemsalazar/Desktop/test.tsv",sep="\t",header=F)
a<-dat$V1
length(a)

# metaG
pos<-match(a,mapping.correct.metaG$PANGAEA_ID)
length(pos)
length(na.exclude(pos))
a[which(is.na(pos))]
correct.incorrect.map[match(a[which(is.na(pos))],correct.incorrect.map$PANGAEA_ID_incorrect),]
mapping.correct.metaG[pos,]
unique(sapply(strsplit(mapping.correct.metaG$old_name[pos],"_"),"[[",5))

# metaT
pos<-sapply(a,function(x){grep(x,mapping.correct.metaT$PANGAEA_ID)})
length(pos)
length(pos[pos>0])
mapping.correct.metaT[pos,]
unique(sapply(strsplit(mapping.correct.metaT$old_name[pos],"_"),"[[",5))
