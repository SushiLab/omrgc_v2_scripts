library(data.table)
to.integer<-function(taula){taula<-round(taula/(max(taula)/1e9))}
##############################################################
# Load the KO data ###########################################
##############################################################
metaG<-fread("../data/processed/KO_metaG.txt",header=T,sep="\t",data.table = F)
rownames(metaG)<-metaG[,1]
metaG<-metaG[,-1]

metaT<-fread("../data/processed/KO_metaT.txt",header=T,sep="\t",data.table = F)
rownames(metaT)<-metaT[,1]
metaT<-metaT[,-1]

metaG.match<-fread("../data/processed/KO_metaG.match.txt",header=T,sep="\t",data.table = F)
rownames(metaG.match)<-metaG.match[,1]
metaG.match<-metaG.match[,-1]

metaT.match<-fread("../data/processed/KO_metaT.match.txt",header=T,sep="\t",data.table = F)
rownames(metaT.match)<-metaT.match[,1]
metaT.match<-metaT.match[,-1]

#######################
# Round and transpose #
#######################
metaG<-t(to.integer(metaG))
metaT<-t(to.integer(metaT))
metaG.match<-t(to.integer(metaG.match))
metaT.match<-t(to.integer(metaT.match))
metaG<-metaG[which(rowSums(metaG)>0),]
metaT<-metaT[which(rowSums(metaT)>0),]
metaG.match<-metaG.match[which(rowSums(metaG.match)>0),]
metaT.match<-metaT.match[which(rowSums(metaT.match)>0),]

cat("Rarefaction depth for metaG.scaled and metaT.scaled: ",min(c(colSums(metaG),colSums(metaT),colSums(metaG.match),colSums(metaT.match))))

###############
## Save data ##
###############
write.table(metaG,file = "../data/processed/KO_metaG.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.match,file = "../data/processed/KO_metaG.match.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT,file = "../data/processed/KO_metaT.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.match,file = "../data/processed/KO_metaT.match.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)


###############################################################
# Load the NOG data ###########################################
###############################################################
metaG<-fread("../data/processed/NOG_metaG.txt",header=T,sep="\t",data.table = F)
rownames(metaG)<-metaG[,1]
metaG<-metaG[,-1]

metaT<-fread("../data/processed/NOG_metaT.txt",header=T,sep="\t",data.table = F)
rownames(metaT)<-metaT[,1]
metaT<-metaT[,-1]

metaG.match<-fread("../data/processed/NOG_metaG.match.txt",header=T,sep="\t",data.table = F)
rownames(metaG.match)<-metaG.match[,1]
metaG.match<-metaG.match[,-1]

metaT.match<-fread("../data/processed/NOG_metaT.match.txt",header=T,sep="\t",data.table = F)
rownames(metaT.match)<-metaT.match[,1]
metaT.match<-metaT.match[,-1]

#######################
# Round and transpose #
#######################
metaG<-t(to.integer(metaG))
metaT<-t(to.integer(metaT))
metaG.match<-t(to.integer(metaG.match))
metaT.match<-t(to.integer(metaT.match))
metaG<-metaG[which(rowSums(metaG)>0),]
metaT<-metaT[which(rowSums(metaT)>0),]
metaG.match<-metaG.match[which(rowSums(metaG.match)>0),]
metaT.match<-metaT.match[which(rowSums(metaT.match)>0),]

cat("Rarefaction depth for metaG.scaled and metaT.scaled: ",min(c(colSums(metaG),colSums(metaT),colSums(metaG.match),colSums(metaT.match))))

###############
## Save data ##
###############
write.table(metaG,file = "../data/processed/NOG_metaG.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.match,file = "../data/processed/NOG_metaG.match.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT,file = "../data/processed/NOG_metaT.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.match,file = "../data/processed/NOG_metaT.match.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)

###############################################################
# Load the NOGplusGF data ###########################################
###############################################################
metaG<-fread("../data/processed/NOGplusGF_metaG.txt",header=T,sep="\t",data.table = F)
rownames(metaG)<-metaG[,1]
metaG<-metaG[,-1]

metaT<-fread("../data/processed/NOGplusGF_metaT.txt",header=T,sep="\t",data.table = F)
rownames(metaT)<-metaT[,1]
metaT<-metaT[,-1]

metaG.match<-fread("../data/processed/NOGplusGF_metaG.match.txt",header=T,sep="\t",data.table = F)
rownames(metaG.match)<-metaG.match[,1]
metaG.match<-metaG.match[,-1]

metaT.match<-fread("../data/processed/NOGplusGF_metaT.match.txt",header=T,sep="\t",data.table = F)
rownames(metaT.match)<-metaT.match[,1]
metaT.match<-metaT.match[,-1]

#######################
# Round and transpose #
#######################
metaG<-t(to.integer(metaG))
metaT<-t(to.integer(metaT))
metaG.match<-t(to.integer(metaG.match))
metaT.match<-t(to.integer(metaT.match))
metaG<-metaG[which(rowSums(metaG)>0),]
metaT<-metaT[which(rowSums(metaT)>0),]
metaG.match<-metaG.match[which(rowSums(metaG.match)>0),]
metaT.match<-metaT.match[which(rowSums(metaT.match)>0),]

cat("Rarefaction depth for metaG.scaled and metaT.scaled: ",min(c(colSums(metaG),colSums(metaT),colSums(metaG.match),colSums(metaT.match))))

###############
## Save data ##
###############
write.table(metaG,file = "../data/processed/NOGplusGF_metaG.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaG.match,file = "../data/processed/NOGplusGF_metaG.match.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT,file = "../data/processed/NOGplusGF_metaT.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(metaT.match,file = "../data/processed/NOGplusGF_metaT.match.scaled.txt",sep="\t",row.names=T,col.names=NA,quote=F)
