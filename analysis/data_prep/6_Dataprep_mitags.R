library(data.table)
library(vegan)

# Load the environmental table and reorder the OTUtab
env.mat<-fread("../data/processed/KO_env.mat.metaG.txt",header=T,data.table = F)

##########
# Domain #
##########

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_domain.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))
colnames(OTU.tab)[colnames(OTU.tab)=="Unknown domain"]<-"unclassified"

# Get the total number of reads
tot<-apply(OTU.tab,1,sum)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_domain.tsv",sep="\t",quote=F,row.names = T)

##########
# Phylum #
##########

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_phylum.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))

# Get the total number of reads
unclassified<-tot-apply(OTU.tab,1,sum)
OTU.tab<-cbind(OTU.tab,unclassified)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_phylum.tsv",sep="\t",quote=F,row.names = T)

#########
# Class #
#########

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_class.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))

# Get the total number of reads
unclassified<-tot-apply(OTU.tab,1,sum)
OTU.tab<-cbind(OTU.tab,unclassified)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_class.tsv",sep="\t",quote=F,row.names = T)

#########
# Order #
#########

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_order.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))

# Get the total number of reads
unclassified<-tot-apply(OTU.tab,1,sum)
OTU.tab<-cbind(OTU.tab,unclassified)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_order.tsv",sep="\t",quote=F,row.names = T)


##########
# Family #
##########

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_family.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))

# Get the total number of reads
unclassified<-tot-apply(OTU.tab,1,sum)
OTU.tab<-cbind(OTU.tab,unclassified)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_family.tsv",sep="\t",quote=F,row.names = T)

##########
# Genus #
##########

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_genus.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))

# Get the total number of reads
unclassified<-tot-apply(OTU.tab,1,sum)
OTU.tab<-cbind(OTU.tab,unclassified)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_genus.tsv",sep="\t",quote=F,row.names = T)


#######
# OTU #
#######

# Load the miTags table and remove Eukaryotes
OTU.tab<-fread("/Users/guillemsalazar/polybox/ETH/PROJECTS/mtags_multistudy_GitLab/Results/mtags_TARA180/Merged_table_OTU.TARA180.merged16S18S.txt",sep="\t",header=T,data.table = F)
rownames(OTU.tab)<-OTU.tab[,1]
OTU.tab<-OTU.tab[,-c(1)]
#OTU.tab<-OTU.tab[-c(grep("Eukaryota;",rownames(OTU.tab))),]
OTU.tab<-t(OTU.tab)
OTU.tab<-OTU.tab[match(rownames(OTU.tab),gsub("-",".",env.mat$Sample_name)),]
rownames(OTU.tab)<-env.mat$Barcode
colnames(OTU.tab)<-gsub(";$","",colnames(OTU.tab))

# Get the total number of reads
unclassified<-tot-apply(OTU.tab,1,sum)
OTU.tab<-cbind(OTU.tab,unclassified)

# Save
fwrite(as.data.frame(OTU.tab),file="../data/processed/miTags/to_release/mitags_tab_otu.tsv",sep="\t",quote=F,row.names = T)
