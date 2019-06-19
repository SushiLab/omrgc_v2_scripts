library(tidyverse)
library(data.table)
tabs1<-fread("Table_S1.tsv",sep="\t",data.table = F,header=T)

correct.map<-fread("correcct_incorrect_PANGAEA_ID_map.tsv",data.table = F,header=T)
pos<-match(tabs1$`PANGAEA sample id`,correct.map$PANGAEA_ID_incorrect)
tabs1$`PANGAEA sample id`<-correct.map$PANGAEA_ID[pos]

allinfo<-fread("/Users/guillemsalazar/GoogleDrive/Tara Metadata/TARA_CONTEXT_230_ALL_STATIONS.txt",sep="\t",header=T,data.table = F)
pos<-match(tabs1$`PANGAEA sample id`,allinfo$`Sample ID : 16065@TARA_barcode# : registered at PANGAEA, Data Publisher for Earth and Environmental Science (www.pangaea.de)`)
tabs1$BioSamples_ID<-allinfo$`Sample ID : 16065@BioSamples accession number (SAMEA#) : registered at the BioSamples database (http://www.ebi.ac.uk/biosamples/)`[pos]
tabs1$ENA_ID<-allinfo$`Sample ID : 16065@ENA sample accession number (ERS#) : registered at the European Nucleotides Archive (http://www.ebi.ac.uk/ENA/)`[pos]

fwrite(tabs1,"Table_S1_namesfixed.tsv",sep="\t",quote = F)
