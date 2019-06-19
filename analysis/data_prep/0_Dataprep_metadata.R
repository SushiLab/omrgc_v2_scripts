library(data.table)

# Load complete metadata file
metadata<-fread("/Users/guillemsalazar/GoogleDrive/Tara Metadata/TARA_CONTEXT_230_ALL_STATIONS.txt",sep="\t",header=T,data.table = F)

# Load the metaG and metaT map
metag.map<-fread("../data/lib/metaG.map.txt",sep="\t",data.table = F,col.names=c("PANGAEA_ID","Sample_name"),header = F)
metat.map<-fread("../data/lib/metaT.map.txt",sep="\t",data.table = F,col.names=c("PANGAEA_ID","Sample_name"),header = F)
metat.map$PANGAEA_ID<-sapply(strsplit(metat.map$PANGAEA_ID,","),"[[",1)

# Select the samples
metag.pos<-match(metag.map$PANGAEA_ID,metadata$`Sample ID : 16065@TARA_barcode# : registered at PANGAEA, Data Publisher for Earth and Environmental Science (www.pangaea.de)`)
metat.pos<-match(metat.map$PANGAEA_ID,metadata$`Sample ID : 16065@TARA_barcode# : registered at PANGAEA, Data Publisher for Earth and Environmental Science (www.pangaea.de)`)

# Select variable and change names
var.pos<-c(1:11,37:43,45:46,48,51,57,60,62,66:68,71,84:87,101,105:107,140,142:145,216:235)
oldnames<-colnames(metadata)[var.pos]
newnames<-c('Barcode','pH','CO2','CO2.partial.pressure','CO2.fugacity','HCO3','CO3','Carbon.total',
            'Alkalinity.total','Calcite.saturation.state', 'Aragonite.saturation.state', 'NO2', 'PO4', 'NO2NO3','Si',
            'Temperature','Conductivity','Salinity','Density','Oxygen','NO3','ChlorophyllA','Fluorescence','PAR.PC', 'PAR.TO','Marine.biome','Ocean.region',
            'Biogeographical.province', 'Ice.free.period', 'Iron.5m', 'Ammonium.5m', 'NO2.5m', 'NO3.5m', 'Gradient.Surface.temp(SST)',
            'Okubo.Weiss','Lyapunov', 'Residence.time', 'Depth.Mixed.Layer', 'Brunt.Väisälä','Depth.Max.O2', 'Depth.Min.O2','Nitracline','row.ref','BioSamples','ENA','Station.label','Campaign.label','Event.label',
            'Event.device','Event.comment','Event.date','Latitude','Longitude','Depth','Env.feature','Depth.nominal','Depth.top/min',
            'Depth.bottom/max','lower.size.fraction','upper.size.fraction','Sample','Method')

# Select samples and variables to save
metag.metadata<-cbind(metag.map$Sample_name,metadata[metag.pos,var.pos])
colnames(metag.metadata)<-c("Sample_name",newnames)
metat.metadata<-cbind(metat.map$Sample_name,metadata[metat.pos,var.pos])
colnames(metat.metadata)<-c("Sample_name",newnames)

# Get the correct layer and save a new variable with the correct name
metag.metadata$Sample_name_goodlayer<-paste(metag.metadata$Station.label,substr(metag.metadata$Depth,2,4),sep="_")
metat.metadata$Sample_name_goodlayer<-paste(metat.metadata$Station.label,substr(metat.metadata$Depth,2,4),sep="_")
metag.metadata$Layer<-substr(metag.metadata$Depth,2,4)
metat.metadata$Layer<-substr(metat.metadata$Depth,2,4)

# Select the sample which are present in both metaG and metaT
both_gt<-names(which(table(c(as.character(metag.metadata$Sample_name_goodlayer),as.character(metat.metadata$Sample_name_goodlayer)))==2))
both_gt.pos.metag<-match(both_gt,metag.metadata$Sample_name_goodlayer)
both_gt.pos.metat<-match(both_gt,metat.metadata$Sample_name_goodlayer)
match.metadata.metaG<-metag.metadata[both_gt.pos.metag,]
match.metadata.metaT<-metat.metadata[both_gt.pos.metat,]
match.metadata.metaG$Barcode<-paste(match.metadata.metaG$Barcode,match.metadata.metaT$Barcode,sep="-")

#for (i in 1:ncol(match.metadata.metaG)){
#  if (is.factor(match.metadata.metaG[,i])) res<-all(as.character(match.metadata.metaG[,i])==as.character(match.metadata.metaT[,i]),na.rm = T) else res<-all(match.metadata.metaG[,i]==match.metadata.metaT[,i],na.rm = T)
#  cat(colnames(match.metadata.metaG)[i],res,"\n")
#}

# Save the tables
fwrite(metag.metadata,file="../data/raw/metaG_metadata.tsv",quote = F,sep = "\t",row.names = F)
fwrite(metat.metadata,file="../data/raw/metaT_metadata.tsv",quote = F,sep = "\t",row.names = F)
fwrite(match.metadata.metaG,file="../data/raw/match_metadata.tsv",quote = F,sep = "\t",row.names = F)
fwrite(data.frame(short.name=newnames,complete.name=oldnames),file="../data/raw/metadata_name_map.tsv",quote = F,sep = "\t",row.names = F)


