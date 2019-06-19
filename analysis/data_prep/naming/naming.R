library(data.table)
library(tidyverse)

# Load sample names sent to Genoscope
orig.names<-fread("read_filenames_forGenoscope_AA.tsv",sep="\t",header=T)

# Select only the variables needed
orig.names <- orig.names %>%
  select(old_name=sample,read_filename=`read filename`,metaG_T=`metaG/metaT`)

# Remove the old metaG samples (243 from TO 2015)
to.remove.pos<-grep("TARA_",orig.names$read_filename)
length(unique(orig.names$old_name[to.remove.pos])) # check that this correspond to 243 samples
orig.names.red<-orig.names[-to.remove.pos,]

# Get the portion of the string in the read filename useful for mapping to the PANGAEA ID
orig.names.red$mapping_id<-paste(sapply(strsplit(orig.names.red$read_filename,"_"),"[[",1),sapply(strsplit(orig.names.red$read_filename,"_"),"[[",2),sapply(strsplit(orig.names.red$read_filename,"_"),"[[",3),sapply(strsplit(orig.names.red$read_filename,"_"),"[[",4),sep="_")

# Load the mapping file from Genoscope
genoscope.map<-fread("Barcodes_Guillem.tsv",sep="\t",header=F)

# Rename genoscope map
genoscope.map <- genoscope.map %>%
  select(PANGAEA_ID=V1,mapping_id=V2) %>%
  mutate(PANGAEA_ID=paste("TARA_",PANGAEA_ID,sep=""))

# Join both tables
merged<-left_join(orig.names.red,genoscope.map,by="mapping_id") %>%
  na.exclude()

# Check that there is 127 metaG samples and 187 metaT samples
merged %>%
  select(old_name,metaG_T) %>%
  unique() %>%
  group_by(metaG_T) %>%
  summarise(n=n())

# Load the info for the 243 metaG samples from Sunagawa 2015 and get the PANGAEA_ID
metag243<-fread("OM.CompanionTables.tsv",sep="\t",header=T)
metag243$`Sample label [TARA_station#_environmental-feature_size-fraction]`<-paste(metag243$`Sample label [TARA_station#_environmental-feature_size-fraction]`,"_G",sep="")
metag243$`Sample label [TARA_station#_environmental-feature_size-fraction]`<-gsub("<","lt",metag243$`Sample label [TARA_station#_environmental-feature_size-fraction]`)

orig.names.243<-orig.names[to.remove.pos,]
orig.names.243$mapping_id<-NA
orig.names.243$PANGAEA_ID<-metag243$`PANGAEA sample identifier`[match(orig.names.243$old_name,metag243$`Sample label [TARA_station#_environmental-feature_size-fraction]`)]

# Combine the old an new data
merged<-rbind(merged,orig.names.243)

# Check that there is 370 metaG samples and 187 metaT samples
merged %>%
  select(old_name,metaG_T) %>%
  unique() %>%
  group_by(metaG_T) %>%
  summarise(n=n())

# Save this data
fwrite(merged,file="MAIN_MAPPING_byfiles.tsv",quote = F,sep = "\t")

# Reshape the mapping file with one old_name per raw
merged.bysample<-merged %>%
  select(old_name,metaG_T,PANGAEA_ID) %>%
  group_by(old_name) %>%
  summarise(PANGAEA_ID=paste(unique(PANGAEA_ID),collapse = ",")) %>%
  as.data.frame()

# Save this data
fwrite(merged.bysample,file="MAIN_MAPPING_bysample.tsv",quote = F,sep = "\t")