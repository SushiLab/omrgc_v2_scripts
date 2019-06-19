# Load palette
library(patchwork)
library(ggplot2)
library(vegan)
library(data.table)
library(tidyverse)
source("../data/lib/sushipal.R")
palette(sushi.palette(alpha=0.7)[c(2,3,1,14)])

dat<-fread("/Users/guillemsalazar/GoogleDrive/Tara Metadata/TARA_CONTEXT_230_ALL_STATIONS.txt",sep="\t",header=T,data.table = F)

tmp<-dat %>%
  select(Latitude="Latitude : 1600@of the closest geographic coordinate on a coast : NA",Temperature="Temperature (degC) : 717@median value (50th percentile, second quartile, Q2) in the selected environmental feature around the sampling date and location : calculated from in situ sensor data calibrated using factory settings",Layer="Environmental feature : 139620@[abbreviation], full name (ENVO:ID) from which this sample was collected : terms registered at EnvO, the Environmental Ontology (http://environmentontology.org/)") %>%
  filter(grepl(paste(c("SRF","DCM"), collapse="|"),Layer)) %>%
  mutate(Layer= fct_recode(Layer, "SRF" = "[SRF] surface water layer (ENVO:00010504)","DCM"="[DCM] deep chlorophyll maximum layer (ENVO:01000326)")) %>%
  filter(!is.na(Temperature) & !is.na(Latitude)) %>%
  distinct()

p1<-ggplot(data=tmp,aes(x=Latitude,y=Temperature,col=Layer)) +
  geom_point(alpha=0.5,size=1) +
  #geom_line() +
  geom_smooth(se=F,span=0.4,size=0.5) +
  xlab("Latitude (deg)") +
  ylab("Temperature (ºC)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values=c("#0cb14c","#b82327","#136fba"))

p2<-ggplot(data=tmp,aes(x=cut(abs(Latitude),c(0,10,30,50,70,90)),y=Temperature,fill=Layer)) +
  geom_boxplot() +
  #geom_violin(scale = "width",draw_quantiles = 0.5) +
  #geom_point(alpha=0.5,size=1) +
  #geom_line() +
  #geom_smooth(se=F,span=0.4,size=0.5) +
  xlab("Absolute latitude (deg)") +
  ylab("Temperature (ºC)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#0cb14c","#b82327","#136fba"))
p1 | p2
ggsave("../results/FigSX_Lat_Temp_fct_ALL_v2.pdf",width=unit(10,"cm"),height=unit(5,"cm"))
