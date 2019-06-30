# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure 1: Map ==========================================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(rgdal)      # for spTransform() & project()
library(ggplot2)    # for ggplot()
library(ggrepel)    # for geom_text_repel() - repel overlapping text labels
library(data.table)
library(tidyverse)

# Variables ------------------------------------------------------------------------------

# Using the Winkel tripel projection
PROJ <- "+proj=wintri +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 

# Load data ------------------------------------------------------------------------------

# Load environmental table
env.mat.metaG<-fread("zcat < ../data/processed/NOG_env.mat.metaG.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
env.mat.metaT<-fread("zcat < ../data/processed/NOG_env.mat.metaT.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)

# Transform data -------------------------------------------------------------------------

env.mat.metaG$Station.label[env.mat.metaG$Station.label=="TARA_148b"]<-"TARA_148"
env.mat.metaT$Station.label[env.mat.metaT$Station.label=="TARA_148b"]<-"TARA_148"
env.mat<-rbind(env.mat.metaG,env.mat.metaT)
env.mat<-env.mat[match(unique(env.mat$Station.label),env.mat$Station.label),]

sequenced<-NULL
for (i in 1:nrow(env.mat)){
  if ((env.mat$Station.label[i] %in% env.mat.metaG$Station.label) & (env.mat$Station.label[i] %in% env.mat.metaT$Station.label)) sequenced<-c(sequenced,"Both")
  if ((env.mat$Station.label[i] %in% env.mat.metaG$Station.label) & (env.mat$Station.label[i] %in% env.mat.metaT$Station.label)==F) sequenced<-c(sequenced,"metagenome")
  if ((env.mat$Station.label[i] %in% env.mat.metaG$Station.label)==F & (env.mat$Station.label[i] %in% env.mat.metaT$Station.label)) sequenced<-c(sequenced,"metatranscriptome")
}
env.mat$sequenced<-sequenced
env.mat <- cbind(env.mat,project(cbind(env.mat$Longitude,env.mat$Latitude), proj = PROJ))
colnames(env.mat)[(ncol(env.mat)-1):ncol(env.mat)]<-c("Lon.proj","Lat.proj")
env.mat.red<-as.data.frame(env.mat)

## WORLD MAP ##
# ======================================================================================
# Create a simple world map in Eckert IV projection with labeled graticules using ggplot.
# Repel overlapping text labels with ggrepel.
# Eckert IV is a nice looking equal-area projection :)
# https://upload.wikimedia.org/wikipedia/commons/c/c5/Ecker_IV_projection_SW.jpg
# ======================================================================================

# ~~~~~~~~~~~ Load ready to use data from GitHub ~~~~~~~~~~~ #
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
# This will load 6 objects:
#   xbl.X & lbl.Y are two data.frames that contain labels for graticule lines
#       They can be created with the code at this link: 
#       https://gist.github.com/valentinitnelav/8992f09b4c7e206d39d00e813d2bddb1
#   NE_box is a SpatialPolygonsDataFrame object and represents a bounding box for Earth 
#   NE_countries is a SpatialPolygonsDataFrame object representing countries 
#   NE_graticules is a SpatialLinesDataFrame object that represents 10 dg latitude lines and 20 dg longitude lines
#           (for creating graticules check also the graticule package or gridlines fun. from sp package)
#           (or check this gist: https://gist.github.com/valentinitnelav/a7871128d58097e9d227f7a04e00134f)
#   NE_places - SpatialPointsDataFrame with city and town points
#   NOTE: data downloaded from http://www.naturalearthdata.com/
#         here is a sample script how to download, unzip and read such shapefiles:
#         https://gist.github.com/valentinitnelav/a415f3fbfd90f72ea06b5411fb16df16

# ~~~~~~~~~~~ Project from long-lat to Eckert IV projection ~~~~~~~~~~~ #
# spTransform() is used for shapefiles and project() in the case of data frame
# for more PROJ.4 strings check the followings
#   http://proj4.org/projections/index.html
#   https://epsg.io/

# or use the short form "+proj=eck4"

# __ project the shapefiles
NE_countries.prj  <- spTransform(NE_countries, CRSobj = PROJ)
NE_graticules.prj <- spTransform(NE_graticules, CRSobj = PROJ)
NE_box.prj        <- spTransform(NE_box, CRSobj = PROJ)

# __ project long-lat coordinates columns for data frames 
# (two extra columns with projected XY are created)

# before projecting, transform NE_places to data frame to use it inside ggplot()
#NE_places.df     <- cbind(NE_places@coords, NE_places@data)
#names(NE_places.df)[1:2] <- c("lon", "lat")
#prj.coord        <- project(cbind(NE_places.df$lon, NE_places.df$lat), proj = PROJ)
#NE_places.df.prj <- cbind(prj.coord, NE_places.df)
#names(NE_places.df.prj)[1:2] <- c("X.prj","Y.prj")


ggplot() +
  geom_polygon(data = NE_countries.prj, 
               aes(long,lat, group = group), 
               colour = "gray70", fill = "gray90", size = .25) +
  geom_polygon(data = NE_box.prj, 
               aes(x = long, y = lat), 
               colour = "black", fill = "transparent", size = .25) +
  geom_point(data = env.mat.red, aes(x = Lon.proj, y = Lat.proj),colour="black",alpha = .5,size=3,pch=19) +
  geom_label_repel(data = env.mat.red, 
                   aes(x = Lon.proj, y = Lat.proj, label = sub("TARA_","",Station.label),fill=sequenced),
                   size = 5,
                   segment.colour = "black",
                   segment.alpha = .8,
                   segment.size = .45,
                   force = 1,
                   max.iter = 10e3,
                   show.legend = TRUE) +
  geom_path(data = NE_graticules.prj, 
            aes(long, lat, group = group), 
            linetype = "dotted", colour = "grey50", size = .25) +
  coord_fixed(ratio = 1) +
  theme_void() + #+ # remove the default background, gridlines & default gray color around legend's symbols
  scale_fill_manual(values=c("#c6e9af","#ed9406","#33ccff"))
ggsave("../results/figures/Figure_Map.pdf", width=40, height=20, units="cm")

## POLAR MAP ##
world<-map_data("world") 
env.mat.red.pol<-env.mat.red[env.mat.red$Latitude>50,]

p<-ggplot(data=world,aes(x=long,y=lat,group=group)) +
  geom_polygon(fill='gray60') +
  geom_point(data = env.mat.red.pol, aes(x=Longitude,y=Latitude,group=NULL),alpha=.5,size=3,color="black") +
  geom_label_repel(data = env.mat.red.pol, aes(x = Longitude, y = Latitude, label = sub("TARA_","",Station.label),group=NULL),size = 2,segment.colour = "black",segment.alpha = .5,segment.size = .25,force = 2,max.iter = 10e3,show.legend = FALSE,fill="#c6e9af") +
  geom_path(color='gray60',aes(x=long,y=lat),size=0.1) +
  geom_hline(yintercept=seq(-180,180,by=20),colour = "grey50", size = .25,alpha=0.5,linetype = "dotted") +
  geom_vline(xintercept=seq(-90,90, by=20),colour = "grey50", size = .25,alpha=0.5,linetype = "dotted") +
  coord_map("harrison",parameters=c(3,0)) +
  theme_void()
ggsave("../results/figures/Figure_Polar_map.pdf",width=15, height=7.5, units="cm")


# Table of stations
numbers = as_tibble(table(sequenced))  %>%
  dplyr::rename(stations = n)

# Table of samples
numbers$samples = c(
  nrow(env.mat.match), # Both
  nrow(env.mat.metaG)-length(which(env.mat.metaG$Barcode %in% sapply(strsplit(as.character(env.mat.match$Barcode),"-"),"[[",1))), # metaG
  nrow(env.mat.metaT)-length(which(env.mat.metaT$Barcode %in% sapply(strsplit(as.character(env.mat.match$Barcode),"-"),"[[",2))) # metaT
)

write_tsv(numbers, '../results/tables/Table_numbers_for_map.tsv')
