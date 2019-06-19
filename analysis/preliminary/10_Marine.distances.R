mar.dist<-function(stations.positions=NULL,plotting=T,open.straits=T,out.in.km=T,z.cutoff=100){
  require(maptools)
  require(raster)
  require(rgdal)
  require(gdistance)
  stations.positions<-as.matrix(stations.positions)
  #stations.positions[,2]<--stations.positions[,2]
  ##################
  # A - Bathymetry #
  ##################
  
  # Creating a raster for the distance calculation :
  cat("Creating raster for distance calculation...\n")
  # Setting the projection (WGS84, official one) : 
  crs <- CRS("+init=epsg:4326")
  
  # make a raster with all land having "NA" value
  # use own shapefile or raster if possible
  # the wrld_simpl data set is from maptools package
  data(wrld_simpl)
  world <- wrld_simpl
  world.shp <- spTransform(world, crs)
  ras <- raster(nrow=300, ncol=300)
  crs(ras) <- crs(world.shp)
  extent(ras) <- extent(world.shp)
  # rasterize will set ocean to NA so we just inverse it
  # and set water to "1"
  # land is equal to zero because it is "NOT" NA
  world.mask <- rasterize(world.shp, ras)
  world.ras <- is.na(world.mask)
  
  # Setting land to NA
  world.ras[world.ras==z.cutoff] <- NA
  
  if (open.straits==T){
    cat("Opening Gibraltar, Red Sea and Arctic Bay straits...\n")
    #Opening Gibraltar
    world.ras[83,146]=1
    world.ras[84,146]=1
    #Opening Red Sea
    world.ras[123,187]=1
    world.ras[123,186]=1
    #Open Arctic Bay (St206)
    world.ras[22,106]=1}
  
  # create a Transition object from the raster
  cat("Creating transition matrix...\n")
  tr <- transition(world.ras, mean, 8)
  tr = geoCorrection(tr, scl=FALSE) # Correcting the transition object.
  
  # distance matrix excluding the land    
  A <- accCost(tr, stations.positions)
  cost.distances=costDistance(tr, stations.positions)
  
  distance<-as.matrix(cost.distances)
  A <- mask(A, world.mask, inverse=TRUE)
  
  # Plotting the distance map
  if (plotting==T){
    colfunc <- colorRampPalette(c("beige", rgb(0.04, 0.34, 0.62)))
    plot(A,col=colfunc(100))
    points(stations.positions,pch=19,col=rgb(0,0,0,0.5))}
  if (out.in.km==T) {
    cat("Marine distance matrix computed. Units in 'km'.\n")
    distance<-distance/1000}
  else {
    cat("Marine distance matrix computed. Units in 'm'.\n")}
  if (length(which(is.infinite(distance)))){
    warning("Infinite distances found. Looks like some coordinates are too close to a land mass.\nTry to increase the 'z.cutoff' value.")}
  colnames(distance)<-rownames(stations.positions)
  rownames(distance)<-rownames(stations.positions)
  distance
}


env.mat.match<-fread("../data/processed/KO_env.mat.match.txt",header=T,sep="\t",data.table = F)
rownames(env.mat.match)<-env.mat.match$Barcode
env.mat.metaG<-fread("../data/processed/KO_env.mat.metaG.txt",header=T,sep="\t",data.table = F)
rownames(env.mat.metaG)<-env.mat.metaG$Barcode
env.mat.metaT<-fread("../data/processed/KO_env.mat.metaT.txt",header=T,sep="\t",data.table = F)
rownames(env.mat.metaT)<-env.mat.metaT$Barcode


mar.dist.match<-mar.dist(env.mat.match[,c("Longitude","Latitude")])
mar.dist.metaG<-mar.dist(env.mat.metaG[,c("Longitude","Latitude")])
mar.dist.metaT<-mar.dist(env.mat.metaT[,c("Longitude","Latitude")])
mar.dist.all<-mar.dist(unique(rbind(env.mat.metaT[,c("Longitude","Latitude")],env.mat.metaG[,c("Longitude","Latitude")])))

write.table(mar.dist.match,file="../data/processed/mardist.match.txt",sep="\t",quote=F,row.names = T,col.names = NA)
write.table(mar.dist.metaG,file="../data/processed/mardist.metaG.txt",sep="\t",quote=F,row.names = T,col.names = NA)
write.table(mar.dist.metaT,file="../data/processed/mardist.metaT.txt",sep="\t",quote=F,row.names = T,col.names = NA)
write.table(mar.dist.all,file="../data/processed/mardist.all.txt",sep="\t",quote=F,row.names = T,col.names = NA)
