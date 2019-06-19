# Distance functions
sqr.euc.dist<-function(taula,dim.div=T){
  taula<-as.matrix(taula)
  if (dim.div==T) n.dim<-ncol(taula) else n.dim<-1
  res<-sapply(1:nrow(taula),function(y){sapply(1:nrow(taula),function(x){sum((taula[y,]-taula[x,])^2)})})/n.dim
  rownames(res)<-rownames(taula)
  colnames(res)<-rownames(taula)
  res
}
sqr.euc.int<-function(taula1,taula2,dim.div=T){
  taula1<-as.matrix(taula1)
  taula2<-as.matrix(taula2)
  if (dim.div==T) n.dim<-ncol(taula1) else n.dim<-1
  res<-sapply(1:nrow(taula1),function(y){sapply(1:nrow(taula1),function(x){sum(2*(taula1[y,]-taula1[x,])*(taula2[y,]-taula2[x,]))})})/n.dim
  rownames(res)<-rownames(taula1)
  colnames(res)<-rownames(taula1)
  res
}

sqr.euc.int.aprox<-function(taula1,taula2,dim.div=T){
  taula1<-as.matrix(taula1)
  taula2<-as.matrix(taula2)
  if (dim.div==T) n.dim<-ncol(taula1) else n.dim<-1
  tmp<-2*as.matrix(dist(taula1))*as.matrix(dist(taula1))
  tmp2<-sapply(1:nrow(taula1),function(y){sapply(1:nrow(taula1),function(x){cor((taula1[y,]-taula1[x,]),(taula2[y,]-taula2[x,]))})})/n.dim
  res<-tmp*tmp2
  rownames(res)<-rownames(taula1)
  colnames(res)<-rownames(taula1)
  res
}

varpart.sqr.euc.mean<-function(mat.T,mat.G,mat.E,tol=1E-09){
  # Distance functions
  sqr.euc.dist<-function(taula,dim.div=T){
    taula<-as.matrix(taula)
    if (dim.div==T) n.dim<-ncol(taula) else n.dim<-1
    res<-sapply(1:nrow(taula),function(y){sapply(1:nrow(taula),function(x){sum((taula[y,]-taula[x,])^2)})})/n.dim
    rownames(res)<-rownames(taula)
    colnames(res)<-rownames(taula)
    res
  }
  sqr.euc.int<-function(taula1,taula2,dim.div=T){
    taula1<-as.matrix(taula1)
    taula2<-as.matrix(taula2)
    if (dim.div==T) n.dim<-ncol(taula1) else n.dim<-1
    res<-sapply(1:nrow(taula1),function(y){sapply(1:nrow(taula1),function(x){sum(2*(taula1[y,]-taula1[x,])*(taula2[y,]-taula2[x,]))})})/n.dim
    rownames(res)<-rownames(taula1)
    colnames(res)<-rownames(taula1)
    res
  }
  
  # Compute distance components
  cat("Computing square euclidean distance of mat.T \n")
  res.T<-as.dist(sqr.euc.dist(mat.T))
  cat("Computing square euclidean distance of mat.G \n")
  res.G<-as.dist(sqr.euc.dist(mat.G))
  cat("Computing square euclidean distance of mat.E \n")
  res.E<-as.dist(sqr.euc.dist(mat.E))
  cat("Computing interaction component\n")
  res.int<-as.dist(sqr.euc.int(mat.G,mat.E))
  
  # Check that res.T=res.G+res.E+res.int
  if (max(c(res.T)-c(c(res.G)+c(res.E)+c(res.int)))>1E-09) {
    plot(res.T,res.G+res.E+res.int)
    stop("The equality 'metaT=Abundance+Expression+Interaction' is not met with tolerance = ",tol,"\n Make sure that the input matrices are correct.")}
  
  # Compute means
  res<-c(mean(res.G),mean(res.E),mean(res.int))
  names(res)<-c("Abundance","Expression","interaction")
  #res.norm<-res/sum(abs(res))
  res.norm<-res/sum(res)
  list(components=res,components.norm=res.norm)
}

varpart.sqr.euc.all<-function(mat.T,mat.G,mat.E,tol=1E-09){
  # Distance functions
  sqr.euc.dist<-function(taula,dim.div=T){
    taula<-as.matrix(taula)
    if (dim.div==T) n.dim<-ncol(taula) else n.dim<-1
    res<-sapply(1:nrow(taula),function(y){sapply(1:nrow(taula),function(x){sum((taula[y,]-taula[x,])^2)})})/n.dim
    rownames(res)<-rownames(taula)
    colnames(res)<-rownames(taula)
    res
  }
  sqr.euc.int<-function(taula1,taula2,dim.div=T){
    taula1<-as.matrix(taula1)
    taula2<-as.matrix(taula2)
    if (dim.div==T) n.dim<-ncol(taula1) else n.dim<-1
    res<-sapply(1:nrow(taula1),function(y){sapply(1:nrow(taula1),function(x){sum(2*(taula1[y,]-taula1[x,])*(taula2[y,]-taula2[x,]))})})/n.dim
    rownames(res)<-rownames(taula1)
    colnames(res)<-rownames(taula1)
    res
  }
  
  # Compute distance components
  cat("Computing square euclidean distance of mat.T \n")
  res.T<-sqr.euc.dist(mat.T)
  cat("Computing square euclidean distance of mat.G \n")
  res.G<-sqr.euc.dist(mat.G)
  cat("Computing square euclidean distance of mat.E \n")
  res.E<-sqr.euc.dist(mat.E)
  cat("Computing interaction component\n")
  res.int<-sqr.euc.int(mat.G,mat.E)
  
  diag(res.G)<-NA
  res.G[upper.tri(res.G)]<-NA
  diag(res.T)<-NA
  res.T[upper.tri(res.T)]<-NA
  diag(res.E)<-NA
  res.E[upper.tri(res.E)]<-NA
  diag(res.int)<-NA
  res.int[upper.tri(res.int)]<-NA
  
  res.T<-mutate(as.data.frame(res.T),sample1=rownames(as.data.frame(res.T))) %>% gather("sample2","metaT",-sample1) %>% filter(!is.na(metaT))
  res.G<-mutate(as.data.frame(res.G),sample1=rownames(as.data.frame(res.G))) %>% gather("sample2","Abundance",-sample1) %>% filter(!is.na(Abundance))
  res.E<-mutate(as.data.frame(res.E),sample1=rownames(as.data.frame(res.E))) %>% gather("sample2","Expression",-sample1) %>% filter(!is.na(Expression))
  res.int<-mutate(as.data.frame(res.int),sample1=rownames(as.data.frame(res.int))) %>% gather("sample2","interaction",-sample1) %>% filter(!is.na(interaction))
  
  # Check that res.T=res.G+res.E+res.int
  if (max(res.T$metaT-c(res.G$Abundance+res.E$Expression+res.int$interaction))>1E-09) {
    plot(res.T$metaT,res.G$Abundance+res.E$Expression+res.int$interaction)
    stop("The equality 'metaT=Abundance+Expression+Interaction' is not met with tolerance = ",tol,"\n Make sure that the input matrices are correct.")}
  
  # Compute means
  res<-cbind(res.G,Expression=res.E$Expression,interaction=res.int$interaction)
  #res.norm<-cbind(res[,1:2],t(apply(res[,3:5],1,function(x){x/sum(abs(x))})))
  res.norm<-cbind(res[,1:2],t(apply(res[,3:5],1,function(x){x/sum(x)})))
  list(components=res,components.norm=res.norm)
}
