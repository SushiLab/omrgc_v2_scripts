cor.cutoff<-function(x,r.min=0.8,blocksize=1000,file=NULL){
  ncols<-ncol(x)
  n<-ceiling(ncols/blocksize)
  ini<-1
  INI<-NULL
  FIN<-NULL
  for (i in 1:n){
    fin<-ini+(blocksize-1)
    if (i==n) fin<-ncols
    INI<-c(INI,ini)
    FIN<-c(FIN,fin)
    ini<-fin+1}
  i<-NULL
  j<-NULL
  RES<-NULL
  for (i in 1:(length(INI)-1)){
    for (j in (i+1):length(INI)){
      cat("Computing block ",i," against block ",j,"\n")
      tmp<-x[,c(INI[i]:FIN[i],INI[j]:FIN[j])] %>%
        cor()
      tmp[upper.tri(tmp)]<-NA
      diag(tmp)<-NA
      res<-tmp %>%
        as.data.frame() %>%
        rownames_to_column(var="x") %>%
        gather(key="y",value="r",-x) %>%
        filter(!is.na(r)) %>%
        filter(r>r.min)
      RES<-bind_rows(RES,res)}
    RES<-RES %>%
      distinct()
  }
  RES
}

cor2.cutoff<-function(x,y,r.min=0.8,blocksize=1000,file=NULL){
  ncols<-ncol(x)
  n<-ceiling(ncols/blocksize)
  ini<-1
  INI<-NULL
  FIN<-NULL
  for (i in 1:n){
    fin<-ini+(blocksize-1)
    if (i==n) fin<-ncols
    INI<-c(INI,ini)
    FIN<-c(FIN,fin)
    ini<-fin+1}
  i<-NULL
  j<-NULL
  RES<-NULL
  for (i in 1:(length(INI)-1)){
    for (j in (i+1):length(INI)){
      cat("Computing block ",i," against block ",j,"\n")
      
      # First matrix
      tmp.x<-x[,c(INI[i]:FIN[i],INI[j]:FIN[j])] %>%
        cor()
      tmp.x[upper.tri(tmp.x)]<-NA
      diag(tmp.x)<-NA
      res.x<-tmp.x %>%
        as.data.frame() %>%
        rownames_to_column(var="x") %>%
        gather(key="y",value="r",-x) %>%
        filter(!is.na(r))
      # Second matrix
      tmp.y<-y[,c(INI[i]:FIN[i],INI[j]:FIN[j])] %>%
        cor()
      tmp.y[upper.tri(tmp.y)]<-NA
      diag(tmp.y)<-NA
      res.y<-tmp.y %>%
        as.data.frame() %>%
        rownames_to_column(var="x") %>%
        gather(key="y",value="r",-x) %>%
        filter(!is.na(r))
      res.y$r<-res.x$r*res.y$r
      res<-res.y %>%
        filter(r>r.min)
      RES<-bind_rows(RES,res)}
    RES<-RES %>%
      distinct()
  }
  RES
}