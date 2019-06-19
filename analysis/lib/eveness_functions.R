SimpE<-function(x,zeros=T){
  if (zeros==F) x <- as.vector(x[x>0])
  S <- length(x)
  x = as.data.frame(x)
  D <- diversity(x, "inv")
  E <- (D)/S
  E}

PielouE<-function(x,zeros=T){
  if (zeros==F) x <- as.vector(x[x>0])
  H <- diversity(x)
  S <- length(x)
  J <- H/log(S) 
  J
}

evenness<-function(x){
  res<-rbind(apply(x,1,SimpE),apply(x,1,PielouE))
  colnames(res)<-rownames(x)
  rownames(res)<-c("SimpE","PielouE")
  res
}
