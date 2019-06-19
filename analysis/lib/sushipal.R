sushi.palette<-function(n=15,alpha=1){
  if (n>15) stop("Maximum number of colors is 15. Set an 'n' value between 1 and 15")
  if (alpha>1) stop("Alpha value must be between 0 1n 1")
  red<-c(0,192,0,141,188,217,255,252,251,204,253,128,190,128,0)/255
  green<-c(112,0,176,211,128,217,255,205,128,235,180,177,186,128,0)/255
  blue<-c(192,0,80,199,189,217,179,229,114,197,98,211,218,128,0)/255
  rgb(red[1:n],green[1:n],blue[1:n],rep(alpha,n))
}

cat("\n The Sushi palette was loaded. \nUsage: sushi.palette(n=15,alpha=1)\n   n       number of colors [1-15]\n   alpha   transparency value [0-1]")
