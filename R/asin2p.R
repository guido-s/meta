asin2p <- function(x, n=NULL){
  if (is.null(n)){
    ## Back transformation of arcsine transformation:
    res <- sin(x)^2
  }
  else{
    ## Back transformation of Freeman-Tukey arcsine transformation:
    res <- 0.5*(1-sign(cos(x))*sqrt(1-(sin(x)+(sin(x)-1/sin(x))/n)**2))
  }
  res
}
