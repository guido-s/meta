nnt.metabias <- function(x, p.c){
  
  res <- list()
  
  if (x$sm=="RD")
    res$fixed <- meta::ci(x$TE.fixed, x$seTE.fixed)
  ##
  if (x$sm=="RR")
    res$fixed <- p.c*meta::ci(x$TE.fixed, x$seTE.fixed)
  
  res  
}
