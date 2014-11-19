summary.rm5 <- function(object, comp.no, outcome.no, ...){
  
  
  ##
  ##
  ## (1) Check for rm5 object
  ##
  ##
  chkclass(object, "rm5")
  
  
  if (missing(comp.no))
    comp.no <- unique(object$comp.no)
  ##
  res <- list()
  ##
  n <- 1
  ##
  for (i in comp.no){
    if (missing(outcome.no))
      jj <- unique(object$outcome.no[object$comp.no==i])
    else
      jj <- outcome.no
    for (j in jj){
      ##
      res[[n]] <- summary(metacr(object, i, j, ...))
      ##
      n <- n+1
      ##
    }
  }
  ##
  class(res) <- "summary.rm5"
  
  
  res
}
