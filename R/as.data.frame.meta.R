as.data.frame.meta <- function(x, row.names=NULL,
                               optional=FALSE, ...){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  
  
  ## Remove element 'call' from object of class meta to get rid
  ## of an error message in meta-analyses with six studies:
  ## 'Error: evaluation nested too deeply: infinite recursion ...'
  ##
  ## NB: Element 'call' which is of length six contains information
  ##     on the function call.
  ##
  x$call <- NULL
  
  sel <- as.vector(lapply(x, length) == length(x$TE))
  
  res <- as.data.frame(x[names(x)[sel]], ...)
  
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}
