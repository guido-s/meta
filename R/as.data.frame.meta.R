as.data.frame.meta <- function(x, row.names=NULL,
                               optional=FALSE, ...){
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  sel <- as.vector(lapply(x, length) == length(x$TE))
  
  res <- as.data.frame(x[names(x)[sel]])
  ##
  res
}
