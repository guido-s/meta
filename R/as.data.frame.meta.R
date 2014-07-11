as.data.frame.meta <- function(x, row.names=NULL,
                               optional=FALSE, ...){
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  ## Upgrade meta objects created with older versions of meta
  ##
  if (!(!is.null(x$version) &&
        as.numeric(unlist(strsplit(x$version, "-"))[1]) >= 3.7))
    x <- update(x, warn=FALSE)
  
  
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
