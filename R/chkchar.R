chkchar <- function(x, single=TRUE, name=NULL){
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (single){
    if (length(x)!=1 || !is.character(x))
      stop(paste("Argument '", name,
                 "' must be a character string.", sep=""))
  }
  else
    if (!is.character(x))
      stop(paste("Argument '", name,
                 "' must be a character vector.", sep=""))
}
