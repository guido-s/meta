chkchar <- function(x, single=TRUE){
  if (single){
    if (length(x)!=1 || !is.character(x))
      stop(paste("Argument '", deparse(substitute(x)),
                 "' must be a character string.", sep=""))
  }
  else
    if (!is.character(x))
      stop(paste("Argument '", deparse(substitute(x)),
                 "' must be a character vector.", sep=""))
}
