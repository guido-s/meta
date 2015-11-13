chknull <- function(x, name=NULL){
  ##
  ## Check whether argument is NULL
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (is.null(x))
    stop("Argument '", name, "' is NULL.", call.=FALSE)
  ##
  invisible(NULL)
}
