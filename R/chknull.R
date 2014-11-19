chknull <- function(x){
  ##
  ## Check whether argument is NULL
  ##
  name <- deparse(substitute(x))
  ##
  if (is.null(x))
    stop("Argument '", name, "' is NULL.", call.=FALSE)
  ##
  invisible(NULL)
}
