chklogical <- function(x){
  ##
  ## Check whether argument is logical
  ##
  name <- deparse(substitute(x))
  ##
  if (length(x)!= 1 || !is.logical(x))
    stop("Argument '", name, "' must be a logical.", call.=FALSE)
  ##
  invisible(NULL)
}
