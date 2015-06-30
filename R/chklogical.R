chklogical <- function(x){
  ##
  ## Check whether argument is logical
  ##
  name <- deparse(substitute(x))
  ##
  if (is.numeric(x))
    x <- as.logical(x)
  ##
  if (length(x)!= 1 || !is.logical(x) || is.na(x))
    stop("Argument '", name, "' must be a logical.", call.=FALSE)
  ##
  invisible(NULL)
}
