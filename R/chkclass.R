chkclass <- function(x, class){
  ##
  ## Check class of R object
  ##
  name <- deparse(substitute(x))
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         class, "\".", call.=FALSE)
  ##
  invisible(NULL)
}
