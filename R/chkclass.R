chkclass <- function(x, class, name=NULL){
  ##
  ## Check class of R object
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         class, "\".", call.=FALSE)
  ##
  invisible(NULL)
}
