chkmiss <- function(x, name=NULL){
  ##
  ## Check for missing values
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (any(is.na(x)))
    stop("Missing values in argument '", name, "'.",
         call.=FALSE)
  ##
  invisible(NULL)
}
