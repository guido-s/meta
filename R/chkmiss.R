chkmiss <- function(x){
  ##
  ## Check for missing values
  ##
  name <- deparse(substitute(x))
  ##
  if (any(is.na(x)))
    stop("Missing values in argument '", name, "'.",
         call.=FALSE)
  ##
  invisible(NULL)
}
