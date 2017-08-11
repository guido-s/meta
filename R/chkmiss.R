chkmiss <- function(x, name = NULL) {
  ##
  ## Check for missing values
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (anyNA(x))
    stop("Missing values in argument '", name, "'.",
         call. = FALSE)
  ##
  invisible(NULL)
}
