chkchar <- function(x, single = TRUE, name = NULL, nchar = NULL) {
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (single) {
    if (length(x) != 1 || !is.character(x))
      stop("Argument '", name, "' must be a character string.")
    ##
    if (!is.null(nchar) && nchar(x) != nchar)
      if (nchar == 1)
        stop("Argument '", name, "' must be a single character.")
      else
        stop("Argument '", name, "' must be a character string of length ",
             nchar, ".")
  }
  else if (!is.character(x))
    stop("Argument '", name, "' must be a character vector.")
  else {
    if (!is.null(nchar) & any(nchar(x) != nchar))
      if (nchar == 1)
        stop("Argument '", name, "' must be a vector of single characters.")
      else
        stop("Argument '", name, "' must be a character vector where ",
             "each element has ", nchar, " characters.")
  }
}
