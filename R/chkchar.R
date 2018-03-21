chkchar <- function(x, single = TRUE, name = NULL, nchar = NULL) {
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (single) {
    if (length(x) != 1 || !is.character(x))
      stop(paste("Argument '", name,
                 "' must be a character string.", sep = ""))
    ##
    if (!is.null(nchar) && nchar(x) != nchar)
      if (nchar == 1)
        stop(paste("Argument '", name,
                   "' must be a single character.", sep = ""))
      else
        stop(paste("Argument '", name,
                   "' must be a character string of length ",
                   nchar, ".", sep = ""))
  }
  else if (!is.character(x))
    stop(paste("Argument '", name,
               "' must be a character vector.", sep = ""))
  else {
    if (!is.null(nchar) & any(nchar(x) != nchar))
      if (nchar == 1)
        stop(paste("Argument '", name,
                   "' must be a vector of single characters.", sep = ""))
      else
        stop(paste("Argument '", name,
                   "' must be a character vector where each element has ", nchar,
                   " characters.", sep = "")) 
  }
}
