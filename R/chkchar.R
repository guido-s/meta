chkchar <- function(x, length = 0, name = NULL, nchar = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a character vector of length ",
         length, ".",
         call. = FALSE)
  ##
  if (length == 1) {
    if (!is.null(nchar) && nchar(x) != nchar)
      if (nchar == 1)
        stop("Argument '", name, "' must be a single character.",
             call. = FALSE)
      else
        stop("Argument '", name, "' must be a character string of length ",
             nchar, ".",
             call. = FALSE)
  }
  ##
  if (!is.character(x))
    stop("Argument '", name, "' must be a character vector.")
  else {
    if (!is.null(nchar) & any(nchar(x) != nchar))
      if (nchar == 1)
        stop("Argument '", name, "' must be a vector of single characters.",
             call. = FALSE)
      else
        stop("Argument '", name, "' must be a character vector where ",
             "each element has ", nchar, " characters.",
             call. = FALSE)
  }
}
