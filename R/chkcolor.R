chkcolor <- function(x, length = 0, name = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must must be a character or ",
         "numeric vector of length ", length, ".",
         call. = FALSE)
  else if (!(is.character(x) || is.numeric(x)))
    stop("Argument '", name, "' must be a character or numeric vector.",
         call. = FALSE)
}
