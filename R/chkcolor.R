chkcolor <- function(x, single = TRUE, name = NULL) {
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (single) {
    if (length(x) != 1)
      stop("Argument '", name, "' must have a single value.")
  }
  else if (!(is.character(x) || is.numeric(x)))
    stop("Argument '", name, "' must be a character or numeric vector.")
}
