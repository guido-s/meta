chklevel <- function(x, length = 0, ci = TRUE, name = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  ##
  ## Check for levels of confidence interval / contour level
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  if (ci)
    "level for confidence interval (range: 0-1)"
  else
    "contour levels (range: 0-1)"
  ##
  if (!is.numeric(x))
    if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
    else
      stop("Argument '", name, "' must be numeric.",
           call. = FALSE)
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  ##
  if (any(x <= 0, na.rm = TRUE) | any(x >= 1, na.rm = TRUE))
    stop("Argument '", name, "' must be a numeric between 0 and 1.",
         call. = FALSE)
  ##
  invisible(NULL)
}
