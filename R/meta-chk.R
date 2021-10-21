## Auxiliary functions to check function arguments
##
## Package: meta
## Author: Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
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
chkclass <- function(x, class, name = NULL) {
  ##
  ## Check class of R object
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         class, "\".", call. = FALSE)
  ##
  invisible(NULL)
}
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
chklength <- function(x, k.all, fun = "", text, name = NULL) {
  ##
  ## Check length of vector
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length(x) != k.all) {
    funcs <- c("metabin", "metacont", "metacor",
               "metagen", "metainc", "metamean",
               "metaprop", "metarate",
               "funnel", "forest.meta")
    args <- c("event.e", "n.e", "cor",
              "TE", "event.e", "n",
              "event", "event",
              "TE", "TE")
    ##
    idx <- charmatch(fun, funcs, nomatch = NA)
    if (!is.na(idx))
      argname <- args[idx]
    else
      argname <- fun
    ##
    if (missing(text))
      stop("Arguments '", argname, "' and '", name,
           "' must have the same length.",
           call. = FALSE)
    else
      stop(text, call. = FALSE)
  }
  ##
  invisible(NULL)
}
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
chklogical <- function(x, name = NULL) {
  ##
  ## Check whether argument is logical
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (is.numeric(x))
    x <- as.logical(x)
  ##
  if (length(x) !=  1 || !is.logical(x) || is.na(x))
    stop("Argument '", name, "' must be a logical.", call. = FALSE)
  ##
  invisible(NULL)
}
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
chknull <- function(x, name = NULL) {
  ##
  ## Check whether argument is NULL
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (is.null(x))
    stop("Argument '", name, "' is NULL.", call. = FALSE)
  ##
  invisible(NULL)
}
chknumeric <- function(x, min, max, zero = FALSE, length = 0,
                       name = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  ##
  ## Check numeric variable
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(NULL)
  ##
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  ##
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  ##
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      stop("Argument '", name, "' must be negative.",
           call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  ##
  invisible(NULL)
}
argid <- function(x, value) {
  if (any(x == value))
    res <- seq(along = x)[x == value]
  else
    res <- NA
  res
}
chkdeprecated <- function(x, new, old) {
  if (!is.na(argid(x, old)))
    warning("Deprecated argument '", old, "' ignored. ",
            "Use argument '", new, "' instead.",
            call. = FALSE)
  invisible(NULL)
}
