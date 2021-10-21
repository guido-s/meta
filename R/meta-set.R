## Auxiliary functions to set arguments
##
## Package: meta
## Author: Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "") {
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  ##
  if (is.numeric(x)) {
    numeric.x <- TRUE
    idx <- x
    idx[idx < 1] <- NA
    idx[idx >= nval + 1] <- NA
  }
  else {
    numeric.x <- FALSE
    ##
    if (length(unique(tolower(x))) != length(unique(x)) |
        length(unique(tolower(val))) != length(unique(val)))
      idx <- charmatch(x, val, nomatch = NA)
    else
      idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  }
  ##
  if (anyNA(idx) || any(idx == 0)) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    ##
    if (missing(text)) {
      if (numeric.x) {
        if (nval == 1)
          vlist <- "1"
        else if (nval == 2)
          vlist <- "1 or 2"
        else
          vlist <- paste("between 1 and", nval)
      }
      else {
        if (nval == 1)
          vlist <- paste0('"', val, '"')
        else if (nval == 2)
          vlist <- paste0('"', val, '"', collapse = " or ")
        else
          vlist <- paste0(paste0('"', val[-nval], '"', collapse = ", "),
                          ', or ', '"', val[nval], '"')
      }
      ##
      if (stop.at.error)
        stop(first, name, "' must be ", vlist, addtext, ".", call. = FALSE)
      else
        return(NULL)
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else
        return(NULL)
    }
  }
  ##
  val[idx]
}
setstudlab <- function(x, k) {
  ##
  ## Set study labels
  ##
  if (is.null(x)) {
    if (k == 1)
      x <- ""
    else
      x <- seq(along = rep(1, k))
  }
  ##
  if (is.factor(x))
    x <- as.character(x)
  ##
  x
}
setunit <- function(x) {
  xname <- deparse(substitute(x))
  
  if (is.character(x)) {
    if (length(x) != 1)
      stop("Argument '", xname, "' must be a character string.")
    ##
    plotunit <- ""
    if (length(grep("cm", tolower(x))) == 1)
      plotunit <- "cm"
    else if (length(grep("inch", tolower(x))) == 1)
      plotunit <- "inch"
    else if (length(grep("mm", tolower(x))) == 1)
      plotunit <- "mm"
    ##
    if (plotunit == "")
      stop("Argument '", xname, "' must contain 'cm', 'inch', or 'mm'.")
    ##
    if (plotunit == "cm") {
      plotval <- substring(x, 1, nchar(x) - 2)
      if (length(plotval) == 0 | is.na(as.numeric(plotval)))
        stop("Argument '", xname, "' must contain at least one number.")
      ##
      res <- grid::unit(as.numeric(plotval), "cm")
    }
    else if (plotunit == "inch") {
      plotval <- substring(x, 1, nchar(x) - 4)
      if (length(plotval) == 0 | is.na(as.numeric(plotval)))
        stop("Argument '", xname, "' must contain at least one number.")
      ##
      res <- grid::unit(as.numeric(plotval), "inch")
    }
    else if (plotunit == "mm") {
      plotval <- substring(x, 1, nchar(x) - 2)
      if (length(plotval) == 0 | is.na(as.numeric(plotval)))
        stop("Argument '", xname, "' must contain at least one number.")
      ##
      res <- grid::unit(as.numeric(plotval), "mm")
    }
  }
  else
    res <- x
  
  res
}
