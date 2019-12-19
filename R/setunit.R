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
