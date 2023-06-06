## Auxiliary functions to set arguments
##
## Package: meta
## Author: Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
## License: GPL (>= 2)
##

setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "",
                    return.NULL = TRUE, nchar.equal = FALSE) {
  val <- unique(val)
  ##
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
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
  }
  ##
  if (is.null(x))
    return(NULL)
  else
    res <- val[idx]
  ##
  if (nchar.equal && nchar(res) != nchar(x))
    res <- x
  ##
  res
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

setmethodbias <- function(x, subset) {
  oldmethod <- setchar(x, gs("meth4bias.old"),
                       stop.at.error = FALSE)
  ##
  if (is.null(oldmethod))
    if (missing(subset))
      res <- setchar(x, gs("meth4bias"), name = "method.bias")
    else
      res <- setchar(x, gs("meth4bias")[subset], name = "method.bias")
  else
    res <- switch(oldmethod,
                  rank = "Begg",
                  linreg = "Egger",
                  mm = "Thompson",
                  count = "Schwarzer",
                  score = "Harbord")
  ##
  res
}

setmethodtau <- function(method.tau, missing.tau,
                         method.predict, missing.predict,
                         warn = TRUE) {
  if (method.tau != "REML" & any(method.predict == "KR")) {
    if (missing.tau & !missing.predict) {
      if (warn)
        warning("Argument 'method.tau' set to \"REML\" as ",
                "'method.predict' = \"KR\".",
                call. = FALSE)
      ##
      method.tau <- "REML"
    }
  }
  ##
  method.tau
}

setmethodpredict <- function(method.predict, missing.predict,
                             method.tau, missing.tau,
                             warn = TRUE) {
  if (method.tau != "REML" & any(method.predict == "KR")) {
    if (!missing.tau & missing.predict) {
      method.predict[method.predict == "KR"] <- "HTS"
    }
    else if (!missing.tau & !missing.predict) {
      if (warn)
        warning("Argument 'method.predict' set to \"HTS\" instead of ",
                "\"KR\" as 'method.tau' != \"REML\".",
                call. = FALSE)
      method.predict[method.predict == "KR"] <- "HTS"
    }
    else
      method.predict[method.predict == "KR"] <- "HTS"
  }
  ##
  method.predict
}

setVal <- function(data, varname, default = NULL) {
  if (isCol(data, varname))
    return(data[[varname]])
  else
    return(default)
}
