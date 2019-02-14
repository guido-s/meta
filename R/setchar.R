setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE) {
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
          vlist <- paste('"', val, '"', sep = "")
        else if (nval == 2)
          vlist <- paste('"', val, '"', collapse = " or ", sep = "")
        else
          vlist <- paste(paste('"', val[-nval],
                               '"', collapse = ", ", sep = ""),
                         ', or ', '"', val[nval], '"', sep = "")
      }
      ##
      if (stop.at.error)
        stop(first, name, "' must be ", vlist, ".", call. = FALSE)
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
