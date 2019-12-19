warnarg <- function(x, y, fun, cl, otherarg) {
  if (x %in% y)
    if (!missing(cl))
      warning("Argument '", x, "' has been removed from R function ", fun,
              ".\nThis argument can be used in R function ", cl, ".",
              call. = FALSE)
    else if (!missing(otherarg))
      warning("Argument '", x, "' has been replaced by argument '", otherarg,
              "' in R function ", fun, ".\nSee help page of R function ",
              fun, " for information on the use of the new argument.",
              call. = FALSE)
  ##
  invisible(NULL)
}
