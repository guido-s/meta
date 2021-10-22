## Auxiliary functions
##
## Package: meta
## Author: Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
bylevs <- function(x) {
  if (is.factor(x))
    res <- levels(factor(x))
  else
    res <- unique(x)
  res
}
byvarname <- function(x) {
  ##
  ## Determine name of subgroup variable
  ##
  res <- as.character(x)
  if (length(res) > 1 & res[1] == "$")
    res <- res[length(res)]
  if (length(res) == 0 || length(res) > 1)
    res <- "subgroup"
  ##
  res
}
int2num <- function(x) {
  ##
  ## Convert integer to numeric
  ##
  if (is.integer(x))
    res <- as.numeric(x)
  else
    res <- x
  ##
  res
}
npn <- function(x) {
  ##
  ## Check for non-positive values in vector
  ##
  selNA <- is.na(x)
  res <- selNA
  if (sum(!selNA) > 0)
    x[!selNA] <- x[!selNA] <= 0
  ##
  res
}
replaceNULL <- function(x, replace = NA) {
  if (is.null(x))
    return(replace)
  x
}
replaceNA <- function(x, replace = NA) {
  if (is.null(x))
    return(x)
  else
    x[is.na(x)] <- replace
  x
}
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
catchvar <- function(varname, x, mf) {
  res <- NULL
  error <-
    try(res <- eval(mf[[match(varname, names(mf))]],
                    x,
                    enclos = sys.frame(sys.parent())),
        silent = TRUE)
  ##
  if (class(error) == "try-error") {
    res <- eval(mf[[match(varname, names(mf))]],
                x$data, enclos = NULL)
  }
  ##
  res
}
augment <- function(x, len, fun) {
  if (length(x) > 1)
    chklength(x, len, fun)
  else
    x <- rep(x, len)
  x
}
stoponly <- function(arg, val, func)
  stop("Argument ", arg, " =\"", val, "\"",
               " only defined for meta-analysis conducted with ",
               func, ".",
               call. = FALSE)
deprecated <- function(newvar, newmiss, args, old, warn = TRUE) {
  ##
  new <- deparse(substitute(newvar))
  ##
  if (length(args) == 0)
    return(newvar)
  ##
  if (is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (!is.na(charmatch(old, additional.arguments)))
    if (!newmiss) {
      if (warn)
        warning("Deprecated argument '", old, "' ignored as ",
                "'", new, "' is also provided.",
                call. = FALSE)
      return(newvar)
    }
    else {
      if (warn)
        warning("Use argument '", new, "' instead of '",
                old, "' (deprecated).",
                call. = FALSE)
      return(args[[charmatch(old, additional.arguments)]])
    }
  else
    return(newvar)
}
deprecated2 <- function(newvar, newmiss, oldvar, oldmiss, warn = TRUE) {
  ##
  new <- deparse(substitute(newvar))
  old <- deparse(substitute(oldvar))
  ##
  if (newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & !oldmiss) {
    if (warn)
      warning("Deprecated argument '", old, "' ignored as ",
              "'", new, "' is also provided.",
              call. = FALSE)
    return(newvar)
  }
  else if (newmiss & !oldmiss) {
    if (warn)
      warning("Use argument '", new, "' instead of '",
              old, "' (deprecated).",
              call. = FALSE)
    return(oldvar)
  }
}
runNN <- function(func, args) {
  args <- args[!sapply(args, is.null)]
  do.call(func, args)
}
