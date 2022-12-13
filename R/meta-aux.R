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

byvarname <- function(argname, matchcall) {
  ##
  ## Determine name of subgroup variable
  ##
  res <- as.character(matchcall[[match(argname, names(matchcall))]])
  ##
  if (length(res) > 1 & res[1] == "$")
    res <- res[length(res)]
  ##
  if (length(res) == 0 || length(res) > 1)
    res <- "subgroup"
  ##
  res
}

catch <- function(argname, matchcall, data, encl) {
  ##
  ## Catch value for argument
  ##
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
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
  if (inherits(error, "try-error")) {
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

runNN <- function(func, args, warn = TRUE) {
  args <- args[!sapply(args, is.null)]
  if (warn)
    do.call(func, args)
  else
    suppressWarnings(do.call(func, args))
}

setNA3 <- function(x) {
  res <- x
  ##
  res$Q.b.random <- NA
  res$df.Q.b <- NA
  res$df.Q.b.random <- NA
  res$pval.Q.b.random <- NA
  ##
  res
}

mismatch <- function(x, y, var) {
  x <- x$data[[var]]
  y <- y$data[[var]]
  ##
  bothnull <- is.null(x) & is.null(y)
  ##
  if (bothnull)
    return(FALSE)
  else {
    if (!is.null(x) & !is.null(y))
      return(any(x != y))
    else
      return(TRUE)
  }
}

addHet <- function(x, hcc, within = TRUE) {
  
  res <- x
  ##
  if (within) {
    res$Q.w.random <- hcc$Q.resid
    res$df.Q.w.random <- hcc$df.Q.resid
    res$pval.Q.w.random <- hcc$pval.Q.resid
  }
  ##
  res$tau2.resid <- hcc$tau2.resid
  res$lower.tau2.resid <- hcc$lower.tau2.resid
  res$upper.tau2.resid <- hcc$upper.tau2.resid
  ##
  res$tau.resid <- hcc$tau.resid
  res$lower.tau.resid <- hcc$lower.tau.resid
  res$upper.tau.resid <- hcc$upper.tau.resid
  res$sign.lower.tau.resid <- hcc$sign.lower.tau.resid
  res$sign.upper.tau.resid <- hcc$sign.upper.tau.resid
  ##
  res$Q.resid <- hcc$Q.resid
  res$df.Q.resid <- hcc$df.Q.resid
  res$pval.Q.resid <- hcc$pval.Q.resid
  ##
  res$H.resid <- hcc$H.resid
  res$lower.H.resid <- hcc$lower.H.resid
  res$upper.H.resid <- hcc$upper.H.resid
  ##
  res$I2.resid <- hcc$I2.resid
  res$lower.I2.resid <- hcc$lower.I2.resid
  res$upper.I2.resid <- hcc$upper.I2.resid
  ##
  res
}

backward <- function(x) {
  res <- x
  ##
  res$w.fixed <- res$w.common
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  res$zval.fixed <- res$zval.common
  ##
  res$text.fixed <- res$text.common
  res$text.w.fixed <- res$text.w.common
  ##
  res$df.hakn <- res$df.hakn.ci
  res$seTE.hakn <- res$seTE.hakn.ci
  res$seTE.hakn.adhoc <- res$seTE.hakn.adhoc.ci
  ##
  if (!is.null(res$subgroup)) {
    res$byvar <- res$subgroup
    res$bylab <- res$subgroup.name
    res$print.byvar <- res$print.subgroup.name
    res$byseparator <- res$sep.subgroup
    ##
    res$TE.fixed.w <- res$TE.common.w
    res$seTE.fixed.w <- res$seTE.common.w
    res$lower.fixed.w <- res$lower.common.w
    res$upper.fixed.w <- res$upper.common.w
    res$statistic.fixed.w <- res$statistic.common.w
    res$pval.fixed.w <- res$pval.common.w
    res$zval.fixed.w <- res$zval.common.w
    res$w.fixed.w <- res$w.common.w
    ##
    res$Q.w.fixed <- res$Q.w.common
    res$pval.Q.w.fixed <- res$pval.Q.w.common
    res$Q.b.fixed <- res$Q.b.common
    res$pval.Q.b.fixed <- res$pval.Q.b.common
  }
  ##
  if (!is.null(res$three.level) && res$three.level)
    res$id <- res$cluster
  ##
  res
}


ci2meta <- function(x, ci.c = NULL, ci.r = NULL) {
  res <- x
  ##
  if (!is.null(ci.c)) {
    res$TE.common <- ci.c$TE
    res$seTE.common <- ci.c$seTE
    res$lower.common <- ci.c$lower
    res$upper.common <- ci.c$upper
    res$statistic.common <- ci.c$statistic
    res$pval.common <- ci.c$p      
    res$zval.common <- ci.c$statistic
  }
  ##
  if (!is.null(ci.r)) {
    res$TE.random <- ci.r$TE
    res$seTE.random <- ci.r$seTE
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p      
    res$zval.random <- ci.r$statistic
  }
  ##
  res
}


setNAwithin <- function(x, condition) {
  res <- x
  ##
  if (condition) {
    res$Q.w.common <- NA
    res$Q.w.random <- NA
    res$df.Q.w <- NA
    res$pval.Q.w.common <- NA
    res$pval.Q.w.random <- NA
  }
  ##
  res
}
