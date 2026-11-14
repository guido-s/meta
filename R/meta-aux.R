# Auxiliary functions
#
# Package: meta
# Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

allNA <- function(x)
  all(is.na(x))

bylevs <- function(x) {
  if (is.factor(x))
    res <- levels(factor(x))
  else
    res <- unique(x)
  res
}

subgroup_name <- function(x) {
  if (quo_is_missing(x))
    return("subgroup")
  #
  res <- as.character(quo_get_expr(x))
  #
  if (length(res) > 1 && res[1] == "$")
    res <- res[length(res)]
  #
  if (length(res) == 0 || length(res) > 1)
    res <- "subgroup"
  #
  res
}

catch <- function(x, data = NULL) {
  if (quo_is_missing(x))
    return(NULL)
  #
  eval_tidy(x, data = data)
}

catch2 <- function(x, varname, return = NULL, fromobject = FALSE) {
  if (!is.null(x[[varname]]) & fromobject)
    return(x[[varname]])
  else if (isCol(x$data, paste0(".", varname)))
    return(x$data[[paste0(".", varname)]])
  else
    return(return)
}

int2num <- function(x) {
  #
  # Convert integer to numeric
  #
  if (is.integer(x))
    res <- as.numeric(x)
  else
    res <- x
  #
  res
}

# Check for non-positive values in vector
#
npn <- function(x)
  is.na(x) | x <= 0

replaceNULL <- function(x, replace = NA, chklength = FALSE) {
  if (is.null(x) || (chklength && length(x) == 0))
    return(replace)
  #
  x
}

replaceNULL_list <- function(x, defaults) {
  for (i in names(defaults))
    x[i] <- list(replaceNULL(x[[i]], defaults[[i]]))
  #
  x
}

replaceNA <- function(x, replace = NA) {
  if (is.null(x))
    return(x)
  else
    x[is.na(x)] <- replace
  x
}

replaceVal <- function(x, old, new) {
  if (is.null(x))
    return(x)
  else
    x[x == old] <- new
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
  #
  invisible(NULL)
}

catchvar <- function(varname, x, mf) {
  res <- NULL
  error <-
    try(res <- eval(mf[[match(varname, names(mf))]],
                    x,
                    enclos = sys.frame(sys.parent())),
        silent = TRUE)
  #
  if (inherits(error, "try-error")) {
    res <- eval(mf[[match(varname, names(mf))]],
                x$data, enclos = NULL)
  }
  #
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

deprecated <- function(newvar, newmiss, args, old, warn = TRUE,
                       name = deparse(substitute(newvar))) {
  #
  if (length(args) == 0)
    return(newvar)
  #
  if (is.list(args[[1]]))
    args <- args[[1]]
  #
  additional.arguments <- names(args)
  #
  if (!is.na(charmatch(old, additional.arguments))) {
    if (!newmiss) {
      if (warn)
        warning("Deprecated argument '", old, "' ignored as ",
                "'", name, "' is also provided.",
                call. = FALSE)
      return(newvar)
    }
    else {
      if (warn)
        warning("Use argument '", name, "' instead of '",
                old, "' (deprecated).",
                call. = FALSE)
      return(args[[charmatch(old, additional.arguments)]])
    }
  }
  else
    return(newvar)
}

deprecated2 <- function(newvar, newmiss, oldvar, oldmiss, warn = TRUE,
                        oldtxt = NULL) {
  #
  new <- deparse(substitute(newvar))
  if (is.null(oldtxt))
    oldtxt <- deparse(substitute(oldvar))
  #
  if (newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & !oldmiss) {
    if (warn)
      warning("Deprecated argument '", oldtxt, "' ignored as ",
              "'", new, "' is also provided.",
              call. = FALSE)
    return(newvar)
  }
  else if (newmiss & !oldmiss) {
    if (warn)
      warning("Use argument '", new, "' instead of '",
              oldtxt, "' (deprecated).",
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
  #
  res$Q.b.random <- NA
  res$df.Q.b <- NA
  res$df.Q.b.random <- NA
  res$pval.Q.b.random <- NA
  #
  res
}

mismatch <- function(x, y, var) {
  x <- x$data[[var]]
  y <- y$data[[var]]
  #
  bothnull <- is.null(x) & is.null(y)
  #
  if (bothnull)
    return(FALSE)
  else {
    if (!is.null(x) & !is.null(y))
      return(any(x != y, na.rm = TRUE))
    else
      return(TRUE)
  }
}

addHet <- function(x, hcc, within = TRUE) {
  
  res <- x
  #
  if (within) {
    res$Q.w.random <- hcc$Q.resid
    res$df.Q.w.random <- hcc$df.Q.resid
    res$pval.Q.w.random <- hcc$pval.Q.resid
  }
  #
  res$tau2.resid <- hcc$tau2.resid
  res$lower.tau2.resid <- hcc$lower.tau2.resid
  res$upper.tau2.resid <- hcc$upper.tau2.resid
  #
  res$tau.resid <- hcc$tau.resid
  res$lower.tau.resid <- hcc$lower.tau.resid
  res$upper.tau.resid <- hcc$upper.tau.resid
  res$sign.lower.tau.resid <- hcc$sign.lower.tau.resid
  res$sign.upper.tau.resid <- hcc$sign.upper.tau.resid
  #
  res$Q.resid <- hcc$Q.resid
  res$df.Q.resid <- hcc$df.Q.resid
  res$pval.Q.resid <- hcc$pval.Q.resid
  #
  res$H.resid <- hcc$H.resid
  res$lower.H.resid <- hcc$lower.H.resid
  res$upper.H.resid <- hcc$upper.H.resid
  #
  res$I2.resid <- hcc$I2.resid
  res$lower.I2.resid <- hcc$lower.I2.resid
  res$upper.I2.resid <- hcc$upper.I2.resid
  #
  # No confidence interval for tau2.w for three-level model
  #
  if (x$tau.common & x$method.tau.ci %in% c("BJ", "J", "QP")) {
    res$lower.tau2.w[is.na(res$lower.tau2.w)] <- res$lower.tau2.resid
    res$upper.tau2.w[is.na(res$upper.tau2.w)] <- res$upper.tau2.resid
  }
  #
  res
}

backward <- function(x) {
  res <- x
  #
  res$zval <- res$statistic
  res$zval.common <- res$statistic.common
  res$zval.random <- res$statistic.random
  #
  res$w.fixed <- res$w.common
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  res$zval.fixed <- res$zval.common
  #
  res$text.fixed <- res$text.common
  res$text.w.fixed <- res$text.w.common
  #
  res$hakn <- any(res$method.random.ci == "HK")
  res$df.hakn <- res$df.hakn.ci
  res$seTE.hakn <- res$seTE.hakn.ci
  res$seTE.hakn.adhoc <- res$seTE.hakn.adhoc.ci
  #
  if (!is.null(res$subgroup)) {
    res$byvar <- res$subgroup
    res$bylab <- res$subgroup.name
    res$print.byvar <- res$print.subgroup.name
    res$byseparator <- res$sep.subgroup
    #
    res$TE.fixed.w <- res$TE.common.w
    res$seTE.fixed.w <- res$seTE.common.w
    res$lower.fixed.w <- res$lower.common.w
    res$upper.fixed.w <- res$upper.common.w
    res$statistic.fixed.w <- res$statistic.common.w
    res$pval.fixed.w <- res$pval.common.w
    res$zval.fixed.w <- res$zval.common.w
    res$w.fixed.w <- res$w.common.w
    #
    res$Q.w.fixed <- res$Q.w.common
    res$pval.Q.w.fixed <- res$pval.Q.w.common
    res$Q.b.fixed <- res$Q.b.common
    res$pval.Q.b.fixed <- res$pval.Q.b.common
  }
  #
  if (!is.null(res$three.level) && any(res$three.level))
    res$id <- res$cluster
  #
  res
}

drop_backward <- function(x) {
  res <- x
  #
  res$zval <- NULL
  res$zval.common <- NULL
  res$zval.random <- NULL
  #
  res$w.fixed <- NULL
  res$TE.fixed <- NULL
  res$seTE.fixed <- NULL
  res$lower.fixed <- NULL
  res$upper.fixed <- NULL
  res$statistic.fixed <- NULL
  res$pval.fixed <- NULL
  res$zval.fixed <- NULL
  #
  res$text.fixed <- NULL
  res$text.w.fixed <- NULL
  #
  res$hakn <- NULL
  res$df.hakn <- NULL
  res$seTE.hakn <- NULL
  res$seTE.hakn.adhoc <- NULL
  #
  res$byvar <- NULL
  res$bylab <- NULL
  res$print.byvar <- NULL
  res$byseparator <- NULL
  #
  res$TE.fixed.w <- NULL
  res$seTE.fixed.w <- NULL
  res$lower.fixed.w <- NULL
  res$upper.fixed.w <- NULL
  res$statistic.fixed.w <- NULL
  res$pval.fixed.w <- NULL
  res$zval.fixed.w <- NULL
  res$w.fixed.w <- NULL
  #
  res$Q.w.fixed <- NULL
  res$pval.Q.w.fixed <- NULL
  res$Q.b.fixed <- NULL
  res$pval.Q.b.fixed <- NULL
  #
  res$id <- NULL
  #
  res
}

ci2meta <- function(x, ci.c = NULL, ci.r = NULL) {
  res <- x
  #
  if (!is.null(ci.c)) {
    res$TE.common <- ci.c$TE
    res$seTE.common <- ci.c$seTE
    res$lower.common <- ci.c$lower
    res$upper.common <- ci.c$upper
    res$statistic.common <- ci.c$statistic
    res$pval.common <- ci.c$p      
    res$zval.common <- ci.c$statistic
  }
  #
  if (!is.null(ci.r)) {
    res$TE.random <- ci.r$TE
    res$seTE.random <- ci.r$seTE
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p      
    res$zval.random <- ci.r$statistic
  }
  #
  res
}

setNAwithin <- function(x, condition) {
  res <- x
  #
  if (condition) {
    res$Q.w.common <- NA
    res$Q.w.random <- NA
    res$df.Q.w <- NA
    res$pval.Q.w.common <- NA
    res$pval.Q.w.random <- NA
  }
  #
  res
}

drop_allNAs <- function(x, keep = character()) {
  all_NA <- vapply(x, function(column)
    if (is.list(column))
      all(vapply(column, function(value) !is.null(value) && allNA(value), logical(1)))
    else
      allNA(column),
    logical(1))
  all_NA[colnames(x) %in% keep] <- FALSE
  #
  res <- x[, !all_NA, drop = FALSE]
  colnames(res) <- colnames(x)[!all_NA]
  res
}

extrVar <- function(x, name)
  x[[name]]

# Function only used with MLM, GLMM, LRP, or USW
#
calcPI <- function(x) {
  
  res <- x
  #
  k <- x$k
  TE.random <- x$TE.random[1]
  seTE.random <- x$seTE.random[1]
  method.predict <- x$method.predict
  level.predict <- x$level.predict
  #
  res$df.predict <- res$upper.predict <- res$lower.predict <-
    res$seTE.predict <- NA
  #
  tau2.calc <- if (is.na(sum(res$tau2))) 0 else sum(res$tau2)
  seTE.predict <- sqrt(seTE.random^2 + tau2.calc)
  #
  df.predict <- set_df_predict(method.predict, k)
  #
  ci.p <- ci(TE.random, seTE.predict, level.predict, df = df.predict)
  #
  res$seTE.predict <- rep_len(seTE.predict, length(df.predict))
  res$lower.predict <- ci.p$lower
  res$upper.predict <- ci.p$upper
  res$df.predict <- ci.p$df
  #
  if (length(method.predict) > 1)
    names(res$seTE.predict) <-
      names(res$lower.predict) <- names(res$upper.predict) <-
      names(res$df.predict) <- method.predict
  #
  res
}

runCRuni <- function(x, method.tau, method.random.ci, level,
                     warn = TRUE, ...) {

  n.methci <- length(method.random.ci)
  #
  list.b <- list(level = 100 * level, method = method.tau, ...)
  #
  res <- list.r <- vector("list", n.methci)
  #
  for (i in seq_len(n.methci))
    list.r[[i]] <- c(x, list.b)
  #
  for (i in seq_len(n.methci)) {
    res[[i]] <- runNN(rma.uni, list.r[[i]], warn = warn)
    #
    if (method.random.ci[i] %in% c("CR0", "CR1", "CR2")) {
      res[[i]] <-
        suppressPackageStartupMessages(
          robust(res[[i]], cluster = x$data$studlab,
                 adjust = method.random.ci[i] != "CR0",
                 clubSandwich = method.random.ci[i] == "CR2"))
    }
  }
  #
  res
}

extrCRuni <- function(x, k, len, sel,
                      method.random.ci, method.predict,
                      level.ma, level.predict) {

  #
  # Random effects model(s)
  #
  res <- list(tau2 = NA, tau = NA, se.tau2 = NA)
  #
  res$k <- k
  res$method.predict <- method.predict
  res$level.predict <- level.predict
  #
  if (k > 1) {
    res$tau2 <- x[[1]]$tau2
    res$tau <- sqrt(res$tau2)
  }
  #
  res$TE.random <- as.numeric(coef(x[[1]]))
  res$seTE.random <- sapply(x, extrVar, "se")
  #
  res$lower.random <- sapply(x, extrVar, "ci.lb")
  res$upper.random <- sapply(x, extrVar, "ci.ub")
  #
  res$statistic.random <- sapply(x, extrVar, "zval")
  res$df.random <- res$df.hakn <- sapply(x, extrVar, "ddf")
  res$pval.random <- sapply(x, extrVar, "pval")
  #
  res$w.random <- rep_len(NA, len)
  res$w.random[sel] <- weights(x[[1]])
  res$w.random[is.na(res$w.random)] <- 0
  #
  res$zval.random <- res$statistic.random
  #
  if (length(method.random.ci) > 1)
    names(res$seTE.random) <-
    names(res$lower.random) <- names(res$upper.random) <-
    names(res$statistic.random) <- names(res$pval.random) <-
    names(res$df.random) <- names(res$df.hakn) <-
    names(res$zval.random) <- method.random.ci
  #
  # Prediction interval
  #
  res <- calcPI(res)
  #
  res$k <- res$method.predict <- res$level.predict <- NULL
  #
  res
}

runMLM <- function(x, method.tau, method.random.ci, level,
                   warn = TRUE, ...) {
  
  n.methci <- length(method.random.ci)
  #
  list.b <- list(level = 100 * level, method = method.tau,
                 random = as.call(~ 1 | cluster / idx), ...)
  #
  res <- list.r <- vector("list", n.methci)
  #
  for (i in seq_len(n.methci))
    list.r[[i]] <- c(x, list.b,
                     test = ifelse(method.random.ci[i] == "HK", "t", "z"))
  #
  for (i in seq_len(n.methci)) {
    res[[i]] <- runNN(rma.mv, list.r[[i]], warn = warn)
    #
    if (method.random.ci[i] %in% c("CR0", "CR1", "CR2")) {
      res[[i]] <-
        suppressPackageStartupMessages(
          robust(res[[i]], cluster = x$data$cluster,
                 adjust = method.random.ci[i] != "CR0",
                 clubSandwich = method.random.ci[i] == "CR2"))
    }
  }
  #
  res
}

extrMLM <- function(x, k, len, sel,
                    method.random.ci, method.predict,
                    level.ma, level.predict, null.effect) {
  
  #
  # Random effects model(s)
  #
  res <- list(tau2 = NA, tau = NA, se.tau2 = NA)
  #
  res$k <- k
  res$method.predict <- method.predict
  res$level.predict <- level.predict
  #
  if (k > 1) {
    res$tau2 <- x[[1]]$sigma2
    res$tau <- sqrt(res$tau2)
  }
  #
  res$TE.random <- as.numeric(coef(x[[1]]))
  res$seTE.random <- sapply(x, extrVar, "se")
  #
  res$lower.random <- sapply(x, extrVar, "ci.lb")
  res$upper.random <- sapply(x, extrVar, "ci.ub")
  #
  res$statistic.random <- sapply(x, extrVar, "zval")
  res$df.random <- res$df.hakn <- sapply(x, extrVar, "ddf")
  res$pval.random <- sapply(x, extrVar, "pval")
  #
  res$w.random <- rep_len(NA, len)
  res$w.random[sel] <- weights(x[[1]], type = "rowsum")
  res$w.random[is.na(res$w.random)] <- 0
  #
  res$zval.random <- res$statistic.random
  #
  if (length(method.random.ci) > 1)
    names(res$seTE.random) <-
      names(res$lower.random) <- names(res$upper.random) <-
      names(res$statistic.random) <- names(res$pval.random) <-
      names(res$df.random) <- names(res$df.hakn) <-
      names(res$zval.random) <- method.random.ci
  #
  # Prediction interval
  #
  res <- calcPI(res)
  #
  res$k <- res$method.predict <- res$level.predict <- NULL
  #
  res
}

runUSW <- function(x, method.tau, method.random.ci, level,
                   warn = TRUE, ...) {
  
  if (method.tau %in% c("FE", "CE")) {
    n.methci <- 1
    #
    list.b <- list(level = 100 * level, method = method.tau, ...)
    #
    res <- list.c <- vector("list", n.methci)
    #
    for (i in seq_len(n.methci)) {
      list.c[[i]] <- c(x, list.b, test = "z")
      #
      res[[i]] <- runNN(rma.uni, list.c[[i]], warn = warn)
    }
  }
  else {
    n.methci <- length(method.random.ci)
    #
    list.b <- list(level = 100 * level, method = method.tau, ...)
    #
    res <- list.r <- vector("list", n.methci)
    #
    for (i in seq_len(n.methci)) {
      list.r[[i]] <- c(x, list.b,
                       test = ifelse(method.random.ci[i] == "HK", "t", "z"))
      #  
      res[[i]] <- runNN(rma.uni, list.r[[i]], warn = warn)
    }
  }
  #
  res
}

extrUSW <- function(x, method.tau, null.effect,
                    k, len, sel,
                    method.random.ci, method.predict,
                    level.ma, level.predict) {
  if (method.tau %in% c("FE", "CE")) {
    res <- list(
      w.common = x[[1]]$weights,
      TE.common = as.numeric(x[[1]]$b),
      seTE.common = as.numeric(x[[1]]$se),
      lower.common = x[[1]]$ci.lb,
      upper.common = x[[1]]$ci.ub)
    #
    ci.c <- ci(res$TE.common - null.effect, res$seTE.common,
         null.effect = null.effect)
    #
    res$statistic.common <- ci.c$statistic
    res$pval.common <- ci.c$p
    res$zval.common <- res$statistic.common
  }
  else {
    #
    # Random effects model(s)
    #
    res <- list(tau2 = NA, tau = NA, se.tau2 = NA)
    #
    res$k <- k
    res$method.predict <- method.predict
    res$level.predict <- level.predict
    #
    if (k > 1) {
      res$tau2 <- x[[1]]$tau2
      res$tau <- sqrt(res$tau2)
    }
    #
    res$w.random <- x[[1]]$weights
    #
    res$TE.random   <- as.numeric(x[[1]]$b)
    res$seTE.random <- sapply(x, extrVar, "se")
    #
    res$lower.random <- sapply(x, extrVar, "ci.lb")
    res$upper.random <- sapply(x, extrVar, "ci.ub")
    #
    ci.r <- ci(res$TE.random, res$seTE.random,
               level = level.ma,
               null.effect = null.effect,
               df = ifelse(method.random.ci == "HK", k - 1, Inf))
    #
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p
    res$zval.random <- res$statistic.random
    #
    if (length(method.random.ci) > 1)
      names(res$seTE.random) <-
      names(res$lower.random) <- names(res$upper.random) <-
      names(res$statistic.random) <- names(res$pval.random) <-
      names(res$zval.random) <- method.random.ci
    #
    # Prediction interval
    #
    res <- calcPI(res)
  }
  #
  res$k <- res$method.predict <- res$level.predict <- NULL
  #
  res
}

runGLMM <- function(x, method.tau, method.random.ci, level,
                    use.random = TRUE, warn = TRUE, ...) {
  
  n.methci <- length(method.random.ci)
  #
  list.b <- list(level = 100 * level, ...)
  list.c <- c(x, list.b, method = "FE", test = "z")
  #
  res.r <- list.r <- vector("list", n.methci)
  #
  for (i in seq_len(n.methci))
    list.r[[i]] <- c(x, list.b, method = method.tau,
                     test = ifelse(method.random.ci[i] == "HK", "t", "z"))
  #
  res.c <- runNN(rma.glmm, list.c, warn = FALSE)
  #
  for (i in seq_len(n.methci)) {
    if (use.random) {
      res.r[[i]] <-
        try(runNN(rma.glmm, list.r[[i]], warn = warn), silent = TRUE)
      #
      if ("try-error" %in% class(res.r[[i]])) {
        if (grepl(paste0("Number of parameters to be estimated is ",
                         "larger than the number of observations"),
                  res.r[[i]]) |
            grepl("Cannot fit ML model.", res.r[[i]])) {
          warning(paste(res.r[[i]],
                        " Using result of common effect model as fallback."))
          #
          res.r[[i]] <- res.c
        }
        else
          stop(res.r[[i]])
      }
    }
    else {
      # Fallback to common effect model due to small number of
      # studies or zero or all events
      res.r[[i]] <- res.c
    }
  }
  
  res <- list(glmm.common = res.c, glmm.random = res.r)
  #
  res
}

addGLMM <- function(x, glmm, method.I2, null.effect = x$null.effect) {
  
  res <- x
  #
  k <- x$k
  len <- length(x$TE)
  level.ma <- x$level.ma
  level.hetstat <- x$level.hetstat
  level.predict <- x$level.predict
  method.random.ci <- x$method.random.ci
  method.predict <- x$method.predict
  null.effect <- replaceNULL(null.effect, 0)
  #
  # Common effect model
  #
  glmm.common <- glmm$glmm.common
  #
  res$TE.common   <- as.numeric(glmm.common$b)
  res$seTE.common <- as.numeric(glmm.common$se)
  res$w.common <- rep(NA, len)
  res$lower.common <- glmm.common$ci.lb
  res$upper.common <- glmm.common$ci.ub
  #
  ci.c <- ci(res$TE.common, res$seTE.common, level = level.ma,
             null.effect = null.effect)
  res$statistic.common <- res$zval.common <- ci.c$statistic
  res$pval.common <- ci.c$p
  #
  # Random effects model(s)
  #
  glmm.random <- glmm$glmm.random
  #
  res$TE.random   <- as.numeric(glmm.random[[1]]$b)
  res$seTE.random <- sapply(glmm.random, extrVar, "se")
  res$w.random <- rep(NA, len)
  res$lower.random <- sapply(glmm.random, extrVar, "ci.lb")
  res$upper.random <- sapply(glmm.random, extrVar, "ci.ub")
  #
  ci.r <-
    ci(res$TE.random, res$seTE.random,
       level = level.ma,
       null.effect = null.effect,
       df = ifelse(method.random.ci == "HK", k - 1, Inf))
  res$statistic.random <- res$zval.random <- ci.r$statistic
  res$pval.random <- ci.r$p
  #
  if (length(method.random.ci) > 1)
    names(res$seTE.random) <-
      names(res$lower.random) <- names(res$upper.random) <-
      names(res$statistic.random) <- names(res$pval.random) <-
      names(res$zval.random) <- method.random.ci
  #
  res$se.tau2 <- res$tau <- res$tau2 <- NA
  #
  if (k > 1) {
    res$tau2 <- glmm.random[[1]]$tau2
    res$tau <- sqrt(res$tau2)
    res$se.tau2 <- glmm.random[[1]]$se.tau2
  }
  #
  # Prediction interval
  #
  res <- calcPI(res)
  #
  # Heterogeneity measures
  #
  gr1 <- glmm.random[[1]]
  #
  Q <- if (gr1$k > 1) gr1$QE.Wld else 0
  df.Q <- gr1$QE.df
  pval.Q <- pvalQ(Q, df.Q)
  #
  Q.LRT <- if (gr1$k > 1) gr1$QE.LRT else 0
  pval.Q.LRT <- pvalQ(Q.LRT, df.Q)
  #
  res$Q <- c(Q, Q.LRT)
  res$df.Q <- c(df.Q, df.Q)
  res$pval.Q <- c(pval.Q, pval.Q.LRT)
  names(res$Q) <- c("Wald", "LRT")
  #
  res$upper.tau <- res$lower.tau <- res$upper.tau2 <- res$lower.tau2 <- NA
  #
  res$sign.upper.tau <- res$sign.lower.tau <- res$method.tau.ci <- ""
  #
  if (method.I2 == "Q") {
    H <- calcH(Q, df.Q, level.hetstat)
    res$H <- H$TE
    res$lower.H <- H$lower
    res$upper.H <- H$upper
    #
    I2 <- isquared(Q, df.Q, level.hetstat)
    res$I2 <- I2$TE
    res$lower.I2 <- I2$lower
    res$upper.I2 <- I2$upper
  }
  else {
    res$H <- sqrt(gr1$H2)
    res$lower.H <- NA
    res$upper.H <- NA
    #
    res$I2 <- gr1$I2 / 100
    res$lower.I2 <- NA
    res$upper.I2 <- NA
  }
  #
  res$upper.Rb <- res$lower.Rb <- res$Rb <- NA
  #
  res$.glmm.common <- glmm.common
  #
  if (length(glmm.random) == 1)
    res$.glmm.random <- glmm.random[[1]]
  else
    res$.glmm.random <- glmm.random
  #
  res$version.metafor <- packageDescription("metafor")$Version
  #
  res
}

hccGLMM <- function(x, glmm, method.I2) {
  Q.r <- glmm$QE.Wld
  df.Q.r <- glmm$k - glmm$p
  #
  if (method.I2 == "Q") {
    H.r  <- calcH(Q.r, df.Q.r, x$level.hetstat)
    I2.r <- isquared(Q.r, df.Q.r, x$level.hetstat)
  }
  else {
    H.r  <- list(TE = sqrt(glmm$H2), lower = NA, upper = NA)
    I2.r <- list(TE = glmm$I2 / 100, lower = NA, upper = NA)
  }
  #
  list(tau2.resid = glmm$tau2,
       lower.tau2.resid = NA,
       upper.tau2.resid = NA,
       #
       tau.resid = sqrt(glmm$tau2),
       lower.tau.resid = NA,
       upper.tau.resid = NA,
       sign.lower.tau.resid = "",
       sign.upper.tau.resid = "",
       #
       Q.resid = Q.r,
       df.Q.resid = df.Q.r,
       pval.Q.resid = pvalQ(Q.r, df.Q.r),
       #
       H.resid = H.r$TE,
       lower.H.resid = H.r$lower,
       upper.H.resid = H.r$upper,
       #
       I2.resid = I2.r$TE,
       lower.I2.resid = I2.r$lower,
       upper.I2.resid = I2.r$upper
       )
}

runLRP <- function(event1, n1, event2, n2, warn = TRUE, ...) {
  
  dat.bin <- data.frame(event1, n1, event2, n2)
  #
  long.bin <-
    longarm("B", "A", event1 = event1, n1 = n1, event2 = event2, n2 = n2,
            data = dat.bin)
  #
  text.formula <- "cbind(events, nonevents) ~ as.factor(treat)"
  #
  if (length(unique(long.bin$studlab)) > 1)
    text.formula <- paste(text.formula, "+ as.factor(studlab)")
  #
  res.lrp <- glm(as.formula(text.formula), data = long.bin,
                 family = binomial(link = "logit"),
                 method = brglm2::brglmFit, type = "MPL_Jeffreys",
                 ...)
  #
  phi <- 1
  if (length(unique(long.bin$studlab)) > 1)
    phi <- phi(res.lrp)
  #
  sel.trt <- grepl("treat", names(coef(res.lrp)))
  #
  TE.common   <- as.numeric(coef(res.lrp)[sel.trt])
  seTE.common <- as.numeric(sqrt(diag(vcov(res.lrp)))[sel.trt])
  #
  res <- list(TE.common = TE.common, seTE.common = seTE.common,
              TE.random = TE.common, seTE.random = seTE.common * phi,
              phi = phi, fit = res.lrp)
  #
  res
}

# Estimate the heterogeneity parameter phi using the
# modified version of Pearson's statistic.
#
phi <- function(x) {
  # Extract number of trials
  n.trials <- x$prior.weights
  #
  if (identical(unique(n.trials), 1))
    stop("The number of successes must be summarized for valid computation of ",
         "c-hat.")
  
  # Pearson chi-square
  chisq <- sum(residuals(x, type = "pearson")^2)
  
  # Extract raw residuals
  raw.res <- residuals(x, type = "response")
  
  # Extract fitted values
  fit.vals <- fitted(x)
  
  # Estimate s.bar
  s.bar <- mean((1 - 2 * fit.vals) / ((n.trials * fit.vals) * (1 - fit.vals)))
  
  # Calculate estimate based on Fletcher estimator
  phi <- (chisq / x$df.residual) / (1 + s.bar)
  
  # Set phi = 1 if phi < 1 to remain consistent with common effect model
  phi <- max(phi, 1)
  
  phi
}


calcPercent <- function(x)
  100 * x / sum(x, na.rm = TRUE)

cond <- function(x, only.finite = TRUE, digits = 2, big.mark = "") {
  if (is.null(x))
    return(x)
  #
  if (is.character(x))
    return(x)
  #
  if (only.finite)
    x <- x[is.finite(x)]
  #
  paste(formatN(unique(round(x, digits = digits)), digits = digits,
                big.mark = big.mark), collapse = ", ")
}

setNA_ifnot <- function(x, y, unequal) {
  if (is.null(x) | is.null(y))
    return(x)
  else
    return(ifelse(y != unequal, NA, x))
}

list2mat <- function(x, y = NULL) {
  if (is.list(x)) {
    lnames <- names(x)
    #
    if (!is.null(y) &&
        (!is.list(y) || length(y) != length(x)))
      stop("Reference blocks must correspond to the result blocks.",
           call. = FALSE)
    #
    for (i in seq_along(x)) {
      value <- x[[i]]
      #
      if (is.null(value))
        next
      #
      if (!is.matrix(value))
        value <- matrix(value, ncol = 1,
                        dimnames = list(names(value), NULL))
      #
      if (!is.null(y)) {
        reference <- y[[i]]
        #
        if (!is.matrix(reference))
          reference <- matrix(
            reference, ncol = 1,
            dimnames = list(names(reference), NULL))
        #
        if (nrow(value) != nrow(reference))
          stop("Result and reference blocks have different numbers of rows.",
               call. = FALSE)
        #
        # A shared subgroup estimate is repeated across its own CI methods.
        #
        if (ncol(value) == 1L && ncol(reference) > 1L)
          value <- matrix(rep(value, ncol(reference)),
                          nrow = nrow(reference),
                          ncol = ncol(reference))
        #
        if (!identical(dim(value), dim(reference)))
          stop("Result and reference blocks have incompatible dimensions.",
               call. = FALSE)
        #
        dimnames(value) <- dimnames(reference)
      }
      #
      if (!is.null(lnames) &&
          !is.na(lnames[i]) && nzchar(lnames[i])) {
        cnames <- colnames(value)
        #
        if (is.null(cnames)) {
          if (ncol(value) == 1L)
            cnames <- lnames[i]
          else
            cnames <- paste(lnames[i], seq_len(ncol(value)), sep = "-")
        }
        else
          cnames <- gsub(
            "-classic$", "",
            paste(lnames[i], cnames, sep = "-"))
        #
        colnames(value) <- cnames
      }
      #
      x[i] <- list(value)
    }
    #
    do.call("cbind", x)
  }
  else if (is.vector(x))
    matrix(x, ncol = 1)
  else
    x
}

list2vec <- function(x) {
  if (is.list(x))
    return(do.call("c", x))
  else
    return(x)
}

em2sm <- function(x, type = NULL) {
  x[x == "Odds Ratio"] <- "OR"
  x[x == "Odds Ratio (Non-event)"] <- "OR"
  x[x == "Peto Odds Ratio"] <- "OR"
  x[x == "Risk Ratio"] <- "RR"
  x[x == "Risk Difference"] <- "RD"
  x[x == "Mean Difference"] <- "MD"
  x[x == "Standardized Mean Difference"] <- "SMD"
  x[x == "Std. Mean Difference"] <- "SMD"
  x[x == "Hazard Ratio"] <- "HR"
  #
  if (!is.null(type))
    x <- ifelse(is.na(x) & type == "I", "", x)
  #
  x
}

sm2meth <- function(x) {
  x[x == "IV"] <- "Inverse"
  x[x == "EXP_O_E_VAR"] <- "Peto"
  #
  x
}

setsv <- function(x) {
  if (is.null(x))
    res <- "desirable"
  else {
    res <- setchar(x, c("good", "bad"), stop.at.error = FALSE)
    #
    if (!is.null(res))
      res <- switch(res, good = "desirable", bad = "undesirable")
    else
      res <- x
  }
  #
  setchar(res, c("desirable", "undesirable"))
}

expand <- function(x, y) {
  if (length(x) == 1 & length(y) > 1)
    return(rep_len(x, length(y)))
  else
    return(x)  
}

expandvar <- function(x, n, length = NULL) {
  res <- x
  if (!is.null(length))
    lenOK <- length(x) == length
  else
    lenOK <- TRUE
  #
  if (lenOK & length(x) != n)
    res <- rep(x, rep_len(n, length(x)))
  #
  res
}

warn_ignore_input <- function(x, cond = TRUE, text = "") {
  if (cond)
    warning("Argument '", deparse(substitute(x)), "' ignored",
            if (text != "") " ", text, ".",
            call. = FALSE)
  #
  invisible(NULL)
}

warn_set_input <- function(x, cond = TRUE, set, text = "") {
  if (cond)
    warning("Argument '", deparse(substitute(x)), " = ", set, "'",
            if (text != "") " ", text, ".",
            call. = FALSE)
  #
  invisible(NULL)
}

drop_from_dots <- function(x, drop = NULL, replace = list(),
                           deprecated = FALSE, text.deprecated = NULL) {
  for (old in drop) {
    if (!is.null(x[[old]])) {
      new <- replace[[old]]
      if (!is.null(new) && length(new) == 1 && nzchar(new))
        warning("Argument '", old, "' ignored; please use argument '",
                new, "' instead.",
                call. = FALSE)
      else if (deprecated)
        warning("Argument '", old, "' is deprecated and ignored.",
                if (!is.null(text.deprecated))
                  paste0(" ", text.deprecated),
                call. = FALSE)
      else
        warning("Argument '", old, "' ignored as it is set internally.",
                call. = FALSE)
      #
      x[[old]] <- NULL
    }
  }
  #
  x
}

warn_ignore_dots <- function(x, func = "") {
  if (length(x) > 0) {
    nam.args <- names(x)
    #
    if (func != "")
      txt <- paste0(func, "(): ")
    else
      txt <- ""
    #
    if (length(nam.args) == 1)
      txt <- paste0(txt,
                    "The argument '",
                    nam.args,
                    "' provided in '...' is ignored.")
    else
      txt <- paste0(txt,
                    "The following arguments provided in '...' are ignored: ",
                    paste0(nam.args, collapse = ", "))
    #
    warning(txt, call. = FALSE)
  }
  #
  invisible(NULL)
}

merge_legacy <- function(x, args, default,
                         prefix = deparse(substitute(x)),
                         inherit = NULL, warn = FALSE) {
  nam.x <- deparse(substitute(x))
  force(prefix)
  #
  if (is.null(x))
    x <- list()
  #
  if (!is.list(x) || (length(x) > 0 && is.null(names(x))))
    stop("Argument '", nam.x, "' must be a named list or NULL.", call. = FALSE)
  #
  if (!is.list(default) || is.null(names(default)))
    stop("Bug: argument 'default' must be a named list.", call. = FALSE)
  #
  supplied <- names(x)
  #
  res <- modifyList(default, x, keep.null = TRUE)
  #
  miss <- setNames(as.list(!(names(default) %in% supplied)), names(default))
  #
  null <-
    setNames(lapply(names(miss), function(i) is.null(res[[i]])), names(miss))
  #
  old.names <- intersect(names(args), paste0(prefix, ".", names(default)))
  #
  for (old in old.names) {
    name <- sub(paste0("^", prefix, "\\."), "", old)
    #
    if (name %in% supplied)
      stop("Argument '", old, "' and element '",
           prefix, "[['", name, "']] cannot both be specified.",
           call. = FALSE)
    #
    if (warn)
      warning("Use element '", prefix,
              "[['", name, "']]' instead of argument '", old,
              "' (deprecated).",
              call. = FALSE)
    #
    res[name] <- list(args[[old]])
    miss[name] <- FALSE
    null[[name]] <- is.null(res[[name]])
    args[[old]] <- NULL
  }
  #
  if (!is.null(inherit)) {
    if (is.null(names(inherit)))
      stop("Argument 'inherit' must be a named list.",
            call. = FALSE)
    #
    inherit.eligible <- Map(function(m, n) m || n, miss, null)
    #
    for (pass in seq_len(length(inherit) + 1L)) {
      changed <- FALSE
      #
      for (child in names(inherit)) {
        parent <- inherit[[child]]
        #
        if (!(child %in% names(miss)) || !(parent %in% names(miss)))
          stop("Unknown child or parent in argument 'inherit'.",
               call. = FALSE)
        #
        if (inherit.eligible[[child]] && !is.null(res[[parent]])) {
          changed <- changed ||
            !identical(res[[child]], res[[parent]]) ||
            null[[child]]
          res[[child]] <- res[[parent]]
          null[[child]] <- FALSE
        }
      }
      #
      if (!changed)
        break
    }
  }
  #
  list(res = res, args = args, miss = miss, null = null)
}

as_named_list <- function(x) {
  name <- deparse(substitute(x))
  #
  if (is.null(x))
    return(x)
  #
  if (is.list(x)) {
    if (is.null(names(x)))
      stop("Argument '", name, "' must be a named list.", call. = FALSE)
    #
    return(x)
  }
  #
  setNames(list(x), name)
}

collapse_choices <- function(x) {
  if (is.numeric(x) && length(x) > 2 && all(diff(x) == 1)) {
    paste("from", x[1], "to", x[length(x)])
  }
  else if (length(x) == 1) {
    x
  }
  else if (length(x) == 2) {
    paste(x, collapse = " or ")
  }
  else {
    paste0(
      paste(x[-length(x)], collapse = ", "), ", or ", x[length(x)])
  }
}

format_seq <- function(x, exclude = NULL, digits = NULL) {
  if (!is.null(exclude))
    x <- x[x != exclude]
  #
  if (!is.null(digits) && is.numeric(x))
    x <- round(x, digits)
  #
  x <- sort(unique(x))
  #
  if (length(x) == 1)
    as.character(x)
  else if (!is.null(exclude) && length(x) == 0)
    as.character(exclude)
  else
    paste0("{", paste(x, collapse = ", "), "}")
}

sum_or_na <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)))
    return(NA)
  #
  sum(x, na.rm = TRUE)
}

makeunique <- function(x, val = NA, orig = FALSE) {
  res <- unique(x)
  #
  if (length(res) != 1) {
    if (orig)
      return(res)
    else
      res <- val
  }
  #
  res
}

addNAs2var <- function(x, n) {
  if (n == 1)
    res <- x
  else {
    res <- vector("numeric", 0)
    #
    for (i in seq_along(x))
      res <- c(res, c(x[i], rep(NA, n - 1)))
  }
  #
  res
}

fs_defaults <- function(fontsize) {
  list(
    main = gs("fs.main"),
    heading = fontsize,
    common = gs("fs.common"),
    random = gs("fs.random"),
    predict = gs("fs.predict"),
    common.labels = gs("fs.common.labels"),
    random.labels = gs("fs.random.labels"),
    predict.labels = gs("fs.predict.labels"),
    study = fontsize,
    study.labels = fontsize,
    hetstat = gs("fs.hetstat"),
    test.overall = gs("fs.test.overall"),
    test.subgroup = gs("fs.test.subgroup"),
    test.effect.subgroup = gs("fs.test.effect.subgroup"),
    addline = gs("fs.addline"),
    axis = fontsize - 1,
    smlab = fontsize,
    xlab = fontsize,
    lr = fontsize,
    rob = gs("fs.rob"),
    rob.symbols = gs("fs.rob.symbols"),
    details = fontsize - 1
  )
}

ff_defaults <- function() {
  list(
    main = gs("ff.main"),
    heading = "bold",
    common = gs("ff.common"),
    random = gs("ff.random"),
    predict = gs("ff.predict"),
    common.labels = gs("ff.common.labels"),
    random.labels = gs("ff.random.labels"),
    predict.labels = gs("ff.predict.labels"),
    study = "plain",
    study.labels = "plain",
    hetstat = gs("ff.hetstat"),
    test.overall = gs("ff.test.overall"),
    test.subgroup = gs("ff.test.subgroup"),
    test.effect.subgroup = gs("ff.test.effect.subgroup"),
    addline = gs("ff.addline"),
    axis = gs("ff.axis"),
    smlab = gs("ff.smlab"),
    xlab = gs("ff.xlab"),
    lr = gs("ff.lr"),
    rob = "plain",
    rob.symbols = "bold",
    details = "plain"
  )
}

col_defaults <- function(x, class = "meta") {
  if (class == "meta") {
    col.inside <- gs("col.inside")
    col.diamond <- gs("col.diamond")
    col.diamond.lines <- gs("col.diamond.lines")
    col.predict <- gs("col.predict")
    col.predict.lines <- gs("col.predict.lines")
    col.subgroup <- gs("col.subgroup")
    col.lines <- gs("col.lines")
    col.circle <- gs("col.circle")
    #
    res <- list(
      main = gs("col.main"),
      common = gs("col.common"),
      random = gs("col.random"),
      label = gs("col.label"),
      study = gs("col.study"),
      study.text = gs("col.study.text"),
      square = gs("col.square"),
      square.lines = gs("col.square.lines"),
      circle = col.circle,
      circle.lines = col.circle,
      inside = col.inside,
      inside.common = col.inside,
      inside.random = col.inside,
      diamond = col.diamond,
      diamond.common = col.diamond,
      diamond.common.subgroup = col.diamond,
      diamond.random = col.diamond,
      diamond.random.subgroup = col.diamond,
      diamond.lines = col.diamond.lines,
      diamond.lines.common = col.diamond.lines,
      diamond.lines.common.subgroup = col.diamond.lines,
      diamond.lines.random = col.diamond.lines,
      diamond.lines.random.subgroup = col.diamond.lines,
      predict = col.predict,
      predict.subgroup = col.predict,
      predict.lines = col.predict.lines,
      predict.lines.subgroup = col.predict.lines,
      subgroup = col.subgroup,
      label.left = x$col.label.left,
      label.right = x$col.label.right,
      lines = col.lines,
      header.line = col.lines,
      jama.line = col.subgroup,
      cid = gs("col.cid")
    )
  }
  else if (class == "metacum" || class == "metainf") {
    res <- list(
      col = NULL,
      bg = NULL,
      border = NULL,
      bg.predict = NULL,
      border.predict = NULL)
  }
  else if (class == "metabind") {
    res <- list(
      square = NULL,
      square.lines = NULL,
      circle = NULL,
      circle.lines = NULL,
      diamond = NULL,
      diamond.common = NULL,
      diamond.random = NULL,
      diamond.lines = NULL,
      diamond.lines.common = NULL,
      diamond.lines.random = NULL,
      predict = NULL,
      predict.lines = NULL)
  }
  #
  res
}

calcwidth_defaults <- function(common, random, overall, subgroup) {
  calcwidth.pooled <-
    (common | random) & (overall | !is.null(subgroup))
  #
  list(
    pooled = calcwidth.pooled,
    common = calcwidth.pooled,
    random = calcwidth.pooled,
    predict = gs("calcwidth.predict"),
    hetstat = gs("calcwidth.hetstat"),
    tests = gs("calcwidth.tests"),
    subgroup = gs("calcwidth.subgroup"),
    addline = gs("calcwidth.addline"),
    details = gs("calcwidth.details")
  )
}

print_defaults <- function(x, class = "meta", func = "forest") {
  if (class %in% c("meta", "summary.meta") && func == "print") {
    res <- list(
      subgroup.name = x$print.subgroup.name,
      CMH = x$print.CMH,
      tau2 = gs("print.tau2"),
      tau2.ci = gs("print.tau2.ci"),
      tau = gs("print.tau"),
      tau.ci = gs("print.tau.ci"),
      Q = gs("print.Q"),
      I2 = gs("print.I2"),
      I2.ci = gs("print.I2.ci"),
      H = gs("print.H"),
      Rb = gs("print.Rb"))
  }
  else if (class == "meta") {
    res <- list(
      subgroup.labels = TRUE,
      subgroup.name = x$print.subgroup.name,
      I2 = gs("forest.I2"),
      I2.ci = gs("forest.I2.ci"),
      tau2 = gs("forest.tau2"),
      tau2.ci = gs("forest.tau2.ci"),
      tau = gs("forest.tau"),
      tau.ci = gs("forest.tau.ci"),
      Q = gs("forest.Q"),
      pval.Q = gs("forest.pval.Q"),
      Rb = gs("forest.Rb"),
      Rb.ci = gs("forest.Rb.ci"),
      stat = gs("forest.stat"),
      tau2.ci.subgroup = gs("print.tau2.ci.subgroup"),
      tau.ci.subgroup = gs("print.tau.ci.subgroup"),
      Q.subgroup = gs("forest.Q.subgroup"))
  }
  else if (class == "metacum" || class == "metainf") {
    if (func == "print") {
      res <- list(
        stat = FALSE,
        tau2 = TRUE,
        tau2.ci = FALSE,
        tau = TRUE,
        tau.ci = FALSE,
        I2 = TRUE,
        I2.ci = FALSE,
        prob = TRUE)
    }
    else
      res <- list(
        stat = FALSE,
        pval = TRUE,
        I2 = TRUE,
        I2.ci = FALSE,
        tau2 = TRUE,
        tau2.ci = FALSE,
        tau = TRUE,
        tau.ci = FALSE,
        Q = FALSE,
        pval.Q = FALSE,
        prob = FALSE)
  }
  else if (class == "metabind") {
    func <- setchar(func, c("forest", "print"), length = 1)
    #
    if (func == "forest") {
      res <- list(
        studlab = TRUE,
        type = NULL,
        func = NULL,
        method = NULL,
        method.interval = NULL,
        method.tau = NULL,
        backtransf = NULL,
        args.backtransf = NULL,
        k = TRUE,
        sm = NULL,
        se = FALSE,
        stat = NULL,
        pval = NULL,
        level.ma = NULL,
        level.predict = NULL,
        level.hetstat = NULL,
        null.effect = FALSE,
        I2 = FALSE,
        I2.ci = FALSE,
        tau2 = gs("forest.tau2"),
        tau2.ci = gs("forest.tau2.ci"),
        tau = gs("forest.tau"),
        tau.ci = gs("forest.tau.ci"),
        Q = FALSE,
        pval.Q = NULL,
        H = FALSE,
        H.ci = FALSE,
        Q.b = NULL,
        pval.Q.b = NULL,
        subgroup.name = FALSE)
    }
    else {
      res <- list(
        func = NULL,
        method = NULL,
        method.interval = NULL,
        backtransf = NULL,
        args.backtransf = NULL,
        k = TRUE,
        sm = NULL,
        se = FALSE,
        stat = TRUE,
        pval = TRUE,
        null.effect = FALSE,
        level.ma = NULL,
        level.predict = NULL,
        level.hetstat = NULL,
        I2 = gs("print.I2"),
        I2.ci = FALSE,
        tau2 = gs("print.tau2"),
        tau2.ci = FALSE,
        tau = gs("print.tau"),
        tau.ci = FALSE,
        Q = gs("print.Q"),
        pval.Q = NULL,
        H = FALSE,
        H.ci = FALSE,
        Q.b = NULL,
        pval.Q.b = NULL,
        subgroup.name = NULL)
    }
  }
  #
  res
}

label_defaults <- function(x) {
  list(
    e = x$label.e,
    c = x$label.c,
    e.attach = gs("label.e.attach"),
    c.attach = gs("label.c.attach"),
    left = x$label.left,
    right = x$label.right,
    n = gs("label.n"),
    events = gs("label.events"),
    mean = gs("label.mean"),
    sd = gs("label.sd"),
    cor = gs("label.cor"),
    time = gs("label.time"),
    pval = gs("label.pval"),
    tau2 = gs("label.tau2"),
    tau = gs("label.tau"),
    I2 = gs("label.I2"),
    cluster = gs("label.cluster"),
    cycles = gs("label.cycles"),
    test.overall.common = NULL,
    test.overall.random = NULL,
    test.subgroup.common = NULL,
    test.subgroup.random = NULL,
    test.effect.subgroup.common = NULL,
    test.effect.subgroup.random = NULL
  )
}

colgap_defaults <- function() {
  colgap <- gs("colgap")
  colgap.forest <- gs("colgap.forest")
  #
  list(
    colgap = colgap,
    left = colgap,
    right = colgap,
    studlab = colgap,
    forest = colgap.forest,
    forest.left = colgap.forest,
    forest.right = colgap.forest,
    rob = "1mm",
    rob.overall = "2mm"
  )
}

just_defaults <- function(layout) {
  just <- if (layout == "JAMA") "left" else "right"
  just.addcols <- gs("just.addcols")
  #
  list(
    just = just,
    main = gs("just.main"),
    studlab = gs("just.studlab"),
    addcols = just.addcols,
    addcols.left = just.addcols,
    addcols.right = just.addcols,
    label.e = just,
    label.c = just,
    rob = "left"
  )
}

digits_defaults <- function(class = "meta", func = "forest") {
  if (class %in% c("meta", "summary.meta") && func == "print") {
    res <- list(
      digits = gs("digits"),
      stat = gs("digits.stat"),
      pval = max(gs("digits.pval"), 2),
      tau2 = gs("digits.tau2"),
      tau = gs("digits.tau"),
      Q = gs("digits.Q"),
      df = gs("digits.df"),
      pval.Q = max(gs("digits.pval.Q"), 2),
      H = gs("digits.H"),
      I2 = gs("digits.I2"))
    #
    if (class == "summary.meta")
      res <- c(res, list(se = gs("digits.se"),
                         prop = gs("digits.prop"),
                         weight = gs("digits.weight")))
  }
  else if (class == "meta") {
    digits <- gs("digits.forest")
    #
    res <- list(
      digits = digits,
      se = gs("digits.se"),
      stat = gs("digits.stat"),
      pval = gs("digits.pval"),
      pval.Q = gs("digits.pval.Q"),
      Q = gs("digits.Q"),
      tau2 = gs("digits.tau2"),
      tau = gs("digits.tau"),
      I2 = gs("digits.I2"),
      weight = gs("digits.weight"),
      mean = gs("digits.mean"),
      sd = gs("digits.sd"),
      cor = digits,
      time = digits,
      n = 0,
      event = 0,
      TE = gs("digits.TE.forest"),
      addcols = digits,
      addcols.right = digits,
      addcols.left = digits)
  }
  else if (class == "metacum" || class == "metainf") {
    res <- list(
      digits = if (func == "print") gs("digits") else gs("digits.forest"),
      stat = gs("digits.stat"),
      pval = gs("digits.pval"),
      tau2 = gs("digits.tau2"),
      tau = gs("digits.tau"),
      I2 = gs("digits.I2"),
      cid = gs("digits.cid"),
      percent = 1)
  }
  else if (class == "metabind") {
    func <- setchar(func, c("forest", "print"), length = 1)
    #
    if (func == "forest") {
      res <- list(
        digits = gs("digits.forest"),
        se = gs("digits.se"),
        stat = gs("digits.stat"),
        pval = max(gs("digits.pval") - 2, 2),
        pval.Q = max(gs("digits.pval.Q") - 2, 2),
        Q = gs("digits.Q"),
        tau2 = gs("digits.tau2"),
        tau = gs("digits.tau"),
        I2 = gs("digits.I2"))
    }
    else {
      res <- list(
        digits = gs("digits.forest"),
        se = gs("digits.se"),
        stat = gs("digits.stat"),
        pval = gs("digits.pval"),
        pval.Q = gs("digits.pval.Q"),
        Q = gs("digits.Q"),
        tau2 = gs("digits.tau2"),
        tau = gs("digits.tau"),
        I2 = gs("digits.I2"))
    }
  }
  #
  res
}

test_defaults <- function(x, common, random, overall,
                          test.overall.default = TRUE, func = "forest") {
  if (func == "print")
    return(list(subgroup = x$test.subgroup,
                subgroup.common = NULL,
                subgroup.random = NULL))
  #
  test.overall <- replaceNULL(gs("test.overall"), test.overall.default)
  test.subgroup <- replaceNULL(x$test.subgroup, TRUE)
  #
  list(
    overall = test.overall,
    overall.common = common & overall & test.overall,
    overall.random = random & overall & test.overall,
    subgroup = test.subgroup,
    subgroup.common = test.subgroup & common,
    subgroup.random = test.subgroup & random,
    effect.subgroup = gs("test.effect.subgroup"),
    effect.subgroup.common = NULL,
    effect.subgroup.random = NULL
  )
}

type_defaults <- function(study.results, class = "meta") {
  if (class == "meta") {
    type.subgroup <- ifelse(study.results, "diamond", "square")
    #
    res <- list(
      study = gs("type.study"),
      common = gs("type.common"),
      random = gs("type.common"),
      subgroup = type.subgroup,
      subgroup.common = type.subgroup,
      subgroup.random = type.subgroup
    )
  }
  else if (class == "metabind") {
    res <- list(
      type = NULL,
      common = NULL,
      random = NULL,
      predict = NULL)
  }
  #
  res
}

fill_defaults <- function() {
  fill.cid <- gs("fill.cid")
  fill.equi <- gs("fill.equi")
  #
  list(
    fill = gs("fill"),
    cid = fill.cid,
    cid.below.null = fill.cid,
    cid.above.null = rev(fill.cid),
    equi = fill.equi,
    lower.equi = fill.equi,
    upper.equi = rev(fill.equi)
  )
}

cid_defaults <- function(x = NULL) {
  list(
    cid = gs("cid"),
    below.null = gs("cid.below.null"),
    above.null = gs("cid.above.null"),
    pooled.only = gs("cid.pooled.only")
  )
}

addrow_defaults <- function() {
  list(
    addrow = gs("addrow"),
    overall = gs("addrow.overall"),
    subgroup = gs("addrow.subgroups")
  )
}

pooled_defaults <- function(common, random) {
  list(
    totals = common | random,
    events = gs("pooled.events"),
    times = gs("pooled.times")
  )
}

text_defaults <- function(x) {
  list(
    common = x$text.common,
    random = x$text.random,
    w.common = x$text.w.common,
    w.random = x$text.w.random,
    predict = x$text.predict,
    subgroup = NULL,
    common.w = x$text.common,
    random.w = x$text.random,
    predict.w = x$text.predict,
    subgroup.nohet = gs("text.subgroup.nohet"),
    addline1 = NULL,
    addline2 = NULL
  )
}

lty_defaults <- function() {
  list(
    common = gs("lty.common"),
    random = gs("lty.random"),
    cid = gs("lty.cid")
  )
}

lwd_defaults <- function() {
  list(
    lwd = gs("lwd"),
    square = gs("lwd.square"),
    diamond = gs("lwd.diamond")
  )
}

rob_defaults <- function(x) {
  list(
    rob = x$rob,
    text = "Risk of Bias",
    legend = TRUE,
    only = FALSE
  )
}
