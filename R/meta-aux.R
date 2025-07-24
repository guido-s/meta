## Auxiliary functions
##
## Package: meta
## Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
## License: GPL (>= 2)
##

allNA <- function(x)
  all(is.na(x))

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

catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)

catch2 <- function(x, varname, return = NULL, fromobject = FALSE) {
  if (!is.null(x[[varname]]) & fromobject)
    return(x[[varname]])
  else if (isCol(x$data, paste0(".", varname)))
    return(x$data[[paste0(".", varname)]])
  else
    return(return)
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

deprecated2 <- function(newvar, newmiss, oldvar, oldmiss, warn = TRUE,
                        oldtxt = NULL) {
  ##
  new <- deparse(substitute(newvar))
  if (is.null(oldtxt))
    oldtxt <- deparse(substitute(oldvar))
  ##
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
      return(any(x != y, na.rm = TRUE))
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
  if (!is.null(res$three.level) && any(res$three.level))
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

extrVar <- function(x, name)
  x[[name]]

# Function only used with MLM, GLMM, LRP, or USW
#
calcPI <- function(x) {
  
  res <- x
  ##
  k <- x$k
  TE.random <- x$TE.random[1]
  seTE.random <- x$seTE.random[1]
  method.predict <- x$method.predict
  level.predict <- x$level.predict
  ##
  res$df.predict <- res$upper.predict <- res$lower.predict <-
    res$seTE.predict <- NA
  ##
  tau2.calc <- if (is.na(sum(res$tau2))) 0 else sum(res$tau2)
  seTE.predict <- sqrt(seTE.random^2 + tau2.calc)
  #
  df.predict <- set_df_predict(method.predict, k)
  #
  ci.p <- ci(TE.random, seTE.predict, level.predict, df = df.predict)
  ##
  res$seTE.predict <- rep_len(seTE.predict, length(df.predict))
  res$lower.predict <- ci.p$lower
  res$upper.predict <- ci.p$upper
  res$df.predict <- ci.p$df
  ##
  if (length(method.predict) > 1)
    names(res$seTE.predict) <-
      names(res$lower.predict) <- names(res$upper.predict) <-
      names(res$df.predict) <- method.predict
  ##
  res
}

runMLM <- function(x, method.tau, method.random.ci, level,
                   warn = TRUE, ...) {
  
  n.methci <- length(method.random.ci)
  ##
  list.b <- list(level = 100 * level, method = method.tau,
                 random = as.call(~ 1 | cluster / idx), ...)
  ##
  res <- list.r <- vector("list", n.methci)
  ##
  for (i in seq_len(n.methci))
    list.r[[i]] <- c(x, list.b,
                     test = ifelse(method.random.ci[i] == "HK", "t", "z"))
  ##
  for (i in seq_len(n.methci))
    res[[i]] <- runNN(rma.mv, list.r[[i]], warn = warn)
  ##
  res
}

extrMLM <- function(x, k, len, sel,
                    method.random.ci, method.predict,
                    level.ma, level.predict, null.effect) {
  
  ##
  ## Random effects model(s)
  ##
  res <- list(tau2 = NA, tau = NA, se.tau2 = NA)
  ##
  res$k <- k
  res$method.predict <- method.predict
  res$level.predict <- level.predict
  ##
  if (k > 1) {
    res$tau2 <- x[[1]]$sigma2
    res$tau <- sqrt(res$tau2)
  }
  ##
  res$TE.random   <- as.numeric(x[[1]]$b)
  res$seTE.random <- sapply(x, extrVar, "se")
  ##
  res$w.random <- rep_len(NA, len)
  res$w.random[sel] <- weights(x[[1]], type = "rowsum")
  res$w.random[is.na(res$w.random)] <- 0
  ##
  res$lower.random <- sapply(x, extrVar, "ci.lb")
  res$upper.random <- sapply(x, extrVar, "ci.ub")
  res$statistic.random <- sapply(x, extrVar, "zval")
  res$pval.random <-
    ci(res$TE.random, res$seTE.random,
       level = level.ma,
       null.effect = null.effect,
       df = ifelse(method.random.ci == "HK", k - 1, Inf))$p
  res$zval.random <- res$statistic.random
  ##
  if (length(method.random.ci) > 1)
    names(res$seTE.random) <-
      names(res$lower.random) <- names(res$upper.random) <-
      names(res$statistic.random) <- names(res$pval.random) <-
      names(res$zval.random) <- method.random.ci
  ##
  ## Prediction interval
  ##
  res <- calcPI(res)
  ##
  res$k <- res$method.predict <- res$level.predict <- NULL
  ##
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
  ##
  list.b <- list(level = 100 * level, ...)
  list.c <- c(x, list.b, method = "FE", test = "z")
  ##
  res.r <- list.r <- vector("list", n.methci)
  ##
  for (i in seq_len(n.methci))
    list.r[[i]] <- c(x, list.b, method = method.tau,
                     test = ifelse(method.random.ci[i] == "HK", "t", "z"))
  ##
  res.c <- runNN(rma.glmm, list.c, warn = FALSE)
  ##
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
      ## Fallback to common effect model due to small number of
      ## studies or zero or all events
      res.r[[i]] <- res.c
    }
  }
  
  res <- list(glmm.common = res.c, glmm.random = res.r)
  ##
  res
}

addGLMM <- function(x, glmm, method.I2, null.effect = x$null.effect) {
  
  res <- x
  ##
  k <- x$k
  len <- length(x$TE)
  level.ma <- x$level.ma
  level.predict <- x$level.predict
  method.random.ci <- x$method.random.ci
  method.predict <- x$method.predict
  null.effect <- replaceNULL(null.effect, 0)
  ##
  ## Common effect model
  ##
  glmm.common <- glmm$glmm.common
  ##
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
  ##
  ## Random effects model(s)
  ##
  glmm.random <- glmm$glmm.random
  ##
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
  ##
  if (length(method.random.ci) > 1)
    names(res$seTE.random) <-
      names(res$lower.random) <- names(res$upper.random) <-
      names(res$statistic.random) <- names(res$pval.random) <-
      names(res$zval.random) <- method.random.ci
  ##
  res$se.tau2 <- res$tau <- res$tau2 <- NA
  ##
  if (k > 1) {
    res$tau2 <- glmm.random[[1]]$tau2
    res$tau <- sqrt(res$tau2)
    res$se.tau2 <- glmm.random[[1]]$se.tau2
  }
  ##
  ## Prediction interval
  ##
  res <- calcPI(res)
  ##
  ## Heterogeneity measures
  ##
  gr1 <- glmm.random[[1]]
  ##
  Q <- if (gr1$k > 1) gr1$QE.Wld else 0
  df.Q <- gr1$QE.df
  pval.Q <- pvalQ(Q, df.Q)
  ##
  Q.LRT <- if (gr1$k > 1) gr1$QE.LRT else 0
  pval.Q.LRT <- pvalQ(Q.LRT, df.Q)
  ##
  res$Q <- c(Q, Q.LRT)
  res$df.Q <- c(df.Q, df.Q)
  res$pval.Q <- c(pval.Q, pval.Q.LRT)
  names(res$Q) <- c("Wald", "LRT")
  ##
  res$upper.tau <- res$lower.tau <- res$upper.tau2 <- res$lower.tau2 <- NA
  ##
  res$sign.upper.tau <- res$sign.lower.tau <- res$method.tau.ci <- ""
  ##
  if (method.I2 == "Q") {
    H <- calcH(Q, df.Q, level.ma)
    res$H <- H$TE
    res$lower.H <- H$lower
    res$upper.H <- H$upper
    #
    I2 <- isquared(Q, df.Q, level.ma)
    res$I2 <- I2$TE
    res$lower.I2 <- I2$lower
    res$upper.I2 <- I2$upper
  }
  else {
    res$H <- sqrt(glmm$H2)
    res$lower.H <- NA
    res$upper.H <- NA
    #
    res$I2 <- glmm$I2 / 100
    res$lower.I2 <- NA
    res$upper.I2 <- NA
  }
  ##
  res$upper.Rb <- res$lower.Rb <- res$Rb <- NA
  ##
  res$.glmm.common <- glmm.common
  ##
  if (length(glmm.random) == 1)
    res$.glmm.random <- glmm.random[[1]]
  else
    res$.glmm.random <- glmm.random
  ##
  res$version.metafor <- packageDescription("metafor")$Version
  ##
  res
}

hccGLMM <- function(x, glmm, method.I2) {
  Q.r <- glmm$QE.Wld
  df.Q.r <- glmm$k - glmm$p
  ##
  if (method.I2 == "Q") {
    H.r  <- calcH(Q.r, df.Q.r, x$level.ma)
    I2.r <- isquared(Q.r, df.Q.r, x$level.ma)
  }
  else {
    H.r  <- list(TE = sqrt(glmm$H2), lower = NA, upper = NA)
    I2.r <- list(TE = glmm$I2 / 100, lower = NA, upper = NA)
  }
  ##
  list(tau2.resid = glmm$tau2,
       lower.tau2.resid = NA,
       upper.tau2.resid = NA,
       ##
       tau.resid = sqrt(glmm$tau2),
       lower.tau.resid = NA,
       upper.tau.resid = NA,
       sign.lower.tau.resid = "",
       sign.upper.tau.resid = "",
       ##
       Q.resid = Q.r,
       df.Q.resid = df.Q.r,
       pval.Q.resid = pvalQ(Q.r, df.Q.r),
       ##
       H.resid = H.r$TE,
       lower.H.resid = H.r$lower,
       upper.H.resid = H.r$upper,
       ##
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
  ##
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
    ##
    if (!is.null(y)) {
      for (i in seq_len(length(x)))
        x[[i]] <- matrix(x[[i]], nrow = length(x[[i]]),
                         ncol = ncol(y) / length(x))
    }
    ##
    lnames <- names(x)
    ##    
    if (!is.null(dimnames(x[[1]])[[2]])) {
      if (!is.null(y)) {
        for (i in seq_along(lnames)) {
          dimnames(x[[i]]) <- dimnames(y)
        }
      }
      else {
        for (i in seq_along(lnames)) {
          dimnames(x[[i]])[[2]] <-
            paste(lnames[i], dimnames(x[[i]])[[2]], sep = "-")
          ##
          dimnames(x[[i]])[[2]] <-
            gsub("-classic$", "", dimnames(x[[i]])[[2]])
        }
      }
    }
    ##
    return(do.call("cbind", x))
  }
  else if (is.vector(x))
    return(matrix(x, ncol = 1))
  else
    return(x)
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
  ##
  if (!is.null(type))
    x <- ifelse(is.na(x) & type == "I", "", x)
  ##
  x
}

sm2meth <- function(x) {
  x[x == "IV"] <- "Inverse"
  x[x == "EXP_O_E_VAR"] <- "Peto"
  ##
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
  ##
  if (lenOK & length(x) != n)
    res <- rep(x, rep_len(n, length(x)))
  ##
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

drop_from_dots <- function(x, old, new) {
  for (i in seq_along(old)) {
    if (!is.null(x[[old[i]]])) {
      if (new[i] != "")
        warning("Argument '", old[i],
                "' ignored; please use argument '",
                new[i], "' instead.",
                call. = FALSE)
      else
        warning("Argument '", old[i],
                "' ignored as it is used internally.",
                call. = FALSE)
      #
      x[[old[i]]] <- NULL
    }
  }
  x
}
