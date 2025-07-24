.onAttach <- function(libname, pkgname) {
  msg <-
    paste0("Loading 'meta' package (version ",
           utils::packageDescription("meta")$Version,
           ").",
           "\nType 'help(meta)' for a brief overview.")
  packageStartupMessage(msg)
}


updateversion <- function(x) {
  ##
  ## Update older meta objects, see gs("major.update") and gs("minor.update")
  ##
  if (update_needed(x$version))
    x <- update(x, warn = FALSE, warn.deprecated = FALSE)
  ##
  x
}


update_needed <- function(version,
                          major = gs("major.update"),
                          minor = gs("minor.update"),
                          verbose = FALSE) {
  if (is.null(version)) {
    version <- 0.1
    major.cur <- 0
    minor.cur <- 1
  }
  else {
    version <- unlist(strsplit(version, "-"))
    divide <- unlist(strsplit(version[1], ".", fixed = TRUE))
    major.cur <- as.numeric(divide[1])
    minor.cur <- as.numeric(paste(divide[2], version[2], sep = "."))
  }
  ##
  res <-
    ifelse(major.cur < major,
           TRUE, ifelse(major.cur > major,
                        FALSE, minor.cur < as.numeric(gsub("-", ".", minor))))
  if (res & verbose)
    message(paste0("Update to meta, version ", major, ".", minor))
  ##
  res
}


taudat <- function(method.tau, detail.tau,
                   tau, lower.tau, upper.tau, print.tau.ci, digits.tau,
                   tau2, lower.tau2, upper.tau2, print.tau2.ci, digits.tau2,
                   sign.lower, sign.upper) {
  dat <- data.frame(method.tau, names = detail.tau, tau, tau2)
  ##
  ## In order to use duplicated()
  ##
  if (!is.null(dat$tau))
    dat$tau <- round(dat$tau, 14)
  if (!is.null(dat$tau2))
    dat$tau2 <- round(dat$tau2, 14)
  ##
  if (!is.null(lower.tau))
    dat$lower.tau <- round(lower.tau, 14)
  if (!is.null(lower.tau2))
    dat$lower.tau2 <- round(lower.tau2, 14)
  ##
  if (!is.null(upper.tau))
    dat$upper.tau <- round(upper.tau, 14)
  if (!is.null(upper.tau2))
    dat$upper.tau2 <- round(upper.tau2, 14)
  ##
  if (!is.null(sign.lower))
    dat$sign.lower <- sign.lower
  if (!is.null(sign.upper))
    dat$sign.upper <- sign.upper
  ##
  if (!is.null(names(tau)))
    dat$names <- names(tau)
  ##
  sel1 <- !duplicated(dat[, c("method.tau", "tau", "tau2")])
  ##
  if (!is.null(dat$lower.tau) & !is.null(dat$lower.tau2) &
      !is.null(dat$upper.tau) & !is.null(dat$upper.tau2))
    sel2 <- !duplicated(dat[, c("method.tau", "tau", "tau2",
                                "lower.tau", "upper.tau",
                                "lower.tau2", "upper.tau2")])
  else
    sel2 <- sel1
  ##
  sel <- sel1 & sel2
  ##
  dat <- dat[sel, , drop = FALSE]
  ##
  if (print.tau.ci) {
    dat$lower.tau <- round(dat$lower.tau, digits.tau)
    dat$upper.tau <- round(dat$upper.tau, digits.tau)
  }
  ##
  if (print.tau2.ci) {
    dat$lower.tau2 <- round(dat$lower.tau2, digits.tau2)
    dat$upper.tau2 <- round(dat$upper.tau2, digits.tau2)
  }
  
  dat
}


hetdat <- function(I2, lowI2, uppI2, print.I2, print.I2.ci, digits.I2,
                   H, lowH, uppH, print.H, digits.H,
                   Rb, lowRb, uppRb, print.Rb) {
  dat <- data.frame()
  ##
  if (print.I2) {
    if (nrow(dat) == 0)
      dat <- data.frame(I2 = round(I2, digits.I2))
    else
      dat$I2 <- round(I2, digits.I2)
    ##
    if (print.I2.ci) {
      dat$lowI2 <- round(lowI2, digits.I2)
      dat$uppI2 <- round(uppI2, digits.I2)
    }
  }
  ##
  if (print.H) {
    if (nrow(dat) == 0)
      dat <- data.frame(H = round(H, digits.H))
    else
      dat$H <- round(H, digits.H)
    ##
    if (print.I2.ci) {
      dat$lowH <- round(lowH, digits.H)
      dat$uppH <- round(uppH, digits.H)
    }
  }
  ##
  if (print.Rb) {
    if (nrow(dat) == 0)
      dat <- data.frame(Rb = round(Rb, digits.I2))
    else
      dat$Rb <- round(Rb, digits.I2)
    ##
    if (print.I2.ci) {
      dat$lowRb <- round(lowRb, digits.I2)
      dat$uppRb <- round(uppRb, digits.I2)
    }
  }
  ##
  if (!is.null(names(I2)))
    dat$names <- names(I2)
  else if (!is.null(names(H)))
    dat$names <- names(H)
  else if (!is.null(names(Rb)))
    dat$names <- names(Rb)
  ##
  dat <- unique(dat)
  dat
}


qdat <- function(Q, df.Q, pval.Q, hetlabel, text.common) {
  if (length(Q) > 1) {
    if (!is.null(names(Q)))
      rownames.Q <- names(Q)
    else if (!is.null(hetlabel) && length(hetlabel) == length(Q))
      rownames.Q <- hetlabel
    else if (length(text.common) == length(Q))
      rownames.Q <- text.common
    else
      rownames.Q <- rep("", length(Q))
  }
  else
    rownames.Q <- rep("", length(Q))
  ##
  sel <- rownames.Q != ""
  rownames.Q[sel] <- paste0(" ", rownames.Q[sel])
  ##
  dat <- data.frame(Q, df.Q, pval.Q, names = rownames.Q)
  dat <- unique(dat)
  dat
}


collapse <- function(x, quote = '"', collapse = ", ", sort = FALSE)
  paste0(paste0(quote, if (sort) sort(x) else x, quote, collapse = collapse))

collapse2 <- function(x, quote = "", collapse = ", ", br1 = "", br2 = "",
                      sort = FALSE, max.len = 5) {
  if (sort)
    x <- sort(x)
  #
  if (length(x) == 0)
    res <- ""
  else if (length(x) == 1)
    res <- x
  else {
    if (is.numeric(x) && (all(diff(x) == 1) & length(x) > max.len))
      res <-
        paste0(br1, min(x, na.rm = TRUE), ", ..., ", max(x, na.rm = TRUE), br2)
    else
      res <- paste0(br1, paste(x, collapse = collapse), br2)
  }
  #
  res
}

condense <- function(x, var)
  unlist(lapply(x, "[[" , var))


cathet <- function(k,
                   method.tau, detail.tau,
                   tau2, lower.tau2, upper.tau2,
                   print.tau2, print.tau2.ci, text.tau2, digits.tau2,
                   tau, lower.tau, upper.tau,
                   print.tau, print.tau.ci, text.tau, digits.tau,
                   sign.lower.tau, sign.upper.tau,
                   I2, lowI2, uppI2, 
                   print.I2, print.I2.ci, text.I2, digits.I2,
                   H, lowH, uppH,
                   print.H, digits.H,
                   Rb, lowRb, uppRb,
                   print.Rb, text.Rb,
                   big.mark,
                   sort.tau = NULL, sort.het = NULL) {
  
  
  if (is.null(lower.tau2))
    lower.tau2 <- NA
  if (is.null(upper.tau2))
    upper.tau2 <- NA
  if (is.null(lower.tau))
    lower.tau <- NA
  if (is.null(upper.tau))
    upper.tau <- NA
  ##
  if (all(is.na(lower.tau2) & is.na(upper.tau2)))
    print.tau2.ci <- FALSE
  if (all(is.na(lower.tau) & all(is.na(upper.tau))))
    print.tau.ci <- FALSE
  
  
  ## Text for tau2 and tau
  ##
  dtau <- taudat(method.tau, detail.tau,
                 tau, lower.tau, upper.tau, print.tau.ci, digits.tau,
                 tau2, lower.tau2, upper.tau2, print.tau2.ci, digits.tau2,
                 sign.lower.tau, sign.upper.tau)
  ##
  ntau <- nrow(dtau)
  sort.tau <- setsort(sort.tau, ntau, "tau2 estimates")
  dtau <- dtau[sort.tau, , drop = FALSE]
  ##
  if (print.tau2 | print.tau) {
    if (ntau > 1) {
      text.tau2 <- paste(text.tau2, seq_len(ntau), sep = ".")
      text.tau <- paste(text.tau, seq_len(ntau), sep = ".")
    }
  }
  ##
  label.tau <-
    if (ntau > 1)
      ifelse(dtau$names == "", "", paste0(" (", dtau$names, ")"))
    else
      ""
  ##
  if (print.tau2.ci) {
    text.tau2.ci <-
      pasteCI(dtau$lower.tau2, dtau$upper.tau2,
              digits.tau2, big.mark,
              dtau$sign.lower, dtau$sign.upper)
    text.tau2.ci[text.tau2.ci == " "] <- ""
  }
  else
    text.tau2.ci <- ""
  ##
  if (print.tau.ci) {
    text.tau.ci <-
      pasteCI(dtau$lower.tau, dtau$upper.tau,
              digits.tau, big.mark,
              dtau$sign.lower, dtau$sign.upper)
    text.tau.ci[text.tau.ci == " "] <- ""
  }
  else
    text.tau.ci <- ""
  ##
  hettxt <- 
    paste0(
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        " ",
      if (print.tau2)
        paste0(formatPT(dtau$tau2,
                        lab = TRUE, labval = text.tau2,
                        digits = digits.tau2,
                        lab.NA = "NA",
                        big.mark = big.mark),
               text.tau2.ci,
               if (!print.tau) label.tau),
      ##
      if (print.tau)
             paste0(
               if (print.tau2) "; " else "",
               formatPT(dtau$tau,
                        lab = TRUE, labval = text.tau,
                        digits = digits.tau,
                        lab.NA = "NA",
                   big.mark = big.mark),
               text.tau.ci,
               label.tau),
      collapse = "\n")
  
  
  ## Text for I2, H and Rb
  ##
  dhet <- hetdat(I2, lowI2, uppI2, print.I2, print.I2.ci, digits.I2,
                 H, lowH, uppH, print.H, digits.H,
                 Rb, lowRb, uppRb, print.Rb)
  ##
  nhet <- nrow(dhet)
  sort.het <- setsort(sort.het, nhet, "heterogeneity estimates")
  dhet <- dhet[sort.het, , drop = FALSE]
  ##
  label.het <-
    if (nhet > 1)
      ifelse(dhet$names == "", "", paste0(" (", dhet$names, ")"))
    else
      ""
  ##
  if (print.I2) {
    if (nhet > 1)
      text.I2 <- paste(text.I2, seq_len(nhet), sep = ".")
  }
  ##
  if (print.H) {
    text.H <- "H"
    if (nhet > 1)
      text.H <- paste(text.H, seq_len(nhet), sep = ".")
  }
  ##  
  if (print.Rb) {
    if (nhet > 1)
      text.Rb <- paste(text.Rb, seq_len(nhet), sep = ".")
  }
  ##
  if (print.I2)
    hettxt <-
      paste0(hettxt,
             ifelse(print.tau2 | print.tau,
             ifelse(ntau > 1 | print.tau2.ci | print.tau.ci |
                    (options()$width < 70 & print.I2.ci),
                    "\n", ";"),
             ""))
  ##
  hettxt <-
    paste0(hettxt,
           paste0(
             if (print.I2)
               paste0(
                 if (print.tau2 | print.tau)
                   " ",
                 text.I2, " = ",
                 ifelse(is.na(dhet$I2), "NA",
                        paste0(formatN(dhet$I2, digits.I2), "%")),
                 if (print.I2.ci)
                   pasteCI(dhet$lowI2, dhet$uppI2, digits.I2,
                           big.mark, unit = "%"),
                 if (nhet > 1 & !print.H)
                   label.het
               ),
             ##
             if (print.H)
               paste0(
                 if (print.tau2 | print.tau | print.I2)
                   "; " else " ",
                 text.H, " = ",
                 ifelse(is.na(dhet$H), "NA",
                        formatN(dhet$H, digits.H, "NA", big.mark = big.mark)),
                 if (print.I2.ci & any(!is.na(lowH) & !is.na(uppH)))
                   pasteCI(dhet$lowH, dhet$uppH, digits.H, big.mark),
                 if (nhet > 1)
                   label.het),
             collapse = "\n")
           )
  ##
  if (print.Rb)
    hettxt <-
      paste0(hettxt,
             ifelse(print.tau2 | print.tau | print.I2 | print.H,
             ifelse(ntau > 1 | nhet > 1 |
                    print.tau2.ci | print.tau.ci |
                    print.I2.ci |
                    (options()$width < 70 & print.I2.ci),
                    "\n", ";"),
             ""))
  ##
  if (print.Rb)
    hettxt <-
      paste0(hettxt,
             paste0(
               " ",
               text.Rb, " = ",
               ifelse(is.na(dhet$Rb), "NA",
                      paste0(formatN(dhet$Rb, digits.I2,
                                     big.mark = big.mark), "%")),
               if (any(!is.na(dhet$lowRb) & !is.na(dhet$uppRb)))
                 pasteCI(dhet$lowRb, dhet$uppRb, digits.I2,
                         big.mark, unit = "%"),
               label.het,
               collapse = "\n")
             )
  
  
  ## Remove spaces before ";", " (" and a line break
  ##
  hettxt <- gsub("(*UCP)\\s+(;)", "\\1", hettxt, perl = TRUE)
  hettxt <- gsub("(*UCP)\\s+( \\()", "\\1", hettxt, perl = TRUE)
  hettxt <- gsub("(*UCP)\\s+(\\\n)", "\\1", hettxt, perl = TRUE)
  
  
  ## Print information on heterogeneity
  ##
  cat(hettxt)
  
  
  ## Empty row(s)
  ##
  if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
    cat("\n")
  
  invisible(NULL)
}


catobsev <- function(var1, var2 = NULL, type = "n", addrow = FALSE,
                     big.mark = gs("big.mark")) {
  if (type == "n") {
    txt <- "observations"
    idx <- "o"
  }
  else if (type == "e") {
    txt <- "events"
    idx <- "e"
  }
  ##
  if (!is.null(var1) & !is.null(var2)) {
    if (!(all(is.na(var1)) | all(is.na(var2)))) {
      sum1 <- sum(var1, na.rm = TRUE)
      sum2 <- sum(var2, na.rm = TRUE)
      ##
      cat(paste0("Number of ", txt, ": ", idx, " = ",
                 format(sum1 + sum2, big.mark = big.mark),
                 if (type == "n")
                   paste0(" (", idx, ".e = ",
                          format(sum1, big.mark = big.mark),
                          ", ", idx, ".c = ",
                          format(sum2, big.mark = big.mark),
                          ")"),
                 "\n"))
    }
  }
  else if (!is.null(var1)) {
    if (!all(is.na(var1))) {
      cat(paste0("Number of ", txt, ": ", idx, " = ",
                 format(sum(var1, na.rm = TRUE),
                        big.mark = big.mark),
                 "\n"))
    }
  }
  else if (!is.null(var2)) {
    if (!all(is.na(var2))) {
      cat(paste0("Number of ", txt, ": ", idx, " = ",
                 format(sum(var2, na.rm = TRUE), big.mark = big.mark),
                 "\n"))
    }
  }
  ##
  if (addrow)
    cat("\n")
  ##
  invisible(NULL)
}


## The following R code is based on the file snowfall-internal.R from
## R package snowfall (Maintainer: Jochen Knaus <jo@imbi.uni-freiburg.de>)
##
## Helpers for managing the internal variables in the package namespace without
## awake the R CMD check for later R versions (which basically blaims many
## global assignments).
##
## The given solution has an advantage: only writing is affected. Reading of the
## objects can remain the same (thanks to Uwe Ligges for the tipp):
##   reading:  gs("CIbracket")
##   writing:  setOption("CIbracket", "(")
##

##
## Set an option in the meta option list.
## (Basically this is the setting of a list entry).
## key - character: object name
## val - object (everything is allowed, even NULL)
##
setOption <- function(key = NULL, val = NULL) {
  if(!is.null(key) && is.character(key)) {
    option <- getVar(".settings") # Get from NS
    option[[key]] <- val
    setVar(".settings", option) # Write to NS
    ##
    return(invisible(TRUE))
  }
  ##
  stop("Argument 'key' or 'val' is NULL or 'key' is no string.")
}


##
## Get a specific variable from the meta namespace.
## var - character: object name
##
getVar <- function(var = NULL) {
  if(!is.null(var) && is.character(var)) {
    tmp <- try(getFromNamespace(var, "meta"))
    ##
    if(inherits(tmp, "try-error"))
      stop("Object", var, "not found in meta package.")
    ##
    return(tmp)
  }

  stop("Argument 'var' is NULL or not a string.")
}


##
## Write a specific variable to the meta namespace.
## var - character: object name
## arg - object (NULL allowed)
##
setVar <- function(var = NULL, arg = NULL) {
  if(!is.null(var) && is.character(var)) {
    assignInNamespace(var, arg, "meta")

    return(invisible(TRUE))
  }

  stop("var is NULL or no character");
}

second <- function(x) x[2]

setOptionDepr <- function(x, new, old, func, ...) {
  newval <- do.call(func, list(argname = new, args = x, ...))
  oldval <- do.call(func, list(argname = old, args = x, ...))
  #
  if (is.na(newval) & !is.na(oldval))
    setOption(new, x[[oldval]])
  #
  invisible(NULL)
}





.settings <- list()
##
## List of internal settings
##
argslist.internal <-
  c("Wan2014.Table1", "Wan2014.Table2",
    "sm4bin", "sm4cont", "sm4cor", "sm4inc", "sm4mean", "sm4prop", "sm4rate",
    "ci4cont", "ci4prop", "ci4rate",
    "meth4bin", "meth4inc", "meth4prop", "meth4rate",
    "meth4tau", "meth4tau.ci", "meth4i2",
    "meth4common.ci",
    "meth4random.ci", "meth4pi",
    "adhoc4hakn.ci", "adhoc4hakn.pi",
    "meth4bias", "meth4bias.old",
    "tool4rob",
    "meth4incr",
    "text.fixed", "text.w.fixed",
    "major.update", "minor.update",
    #
    "special.characters")
##
setOption("argslist.internal", argslist.internal)
##
## Set defaults (for internal options)
##
setOption("sm4bin", c("OR", "RD", "RR", "ASD", "DOR", "VE"))
setOption("sm4cont", c("MD", "SMD", "ROM"))
setOption("sm4cor", c("ZCOR", "COR"))
setOption("sm4inc", c("IRR", "IRD", "IRSD", "VE"))
setOption("sm4mean", c("MRAW", "MLN"))
setOption("sm4prop", c("PLOGIT", "PLN", "PRAW", "PAS", "PFT"))
setOption("sm4rate", c("IR", "IRLN", "IRS", "IRFT"))
##
setOption("ci4cont", c("z", "t"))
setOption("ci4prop", c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
setOption("ci4rate", c("NAsm", "Poisson"))
##
setOption("meth4bin", c("Inverse", "MH", "Peto", "GLMM", "LRP", "SSW"))
setOption("meth4inc", c("Inverse", "MH", "Cochran", "GLMM"))
setOption("meth4prop", c("Inverse", "GLMM"))
setOption("meth4rate", c("Inverse", "GLMM"))
##
setOption("meth4tau", c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
setOption("meth4tau.ci", c("QP", "BJ", "J", "PL", ""))
setOption("meth4i2", c("Q", "tau2"))
setOption("meth4common.ci", c("classic", "IVhet"))
setOption("meth4random.ci", c("classic", "HK", "KR"))
setOption("meth4pi",
          c("V", "HTS", "HK", "HK-PR", "KR", "KR-PR", "NNF", "S", ""))
setOption("adhoc4hakn.ci", c("", "se", "ci", "IQWiG6"))
setOption("adhoc4hakn.pi", c("", "se"))
##
setOption("meth4bias.old", c("rank", "linreg", "mm", "count", "score"))
setOption("meth4bias", c("Begg", "Egger", "Thompson", "Schwarzer",
                         "Harbord", "Peters", "Deeks",
                         "Pustejovsky", "Macaskill"))
##
setOption("tool4rob",
          c("RoB1", "RoB2", "RoB2-cluster", "RoB2-crossover",
            "ROBINS-I", "ROBINS-E"))
##
setOption("meth4incr", c("only0", "if0all", "all", "user"))
#
setOption("special.characters", c("+", ".", "&", "$", "#", "|", "*", "^"))
#
setOption("major.update", 5)
setOption("minor.update", 6)
##
## List of arguments that can be changed by user
##
argslist <-
  c("level", "level.ma", "common", "random",
    "method.common.ci", "method.random.ci", "method.predict",
    "adhoc.hakn.ci", "adhoc.hakn.pi",
    "method.tau", "method.tau.ci", "level.hetstat", "tau.common",
    "method.I2",
    "prediction", "level.predict",
    "method.bias",
    "tool.rob",
    "overall.hetstat",
    "text.common", "text.random", "text.predict",
    "text.w.common", "text.w.random",
    "title", "complab",
    "CIbracket", "CIseparator", "CIlower.blank", "CIupper.blank",
    "print.subgroup.name", "sep.subgroup",
    "keepdata", "keeprma", "warn", "warn.deprecated",
    "transf", "backtransf",
    "smbin", "smcont", "smcor", "sminc", "smmean", "smprop", "smrate",
    "incr", "method.incr",
    "method", "allstudies", "MH.exact",
    "RR.Cochrane", "Q.Cochrane", "model.glmm", "print.CMH",
    "pooledvar", "method.smd", "sd.glass", "exact.smd",
    "method.ci.cont", "method.ci.prop", "method.ci.rate",
    "label.e", "label.c", "label.left", "label.right",
    "layout", "forest.details",
    "test.overall", "test.subgroup", "prediction.subgroup",
    "test.effect.subgroup",
    "digits", "digits.mean", "digits.sd", "digits.se", "digits.stat",
    "digits.Q", "digits.tau2", "digits.tau", "digits.H", "digits.I2",
    "digits.prop", "digits.weight",
    "digits.pval", "digits.pval.Q",
    "digits.forest", "digits.TE.forest",
    "digits.df", "digits.cid",
    "scientific.pval", "big.mark", "zero.pval", "JAMA.pval",
    "details",
    "print.tau2", "print.tau2.ci", "print.tau", "print.tau.ci",
    "print.I2", "print.I2.ci", "print.H", "print.Rb",
    "text.tau2", "text.tau", "text.I2", "text.Rb",
    "print.Q",
    ##
    "lty.common", "lty.random", "col.common", "col.random",
    "sort.subgroup",
    "pooled.events", "pooled.times", "study.results",
    "cid", "cid.below.null", "cid.above.null", "lty.cid", "col.cid", "fill.cid",
    "cid.pooled.only",
    "fill", "fill.equi",
    "leftcols", "rightcols", "leftlabs", "rightlabs", 
    "label.e.attach", "label.c.attach",
    "bottom.lr",
    "lab.NA", "lab.NA.effect", "lab.NA.weight",
    "lwd", "lwd.square", "lwd.diamond",
    "arrow.type", "arrow.length",
    "type.study", "type.common",
    "col.study", "col.square", "col.square.lines", "col.circle", "col.inside",
    "col.diamond", "col.diamond.lines",
    "col.predict", "col.predict.lines",
    "col.subgroup",
    "col.label.right", "col.label.left",
    "col.lines", "col.label",
    "hetlab", "resid.hetstat", "resid.hetlab",
    "forest.I2", "forest.I2.ci", "forest.tau2", "forest.tau2.ci",
    "forest.tau", "forest.tau.ci", "forest.Q", "forest.pval.Q",
    "forest.Rb", "forest.Rb.ci",
    "text.subgroup.nohet",
    "LRT",
    "forest.stat", "forest.Q.subgroup",
    "header.line",
    "fontsize", "fontfamily",
    "fs.common", "fs.random", "fs.predict",
    "fs.common.labels", "fs.random.labels", "fs.predict.labels",
    "fs.hetstat", "fs.test.overall", "fs.test.subgroup",
    "fs.test.effect.subgroup", "fs.addline",
    "ff.heading", "ff.common", "ff.random", "ff.predict",
    "ff.common.labels", "ff.random.labels", "ff.predict.labels",
    "ff.study", "ff.hetstat", "ff.test.overall", "ff.test.subgroup",
    "ff.test.effect.subgroup", "ff.addline",
    "ff.axis", "ff.smlab", "ff.xlab", "ff.lr",
    "colgap", "colgap.forest",
    "width",
    "calcwidth.predict", "calcwidth.hetstat",
    "calcwidth.tests", "calcwidth.subgroup", "calcwidth.addline",
    "just.studlab", "just.addcols",
    "spacing",
    "addrow", "addrow.overall", "addrow.subgroups", "addrows.below.overall"
    )
#
args.depr <- c("fixed", "comb.fixed", "comb.random", "level.comb",
               "hakn", "adhoc.hakn",
               "digits.zval", "print.byvar", "byseparator",
               "addincr", "allincr",
               "lower.equi", "upper.equi", "lty.equi", "col.equi")
#
setOption("argslist", c(argslist, args.depr))
#
setOption("argslist.meta", c(argslist, args.depr))
##
## General settings
##
setOption("level", 0.95)
setOption("level.ma", 0.95)
setOption("level.comb", 0.95)
##
setOption("common", TRUE)
setOption("fixed", TRUE)
setOption("comb.fixed", TRUE)
setOption("random", TRUE)
setOption("comb.random", TRUE)
setOption("method.common.ci", "classic")
setOption("method.random.ci", "classic")
setOption("hakn", FALSE)
setOption("adhoc.hakn", "")
setOption("adhoc.hakn.ci", "")
setOption("adhoc.hakn.pi", "")
setOption("prediction", FALSE)
setOption("level.predict", 0.95)
setOption("method.predict", "V")
setOption("method.tau", "REML")
setOption("method.tau.ci", NULL)
setOption("level.hetstat", 0.95)
setOption("tau.common", FALSE)
setOption("method.I2", "Q")
setOption("method.bias", "Egger")
setOption("tool.rob", NULL)
setOption("overall.hetstat", NULL)
setOption("text.common", "Common effect model")
setOption("text.fixed", "Common effect model")
setOption("text.random", "Random effects model")
setOption("text.predict", "Prediction interval")
setOption("text.w.common", "common")
setOption("text.w.fixed", "common")
setOption("text.w.random", "random")
setOption("title", "")
setOption("complab", "")
setOption("CIbracket", "[")
setOption("CIseparator", "; ")
setOption("CIlower.blank", TRUE)
setOption("CIupper.blank", TRUE)
setOption("print.subgroup.name", TRUE)
setOption("print.byvar", TRUE)
setOption("sep.subgroup", " = ")
setOption("byseparator", " = ")
setOption("test.subgroup", TRUE)
setOption("prediction.subgroup", FALSE)
setOption("keepdata", TRUE)
setOption("keeprma", FALSE)
setOption("warn", TRUE)
setOption("warn.deprecated", TRUE)
setOption("transf", TRUE)
setOption("backtransf", TRUE)
setOption("digits", 4)
setOption("digits.mean", 2)
setOption("digits.sd", 4)
setOption("digits.se", 4)
setOption("digits.stat", 2)
setOption("digits.zval", 2)
setOption("digits.Q", 2)
setOption("digits.tau2", 4)
setOption("digits.tau", 4)
setOption("digits.H", 2)
setOption("digits.I2", 1)
setOption("digits.prop", 4)
setOption("digits.weight", 1)
setOption("digits.pval", 4)
setOption("digits.pval.Q", 4)
setOption("digits.df", 4)
setOption("digits.cid", 4)
setOption("scientific.pval", FALSE)
setOption("big.mark", "")
setOption("zero.pval", TRUE)
setOption("JAMA.pval", FALSE)
##
setOption("details", TRUE)
##
setOption("print.tau2", TRUE)
setOption("print.tau2.ci", TRUE)
setOption("print.tau", TRUE)
setOption("print.tau.ci", TRUE)
setOption("print.I2", TRUE)
setOption("print.I2.ci", TRUE)
setOption("print.H", TRUE)
setOption("print.Rb", FALSE)
setOption("text.tau2", "tau^2")
setOption("text.tau", "tau")
setOption("text.I2", "I^2")
setOption("text.Rb", "Rb")
##
setOption("print.Q", TRUE)
##
## Default summary measure
##
setOption("smbin", "RR")
setOption("smcont", "MD")
setOption("smcor", "ZCOR")
setOption("sminc", "IRR")
setOption("smmean", "MRAW")
setOption("smprop", "PLOGIT")
setOption("smrate", "IRLN")
##
## Settings for R functions metabin, metainc, metaprop
##
setOption("incr", 0.5)
setOption("method.incr", "only0")
setOption("allincr", FALSE)
setOption("addincr", FALSE)
##
## Additional settings for R function metabin
##
setOption("method", "MH")
setOption("allstudies", FALSE)
setOption("MH.exact", FALSE)
setOption("RR.Cochrane", FALSE)
setOption("Q.Cochrane", TRUE)
setOption("model.glmm", "UM.FS")
setOption("print.CMH", FALSE)
##
## Additional setting for R function metacont
##
setOption("pooledvar", FALSE)
setOption("method.smd", "Hedges")
setOption("sd.glass", "control")
setOption("exact.smd", TRUE)
setOption("method.ci.cont", "z")
##
## Additional setting for R function metaprop
##
setOption("method.ci.prop", "CP")
##
## Additional setting for R function metarate
##
setOption("method.ci.rate", "NAsm")
##
## Settings for R functions comparing two treatments
##
setOption("label.e", "Experimental")
setOption("label.c", "Control")
setOption("label.left", "")
setOption("label.right", "")
##
## Settings for R function forest.meta
##
setOption("layout", "meta")
setOption("forest.details", FALSE)
setOption("test.overall", FALSE)
setOption("test.effect.subgroup", FALSE)
setOption("digits.forest", 2)
setOption("digits.TE.forest", 4)
##
setOption("lty.common", 2)
setOption("lty.random", 3)
setOption("col.common", "black")
setOption("col.random", "black")
##
setOption("sort.subgroup", FALSE)
##
setOption("pooled.events", FALSE)
setOption("pooled.times", FALSE)
setOption("study.results", TRUE)
#
setOption("cid", NA)
setOption("cid.below.null", NA)
setOption("cid.above.null", NA)
setOption("lty.cid", 1)
setOption("col.cid", "blue")
setOption("fill.cid", "transparent")
setOption("cid.pooled.only", FALSE)
#
setOption("lower.equi", NA)
setOption("upper.equi", NA)
setOption("lty.equi", 1)
setOption("col.equi", "blue")
setOption("fill.equi", "transparent")
#
setOption("fill", "transparent")
##
setOption("leftcols", NULL)
setOption("rightcols", NULL)
setOption("leftlabs", NULL)
setOption("rightlabs", NULL)
##
setOption("label.e.attach", NULL)
setOption("label.c.attach", NULL)
##
setOption("bottom.lr", TRUE)
##
setOption("lab.NA", ".")
setOption("lab.NA.effect", NULL)
setOption("lab.NA.weight", ".")
##
setOption("lwd", 1)
setOption("lwd.square", 1)
setOption("lwd.diamond", 1)
##
setOption("arrow.type", "open")
setOption("arrow.length", 0.05)
##
setOption("type.study", "square")
setOption("type.common", "diamond")
##
setOption("col.study", "black")
setOption("col.square", "gray")
setOption("col.square.lines", "gray")
setOption("col.circle", "royalblue")
setOption("col.inside", "white")
setOption("col.diamond", "gray")
setOption("col.diamond.lines", "black")
setOption("col.predict", "red")
setOption("col.predict.lines", "black")
setOption("col.subgroup", "black")
setOption("col.label.right", "black")
setOption("col.label.left", "black")
##
setOption("col.lines", "black")
setOption("col.label", "black")
##
setOption("hetlab", "Heterogeneity: ")
setOption("resid.hetstat", NULL)
setOption("resid.hetlab", "Residual heterogeneity: ")
##
setOption("forest.I2", NULL)
setOption("forest.I2.ci", FALSE)
setOption("forest.tau2", NULL)
setOption("forest.tau2.ci", FALSE)
setOption("forest.tau", FALSE)
setOption("forest.tau.ci", FALSE)
setOption("forest.Q", FALSE)
setOption("forest.pval.Q", NULL)
setOption("forest.Rb", FALSE)
setOption("forest.Rb.ci", FALSE)
##
setOption("text.subgroup.nohet", "not applicable")
##
setOption("LRT", FALSE)
##
setOption("forest.stat", TRUE)
setOption("forest.Q.subgroup", TRUE)
##
setOption("header.line", FALSE)
##
setOption("fontsize", 12)
setOption("fontfamily", NULL)
setOption("fs.common", NULL)
setOption("fs.random", NULL)
setOption("fs.predict", NULL)
setOption("fs.common.labels", NULL)
setOption("fs.random.labels", NULL)
setOption("fs.predict.labels", NULL)
setOption("fs.hetstat", NULL)
setOption("fs.test.overall", NULL)
setOption("fs.test.subgroup", NULL)
setOption("fs.test.effect.subgroup", NULL)
setOption("fs.addline", NULL)
##
setOption("ff.heading", "bold")
setOption("ff.common", NULL)
setOption("ff.random", NULL)
setOption("ff.predict", NULL)
setOption("ff.common.labels", NULL)
setOption("ff.random.labels", NULL)
setOption("ff.predict.labels", NULL)
setOption("ff.study", "plain")
setOption("ff.hetstat", NULL)
setOption("ff.test.overall", NULL)
setOption("ff.test.subgroup", NULL)
setOption("ff.test.effect.subgroup", NULL)
setOption("ff.addline", NULL)
setOption("ff.axis", "plain")
setOption("ff.smlab", "bold")
setOption("ff.xlab", "plain")
setOption("ff.lr", "plain")
##
setOption("colgap", "2mm")
setOption("colgap.forest", "2mm")
##
setOption("width", NULL)
##
setOption("calcwidth.predict", FALSE)
setOption("calcwidth.hetstat", FALSE)
setOption("calcwidth.tests", FALSE)
setOption("calcwidth.subgroup", FALSE)
setOption("calcwidth.addline", FALSE)
##
setOption("just.studlab", "left")
setOption("just.addcols", "center")
##
setOption("spacing", 1)
setOption("addrow", NULL)
setOption("addrow.overall", NULL)
setOption("addrow.subgroups", NULL)
setOption("addrows.below.overall", NULL)


setOption("Wan2014.Table1",
          c(0.000, 1.128, 1.693, 2.059, 2.326,
            2.534, 2.704, 2.847, 2.970, 3.078,
            3.173, 3.259, 3.336, 3.407, 3.472,
            3.532, 3.588, 3.640, 3.689, 3.735,
            3.778, 3.819, 3.858, 3.895, 3.931,
            3.964, 3.997, 4.027, 4.057, 4.086,
            4.113, 4.139, 4.165, 4.189, 4.213,
            4.236, 4.259, 4.280, 4.301, 4.322,
            4.341, 4.361, 4.379, 4.398, 4.415,
            4.433, 4.450, 4.466, 4.482, 4.498))
##
setOption("Wan2014.Table2",
          c(0.990, 1.144, 1.206, 1.239, 1.260,
            1.274, 1.284, 1.292, 1.298, 1.303,
            1.307, 1.311, 1.313, 1.316, 1.318,
            1.320, 1.322, 1.323, 1.324, 1.326,
            1.327, 1.328, 1.329, 1.330, 1.330,
            1.331, 1.332, 1.332, 1.333, 1.333,
            1.334, 1.334, 1.335, 1.335, 1.336,
            1.336, 1.336, 1.337, 1.337, 1.337,
            1.338, 1.338, 1.338, 1.338, 1.339,
            1.339, 1.339, 1.339, 1.339, 1.340))
