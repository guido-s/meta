.onAttach <- function(libname, pkgname) {
  msg <-
    paste0("Loading 'meta' package (version ",
           utils::packageDescription("meta")$Version,
           ").",
           "\nType 'help(meta)' for a brief overview.",
           "\nReaders of 'Meta-Analysis with R (Use R!)' should install",
           "\nolder version of 'meta' package: ",
           "https://tinyurl.com/dt4y5drs")
  packageStartupMessage(msg)
}


updateversion <- function(x) {
  ##
  ## Update older meta objects - see gs("major.update"), gs("minor.update")
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
    version <- unlist(strsplit(version, "-")[1])
    major.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[1])
    minor.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[2])
  }
  ##
  res <-
    ifelse(major.cur < major,
           TRUE, ifelse(major.cur > major,
                        FALSE, minor.cur < minor))
  if (res & verbose)
    message(paste0("Update to meta, version ", major, ".", minor))
  ##
  res
}


cathet <- function(k,
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
                   detail.tau = "") {
  
  
  if (is.null(lower.tau2))
    lower.tau2 <- NA
  if (is.null(upper.tau2))
    upper.tau2 <- NA
  if (is.null(lower.tau))
    lower.tau <- NA
  if (is.null(upper.tau))
    upper.tau <- NA
  ##
  if (all(is.na(lower.tau2)) && all(is.na(upper.tau2)))
    print.tau2.ci <- FALSE
  if (all(is.na(lower.tau)) && all(is.na(upper.tau)))
    print.tau.ci <- FALSE
  
  
  stau <- length(tau) == 1
  ##
  if (!stau) {
    text.tau2 <- paste(text.tau2, seq_along(tau), sep = ".")
    text.tau <- paste(text.tau, seq_along(tau), sep = ".")
  }
  ##
  detail.tau <- ifelse(detail.tau != "", paste0(" (", detail.tau, ")"), "")
  
  
  cat(
    paste(
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        " ",
      if (print.tau2)
        paste0(formatPT(tau^2,
                        lab = TRUE, labval = text.tau2,
                        digits = digits.tau2,
                        lab.NA = "NA",
                        big.mark = big.mark),
               if (print.tau2.ci)
                 pasteCI(lower.tau2, upper.tau2, digits.tau2, big.mark,
                         sign.lower.tau, sign.upper.tau),
               if (!print.tau) detail.tau),
      ##
      if (print.tau)
        paste0(
          if (print.tau2) "; " else "",
          formatPT(tau,
                   lab = TRUE, labval = text.tau,
                   digits = digits.tau,
                   lab.NA = "NA",
                   big.mark = big.mark),
          if (print.tau.ci)
            pasteCI(lower.tau, upper.tau, digits.tau, big.mark,
                    sign.lower.tau, sign.upper.tau),
          detail.tau),
      sep = "", collapse = "\n")
  )
  ##
  cat(
    paste0(
      if (print.I2)
        paste0(
          ifelse(
            print.tau2 | print.tau,
          ifelse(!stau | print.tau2.ci | print.tau.ci |
                 (options()$width < 70 & print.I2.ci),
                 "\n", ";"),
          ""),
          if (print.tau2 | print.tau)
            " ",
          text.I2, " = ",
          if (is.na(I2))
            "NA"
          else
            paste0(formatN(I2, digits.I2), "%"),
          if (print.I2.ci)
            pasteCI(lowI2, uppI2, digits.I2, big.mark, unit = "%")
        ),
      ##
      if (print.H)
        paste0(
          if (print.tau2 | print.tau | print.I2)
            "; ",
          "H = ",
          if (is.na(H))
            "NA"
          else
            formatN(H, digits.H, "NA", big.mark = big.mark),
          if (!(is.na(lowH) | is.na(uppH)))
            pasteCI(lowH, uppH, digits.H, big.mark)
        ),
      ##
      if (print.Rb)
        paste0(
          if (print.tau2 | print.tau | print.I2 | print.H)
            ";\n",
          text.Rb, " = ",
          if (is.na(Rb))
            "NA"
          else
            paste0(formatN(Rb, digits.I2, big.mark = big.mark), "%"),
          if (!(is.na(lowRb) | is.na(uppRb)))
            pasteCI(lowRb, uppRb, digits.I2, big.mark, unit = "%")
        ),
      ##
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        "\n"
    )
  )
  
  
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
                 ##" (", idx, ".e = ",
                 ##format(sum1, big.mark = big.mark),
                 ##", ", idx, ".c = ",
                 ##format(sum2, big.mark = big.mark),
                 ##")",
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






.settings <- list()
##
## List of internal settings
##
argslist.internal <-
  c("fixed", "comb.fixed", "comb.random", "level.comb", "digits.zval",
    "print.byvar", "byseparator",
    "Wan2014.Table1", "Wan2014.Table2",
    "sm4bin", "sm4cont", "sm4cor", "sm4inc", "sm4mean", "sm4prop", "sm4rate",
    "ci4cont", "ci4prop", "ci4rate",
    "meth4bin", "meth4inc", "meth4prop", "meth4rate",
    "meth4tau", "meth4tau.ci",
    "adhoc4hakn",
    "meth4bias", "meth4bias.old",
    "meth4incr",
    "text.fixed", "text.w.fixed",
    "major.update", "minor.update")
##
setOption("argslist.internal", argslist.internal)
##
## Set defaults (for internal options)
##
setOption("sm4bin", c("OR", "RD", "RR", "ASD", "DOR"))
setOption("sm4cont", c("MD", "SMD", "ROM"))
setOption("sm4cor", c("ZCOR", "COR"))
setOption("sm4inc", c("IRR", "IRD", "IRSD"))
setOption("sm4mean", c("MRAW", "MLN"))
setOption("sm4prop", c("PLOGIT", "PLN", "PRAW", "PAS", "PFT"))
setOption("sm4rate", c("IR", "IRLN", "IRS", "IRFT"))
##
setOption("ci4cont", c("z", "t"))
setOption("ci4prop", c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
setOption("ci4rate", c("NAsm", "Poisson"))
##
setOption("meth4bin", c("Inverse", "MH", "Peto", "GLMM", "SSW"))
setOption("meth4inc", c("Inverse", "MH", "Cochran", "GLMM"))
setOption("meth4prop", c("Inverse", "GLMM"))
setOption("meth4rate", c("Inverse", "GLMM"))
##
setOption("meth4tau", c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
setOption("meth4tau.ci", c("QP", "BJ", "J", "PL", ""))
setOption("adhoc4hakn", c("", "se", "ci", "iqwig6"))
##
setOption("meth4bias.old", c("rank", "linreg", "mm", "count", "score"))
setOption("meth4bias", c("Begg", "Egger", "Thompson", "Schwarzer",
                         "Harbord", "Peters", "Deeks",
                         "Pustejovsky", "Macaskill"))
##
setOption("meth4incr", c("only0", "if0all", "all"))
##
setOption("major.update", 5)
setOption("minor.update", 5)
##
## List of arguments that can be changed by user
##
argslist <-
  c("level", "level.ma", "common", "random",
    "hakn", "adhoc.hakn", "method.tau", "method.tau.ci", "tau.common",
    "prediction", "level.predict",
    "method.bias",
    "text.common", "text.random", "text.predict",
    "text.w.common", "text.w.random",
    "title", "complab",
    "CIbracket", "CIseparator", "CIlower.blank", "CIupper.blank",
    "print.subgroup.name", "sep.subgroup",
    "keepdata", "warn", "warn.deprecated",
    "backtransf",
    "smbin", "smcont", "smcor", "sminc", "smmean", "smprop", "smrate",
    "incr", "method.incr", "allincr", "addincr",
    "method", "allstudies", "MH.exact",
    "RR.Cochrane", "Q.Cochrane", "model.glmm", "print.CMH",
    "pooledvar", "method.smd", "sd.glass", "exact.smd",
    "method.ci.cont", "method.ci.prop", "method.ci.rate",
    "label.e", "label.c", "label.left", "label.right",
    "layout",
    "test.overall", "test.subgroup", "prediction.subgroup",
    "test.effect.subgroup",
    "digits", "digits.se", "digits.zval", "digits.stat",
    "digits.Q", "digits.tau2", "digits.tau", "digits.H", "digits.I2",
    "digits.prop", "digits.weight",
    "digits.pval", "digits.pval.Q", "digits.forest",
    "scientific.pval", "big.mark", "zero.pval", "JAMA.pval",
    "print.I2", "print.H", "print.Rb",
    "text.tau2", "text.tau", "text.I2", "text.Rb"
    )
##
setOption("argslist", argslist)
##
## General settings
##
setOption("level", 0.95)
setOption("level.ma", 0.95)
setOption("level.comb", 0.95)
setOption("common", TRUE)
setOption("fixed", TRUE)
setOption("comb.fixed", TRUE)
setOption("random", TRUE)
setOption("comb.random", TRUE)
setOption("hakn", FALSE)
setOption("adhoc.hakn", "")
setOption("method.tau", "REML")
setOption("method.tau.ci", NULL)
setOption("tau.common", FALSE)
setOption("prediction", FALSE)
setOption("level.predict", 0.95)
setOption("method.bias", "Egger")
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
setOption("warn", TRUE)
setOption("warn.deprecated", FALSE)
setOption("backtransf", TRUE)
setOption("digits", 4)
setOption("digits.se", 4)
setOption("digits.zval", 2)
setOption("digits.stat", 2)
setOption("digits.Q", 2)
setOption("digits.tau2", 4)
setOption("digits.tau", 4)
setOption("digits.H", 2)
setOption("digits.I2", 1)
setOption("digits.prop", 4)
setOption("digits.weight", 1)
setOption("digits.pval", 4)
setOption("digits.pval.Q", 4)
setOption("scientific.pval", FALSE)
setOption("big.mark", "")
setOption("zero.pval", TRUE)
setOption("JAMA.pval", FALSE)
setOption("print.I2", TRUE)
setOption("print.H", TRUE)
setOption("print.Rb", FALSE)
setOption("text.tau2", "tau^2")
setOption("text.tau", "tau")
setOption("text.I2", "I^2")
setOption("text.Rb", "Rb")
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
setOption("test.overall", FALSE)
setOption("test.effect.subgroup", FALSE)
setOption("digits.forest", 2)


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
