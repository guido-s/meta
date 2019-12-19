.onAttach <- function(libname, pkgname) {
  msg <- paste0("Loading 'meta' package (version ",
                utils::packageDescription("meta")$Version,
                ").",
                "\nType 'help(meta)' for a brief overview.")
  packageStartupMessage(msg)
}


isCol <- function(data, varname) {
  !is.null(data) & varname %in% names(data)
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
##   reading:  .settings$CIbracket
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


pasteCI <- function(lower, upper, digits, big.mark,
                    sign.lower = "", sign.upper = "", text.NA = "NA",
                    unit = "")
  paste0(" ",
         formatCI(paste0(sign.lower,
                         formatN(lower, digits, big.mark = big.mark,
                                 text.NA = text.NA), unit),
                  paste0(sign.upper,
                         formatN(upper, digits, big.mark = big.mark,
                                 text.NA = text.NA), unit)))


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  ##
  res
}


.settings <- list()
##
## Set defaults (for internal options)
##
setOption("metafor", "1.9.9")
##
setOption("sm4bin", c("OR", "RD", "RR", "ASD"))
setOption("sm4cont", c("MD", "SMD", "ROM"))
setOption("sm4cor", c("ZCOR", "COR"))
setOption("sm4inc", c("IRR", "IRD"))
setOption("sm4mean", c("MRAW", "MLN"))
setOption("sm4prop", c("PLOGIT", "PLN", "PRAW", "PAS", "PFT"))
setOption("sm4rate", c("IR", "IRLN", "IRS", "IRFT"))
##
setOption("ci4prop", c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
##
setOption("meth4tau", c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
setOption("meth4tau.ci", c("QP", "BJ", "J", ""))
##
## List of arguments that can be changed by user
##
argslist <- c("level", "level.comb", "comb.fixed", "comb.random",
              "hakn", "method.tau", "tau.common",
              "prediction", "level.predict",
              "method.bias", "title", "complab", "CIbracket", "CIseparator",
              "print.byvar", "byseparator", "keepdata", "warn",
              "backtransf",
              "smbin", "smcont", "smcor", "sminc", "smmean", "smprop", "smrate",
              "incr", "allincr", "addincr",
              "method", "allstudies", "MH.exact",
              "RR.Cochrane", "Q.Cochrane", "model.glmm", "print.CMH",
              "pooledvar", "method.smd", "sd.glass", "exact.smd", "method.ci",
              "label.e", "label.c", "label.left", "label.right",
              "layout",
              "test.overall", "test.subgroup", "test.effect.subgroup",
              "digits", "digits.se", "digits.zval",
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
setOption("level.comb", 0.95)
setOption("comb.fixed", TRUE)
setOption("comb.random", TRUE)
setOption("hakn", FALSE)
setOption("method.tau", "DL")
setOption("tau.common", FALSE)
setOption("prediction", FALSE)
setOption("level.predict", 0.95)
setOption("method.bias", "linreg")
setOption("title", "")
setOption("complab", "")
setOption("CIbracket", "[")
setOption("CIseparator", "; ")
setOption("print.byvar", TRUE)
setOption("byseparator", " = ")
setOption("keepdata", TRUE)
setOption("warn", TRUE)
setOption("backtransf", TRUE)
setOption("digits", 4)
setOption("digits.se", 4)
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
setOption("exact.smd", FALSE)
##
## Additional setting for R function metaprop
##
setOption("method.ci", "CP")
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
setOption("test.subgroup", FALSE)
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
