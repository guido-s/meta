.onLoad <- function(libname, pkgname){
   library.dynam("meta", pkgname, libname)
}

.onUnload <- function(libpath){
  library.dynam.unload("meta", libpath)
}

.onAttach <- function (libname, pkgname){
  msg <- paste("Loading 'meta' package (version ",
               utils::packageDescription("meta")$Version,
               ").", sep="")
  packageStartupMessage(msg)
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
setOption <- function(key=NULL, val=NULL){
  if(!is.null(key) && is.character(key)){
    option <- getVar(".settings")   ## Get from NS
    option[[key]] <- val
    setVar(".settings", option)     ## Write to NS
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
getVar <- function(var=NULL) {
  if(!is.null(var) && is.character(var)){
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
setVar <- function( var=NULL, arg=NULL ) {
  if( !is.null( var ) && is.character( var ) ) {
    assignInNamespace( var, arg, "meta" )

    return( invisible( TRUE ) )
  }

  stop( "var is NULL or no character" );
}


.settings <- list()
##
## Set defaults
##
setOption("CIbracket", "[")
setOption("CIseparator", "; ")
##
argslist <- c("level", "level.comb", "comb.fixed", "comb.random",
              "hakn", "method.tau", "tau.common",
              "prediction", "level.predict",
              "method.bias", "title", "complab",
              "print.byvar", "keepdata", "warn",
              "backtransf",
              "smbin", "smcont", "smcor", "sminc", "smprop",
              "incr", "allincr", "addincr",
              "method", "allstudies", "MH.exact",
              "RR.cochrane", "print.CMH",
              "pooledvar", "method.smd", "sd.glass", "exact.smd", "method.ci",
              "label.e", "label.c", "label.left", "label.right",
              "test.overall", "test.subgroup")
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
setOption("print.byvar", TRUE)
setOption("keepdata", TRUE)
setOption("warn", TRUE)
setOption("backtransf", TRUE)
##
## Default summary measure
##
setOption("smbin", "RR")
setOption("smcont", "MD")
setOption("smcor", "ZCOR")
setOption("sminc", "IRR")
setOption("smprop", "PLOGIT")
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
setOption("RR.cochrane", FALSE)
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
setOption("test.overall", FALSE)
setOption("test.subgroup", FALSE)
