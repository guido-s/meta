.onLoad <- function(libname, pkgname)
{
   library.dynam("meta", pkgname, libname)
}

.onUnload <- function(libpath)
{   
   library.dynam.unload("meta", libpath)
}

.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'meta' package (version ",
               utils::packageDescription("meta")$Version,
               ").", sep="")
  packageStartupMessage(msg)
}


## The following R code is based on the file snowfall-internal.R from
## R package snowfall (Maintainer: Jochen Knaus <jo@imbi.uni-freiburg.de>)


##*****************************************************************************
## Helpers for managing the internal variables in the package namespace without
## awake the R CMD check for later R versions (which basically blaims many
## global assignings).
##
## The given solution has an advantage: only writing is affected. Reading of the
## objects can remain the same (thanks to Uwe Ligges for the tipp):
##   reading:  .metaOptions$CIbracket
##   writing:  setOption("CIbracket", "(")
##*****************************************************************************

##*****************************************************************************
## Set an option in the meta option list.
## (Basically this is the setting of a list entry).
## key - character: object name
## val - object (everything is allowed, even NULL)
##*****************************************************************************
setOption <- function( key=NULL, val=NULL ) {
  if( !is.null(key) && is.character( key ) ) {
    option <- getVar( ".metaOptions" )   ## Get from NS
    option[[key]] <- val
    setVar( ".metaOptions", option )     ## Write to NS

    return( invisible( TRUE ) )
  }

  stop( "key or val is NULL or key no string." )
}

##*****************************************************************************
## Get a specific variable from the meta namespace.
## var - character: object name
##*****************************************************************************
getVar <- function( var=NULL ) {
  if( !is.null( var ) && is.character( var ) ) {
    tmp <- try( getFromNamespace( var, "meta" ) )

    if( inherits( tmp, "try-error" ) )
      stop( paste( "Object", var, "not found in package" ) )

    return( tmp )
  }

  stop( "var is NULL or not a string." )
}

##*****************************************************************************
## Write a specific variable to the meta namespace.
## var - character: object name
## arg - object (NULL allowed)
##*****************************************************************************
setVar <- function( var=NULL, arg=NULL ) {
  if( !is.null( var ) && is.character( var ) ) {
    assignInNamespace( var, arg, "meta" )

    return( invisible( TRUE ) )
  }

  stop( "var is NULL or no character" );
}


.metaOptions <- list()
##
setOption("CIbracket", "[")
setOption("CIseparator", "; ")
