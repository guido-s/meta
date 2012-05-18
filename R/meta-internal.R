#".First.lib" <- function(lib,pkg) {   
#   library.dynam("meta", pkg, lib)
#}

.onLoad <- function(libname, pkgname)
{
   packageStartupMessage("Loading 'meta' package (version ", utils::packageDescription("meta")$Version, ").")
   library.dynam("meta", pkgname, libname)
}

.onUnload <- function(libpath)
{   
   library.dynam.unload("meta", libpath)
}

