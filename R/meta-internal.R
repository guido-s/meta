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
