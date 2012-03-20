is.installed.metafor <- function(x){
  if (!any(as.data.frame(installed.packages())$Package=="metafor"))
    if (missing(x))
      stop("Please install library 'metafor' (R command: 'install.packages(\"metafor\")')")
    else
      stop(paste("Please install library 'metafor' before using function ",
                 x, " (R command: 'install.packages(\"metafor\")')", sep=""))
  ##
  invisible(NULL)
}
