is.installed.package <- function(pkg, text){
  if (!any(as.data.frame(installed.packages())$Package==pkg))
    if (missing(text))
      stop(paste("Please install library '", pkg,
                 "'\n       ",
                 "(R command: 'install.packages(\"", pkg, "\")').",
                 sep=""),
           call.=FALSE)
    else
      stop(paste("Please install library '", pkg,
                 "' before using function ", text,
                 "\n       ",
                 "(R command: 'install.packages(\"", pkg, "\")').",
                 sep=""),
           call.=FALSE)
  ##
  invisible(NULL)
}
