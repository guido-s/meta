is.installed.package <- function(pkg, func, argument, value,
                                 chksettings = FALSE, stop = TRUE) {

  pkginstalled <- any(as.data.frame(installed.packages())$Package == pkg)
  
  if (stop & !pkginstalled) {

    if (chksettings)
      warning(paste("Argument '", argument, "' not changed as necessary R library is missing.",
                   "\n  ",
                    "Please install library '", pkg,
                    "' in order to use argument '", argument, " = \"", value, "\"",
                    "'\n  ",
                    "(R command: 'install.packages(\"", pkg, "\")').",
                    sep = ""))
    else
      if (missing(func))
        stop(paste("Please install library '", pkg,
                   "'\n       ",
                   "(R command: 'install.packages(\"", pkg, "\")').",
                   sep = ""),
             call. = FALSE)
      else
        stop(paste("Please install library '", pkg,
                   "' before using function '", func,
                   if (missing(argument))
                     "'."
                   else
                     paste("' with argument '", argument,
                           ifelse(missing(value), "'.",
                                  paste(" ", value, "'.", sep = "")),
                           sep = ""),
                   "\n       ",
                   "(R command: 'install.packages(\"", pkg, "\")').",
                   sep = ""),
             call. = FALSE)
    
  }
  
  invisible(pkginstalled)
}
