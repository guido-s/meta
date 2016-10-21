is.installed.package <- function(pkg, func, argument, value,
                                 chksettings = FALSE, stop = TRUE,
                                 version = NULL) {
  
  pkginstalled <- requireNamespace(pkg, quietly = TRUE)
  
  oldpkg <- pkginstalled && !is.null(version) && packageVersion(pkg) < version
  
  if (stop & (!pkginstalled) | oldpkg) {
    
    if (oldpkg) {
      oldmsg <- paste("Library '", pkg, "' is too old. ", sep = "")
      oldinst <- "re"
    }
    else {
      oldmsg <- ""
      oldinst <- ""
    }
    
    if (chksettings)
      warning(paste(oldmsg,
                    "Argument '", argument, "' not changed.",
                    "\n  ",
                    "Please ", oldinst,
                    "install library '", pkg,
                    "' in order to use argument '", argument, " = \"", value, "\"",
                    "'\n  ",
                    "(R command: 'install.packages(\"", pkg, "\")').",
                    sep = ""))
    else
      if (missing(func))
        stop(paste(oldmsg,
                   "Please ", oldinst,
                   "install library '", pkg,
                   "'\n       ",
                   "(R command: 'install.packages(\"", pkg, "\")').",
                   sep = ""),
             call. = FALSE)
      else
        stop(paste(oldmsg,
                   "Please ", oldinst,
                   "install library '", pkg,
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
