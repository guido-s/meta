is.installed.package <- function(pkg, func, argument, value,
                                 chksettings = FALSE, stop = TRUE,
                                 version = NULL) {
  
  pkginstalled <- requireNamespace(pkg, quietly = TRUE)
  
  oldpkg <- pkginstalled && !is.null(version) && packageVersion(pkg) < version
  
  if (stop & (!pkginstalled) | oldpkg) {
    
    if (oldpkg) {
      oldmsg <- paste0("Library '", pkg, "' is too old. ")
      oldinst <- "re"
    }
    else {
      oldmsg <- ""
      oldinst <- ""
    }
    
    if (chksettings)
      warning(oldmsg,
              "Argument '", argument, "' not changed.\n  ",
              "Please ", oldinst, "install library '", pkg,
              "' in order to use argument '", argument, " = \"",
              value, "\"'\n  ",
              "(R command: 'install.packages(\"", pkg, "\")').")
    else
      if (missing(func))
        stop(oldmsg,
             "Please ", oldinst, "install library '", pkg,
             "'\n       ",
             "(R command: 'install.packages(\"", pkg, "\")').",
             call. = FALSE)
      else
        stop(oldmsg,
             "Please ", oldinst, "install library '", pkg,
             "' before using function '", func,
             if (missing(argument))
               "'."
             else
               paste0("' with argument '", argument,
                      ifelse(missing(value), "'.",
                             paste0(" ", value, "'."))),
             "\n       ",
             "(R command: 'install.packages(\"", pkg, "\")').",
             call. = FALSE)
  }
  
  invisible(pkginstalled)
}
