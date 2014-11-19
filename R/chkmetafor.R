chkmetafor <- function(x, name){
  ##
  ## Check whether R package metafor is installed
  ##
  ##
  if (!(x == "DL" | x == "PM"))
    is.installed.package("metafor",
                         paste("'", name,
                               "' with\n       argument 'method.tau' ",
                               "unequal to 'DL' and 'PM'.",
                               sep=""))
  ##
  invisible(NULL)
}
