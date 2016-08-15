bylabel <- function(bylab, bylevs, print.byvar, byseparator) {
  if (print.byvar) {
    if (length(bylab) == 0 || bylab == "")
      res <- format(bylevs)
    else
      res <- paste(bylab, byseparator, format(bylevs), sep = "")
  }
  else
    res <- format(bylevs)
  ##
  res
}
