bylabel <- function(bylab, bylevs, print.byvar, byseparator,
                    big.mark = "") {
  if (print.byvar) {
    if (length(bylab) == 0 || bylab == "")
      res <- format(bylevs, big.mark = big.mark)
    else
      res <- paste0(bylab, byseparator,
                    format(bylevs, big.mark = big.mark))
  }
  else
    res <- format(bylevs, big.mark = big.mark)
  ##
  res
}
