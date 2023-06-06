addvars2long <- function(x) {
  levs <- NULL
  fac <- is.factor(x$var1) | is.factor(x$var2)
  ##
  if (!(is.numeric(x$var1) & is.numeric(x$var2))) {
    levs <- unique(c(x$var1, x$var2))
    x$var1 <- as.numeric(factor(x$var1, levels = levs))
    x$var2 <- as.numeric(factor(x$var2, levels = levs))
  }
  ##
  res <-
    to.long("RD", ai = x$var1, n1i = x$var1, ci = x$var2, n2i = x$var2)$out1
  ##
  if (!is.null(levs)) {
    res <- factor(res, levels = seq_along(levs), labels = levs)
    ##
    if (!fac)
      res <- as.character(res)
  }
  ##
  res
}
