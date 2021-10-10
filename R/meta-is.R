isCol <- function(data, varname) {
  !is.null(data) & varname %in% names(data)
}
is.log.effect <- function(x)
  x %in% c("PLN", "IRLN", "MLN")
is.relative.effect <- function(x)
  x %in% c("HR", "OR", "RR", "IRR", "ROM", "DOR")
is.untransformed <- function(x)
  x %in% c("COR", "PRAW", "IR", "RD", "IRD", "MD", "SMD", "MRAW")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  ##
  res
}
is.zero <- function(x, n = 10)
  abs(x) < n * .Machine$double.eps
