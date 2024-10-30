## Auxiliary functions to check inputs
##
## Package: meta
## Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
## License: GPL (>= 2)
##

isCol <- function(data, varname) {
  !is.null(data) & varname %in% names(data)
}

is_log_effect <- function(x)
  x %in% c("PLN", "IRLN", "MLN")

is_relative_effect <- function(x)
  x %in% c("HR", "OR", "RR", "IRR", "ROM", "DOR")

is_untransformed <- function(x)
  x %in% c("COR", "PRAW", "IR", "RD", "IRD", "MD", "SMD", "MRAW")

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  ##
  res
}

is_zero <- function(x, n = 10)
  abs(x) < n * .Machine$double.eps

is_cor <- function(x)
  x %in% gs("sm4cor")

is_mean <- function(x)
  x %in% gs("sm4mean")

is_prop <- function(x)
  x %in% gs("sm4prop")

is_rate <- function(x)
  x %in% gs("sm4rate")
