setstudlab <- function(x, k) {
  ##
  ## Set study labels
  ##
  if (is.null(x)) {
    if (k == 1)
      x <- ""
    else
      x <- seq(along = rep(1, k))
  }
  ##
  if (is.factor(x))
    x <- as.character(x)
  ##
  x
}
