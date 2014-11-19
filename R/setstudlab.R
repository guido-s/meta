setstudlab <- function(x, k){
  ##
  ## Set study labels
  ##
  if (is.null(x))
    x <- seq(along=rep(1, k))
  if (is.factor(x))
    x <- as.character(x)
  ##
  x
}
