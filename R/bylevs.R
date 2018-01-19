bylevs <- function(x) {
  if (is.factor(x))
    res <- levels(factor(x))
  else
    res <- unique(x)
  res
}
