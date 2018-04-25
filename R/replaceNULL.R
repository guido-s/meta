replaceNULL <- function(x, replace) {
  if (is.null(x))
    res <- replace
  else
    res <- x
  res
}
