replaceNULL <- function(x, replace = NA) {
  if (is.null(x))
    res <- replace
  else
    res <- x
  res
}
