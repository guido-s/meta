updateversion <- function(x) {
  ##
  ## Update older meta objects
  ##
  if (is.null(x$version) || x$version != packageDescription("meta")$Version)
    x <- update(x, warn = FALSE)
  ##
  x
}
