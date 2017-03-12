updateversion <- function(x) {
  ##
  ## Update older meta objects
  ##
  if (!(!is.null(x$version) &&
        as.numeric(unlist(strsplit(x$version, "-"))[1]) >= 3.2))
    x <- update(x, warn = FALSE)
  ##
  x
}
