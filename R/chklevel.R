chklevel <- function(x, single=TRUE, ci=TRUE, name=NULL){
  ##
  ## Check for levels of confidence interval / contour level
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  if (ci)
    "level for confidence interval (range: 0-1)"
  else
    "contour levels (range: 0-1)"
  ##
  if (!is.numeric(x))
    if (single)
      stop("Argument '", name, "' must be a numeric of length 1.",
           call.=FALSE)
    else
      stop("Argument '", name, "' must be numeric.",
           call.=FALSE)
  ##
  if (single & length(x)!=1)
    stop("Argument '", name, "' must be a numeric of length 1.",
         call.=FALSE)
  ##
  if (any(x <= 0) | any(x >= 1))
    stop("Argument '", name, "' must be a numeric between 0 and 1.",
         call.=FALSE)
  ##
  invisible(NULL)
}
