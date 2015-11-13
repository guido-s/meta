chknumeric <- function(x, min, max, zero=FALSE, single=FALSE,
                       name=NULL){
  ##
  ## Check numeric variable
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(invisible(NULL))
  ##
  if (!(is.numeric(x)))
    stop("Non-numeric value for argument '", name, "'.",
         call.=FALSE)
  ##
  if (single & length(x)!=1)
    stop("Argument '", name, "' must be a numeric of length 1.",
         call.=FALSE)
  ##
  if (!missing(min) & missing(max)){
    if (zero & min==0 & any(x <= min))
      stop("Argument '", name, "' must be positive.",
           call.=FALSE)
    else if (any(x < min))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call.=FALSE)
  }
  ##
  if (missing(min) & !missing(max)){
    if (zero & max==0 & any(x >= max))
      stop("Argument '", name, "' must be negative.",
           call.=FALSE)
    else if (any(x > max))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call.=FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min) | any(x > max)))
      stop("Argument '", name, "' must be between ",
           min, " and ", max, ".", call.=FALSE)
  ##
  invisible(NULL)
}
