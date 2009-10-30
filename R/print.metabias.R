print.metabias <- function(x, ...){
  
  if (!inherits(x, "metabias"))
    stop("Argument 'x' must be an object of class \"metabias\"")
  
  
  tl <- options()$width-12
  ##
  if (!is.null(x$title))
    if (x$title!="")
      if (nchar(x$title) <= tl)
        cat("Review:     ", x$title, "\n", sep="")
      else
        cat("Review:     ", substring(x$title, 1, tl-4),
            " ...\n", sep="")
  if (!is.null(x$complab))
    if (x$complab!="")
      if (nchar(x$complab) <= tl)
        cat("Comparison: ", x$complab, "\n", sep="")
      else
        cat("Comparison: ", substring(x$complab, 1, tl-4),
            " ...\n", sep="")
  if (!is.null(x$outclab))
    if (x$outclab!="")
      if (nchar(x$outclab) <= tl)
        cat("Outcome:    ", x$outclab, "\n", sep="")
      else
        cat("Outcome:    ", substring(x$outclab, 1, tl-4),
            " ...\n", sep="")

  class(x) <- "htest"

  print(x)

  invisible(NULL)
}
