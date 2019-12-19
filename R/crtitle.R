crtitle <- function(x) {
  tl <- options()$width - 12
  ##
  if (!is.null(x$title))
    if (x$title != "")
      if (nchar(x$title) <= tl)
        cat(paste0("Review:     ", x$title, "\n"))
      else
        cat(paste0("Review:     ", substring(x$title, 1, tl - 4), " ...\n"))
  if (!is.null(x$complab))
    if (x$complab != "")
      if (nchar(x$complab) <= tl)
        cat(paste0("Comparison: ", x$complab, "\n"))
      else
        cat(paste0("Comparison: ", substring(x$complab, 1, tl - 4), " ...\n"))
  if (!is.null(x$outclab))
    if (x$outclab != "")
      if (nchar(x$outclab) <= tl)
        cat(paste0("Outcome:    ", x$outclab, "\n\n"))
      else
        cat(paste0("Outcome:    ", substring(x$outclab, 1, tl - 4), " ...\n\n"))
}
