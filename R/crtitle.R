crtitle <- function(x) {
  tl <- options()$width - 12
  newline <- FALSE
  ##  
  if (!is.null(x$title)) {
    if (x$title != "") {
      newline <- TRUE
      if (nchar(x$title) <= tl)
        cat(paste0("Review:     ", x$title, "\n"))
      else
        cat(paste0("Review:     ", substring(x$title, 1, tl - 4), " ...\n"))
    }
  }
  if (!is.null(x$complab)) {
    if (x$complab != "") {
      newline <- TRUE
      if (nchar(x$complab) <= tl)
        cat(paste0("Comparison: ", x$complab, "\n"))
      else
        cat(paste0("Comparison: ", substring(x$complab, 1, tl - 4), " ...\n"))
    }
  }
  if (!is.null(x$outclab)) {
    if (x$outclab != "") {
      newline <- TRUE
      if (nchar(x$outclab) <= tl)
        cat(paste0("Outcome:    ", x$outclab, "\n"))
      else
        cat(paste0("Outcome:    ", substring(x$outclab, 1, tl - 4), " ...\n"))
    }
  }
  ##
  if (newline)
    cat("\n")
}
