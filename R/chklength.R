chklength <- function(x, k.all, fun, text, name = NULL) {
  ##
  ## Check length of vector
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length(x) != k.all) {
    funcs <- c("metabin", "metacont", "metacor",
               "metagen", "metainc", "metamean",
               "metaprop", "metarate",
               "funnel", "forest.meta")
    args <- c("event.e", "n.e", "cor",
              "TE", "event.e", "n",
              "event", "event",
              "TE", "TE")
    ##
    idx <- charmatch(fun, funcs, nomatch = NA)
    if (!is.na(idx))
      argname <- args[idx]
    else
      argname <- fun
    ##
    if (missing(text))
      stop("Arguments '", argname, "' and '", name,
           "' must have the same length.",
           call. = FALSE)
    else
      stop(text, call. = FALSE)
  }
  ##
  invisible(NULL)
}
