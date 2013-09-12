warnarg <- function(x, y, fun, cl){
  if (x %in% y)
    warning(paste("Argument '", x,
                  "' has been removed from R function ",
                  fun, ".\n  ",
                  "This argument can either be used in R function update.meta or ",
                  cl,
                  ".", sep=""))
}
