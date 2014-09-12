warnarg <- function(x, y, fun, cl, otherarg){
  if (x %in% y){
    if (!missing(cl))
      warning(paste("Argument '", x,
                    "' has been removed from R function ",
                    fun, ".\n",
                    "This argument can either be used in R function update.meta or ",
                    cl,
                    ".", sep=""),
              call.=FALSE)
    else if (!missing(otherarg))
      warning(paste("Argument '", x,
                    "' has been replaced by argument '",
                    otherarg, "' in R function ",
                    fun, ".\nSee help page of R function ",
                    fun, " for information on the use of the new argument.", sep=""),
              call.=FALSE)
  }
}
