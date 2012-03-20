print.rm5 <- function(x, ...){
  
  if (!inherits(x, "rm5"))
    stop("Argument 'x' must be an object of class \"rm5\"")
  
  print.data.frame(x, ...)
  
  invisible(NULL)
}

