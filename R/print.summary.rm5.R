print.summary.rm5 <- function(x, ...){
  
  if (!inherits(x, "summary.rm5"))
    stop("Argument 'object' must be an object of class \"summary.rm5\"")
  
  n <- 1
  for (i in 1:length(x)){
    if (n>1)
      cat("\n*****\n\n")
    print(x[[i]])
    n <- n+1
  }
  
  invisible(NULL)
}
