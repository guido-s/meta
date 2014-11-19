print.rm5 <- function(x, ...){
  
  
  ##
  ##
  ## (1) Check for rm5 object
  ##
  ##
  chkclass(x, "rm5")
  
  
  print.data.frame(x, ...)
  
  
  invisible(NULL)
}
