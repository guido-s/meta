gs <- function(x) {

  if (missing(x))
    return(NULL)
  
  chkchar(x)
  
  val <- settings.meta()[[x]]
  
  val
}
