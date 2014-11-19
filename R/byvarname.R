byvarname <- function(x){
  ##
  ## Determine name of byvar
  ##
  res <- as.character(x)
  if (length(res)>1 & res[1]=="$")
      res <- res[length(res)]
    if (length(res)>1)
      res <- "byvar"
  ##
  res
}
