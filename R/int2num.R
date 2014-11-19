int2num <- function(x){
  ##
  ## Convert integer to numeric
  ##
  if (is.integer(x))
    res <- as.numeric(x)
  else
    res <- x
  ##
  res
}
