npn <- function(x){
  ##
  ## Check for non-positive values in vector
  ##
  selNA <- is.na(x)
  res <- selNA
  if (sum(!selNA) > 0)
    x[!selNA] <- x[!selNA] <= 0
  ##
  res
}
