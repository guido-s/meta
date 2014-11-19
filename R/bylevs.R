bylevs <- function(x){
  if (is.factor(x))
    res <- levels(x)
  else
    res <- unique(x)
  res
}
