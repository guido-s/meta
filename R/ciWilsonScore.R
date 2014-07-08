ciWilsonScore <- function(event, n, level=0.95, correct=FALSE){
  
  ## Reference:
  ## Newcombe RG. Two-sided confidence intervals for the single
  ## proportion: Comparison of seven methods.
  ## Stat Med 1998, Apr 30;17(8): 857-72
  
  if (level <= 0 | level >= 1)
    stop("no valid level for confidence interval")
  
  p <- event/n
  q <- 1-p
  z <- qnorm(1 - (1-level)/2)

  if (!correct){
    lower  <- (2*n*p + z^2 - z*sqrt(z^2 + 4*n*p*q)) / (2*(n+z^2))
    upper  <- (2*n*p + z^2 + z*sqrt(z^2 + 4*n*p*q)) / (2*(n+z^2))
    }
  else{
    lower <- (2*n*p + z^2 - 1 - z*sqrt(z^2 - 2 - 1/n + 4*p*(n*q+1))) / (2*(n+z^2))
    lower[lower<0] <- 0
    upper <- (2*n*p + z^2 + 1 + z*sqrt(z^2 + 2 - 1/n + 4*p*(n*q-1))) / (2*(n+z^2))
    upper[upper>1] <- 1
  }
  
  list(n=n, p=p, lower=lower, upper=upper, level=level)
}
