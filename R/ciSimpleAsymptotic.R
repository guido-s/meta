ciSimpleAsymptotic <- function(event, n, level=0.95, correct=FALSE){
  
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
    lower  <- p - z*sqrt(p*q/n)
    upper  <- p + z*sqrt(p*q/n)
    }
  else{
    lower  <- p - (z*sqrt(p*q/n) + 1/(2*n))
    upper  <- p + (z*sqrt(p*q/n) + 1/(2*n))
  }
  ##
  lower[lower<0] <- 0
  upper[upper>1] <- 1
  
  list(n=n, p=p, lower=lower, upper=upper, level=level)
}
