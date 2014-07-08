ciAgrestiCoull <- function(event, n, level=0.95){
  
  if (level <= 0 | level >= 1)
    stop("no valid level for confidence interval")
  
  z <- qnorm(1 - (1-level)/2)
  ##
  n <- n + z^2
  p <- 1/n * (event + 0.5*z^2)
  
  lower  <- p - z*sqrt(1/n*p*(1-p))
  upper  <- p + z*sqrt(1/n*p*(1-p))
  
  list(n=n, p=p, lower=lower, upper=upper, level=level)
}
