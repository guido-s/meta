ci <- function(TE, seTE, level=0.95){

  if (level <= 0 | level >= 1)
    stop("no valid level for confidence interval")

  alpha <- 1-level
  
  lower  <- TE - qnorm(1-alpha/2)*seTE
  upper  <- TE + qnorm(1-alpha/2)*seTE
  zscore <- TE/seTE
  pval   <- 2*(1-pnorm(abs(zscore)))

  list(TE=TE, seTE=seTE,
       lower=lower, upper=upper,
       z=zscore, p=pval, level=level)
}
