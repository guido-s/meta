estimate.missing <- function(TE, TE.sum, type){
  ##
  ## 1. Centre around mean
  ##
  TE.c <- TE - TE.sum
  n <- length(TE.c)
  ##
  ## 2. Rank absolute values of centred values
  ##
  r.star <- rank(abs(TE.c))*sign(TE.c)
  ##
  if (type=="L"){
    ##
    ## 3. Sum the positive ranks only
    ##
    S.rank <- sum(r.star[r.star>0])
    ##
    ## 4. Estimate for L0
    ##
    res0 <- (4*S.rank - n*(n+1))/(2*n-1)
    res0.plus <- max(0, res0 + 0.5) %/% 1
  }
  if (type=="R"){
    ##
    ## 5. Estimate for R0
    ##
    res0 <- n - abs(min(r.star)) - 1.5
    res0.plus <- max(0, res0 + 0.5) %/% 1
  }
  ##
  res <- list(res0=res0, res0.plus=res0.plus)
  res
}
