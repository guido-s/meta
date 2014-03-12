isquared <- function(Q, df, level){
  ##
  ## Calculate I-Squared
  ## Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  ##
  tres <- calcH(Q, df, level)
  ##
  func.t <- function(x) (x^2-1)/x^2
  ##
  res <- list(TE=func.t(tres$TE),
              lower=func.t(tres$lower),
              upper=func.t(tres$upper))
  
  res
}
