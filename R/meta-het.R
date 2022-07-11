## Auxiliary functions to calculate heterogeneity measures
##
## Package: meta
## Author: Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
Rb <- function(seTE, seTE.random, tau2, Q, df.Q, level) {
  ##
  ## Calculate Rb
  ## Crippa et al. (2016), Statistics in Medicine, 35 3661-75.
  ##
  
  
  if (is.na(df.Q) || df.Q == 0)
    return(list(TE = NaN, lower = NaN, upper = NaN))
  ##
  k <- df.Q + 1
  varTE <- seTE^2
  w.common <- 1 / varTE
  
  
  ##
  ## Equation (4) in Crippa et al. (2016)
  ##
  ## Rb <- 1 / k * sum(tau2 / (varTE + tau2))
  Rb <- tau2 / (k * seTE.random^2)
  
  
  ##
  ## Appendix D in Crippa et al. (2016)
  ##
  S <- function(n, w) sum(w^n)
  ##
  a <- varTE * (S(1, w.common) - S(2, w.common) / S(1, w.common))
  ##
  seQ <- sqrt(2 * (k - 1) +
              4 * (S(1, w.common) - S(2, w.common) / S(1, w.common)) * tau2 +
              2 * (S(2, w.common) - 2 * S(3, w.common) / S(1, w.common) +
                   S(2, w.common)^2 / S(1, w.common)^2) * tau2^2
              )
  ##  
  seRb <- sqrt((1 / k * sum(a / (Q + a - (k - 1))^2))^2 * seQ^2)
  
  cint <- ci(Rb, seRb, level)
  ##
  cint$lower <- ifelse(cint$lower < 0, 0, cint$lower)
  cint$upper <- ifelse(cint$upper > 1, 1, cint$upper)
  
  res <- list(TE = cint$TE, lower = cint$lower, upper = cint$upper)
  ##
  res
}
calcH <- function(Q, df, level) {
  ##
  ## Calculate H
  ## Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  ##
  k <- df + 1
  ##
  if (!is.na(k)) {
    if (k > 1)
      H <- sqrt(Q / (k - 1))
    else
      H <- NA
    ##
    selogH <- ifelse(Q > k,
                     ifelse(k >= 2,
                            0.5 * (log(Q) - log(k - 1)) /
                            (sqrt(2 * Q) - sqrt(2 * k - 3)),
                            NA),
                     ifelse(k > 2,
                            sqrt(1 / (2 * (k - 2)) * (1 - 1 / (3 * (k - 2)^2))),
                            NA))
  }
  else {
    H <- NA
    selogH <- NA
  }
  ##
  tres <- ci(log(max(c(H, 1))), selogH, level)
  ##
  res <- list(TE = exp(tres$TE),
              lower = max(exp(tres$lower), 1),
              upper = max(exp(tres$upper), 1))
  
  res
}
isquared <- function(Q, df, level) {
  ##
  ## Calculate I-Squared
  ## Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  ##
  tres <- calcH(Q, df, level)
  ##
  func.t <- function(x) (x^2 - 1) / x^2
  ##
  res <- list(TE = func.t(tres$TE),
              lower = func.t(tres$lower),
              upper = func.t(tres$upper))
  
  res
}
pvalQ <- function(Q, df, lower.tail = FALSE) {
  ##
  if (length(df) == 1 & length(Q) > 1)
    df <- rep(df, length(Q))
  else if (length(df) > 1 & length(Q) == 1)
    Q <- rep(Q, length(df))
  else if (length(df) > 1 & length(Q) > 1 & length(df) != length(Q))
    stop("Length of arguments 'Q' and 'df' do not match.")
  ##
  res <- ifelse(is.na(df) | df < 1,
                NA,
                pchisq(Q, df, lower.tail = lower.tail))
  ##
  res
}
