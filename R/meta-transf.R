## Auxiliary functions to (back-)transform effect measures
##
## Package: meta
## Author: Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
## License: GPL (>= 2)
##

asin2ir <- function(x, time = NULL, value = "mean") {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  
  ##
  ## Calculate possible minimum for each transformation
  ##
  minimum <- 0.5 * (sqrt(0 / time) + sqrt((0 + 1) / time))
  ##
  sel0 <- x < minimum
  
  
  res <- rep(NA, length(x))
  ##
  sel <- !sel0
  sel <- !is.na(sel) & sel
  ##
  res[sel0] <- 0
  ##
  ## Back transformation of Freeman-Tukey double arcsine transformation:
  ##
  res[sel] <- (1 / time[sel] - 8 * x[sel]^2 + 16 * time[sel] * x[sel]^4) /
    (16 * x[sel]^2 * time[sel])
  ##
  res[res < 0] <- 0
  
  
  res
}


asin2p <- function(x, n = NULL, value = "mean") {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  
  ##
  ## Calculate possible minimum and maximum
  ## for each transformation
  ##
  if (is.null(n)) {
    minimum <- asin(sqrt(0))
    maximum <- asin(sqrt(1))
  }
  else {
    minimum <- 0.5 * (asin(sqrt(0 / (n + 1))) + asin(sqrt((0 + 1) / (n + 1))))
    maximum <- 0.5 * (asin(sqrt(n / (n + 1))) + asin(sqrt((n + 1) / (n + 1))))
  }
  ##
  sel0 <- x < minimum
  sel1 <- x > maximum
  
  
  res <- rep(NA, length(x))
  ##
  sel <- !(sel0 | sel1)
  sel <- !is.na(sel) & sel
  ##
  res[sel0] <- 0
  res[sel1] <- 1
  ##
  if (is.null(n)) {
    ##
    ## Back transformation of arcsine transformation:
    ##
    res[sel] <- sin(x[sel])^2
  }
  else {
    ##
    ## Back transformation of Freeman-Tukey double arcsine transformation:
    ##
    res[sel] <- 0.5 * (1 - sign(cos(2 * x[sel])) *
                       sqrt(1 - (sin(2 * x[sel]) +
                                 (sin(2 * x[sel]) -
                                  1 / sin(2 * x[sel])) / n[sel])^2))
  }
  res
}


logit2p <- function(x)
  1 / (1 + exp(-x))


logVR2VE <- function(x)
  100 * (1 - exp(x))


z2cor <- function(x)
  tanh(x)


backtransf <- function(x, sm, value, n, func = NULL, args = NULL) {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  if (!is.null(func))
    res <- do.call(func, c(list(x), args))
  ##
  else if (is_relative_effect(sm) | is_log_effect(sm))
    res <- exp(x)
  ##
  else if (sm == "ZCOR")
    res <- z2cor(x)
  ##
  else if (sm == "PLOGIT")
    res <- logit2p(x)
  ##
  else if (sm == "PAS")
    res <- asin2p(x, value = value)
  ##
  else if (sm == "PFT")
    res <- asin2p(x, n, value = value)
  ##
  else if (sm == "IRS")
    res <- x^2
  ##
  else if (sm == "IRFT")
    res <- asin2ir(x, n, value = value)
  ##
  else if (sm == "VE")
    res <- logVR2VE(x)
  ##
  else
    res <- x

  if (sm == "PRAW") {
    sel0 <- res[!is.na(res)] < 0 & value == "lower"
    sel1 <- res[!is.na(res)] > 1 & value == "upper"
    ##
    if (any(sel0, na.rm = TRUE) & value == "lower")
      res[sel0] <- 0
    else if (any(sel1, na.rm = TRUE) & value == "upper")
      res[sel1] <- 1
  }
  
  if (sm == "PLN") {
    sel0 <- res[!is.na(res)] < 0 & value == "lower"
    sel1 <- res[!is.na(res)] > 1 & value == "upper"
    ##
    if (any(sel0, na.rm = TRUE) & value == "lower")
      res[sel0] <- 0
    else if (any(sel1, na.rm = TRUE) & value == "upper")
      res[sel1] <- 1
  }
  
  if (sm == "IR") {
    sel0 <- res[!is.na(res)] < 0 & value == "lower"
    ##
    if (any(sel0, na.rm = TRUE) & value == "lower")
      res[sel0] <- 0
  }

  res
}


cor2z <- function(x)
  0.5 * log((1 + x) / (1 - x))


p2asin <- function(x)
  asin(sqrt(x))


p2logit <- function(x)
  qlogis(x)


VE2logVR <- function(x)
  log(1 - x / 100)


transf <- function(x, sm, func = NULL, args = NULL) {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)

  if (!is.null(func))
    res <- do.call(func, c(list(x), args))
  ##
  else if (is_relative_effect(sm) | is_log_effect(sm))
    res <- log(x)
  ##
  else if (sm == "ZCOR")
    res <- cor2z(x)
  ##
  else if (sm == "PLOGIT")
    res <- p2logit(x)
  ##
  else if (sm == "PAS")
    res <- p2asin(x)
  ##
  else if (sm == "IRS")
    res <- sqrt(x)
  ##
  else if (sm == "VE")
    res <- VE2logVR(x)
  ##
  else
    res <- x
  
  res
}

