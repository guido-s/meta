## Auxiliary functions to (back-)transform effect measures
##
## Package: meta
## Author: Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
asin2ir <- function(x, time = NULL, value = "mean", warn = TRUE) {
  
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
  
  
  ##
  ## Check for (impossible) negative values
  ##
  if (any(sel0, na.rm = TRUE)) {
    if (warn)
      warning("Too small value for ",
              if (length(x) > 1)
                "at least one ",
              if (value == "mean")
                paste0("transformed proportion using Freeman-Tukey double ",
                       "arcsine transformation.\n  Rate set to 0."),
              if (value == "lower")
                paste0("lower confidence limit using Freeman-Tukey double ",
                       "arcsine transformation.",
                       "\n  Lower confidence limit set to 0."),
              if (value == "upper")
                paste0("upper confidence limit using Freeman-Tukey double ",
                       "arcsine transformation.",
                       "\n  Upper confidence limit set to 0."))
  }
  
  
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
asin2p <- function(x, n = NULL, value = "mean", warn = TRUE) {
  
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
  
  
  ##
  ## Check for (impossible) negative values
  ##
  if (any(sel0, na.rm = TRUE)) {
    if (is.null(n)) {
      if (warn)
        warning("Negative value for ",
                if (length(x) > 1)
                  "at least one ",
                if (value == "mean")
                  paste0("transformed proportion using arcsine transformation.",
                         "\n  Proportion set to 0."),
                if (value == "lower")
                  paste0("lower confidence limit using arcsine transformation.",
                         "\n  Lower confidence limit set to 0."),
                if (value == "upper")
                  paste0("upper confidence limit using arcsine transformation.",
                         "\n  Upper confidence limit set to 0."))
    }
    else {
      if (warn)
        warning("Too small value for ",
                if (length(x) > 1)
                  "at least one ",
                if (value == "mean")
                  paste0("transformed proportion using Freeman-Tukey double ",
                         "arcsine transformation.\n  Proportion set to 0."),
                if (value == "lower")
                  paste0("lower confidence limit using Freeman-Tukey double ",
                         "arcsine transformation.",
                         "\n  Lower confidence limit set to 0."),
                if (value == "upper")
                  paste0("upper confidence limit using Freeman-Tukey double ",
                         "arcsine transformation.",
                         "\n  Upper confidence limit set to 0."))
    }
  }
  
  ##
  ## Check for (impossible) large values
  ##
  if (any(sel1, na.rm = TRUE)) {
    if (is.null(n)) {
      if (warn)
        warning("Too large value for ",
                if (length(x) > 1)
                  "at least one ",
                if (value == "mean")
                  paste0("transformed proportion using arcsine transformation.",
                         "\n  Proportion set to 1."),
                if (value == "lower")
                  paste0("lower confidence limit using arcsine transformation.",
                         "\n  Lower confidence limit set to 1."),
                if (value == "upper")
                  paste0("upper confidence limit using arcsine transformation.",
                         "\n  Upper confidence limit set to 1."))
    }
    else {
      if (warn)
        warning("Too large value for ",
                if (length(x) > 1)
                  "at least one ",
                if (value == "mean")
                  paste0("transformed proportion using Freeman-Tukey double ",
                         "arcsine transformation.\n  Proportion set to 1."),
                if (value == "lower")
                  paste0("lower confidence limit using Freeman-Tukey double ",
                         "arcsine transformation.",
                         "\n  Lower confidence limit set to 1."),
                if (value == "upper")
                  paste0("upper confidence limit using Freeman-Tukey double ",
                         "arcsine transformation.",
                         "\n  Upper confidence limit set to 1."))
    }
  }
  
  
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
logit2p <- function(x) {
  res <- 1 / (1 + exp(-x))
  res
}
p2logit <- function(x) {
  res <- log(x) - log(1 - x)
  res
}
z2cor <- function(x) {
  res <- (exp(2 * x) - 1) / (exp(2 * x) + 1)
  res
}
backtransf <- function(x, sm, value, n, warn = FALSE) {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  if (is.relative.effect(sm) | is.log.effect(sm))
    res <- exp(x)
  ##
  else if (sm == "ZCOR")
    res <- z2cor(x)
  ##
  else if (sm == "PLOGIT")
    res <- logit2p(x)
  ##
  else if (sm == "PAS")
    res <- asin2p(x, value = value, warn = warn)
  ##
  else if (sm == "PFT")
    res <- asin2p(x, n, value = value, warn = warn)
  ##
  else if (sm == "IRS")
    res <- x^2
  ##
  else if (sm == "IRFT")
    res <- asin2ir(x, n, value = value, warn = warn)
  ##
  else
    res <- x

  if (sm == "PRAW") {
    sel0 <- res[!is.na(res)] < 0 & value == "lower"
    sel1 <- res[!is.na(res)] > 1 & value == "upper"
    ##
    if (warn & any(sel0 | sel1, na.rm = TRUE))
      warning("Negative value for ",
              if (length(x) > 1)
                "at least one ",
              if (value == "lower")
                paste0("lower confidence limit of raw proportions.",
                       "\n  Lower confidence limit set to 0."),
              if (value == "upper")
                paste0("upper confidence limit of raw proportions.",
                       "\n  Upper confidence limit set to 1."))
    if (any(sel0, na.rm = TRUE) & value == "lower")
      res[sel0] <- 0
    else if (any(sel1, na.rm = TRUE) & value == "upper")
      res[sel1] <- 1
  }
  
  if (sm == "PLN") {
    sel0 <- res[!is.na(res)] < 0 & value == "lower"
    sel1 <- res[!is.na(res)] > 1 & value == "upper"
    ##
    if (warn & any(sel0 | sel1, na.rm = TRUE))
      warning("Negative value for ",
              if (length(x) > 1)
                "at least one ",
              if (value == "lower")
                paste0("lower confidence limit using log transformation for ",
                       "proportions.\n  Lower confidence limit set to 0."),
              if (value == "upper")
                paste0("upper confidence limit using log transformation for ",
                       "proportions.\n  Upper confidence limit set to 1."))
    if (any(sel0, na.rm = TRUE) & value == "lower")
      res[sel0] <- 0
    else if (any(sel1, na.rm = TRUE) & value == "upper")
      res[sel1] <- 1
  }
  
  if (sm == "IR") {
    sel0 <- res[!is.na(res)] < 0 & value == "lower"
    ##
    if (warn & any(sel0, na.rm = TRUE))
      warning("Negative value for ",
              if (length(x) > 1)
                "at least one ",
              if (value == "lower")
                paste0("lower confidence limit of incidence rates.",
                       "\n  Lower confidence limit set to 0."))
    if (any(sel0, na.rm = TRUE) & value == "lower")
      res[sel0] <- 0
  }

  res
}
