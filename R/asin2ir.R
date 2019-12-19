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
