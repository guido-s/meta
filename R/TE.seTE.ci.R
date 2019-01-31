TE.seTE.ci <- function(lower, upper, level = 0.95,
                       df = rep_len(NA, length(lower))) {
  
  
  ##
  ## Check arguments
  ##
  if (missing(lower))
    stop("Mandatory argument 'lower' missing.", call. = FALSE)
  if (missing(upper))
    stop("Mandatory argument 'upper' missing.", call. = FALSE)
  ##
  k <- length(lower)
  arg <- "lower"
  chklength(upper, k, arg)
  chklength(df, k, arg)
  ##
  if (any(lower >= upper, na.rm = TRUE))
    stop("Lower limit must be smaller than upper limit.", call. = FALSE)
  ##
  chklevel(level, single = FALSE)
  
  
  ##
  ## Parmar et al. (1998), Stat Med
  ##
  ## Section 4.1 Indirect variance estimation
  ##
  ## Equation (7)
  ##
  varTE <- ifelse(is.na(df),
                  ((upper - lower) /
                   (2 * qnorm((1 - level) / 2, lower.tail = FALSE)))^2,
                  ((upper - lower) /
                   (2 * qt((1 - level) / 2, df = df, lower.tail = FALSE)))^2)
  ##
  seTE <- sqrt(varTE)
  
  
  res <- list(TE = lower + (upper - lower) / 2, seTE = seTE,
              lower = lower, upper = upper,
              level = level, df = df)
  ##
  res
}
