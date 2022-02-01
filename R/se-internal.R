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
  chklevel(level, length = 0)
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
  ##
  res <- list(TE = lower + (upper - lower) / 2, seTE = seTE,
              lower = lower, upper = upper,
              level = level, df = df)
  ##
  res
}
seTE.pval <- function(TE, pval, df = rep_len(NA, length(TE))) {
  ##
  ## Check arguments
  ##
  if (missing(TE))
    stop("Mandatory argument 'TE' missing.", call. = FALSE)
  if (missing(pval))
    stop("Mandatory argument 'pval' missing", call. = FALSE)
  ##
  k <- length(TE)
  arg <- "TE"
  chklength(pval, k, arg)
  chklength(df, k, arg)
  ##
  if (any(pval <= 0, na.rm = TRUE) | any(pval >= 1, na.rm = TRUE))
    stop("No valid value for p-value", call. = FALSE)
  ##
  ## Parmar et al. (1998), Stat Med
  ##
  ## Equation (7)
  ##
  varTE <- ifelse(is.na(df),
  (TE / qnorm(pval / 2, lower.tail = FALSE))^2,
  (TE / qt(pval / 2, df = df, lower.tail = FALSE))^2)
  seTE <- sqrt(varTE)
  ##
  res <- list(TE = TE, seTE = seTE, pval = pval)
  ##
  res
}


NULL
