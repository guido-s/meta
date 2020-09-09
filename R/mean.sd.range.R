mean.sd.range <- function(n, median, min, max, method.mean = "Luo") {
  
  
  ##
  ## Check arguments
  ##
  if (missing(n))
    stop("Mandatory argument 'n' missing.", call. = FALSE)
  if (missing(median))
    stop("Mandatory argument 'median' missing.", call. = FALSE)
  if (missing(min))
    stop("Mandatory argument 'min' missing.", call. = FALSE)
  if (missing(max))
    stop("Mandatory argument 'max' missing.", call. = FALSE)
  ##
  chknumeric(n, min = 0, zero = TRUE)
  ##
  k <- length(n)
  fun <- "TE.seTE.range"
  chklength(median, k, fun,
            text = "Arguments 'n' and 'median' must have the same length.")
  chklength(min, k, fun,
            text = "Arguments 'n' and 'min' must have the same length.")
  chklength(max, k, fun,
            text = "Arguments 'n' and 'max' must have the same length.")
  ##
  if (any(median < min, na.rm = TRUE))
    stop("Median must be larger than minumum.", call. = FALSE)
  if (any(median > max, na.rm = TRUE))
    stop("Median must be smaller than maximum.", call. = FALSE)
  if (any(min >= max, na.rm = TRUE))
    stop("Minimum must be smaller than maximum.", call. = FALSE)
  
  
  ##
  ## Estimation of mean
  ##
  if (tolower(method.mean) == "luo") {
    ## Luo et al. (2018), equation (15)
    mean <-
      4 / (4 + n^0.75) * (min + max) / 2 +
      n^0.75 / (4 + n^0.75) * median
  }
  else if (tolower(method.mean) == "wan") {
    ## Wan et al. (2014), equation (2)
    mean <- (min + 2 * median + max) / 4 +
      ifelse(is.na(n), 0, (min - 2 * median + max) / (4 * n))
  }
  else
    mean <- NA
  
  
  ##
  ## Estimation of standard deviation
  ## Wan et al. (2014), equations (7) and (9)
  ##
  sd <- (max - min) /
    ifelse(n > 50, 2 * qnorm((n - 0.375) / (n + 0.25)),
           .settings$Wan2014.Table1[n])
  
  
  ##
  ## Calculation of standard error
  ##
  se <- sd / sqrt(n)
  
  
  res <- list(mean = mean, sd = sd, se = se,
              median = median, min = min, max = max, n = n,
              method.mean = method.mean)
  ##
  res
}
