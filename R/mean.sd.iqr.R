mean.sd.iqr <- function(n, median, q1, q3, method.mean = "Luo") {
  
  
  ##
  ## Check arguments
  ##
  if (missing(n))
    stop("Mandatory argument 'n' missing.", call. = FALSE)
  if (missing(median))
    stop("Mandatory argument 'median' missing.", call. = FALSE)
  if (missing(q1))
    stop("Mandatory argument 'q1' missing.", call. = FALSE)
  if (missing(q3))
    stop("Mandatory argument 'q3' missing.", call. = FALSE)
  ##
  chknumeric(n, min = 0, zero = TRUE)
  ##
  k <- length(n)
  arg <- "n"
  chklength(median, k, arg)
  chklength(q1, k, arg)
  chklength(q3, k, arg)
  ##
  if (any(median < q1, na.rm = TRUE))
    stop("Median must be larger than first quartile.", call. = FALSE)
  if (any(median > q3, na.rm = TRUE))
    stop("Median must be smaller than third quartile.", call. = FALSE)
  if (any(q1 >= q3, na.rm = TRUE))
    stop("First quartile must be smaller than third quartile.", call. = FALSE)
  
  
  ##
  ## Estimation of mean
  ##
  if (tolower(method.mean) == "luo") {
    ## Luo et al. (2018), equation (15)
    mean <-
      (0.7 + 0.39 / n) * (q1 + q3) / 2 +
      (0.3 - 0.39 / n) * median
  }
  else if (tolower(method.mean) == "wan") {
    ## Wan et al. (2014), equation (14)
    mean <- (q1 + median + q3) / 3
  }
  else
    mean <- NA
  
  
  ##
  ## Estimation of standard deviation
  ## Wan et al. (2014), equations (15) and (16)
  ##
  sd <- (q3 - q1) /
    ifelse(n > 50, 2 * qnorm((0.75 * n - 0.125) / (n + 0.25)),
           .settings$Wan2014.Table2[n])
  
  
  ##
  ## Calculation of standard error
  ##
  se <- sd / sqrt(n)
  
  
  res <- list(mean = mean, sd = sd, se = se,
              median = median, q1 = q1, q3 = q3, n = n,
              method.mean = method.mean)
  ##
  res
}
