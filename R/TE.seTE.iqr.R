TE.seTE.iqr <- function(n, median, q1, q3) {
  
  
  ##
  ## Estimate mean and its standard error from median and
  ## interquartile range using method by Wan et al. (2014), BMC Med
  ## Res Methodol 14 (1): 135
  ##
  
  
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
  ## Equation (14)
  ##
  TE <- (q1 + median + q3) / 3
  ##
  ## Equations (15) and (16)
  ##
  seTE <- (q3 - q1) /
    ifelse(n > 50, 2 * qnorm((0.75 * n - 0.125) / (n + 0.25)),
           .settings$Wan2014.Table2[n])
  
  
  res <- list(TE = TE, seTE = seTE,
              median = median, q1 = q1, q3 = q3, n = n)
  ##
  res
}
