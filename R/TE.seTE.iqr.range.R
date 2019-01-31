TE.seTE.iqr.range <- function(n, median, q1, q3, min, max) {
  
  
  ##
  ## Estimate mean and its standard error from median, interquartile
  ## range and range using method by Wan et al. (2014), BMC Med Res
  ## Methodol 14 (1): 135
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
  if (missing(min))
    stop("Mandatory argument 'min' missing.", call. = FALSE)
  if (missing(max))
    stop("Mandatory argument 'max' missing.", call. = FALSE)
  ##
  chknumeric(n, min = 0, zero = TRUE)
  ##
  k <- length(n)
  arg <- "n"
  chklength(median, k, arg)
  chklength(q1, k, arg)
  chklength(q3, k, arg)
  chklength(min, k, arg)
  chklength(max, k, arg)
  ##
  if (any(median < q1, na.rm = TRUE))
    stop("Median must be larger than first quartile.", call. = FALSE)
  if (any(median > q3, na.rm = TRUE))
    stop("Median must be smaller than third quartile.", call. = FALSE)
  if (any(q1 >= q3, na.rm = TRUE))
    stop("First quartile must be smaller than third quartile.", call. = FALSE)
  ##
  if (any(median < min, na.rm = TRUE))
    stop("Median must be larger than minumum.", call. = FALSE)
  if (any(median > max, na.rm = TRUE))
    stop("Median must be smaller than maximum.", call. = FALSE)
  if (any(min >= max, na.rm = TRUE))
    stop("Minimum must be smaller than maximum.", call. = FALSE)
  ##
  if (any(q1 < min, na.rm = TRUE))
    stop("First quartile must be larger than minumum.", call. = FALSE)
  if (any(q3 > max, na.rm = TRUE))
    stop("Third quartile must be smaller than maximum.", call. = FALSE)
  
  
  ##
  ## Equation (10)
  ##
  TE <- (min + 2 * q1 + 2 * median + 2 * q3 + max) / 8
  ##
  ## Equations (12) and (13)
  ##
  seTE <- 0.5 * (
    (max - min) / ifelse(n > 50, 2 * qnorm((n - 0.375) / (n + 0.25)),
                         .settings$Wan2014.Table1[n]) +
    (q3 - q1) / ifelse(n > 50, 2 * qnorm((0.75 * n - 0.125) / (n + 0.25)),
                       .settings$Wan2014.Table2[n])
  )
  
  
  res <- list(TE = TE, seTE = seTE,
              median = median, q1 = q1, q3 = q3, n = n)
  ##
  res
}
