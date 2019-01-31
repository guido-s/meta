TE.seTE.range <- function(n, median, min, max) {
  
  
  ##
  ## Estimate mean and its standard error from median and range using
  ## method by Wan et al. (2014), BMC Med Res Methodol 14 (1): 135
  ##
  
  
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
  ## Equation (2)
  ##
  TE <- (min + 2 * median + max) / 4 +
    ifelse(is.na(n), 0, (min - 2 * median + max) / (4 * n))
  ##
  ## Equations (7) and (9)
  ##
  seTE <- (max - min) /
    ifelse(n > 50, 2 * qnorm((n - 0.375) / (n + 0.25)),
           .settings$Wan2014.Table1[n])
  
  
  res <- list(TE = TE, seTE = seTE,
              median = median, min = min, max = max, n = n)
  ##
  res
}
