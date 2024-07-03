mean_sd_iqr <- function(n, median, q1, q3, method.mean = "Luo") {
  
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
  if (method.mean == "Luo") {
    ## Luo et al. (2018), equation (15)
    mean <-
      (0.7 + 0.39 / n) * (q1 + q3) / 2 +
      (0.3 - 0.39 / n) * median
  }
  else if (method.mean == "Wan") {
    ## Wan et al. (2014), equation (14)
    mean <- (q1 + median + q3) / 3
  }
  else if (method.mean == "Cai") {
    ## Cai et al. (2021)
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::mln.mean.sd(q1.val = q1[i],
                               med.val = median[i],
                               q3.val = q3[i],
                               n = n[i])$est.mean
  }
  else if (method.mean == "QE-McGrath") {
    ## McGrath et al. (2020), QE method
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::qe.mean.sd(q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              n = n[i])$est.mean
  }
  else if (method.mean == "BC-McGrath") {
    ## McGrath et al. (2020), BC method
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::bc.mean.sd(q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              n = n[i])$est.mean
  }
  else
    mean <- NA
  
  
  #
  # Estimation of standard deviation
  # Wan et al. (2014), equations (15) and (16)
  #
  method.sd <- method.mean
  #
  if (method.sd %in% c("Luo", "Wan")) {
    sd <- (q3 - q1) /
      ifelse(n > 201, 2 * qnorm((0.75 * n - 0.125) / (n + 0.25)),
             gs("Wan2014.Table2")[ceiling(0.25 * (n - 1))])
  }
  else if (method.sd == "Cai") {
    # Cai et al. (2021)
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::mln.mean.sd(q1.val = q1[i],
                               med.val = median[i],
                               q3.val = q3[i],
                               n = n[i])$est.sd
  }
  else if (method.sd == "QE-McGrath") {
    # McGrath et al. (2020), QE method
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::qe.mean.sd(q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              n = n[i])$est.sd
  }
  else if (method.sd == "BC-McGrath") {
    # McGrath et al. (2020), BC method
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::bc.mean.sd(q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              n = n[i])$est.sd
  }

  
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


mean_sd_iqr_range <- function(n, median, q1, q3, min, max,
                              method.mean = "Luo", method.sd = "Shi") {
  
  
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
  ## Estimation of mean
  ##
  if (method.mean == "Luo") {
    ## Luo et al. (2018), equation (15)
    mean <-
      2.2 / (2.2 + n^0.75) * (min + max) / 2 +
      (0.7 - 0.72 / n^0.55) * (q1 + q3) / 2 +
      (0.3 + 0.72 / n^0.55 - 2.2 / (2.2 + n^0.75)) * median
  }
  else if (method.mean == "Wan") {
    ## Wan et al. (2014), equation (10)
    mean <- (min + 2 * q1 + 2 * median + 2 * q3 + max) / 8
  }
  else if (method.mean == "Cai") {
    ## Cai et al. (2021)
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::mln.mean.sd(min.val = min[i],
                               q1.val = q1[i],
                               med.val = median[i],
                               q3.val = q3[i],
                               max.val = max[i],
                               n = n[i])$est.mean
  }
  else if (method.mean == "QE-McGrath") {
    ## McGrath et al. (2020), QE method
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::qe.mean.sd(min.val = min[i],
                              q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              max.val = max[i],
                              n = n[i])$est.mean
  }
  else if (method.mean == "BC-McGrath") {
    ## McGrath et al. (2020), BC method
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::bc.mean.sd(min.val = min[i],
                              q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              max.val = max[i],
                              n = n[i])$est.mean
  }
  else
    mean <- rep_len(NA, k)
  
  
  ##
  ## Estimation of standard deviation
  ##
  if (method.sd == "Shi") {
    ## Shi et al. (2020), equation (11)
    theta1 <- (2 + 0.14 * n^0.6) * qnorm((n - 0.375) / (n + 0.25))
    theta2 <- (2 + 2 / (0.07 * n^0.6)) * qnorm((0.75 * n - 0.125) / (n + 0.25))
    sd <- (max - min) / theta1 + (q3 - q1) / theta2
  }
  else if (method.sd == "Wan") {
    ## Wan et al. (2014), equations (12) and (13)
    ##
    sd <- 0.5 * (
      (max - min) / ifelse(n > 50, 2 * qnorm((n - 0.375) / (n + 0.25)),
                           gs("Wan2014.Table1")[n]) +
      (q3 - q1) / ifelse(n > 201, 2 * qnorm((0.75 * n - 0.125) / (n + 0.25)),
                         gs("Wan2014.Table2")[ceiling(0.25 * (n - 1))])
    )
  }
  else if (method.sd == "Cai") {
    ## Cai et al. (2021)
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::mln.mean.sd(min.val = min[i],
                               q1.val = q1[i],
                               med.val = median[i],
                               q3.val = q3[i],
                               max.val = max[i],
                               n = n[i])$est.sd
  }
  else if (method.sd == "QE-McGrath") {
    ## McGrath et al. (2020), QE method
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::qe.mean.sd(min.val = min[i],
                              q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              max.val = max[i],
                              n = n[i])$est.sd
  }
  else if (method.sd == "BC-McGrath") {
    ## McGrath et al. (2020), BC method
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::bc.mean.sd(min.val = min[i],
                              q1.val = q1[i],
                              med.val = median[i],
                              q3.val = q3[i],
                              max.val = max[i],
                              n = n[i])$est.sd
  }
  else
    sd <- rep_len(NA, k)
  
  
  ##
  ## Calculation of standard error
  ##
  se <- sd / sqrt(n)
  
  
  res <- list(mean = mean, sd = sd, se = se,
              median = median, q1 = q1, q3 = q3, min = min, max = max, n = n,
              method.mean = method.mean)
  ##
  res
}


mean_sd_range <- function(n, median, min, max, method.mean = "Luo") {
  
  
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
  #
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
  if (method.mean == "Luo") {
    ## Luo et al. (2018), equation (15)
    mean <-
      4 / (4 + n^0.75) * (min + max) / 2 +
      n^0.75 / (4 + n^0.75) * median
  }
  else if (method.mean == "Wan") {
    ## Wan et al. (2014), equation (2)
    mean <- (min + 2 * median + max) / 4 +
      ifelse(is.na(n), 0, (min - 2 * median + max) / (4 * n))
  }
  else if (method.mean == "Cai") {
    ## Cai et al. (2021)
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::mln.mean.sd(min.val = min[i],
                               med.val = median[i],
                               max.val = max[i],
                               n = n[i])$est.mean
  }
  else if (method.mean == "QE-McGrath") {
    ## McGrath et al. (2020), QE method
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::qe.mean.sd(min.val = min[i],
                              med.val = median[i],
                              max.val = max[i],
                              n = n[i])$est.mean
  }
  else if (method.mean == "BC-McGrath") {
    ## McGrath et al. (2020), BC method
    mean <- vector("numeric", k)
    for (i in seq_len(k))
      mean[i] <-
        estmeansd::bc.mean.sd(min.val = min[i],
                              med.val = median[i],
                              max.val = max[i],
                              n = n[i])$est.mean
  }
  else
    mean <- rep_len(NA, k)
  
  
  #
  # Estimation of standard deviation
  # Wan et al. (2014), equations (7) and (9)
  #
  method.sd <- method.mean
  #
  if (method.sd %in% c("Luo", "Wan")) {
    sd <- (max - min) /
      ifelse(n > 50, 2 * qnorm((n - 0.375) / (n + 0.25)),
             gs("Wan2014.Table1")[n])
  }
  else if (method.sd == "Cai") {
    # Cai et al. (2021)
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::mln.mean.sd(min.val = min[i],
                               med.val = median[i],
                               max.val = max[i],
                               n = n[i])$est.sd
  }
  else if (method.sd == "QE-McGrath") {
    # McGrath et al. (2020), QE method
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::qe.mean.sd(min.val = min[i],
                              med.val = median[i],
                              max.val = max[i],
                              n = n[i])$est.sd
  }
  else if (method.sd == "BC-McGrath") {
    # McGrath et al. (2020), BC method
    sd <- vector("numeric", k)
    for (i in seq_len(k))
      sd[i] <-
        estmeansd::bc.mean.sd(min.val = min[i],
                              med.val = median[i],
                              max.val = max[i],
                              n = n[i])$est.sd
  }
  
  
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
