ciAgrestiCoull <- function(event, n, level) {
  chknumeric(event, 0)
  chknumeric(n, 0, zero = TRUE)
  chklevel(level)
  ##
  if (length(event) == 1 & length(n) > 1)
    event <- rep(event, length(n))
  else if (length(event) > 1 & length(n) == 1)
    n <- rep(n, length(event))
  else
    if (length(event) != length(n))
      stop("Arguments 'event' and 'n' must be of same length.", call. = FALSE)
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must be smaller equal sample size.", call. = FALSE)
  ##
  z <- qnorm(1 - (1 - level) / 2)
  ##
  n <- n + z^2
  prop <- 1 / n * (event + 0.5 * z^2)
  ##
  lower <- prop - z * sqrt(1 / n * prop * (1 - prop))
  upper <- prop + z * sqrt(1 / n * prop * (1 - prop))
  ##
  lower[lower < 0] <- 0
  upper[upper > 1] <- 1
  ##
  list(event = event, n = n,
       prop = prop, lower = lower, upper = upper,
       statistic = NA, p = NA, level = level,
       df = NA, null.effect = NA)
}


ciClopperPearson <- function(event, n, level, null.effect) {
  chknumeric(event, 0)
  chknumeric(n, 0, zero = TRUE)
  chklevel(level)
  chknumeric(null.effect, min = 0, max = 1, length = 1)
  ##
  if (length(event) == 1 & length(n) > 1)
    event <- rep(event, length(n))
  else if (length(event) > 1 & length(n) == 1)
    n <- rep(n, length(event))
  else
    if (length(event) != length(n))
      stop("Arguments 'event' and 'n' must be of same length.", call. = FALSE)
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must be smaller equal sample size.", call. = FALSE)
  ##
  k <- length(event)
  lower <- upper <- statistic <- pval <- rep(NA, k)
  ##
  for (i in seq_len(k)) {
    if (!is.na(event[i] & !is.na(n[i]))) {
      cint <- binom.test(event[i], n[i], conf.level = level,
                         p = if (!is.na(null.effect)) null.effect else 0.5)
      ##
      lower[i] <- cint$conf.int[[1]]
      upper[i] <- cint$conf.int[[2]]
      if (!is.na(null.effect))
        pval[i] <- cint$p.value
    }
    else {
      lower[i] <- NA
      upper[i] <- NA
      pval[i] <- NA
    }
  }
  ##
  list(event = event, n = n,
       prop = event / n, lower = lower, upper = upper,
       statistic = statistic, p = pval, level = level,
       df = NA, null.effect = null.effect)
}


ciSimpleAsymptotic <- function(event, n, level, correct = FALSE) {
  ##
  ## Newcombe RG. Two-sided confidence intervals for the single
  ## proportion: Comparison of seven methods.
  ## Stat Med 1998, Apr 30;17(8): 857-72
  ##
  chknumeric(event, 0)
  chknumeric(n, 0, zero = TRUE)
  chklevel(level)
  chklogical(correct)
  ##
  if (length(event) == 1 & length(n) > 1)
    event <- rep(event, length(n))
  else if (length(event) > 1 & length(n) == 1)
    n <- rep(n, length(event))
  else
    if (length(event) != length(n))
      stop("Arguments 'event' and 'n' must be of same length.", call. = FALSE)
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must be smaller equal sample size.", call. = FALSE)
  ##
  prop <- event / n
  z <- qnorm(1 - (1 - level) / 2)
  ##
  if (!correct) {
    lower  <- prop - z * sqrt(prop * (1 - prop) / n)
    upper  <- prop + z * sqrt(prop * (1 - prop) / n)
  }
  else {
    lower  <- prop - (z * sqrt(prop * (1 - prop) / n) + 1 / (2 * n))
    upper  <- prop + (z * sqrt(prop * (1 - prop) / n) + 1 / (2 * n))
  }
  ##
  lower[lower < 0] <- 0
  upper[upper > 1] <- 1
  ##
  list(event = event, n = n,
       prop = prop, lower = lower, upper = upper,
       statistic = NA, p = NA, level = level,
       df = NA, null.effect = NA)
}


ciWilsonScore <- function(event, n, level, correct = FALSE) {
  ##
  ## Newcombe RG. Two-sided confidence intervals for the single
  ## proportion: Comparison of seven methods.
  ## Stat Med 1998, Apr 30;17(8): 857-72
  ##
  chknumeric(event, 0)
  chknumeric(n, 0, zero = TRUE)
  chklevel(level)
  ##
  if (length(event) == 1 & length(n) > 1)
    event <- rep(event, length(n))
  else if (length(event) > 1 & length(n) == 1)
    n <- rep(n, length(event))
  else
    if (length(event) != length(n))
      stop("Arguments 'event' and 'n' must be of same length.", call. = FALSE)
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must be smaller equal sample size.", call. = FALSE)
  ##
  prop <- event / n
  z <- qnorm(1 - (1 - level) / 2)
  ##
  if (!correct) {
    lower  <- (2 * n * prop + z^2 -
               z * sqrt(z^2 + 4 * n * prop * (1 - prop))) / (2 * (n + z^2))
    upper  <- (2 * n * prop + z^2 +
               z * sqrt(z^2 + 4 * n * prop * (1 - prop))) / (2 * (n + z^2))
    }
  else {
      lower <- (2 * n * prop + z^2 - 1 -
                z * sqrt(z^2 - 2 - 1/n + 4 * prop * (n * (1 - prop) + 1))) /
        (2 * (n + z^2))
    lower[lower < 0] <- 0
      upper <- (2 * n * prop + z^2 + 1 +
                z * sqrt(z^2 + 2 - 1/n + 4 * prop * (n * (1 - prop) - 1))) /
        (2 * (n + z^2))
    upper[upper > 1] <- 1
  }
  ##
  list(event = event, n = n,
       prop = prop, lower = lower, upper = upper,
       statistic = NA, p = NA, level = level,
       df = NA, null.effect = NA)
}


ciPoisson <- function(event, time, level, null.effect) {
  chknumeric(event, 0)
  chknumeric(time, 0, zero = TRUE)
  chklevel(level)
  chknumeric(null.effect, min = 0, max = 1, length = 1)
  ##
  if (length(event) == 1 & length(time) > 1)
    event <- rep(event, length(time))
  else if (length(event) > 1 & length(time) == 1)
    time <- rep(time, length(event))
  else
    if (length(event) != length(time))
      stop("Arguments 'event' and 'time' must be of same length.",
           call. = FALSE)
  ##
  k <- length(event)
  lower <- upper <- statistic <- pval <- rep(NA, k)
  ##
  for (i in seq_len(k)) {
    if (!is.na(event[i] & !is.na(time[i]))) {
      cint <-
        poisson.test(event[i], time[i], conf.level = level,
                     r = if (!is.na(null.effect)) null.effect else 1)
      ##
      lower[i] <- cint$conf.int[[1]]
      upper[i] <- cint$conf.int[[2]]
      if (!is.na(null.effect))
        pval[i] <- cint$p.value
    }
    else {
      lower[i] <- NA
      upper[i] <- NA
      pval[i] <- NA
    }
  }
  ##
  list(event = event, time = time,
       rate = event / time, lower = lower, upper = upper,
       statistic = statistic, p = pval, level = level,
       df = NA, null.effect = null.effect)
}
