estimate.missing <- function(TE, TE.sum, type) {
  ##
  ## 1. Center around mean
  ##
  TE.c <- TE - TE.sum
  n <- length(TE.c)
  ##
  ## 2. Rank absolute values of centered values
  ##
  r.star <- rank(abs(TE.c)) * sign(TE.c)
  ##
  if (type == "L") {
    ##
    ## 3. Sum the positive ranks only
    ##
    S.rank <- sum(r.star[r.star>0])
    ##
    ## 4. Estimate for L0
    ##
    res0 <- (4 * S.rank - n * (n + 1)) / (2 * n - 1)
    res0.plus <- max(0, res0 + 0.5) %/% 1
  }
  if (type == "R") {
    ##
    ## 5. Estimate for R0
    ##
    res0 <- n - abs(min(r.star)) - 1.5
    res0.plus <- max(0, res0 + 0.5) %/% 1
  }
  ##
  res <- list(res0 = res0, res0.plus = res0.plus)
  res
}
hypergeometric <- function(n1, m1, N, psi) {
  ##
  ## R program for computing the mean, variance, density, cumulative
  ## distribution and generating random deviates.
  ##
  ## Based on
  ## Liao and Rosen (2001): Fast and Stable Algorithms for Computing
  ## and Sampling from the Noncentral Hypergeometric Distribution,
  ## The American Statistician, 55, 236-369.
  ##
  ## this is how to use the function
  ##
  ## n1 <- 100
  ## n2 <- 100
  ## m1 <- 100
  ## N <- n1 + n2
  ## odds.ratio <- 3
  ## obj <- hypergeometric(n1, m1, N, odds.ratio)
  ## obj$mean()
  ## obj$var()
  ## obj$d(40)
  ## obj$p(40)
  ## obj$r()
  ##
  n2 <- N - n1
  ##
  if (n1 < 0 | n2 < 0 | m1 < 0 | m1 > N | psi <= 0)
    stop("Wrong argument in hypergeometric")
  ##
  mode.compute <- function() {
    a <- psi - 1
    b <- -((m1 + n1 + 2) * psi + n2 - m1)
    c <- psi * (n1 + 1) * (m1 + 1)
    q <- b + sign(b) * sqrt(b * b - 4 * a * c)
    q <- -q / 2
    ##
    mode <- trunc(c / q)
    if (uu >= mode && mode >= ll)
      return(mode)
    else
      return(trunc(q / a))
  }
  ##
  r.function <- function(i)
    (n1 - i + 1) * (m1 - i + 1) / i / (n2 - m1 + i) * psi
  ##
  mean <- function()
    sum(prob[(ll:uu) + shift] * (ll:uu))
  ##
  var <- function()
    sum(prob[(ll:uu) + shift] * (ll:uu)^2) - mean()^2
  ##
  d <- function(x)
    return(prob[x + shift])
  ##
  p <- function(x, lower.tail = TRUE) {
    if (lower.tail)
      return(sum(prob[ll:(x + shift)]))
    else
      return(sum(prob[(x + shift):uu]))
  }
  ##
  sample.low.to.high <- function(lower.end, ran) {
    for (i in lower.end:uu) {
      if (ran <= prob[i + shift]) return(i)
      ran <- ran - prob[i + shift]
    }
  }
  ##
  sample.high.to.low <- function(upper.end, ran) {
    for (i in upper.end:ll) {
      if (ran <= prob[i + shift])
        return(i)
      ran <- ran - prob[i + shift]
    }
  }
  ##
  r <- function() {
    ran <- runif(1)
    ##
    if (mode == ll)
      return(sample.low.to.high(ll, ran))
    ##
    if (mode == uu)
      return(sample.high.to.low(uu, ran))
    ##
    if (ran < prob[mode + shift])
      return(mode)
    ##
    ran <- ran - prob[mode + shift]
    ##
    lower <- mode - 1
    upper <- mode + 1
    ##
    repeat {
      if (prob[upper + shift] >= prob[lower + shift]) {
        if (ran < prob[upper + shift])
          return(upper)
        ran <- ran - prob[upper + shift]
        if (upper == uu)
          return(sample.high.to.low(lower, ran))
        upper <- upper + 1
      }
      else {
        if (ran < prob[lower + shift])
          return(lower)
        ran <- ran - prob[lower + shift]
        if (lower == ll)
          return(sample.low.to.high(upper, ran))
        lower <- lower - 1
      }
    }
  }
  ##
  ll <- max(0, m1 - n2)
  uu <- min(n1, m1)
  mode <- mode.compute()
  ##
  prob <- array(1, uu - ll + 1)
  ##
  shift <- 1 - ll
  if (mode < uu) {
    ## note the shift of location
    r1 <- r.function((mode + 1):uu)
    prob[(mode + 1 + shift):(uu + shift)] <- cumprod(r1)
  }
  ##
  if (mode > ll) {
    r1 <- 1 / r.function(mode:(ll + 1))
    prob[(mode - 1 + shift):(ll + shift)] <- cumprod(r1)
  }
  ##
  prob <- prob / sum(prob)
  ##
  return(list(mean = mean, var = var, d = d, p = p, r = r))
}
kentau <- function(x, y, correct = FALSE, keep.data = FALSE) {
  ##
  ## Check:
  ##
  if (length(x) != length(y))
    stop("length of argument x and y must be equal")
  ##
  sel <- !is.na(x) & !is.na(y)
  if (length(x) != sum(sel))
    warning(paste(length(x) - sum(sel),
                  "observation(s) dropped due to missing values"))
  ##
  x <- x[sel]
  y <- y[sel]
  n <- length(x)
  ##
  t <- rle(sort(x))$lengths
  u <- rle(sort(y))$lengths
  ##
  N <- 0.5 * n * (n - 1)
  N1 <- N - sum(t * (t - 1) / 2)
  N2 <- N - sum(u * (u - 1) / 2)
  ##
  ks <- sqrt(N1) * sqrt(N2) * cor(x, y, method = "kendall")
  ##
  ## ks <- .C(kenscore,
  ##          kenscore = as.double(0),
  ##          x = as.double(x),
  ##          y = as.double(y),
  ##          n = as.integer(n))$kenscore
  ##
  ## Calculate S and s.e(S) according to
  ## Stata, release 5, Reference P-Z, p.239-240
  ##
  ## see also Kendall, Gibbons (1990), Rank Correlation Methods
  ## p. 66-68
  ##
  se.ks <- sqrt(1 / 18 * (n * (n - 1) * (2 * n + 5) -
                          sum(t * (t - 1) * (2 * t + 5)) -
                          sum(u * (u - 1) * (2 * u + 5))) +
                1 / (9 * n * (n - 1) * (n - 2))  *
                sum(t * (t - 1) * (t - 2)) *
                sum(u * (u - 1) * (u - 2)) +
                1 / (2 * n * (n - 1))  *
                sum(t * (t - 1))  *
                sum(u * (u - 1)))
  ##
  if (as.logical(correct) &
      any(c(length(unique(x)), length(unique(y))) == 2))
    warning("Continuity corrected statistic may be inappropriate,\n",
            "see Kendall, Gibbons (1990), Rank Correlation Methods, p.67.")
  ##
  statistic <- (ks - sign(ks) * as.logical(correct)) / se.ks
  p.value <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
  ##
  res <- list(tau.a = ks / N,
              tau.b = ks / (sqrt(N1) * sqrt(N2)),
              ks = ks - sign(ks) * as.logical(correct),
              se.ks = se.ks,
              statistic = statistic,
              p.value = p.value,
              correction = as.logical(correct))
  ##
  if (keep.data) {
    res$x <- x
    res$y <- y
  }
  ##
  res
}
linregcore <- function(TE, seTE, covar = NULL,
                       model = "lm", method.tau = "DL",
                       ...) {
  ##
  if (is.null(covar))
    predictor <- "sei"
  else
    predictor <- "ni"
  ##
  rma1 <-
    suppressWarnings(
      runNN(rma.uni,
            list(yi = TE, sei = seTE, ni = covar, method = method.tau,
                 test = "t", ...)))
  ##
  reg <- suppressWarnings(regtest(rma1, predictor = predictor,
                                  model = model))
  ##
  if (inherits(reg$fit, c("lm", "summary.lm"))) {
    if (!inherits(reg$fit, "summary.lm"))
      creg <- suppressWarnings(coef(summary(reg$fit)))
    else
      creg <- suppressWarnings(coef(reg$fit))
    ##
    res <- list(intercept = creg[1, 1], se.intercept = creg[1, 2],
                slope = creg[2, 1], se.slope = creg[2, 2])
  }
  else
    res <- list(intercept = reg$fit$beta[1], se.intercept = reg$fit$se[1],
                slope = reg$fit$beta[2], se.slope = reg$fit$se[2])
  ##
  res$statistic <- res$slope / res$se.slope
  res$df <- reg$dfs
  res$pval <- 2 * pt(abs(res$statistic), df = res$df, lower.tail = FALSE)
  ##
  if (inherits(reg$fit, c("lm", "summary.lm"))) {
    if (!inherits(reg$fit, "summary.lm"))
      res$tau <- suppressWarnings(summary(reg$fit)$sigma)
    else
      res$tau <- reg$fit$sigma
  }
  else
    res$tau <- sqrt(reg$fit$tau2)
  ##
  res$TE <- TE
  res$seTE <- seTE
  res$covar <- covar
  res$predictor <- predictor
  res$model <- model
  res$method.tau <- method.tau
  ##
  names(res$statistic) <- "t"
  ##
  res
}
