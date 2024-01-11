#' Auxiliary functions for (back) transformations
#' 
#' @description
#' Auxiliary functions to (back) transform effect estimates or
#' confidence / prediction interval limit(s).
#' 
#' @details
#' Often in a meta-analysis, effect estimates are transformed before
#' calculating a weighted average. For example, the log odds ratio and
#' its standard error is used instead of the odds ratio in R function
#' \code{\link{metagen}}. To report the results of a meta-analysis,
#' effect estimates are typically back transformed to the original
#' scale. R package \bold{meta} provides some auxiliary functions for
#' (back) transformations.
#'
#' \subsection{Transformations}{
#' The following auxiliary functions are provided by R package
#' \bold{meta} to transform effect estimates or confidence /
#' prediction interval limits.
#'
#' \tabular{ll}{
#' \bold{Function} \tab \bold{Transformation} \cr
#' \code{cor2z} \tab Correlations to Fisher's Z transformed correlations \cr
#' \code{p2logit} \tab Proportions to logit transformed proportions \cr
#' \code{p2asin} \tab Proportions to arcsine transformed proportions \cr
#' \code{VE2logVR} \tab Vaccine efficacy / effectiveness to log vaccine ratio
#' }
#'
#' Note, no function for the Freeman-Tukey arcsine transformation is
#' provided as this transformation is based on the number of events
#' and sample sizes instead of the effect estimates.
#' 
#' R function \code{transf} is a wrapper function for the above and
#' additional transformations, e.g., the log transformation using
#' \code{\link[base]{log}} for odds or risk ratios. Argument \code{sm}
#' is mandatory to specify the requested transformation. It is also
#' possible to specify a different function with arguments \code{func}
#' and \code{args}.
#' }
#'
#' \subsection{Back transformations}{
#' The following auxiliary functions are available to back transform
#' effect estimates or confidence / prediction interval limits.
#'
#' \tabular{ll}{
#' \bold{Function} \tab \bold{Transformation} \cr
#' \code{asin2ir} \tab Freeman-Tukey arcsine transformed rates to rates \cr
#' \code{asin2p} \tab (Freeman-Tukey) arcsine transformed proportions to proportions \cr
#' \code{logit2p} \tab Logit transformed proportions to proportions \cr
#' \code{logVR2VE} \tab Log vaccine ratio to vaccine efficacy / effectiveness \cr
#' \code{z2cor} \tab Fisher's Z transformed correlations  to correlations
#' }
#' 
#' Argument \code{time} is mandatory in R function \code{asin2ir}.
#'
#' If argument \code{n} is provided in R function \code{asin2p},
#' Freeman-Tukey arcsine transformed proportions are
#' back transformed. Otherwise, arcsine transformed proportions are
#' back transformed.
#'
#' R function \code{backtransf} is a wrapper function for the above
#' and additional transformations, e.g., the exponential
#' transformation using \code{\link[base]{exp}} for log odds or log
#' risk ratios. Argument \code{sm} is mandatory to specify the
#' requested transformation. For the Freeman-Tukey transformations,
#' argument \code{n} or \code{time} is mandatory.
#'
#' It is also possible to specify a different function with arguments
#' \code{func} and \code{args}.
#' }
#' 
#' @param x Numerical vector with effect estimates, lower or upper
#'   confidence / prediction interval limit(s).
#' @param sm Summary measure.
#' @param func User-specified function for (back) transformation.
#' @param args Function arguments for user-specified function.
#' @param n Sample size(s) to back transform Freeman-Tukey transformed
#'   proportions.
#' @param time Time(s) to back transform Freeman-Tukey transformed
#'   incidence rates.
#'
#' @name meta-transf
#' 
#' @aliases meta-transf
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-sm}}
#' 
#' @examples
#' logit2p(p2logit(0.5))
#'
#' @rdname meta-transf
#' @export transf


transf <- function(x, sm, func = NULL, args = NULL) {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)

  if (!is.null(func))
    res <- do.call(func, c(list(x), args))
  ##
  else if (is_relative_effect(sm) | is_log_effect(sm))
    res <- log(x)
  ##
  else if (sm == "ZCOR")
    res <- cor2z(x)
  ##
  else if (sm == "PLOGIT")
    res <- p2logit(x)
  ##
  else if (sm == "PAS")
    res <- p2asin(x)
  ##
  else if (sm == "IRS")
    res <- sqrt(x)
  ##
  else if (sm == "VE")
    res <- VE2logVR(x)
  ##
  else
    res <- x
  
  res
}





#' @rdname meta-transf
#' @export cor2z


cor2z <- function(x)
  0.5 * log((1 + x) / (1 - x))





#' @rdname meta-transf
#' @export p2asin


p2asin <- function(x)
  asin(sqrt(x))





#' @rdname meta-transf
#' @export p2logit


p2logit <- function(x)
  qlogis(x)





#' @rdname meta-transf
#' @export VE2logVR


VE2logVR <- function(x)
  log(1 - x / 100)





#' @rdname meta-transf
#' @export backtransf


backtransf <- function(x, sm, n, time, func = NULL, args = NULL) {
  
  chkchar(sm, length = 1)
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  if (!is.null(func))
    res <- do.call(func, c(list(x), args))
  ##
  else if (is_relative_effect(sm) | is_log_effect(sm))
    res <- exp(x)
  ##
  else if (sm == "ZCOR")
    res <- z2cor(x)
  ##
  else if (sm == "PLOGIT")
    res <- logit2p(x)
  ##
  else if (sm == "PAS")
    res <- asin2p(x)
  ##
  else if (sm == "PFT") {
    if (missing(n))
      stop("Argument 'n' must be provided to back transform ",
           "Freeman-Tukey transformed proportions.")
    res <- asin2p(x, n)
  }
  ##
  else if (sm == "IRS")
    res <- x^2
  ##
  else if (sm == "IRFT") {
    if (missing(time))
      stop("Argument 'time' must be provided to back transform ",
           "Freeman-Tukey transformed incidence rates.")
    res <- asin2ir(x, time)
  }
  ##
  else if (sm == "VE")
    res <- logVR2VE(x)
  ##
  else
    res <- x

  if (sm == "PRAW") {
    sel0 <- res[!is.na(res)] < 0
    sel1 <- res[!is.na(res)] > 1
    ##
    if (any(sel0, na.rm = TRUE))
      res[sel0] <- 0
    if (any(sel1, na.rm = TRUE))
      res[sel1] <- 1
  }
  ##
  if (sm == "PLN") {
    sel0 <- res[!is.na(res)] < 0
    sel1 <- res[!is.na(res)] > 1
    ##
    if (any(sel0, na.rm = TRUE))
      res[sel0] <- 0
    else if (any(sel1, na.rm = TRUE))
      res[sel1] <- 1
  }
  
  if (sm == "IR") {
    sel0 <- res[!is.na(res)] < 0
    ##
    if (any(sel0, na.rm = TRUE))
      res[sel0] <- 0
  }

  res
}





#' @rdname meta-transf
#' @export asin2ir


asin2ir <- function(x, time = NULL) {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  
  ##
  ## Calculate possible minimum for each transformation
  ##
  minimum <- 0.5 * (sqrt(0 / time) + sqrt((0 + 1) / time))
  ##
  sel0 <- x < minimum
  
  
  res <- rep(NA, length(x))
  ##
  sel <- !sel0
  sel <- !is.na(sel) & sel
  ##
  res[sel0] <- 0
  ##
  ## Back transformation of Freeman-Tukey double arcsine transformation:
  ##
  res[sel] <- (1 / time[sel] - 8 * x[sel]^2 + 16 * time[sel] * x[sel]^4) /
    (16 * x[sel]^2 * time[sel])
  ##
  res[res < 0] <- 0
  
  
  res
}





#' @rdname meta-transf
#' @export asin2p


asin2p <- function(x, n = NULL) {
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  
  ##
  ## Calculate possible minimum and maximum
  ## for each transformation
  ##
  if (is.null(n)) {
    minimum <- asin(sqrt(0))
    maximum <- asin(sqrt(1))
  }
  else {
    minimum <- 0.5 * (asin(sqrt(0 / (n + 1))) + asin(sqrt((0 + 1) / (n + 1))))
    maximum <- 0.5 * (asin(sqrt(n / (n + 1))) + asin(sqrt((n + 1) / (n + 1))))
  }
  ##
  sel0 <- x < minimum
  sel1 <- x > maximum
  
  
  res <- rep(NA, length(x))
  ##
  sel <- !(sel0 | sel1)
  sel <- !is.na(sel) & sel
  ##
  res[sel0] <- 0
  res[sel1] <- 1
  ##
  if (is.null(n)) {
    ##
    ## Back transformation of arcsine transformation:
    ##
    res[sel] <- sin(x[sel])^2
  }
  else {
    ##
    ## Back transformation of Freeman-Tukey double arcsine transformation:
    ##
    res[sel] <- 0.5 * (1 - sign(cos(2 * x[sel])) *
                       sqrt(1 - (sin(2 * x[sel]) +
                                 (sin(2 * x[sel]) -
                                  1 / sin(2 * x[sel])) / n[sel])^2))
  }
  res
}





#' @rdname meta-transf
#' @export logit2p


logit2p <- function(x)
  1 / (1 + exp(-x))





#' @rdname meta-transf
#' @export logVR2VE


logVR2VE <- function(x)
  100 * (1 - exp(x))





#' @rdname meta-transf
#' @export z2cor


z2cor <- function(x)
  tanh(x)
