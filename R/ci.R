#' Calculation of confidence intervals (based on normal approximation
#' or t-distribution)
#' 
#' @description
#' Calculation of confidence intervals; based on normal approximation
#' or t-distribution.
#' 
#' @param TE Estimated treatment effect.
#' @param seTE Standard error of treatment estimate.
#' @param level The confidence level required.
#' @param df Degrees of freedom (for confidence intervals based on
#'   t-distribution).
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' 
#' @return
#' List with components
#' \item{TE}{Estimated treatment effect}
#' \item{seTE}{Standard error of treatment estimate}
#' \item{lower}{Lower confidence limits}
#' \item{upper}{Upper confidence limits}
#' \item{statistic}{Test statistic (either z-score or t-score)}
#' \item{p}{P-value of text.with null hypothesis \code{TE=0}}
#' \item{level}{The confidence level required}
#' \item{df}{Degrees of freedom (t-distribution)}
#' 
#' @note
#' This function is primarily called from other functions of the
#' library \code{meta}, e.g. \code{forest.meta}, \code{summary.meta}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @examples
#' data.frame(ci(170, 10))
#' data.frame(ci(170, 10, 0.99))
#' data.frame(ci(1.959964, 1))
#' data.frame(ci(2.2621571628, 1, df = 9))
#' 
#' @export ci


ci <- function(TE, seTE, level = 0.95, df = NULL, null.effect = 0) {
  
  chklevel(level)
  alpha <- 1 - level
  
  df.orig <- df
  
  if (is.null(df)) {
    lower <- TE - qnorm(1 - alpha / 2) * seTE
    upper <- TE + qnorm(1 - alpha / 2) * seTE
    statistic <- (TE - null.effect) / seTE
    pval <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
    ##
    df <- Inf
    df.orig <- "NULL"
  }
  else {
    statistic <- (TE - null.effect) / seTE
    ##
    if (length(df) == 1) {
      if (is.na(df) || df <= 0) {
        lower <- TE - qnorm(1 - alpha / 2) * seTE
        upper <- TE + qnorm(1 - alpha / 2) * seTE
        ##
        pval <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
      }
      else {
        lower <- TE - qt(1 - alpha / 2, df = df) * seTE
        upper <- TE + qt(1 - alpha / 2, df = df) * seTE
        ##
        pval <- 2 * pt(abs(statistic), df = df, lower.tail = FALSE)
      }
    }
    else {
      sel0 <- is.na(df) | df <= 0
      df[sel0] <- NA
      lower <- ifelse(!is.na(df),
                      TE - qt(1 - alpha / 2, df = df) * seTE,
                      TE - qnorm(1 - alpha / 2) * seTE)
      upper <- ifelse(!is.na(df),
                      TE + qt(1 - alpha / 2, df = df) * seTE,
                      TE + qnorm(1 - alpha / 2) * seTE)
      ##        
      pval <- ifelse(!is.na(df),
                     2 * pt(abs(statistic), df = df, lower.tail = FALSE),
                     2 * pnorm(abs(statistic), lower.tail = FALSE))
      ##
      df[sel0] <- df.orig[sel0]
    }
    ##
    df[is.na(df)] <- Inf
  }
  
  res <- list(TE = TE, seTE = seTE,
              lower = lower, upper = upper,
              statistic = statistic, p = pval, level = level,
              df = df, null.effect = null.effect)
  
  res
}
