#' Print results of a cumulative meta-analysis
#' 
#' @description
#' Print results of a cumulative meta-analysis
#' 
#' @aliases print.metacum
#' 
#' @param x An object of class \code{\link{metacum}}.
#' @param prediction A logical indicating whether prediction
#'   intervals should be printed.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value of test for overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau2.ci A logical value indicating whether to print
#'   the confidence interval of \eqn{\tau^2}.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param print.tau.ci A logical value indicating whether to print the
#'   confidence interval of \eqn{\tau}.
#' @param print.I2 A logical specifying whether heterogeneity
#'   statistic I\eqn{^2} should be printed.
#' @param print.I2.ci A logical specifying whether confidence interval for
#'   heterogeneity statistic I\eqn{^2} should be printed.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional arguments (ignored).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metacum}}, \code{\link{settings.meta}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' m1
#' metacum(m1)
#' metacum(m1, pooled = "random")
#' metacum(m1, pooled = "random", prediction = TRUE)
#'
#' @method print metacum
#' @export

print.metacum <- function(x, prediction = x$prediction,
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.tau = gs("digits.tau"),
                          digits.I2 = gs("digits.I2"),
                          #
                          big.mark = gs("big.mark"),
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
                          #
                          print.tau2 = TRUE,
                          print.tau2.ci = FALSE,
                          print.tau = TRUE,
                          print.tau.ci = FALSE,
                          #
                          print.I2 = TRUE,
                          print.I2.ci = FALSE,
                          #
                          text.tau2 = gs("text.tau2"),
                          text.tau = gs("text.tau"),
                          text.I2 = gs("text.I2"),
                          #
                          details.methods = gs("details"),
                          ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "metacum")
  x <- updateversion(x)
  #
  chklogical(prediction)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  #
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(print.tau2)
  chklogical(print.tau2.ci)
  chklogical(print.tau)
  chklogical(print.tau.ci)
  #
  if (all(x$method == "LRP") & all(x$method.random == "LRP")) {
    print.tau2 <- FALSE
    print.tau <- FALSE
  }
  #
  chklogical(print.I2)
  chklogical(print.I2.ci)
  #
  chklogical(details.methods)
  #
  print.tau2.ci <- print.tau2 & print.tau2.ci
  print.tau.ci <- print.tau & print.tau.ci
  print.I2.ci <- print.I2 & print.I2.ci
  
  #
  #
  # (3) Create data set
  #
  #
  
  # Get rid of warning 'no visible binding for global variable'
  TE <- lower <- upper <- statistic <- pval <-
    lower.predict <- upper.predict <-
    tau2 <- lower.tau2 <- upper.tau2 <-
    tau <- lower.tau <- upper.tau <-
    I2 <- lower.I2 <- upper.I2 <- NULL
  #
  dat.cum <-
    with(x,
         data.frame(TE, lower, upper,
                    statistic, pval,
                    lower.predict, upper.predict,
                    tau2, lower.tau2, upper.tau2,
                    tau, lower.tau, upper.tau,
                    I2, lower.I2, upper.I2,
                    row.names = studlab))
  #
  space <- data.frame(TE = NA, lower = NA, upper = NA,
                      statistic = NA, pval = NA,
                      lower.predict = NA, upper.predict = NA,
                      tau2 = NA, lower.tau2 = NA, upper.tau2 = NA,
                      tau = NA, lower.tau = NA, upper.tau = NA, 
                      I2 = NA, lower.I2 = NA, upper.I2 = NA,
                      row.names = "")
  #
  dat.pooled <-
    with(x,
         data.frame(TE.pooled, lower.pooled, upper.pooled,
                    statistic.pooled, pval.pooled,
                    lower.predict.pooled, upper.predict.pooled,
                    tau2.pooled, lower.tau2.pooled, upper.tau2.pooled,
                    tau.pooled, lower.tau.pooled, upper.tau.pooled,
                    I2.pooled, lower.I2.pooled, upper.I2.pooled)) %>%
    rename_with(~ gsub(".pooled", "", .x, fixed = TRUE))
  #
  rownames(dat.pooled) <- x$text.pooled
  #
  dat <- rbind(dat.cum, space, dat.pooled) %>%
    mutate(TE = formatN(TE, digits = digits, text.NA = "", big.mark = big.mark),
           lower = 
             if_else(is.na(lower) & is.na(upper), "",
                     formatCI(formatN(lower, digits = digits, text.NA = "",
                                      big.mark = big.mark),
                              formatN(upper, digits = digits, text.NA = "",
                                      big.mark = big.mark))),
           #
           statistic = formatN(statistic, digits = digits.stat, text.NA = ""),
           #
           pval = formatPT(pval, digits = digits.pval, lab.NA = "",
                           scientific = scientific.pval,
                           zero = zero.pval, JAMA = JAMA.pval),
           #
           lower.predict = 
             if_else(is.na(lower.predict) & is.na(upper.predict), "",
                     formatCI(formatN(lower.predict, digits = digits,
                                      text.NA = "", big.mark = big.mark),
                              formatN(upper.predict, digits = digits,
                                      text.NA = "", big.mark = big.mark))),
           #
           tau2 = formatPT(tau2, digits = digits.tau2, lab.NA = "",
                           big.mark = big.mark),
           lower.tau2 = 
             if_else(is.na(lower.tau2) & is.na(upper.tau2), "",
                     formatCI(formatPT(lower.tau2, digits = digits,
                                       lab.NA = "", big.mark = big.mark),
                              formatPT(upper.tau2, digits = digits,
                                       lab.NA = "", big.mark = big.mark))),
           #
           tau = formatPT(tau, digits = digits.tau, lab.NA = "",
                          big.mark = big.mark),
           lower.tau = 
             if_else(is.na(lower.tau) & is.na(upper.tau), "",
                     formatCI(formatPT(lower.tau, digits = digits,
                                       lab.NA = "", big.mark = big.mark),
                              formatPT(upper.tau, digits = digits,
                                       lab.NA = "", big.mark = big.mark))),
           #
           I2 = if_else(is.na(I2), "",
                        paste0(formatN(100 * I2, digits = digits.I2,
                                              text.NA = ""), "%")),
           lower.I2 = 
             if_else(is.na(lower.I2) & is.na(upper.I2), "",
                     formatCI(formatPT(lower.I2, digits = digits,
                                       lab.NA = "", big.mark = big.mark),
                              formatPT(upper.I2, digits = digits,
                                       lab.NA = "", big.mark = big.mark)))
    ) %>%
    select(-upper, -upper.predict, -upper.tau2, -upper.tau, -upper.I2)
  #
  names(dat)[names(dat) == "TE"] <- x$sm
  #
  names(dat)[names(dat) == "lower"] <-
    paste0(round(100 * x$level.ma, 1), "%-CI")
  #
  if (prediction)
    names(dat)[names(dat) == "lower.predict"] <-
    paste0(round(100 * x$level.predict, 1), "%-PI")
  else
    dat$lower.predict <- NULL
  #
  names(dat)[names(dat) == "pval"] <- "p-value"
  #
  if (x$pooled == "random" & x$method.random.ci %in% c("HK", "KR"))
    names(dat)[names(dat) == "statistic"] <- "t"
  else
    names(dat)[names(dat) == "statistic"] <- "z"
  #
  if (print.tau2)
    names(dat)[names(dat) == "tau2"] <- text.tau2
  else {
    dat$tau2 <- NULL
    dat$lower.tau2 <- NULL
  }
  #
  if (print.tau2.ci)
    names(dat)[names(dat) == "lower.tau2"] <- "x%-CI"
  else
    dat$lower.tau2 <- NULL
  #
  if (print.tau)
    names(dat)[names(dat) == "tau"] <- text.tau
  else
    dat$tau <- NULL
  #
  if (print.tau.ci)
    names(dat)[names(dat) == "lower.tau"] <- "x%-CI"
  else
    dat$lower.tau <- NULL
  #
  if (print.I2)
    names(dat)[names(dat) == "I2"] <- text.I2
  else
    dat$I2 <- NULL
  #
  if (print.I2.ci)
    names(dat)[names(dat) == "lower.I2"] <- "x%-CI"
  else
    dat$lower.I2 <- NULL
  #
  dat <- replaceNA(dat, "")
  
  cat(paste0("Cumulative meta-analysis (",
             if (x$pooled == "common") "common effect" else "random effects",
             " model)\n\n"))
  #
  prmatrix(dat, quote = FALSE, right = TRUE, ...)
  
  if (details.methods)
    details <-
    catmeth(x$x,
            x$pooled == "common", x$pooled == "random", prediction,
            TRUE, TRUE,
            #
            func.transf = x$func.transf,
            backtransf = x$backtransf,
            func.backtransf = x$func.backtransf,
            #
            big.mark = big.mark, digits = digits,
            digits.tau = digits.tau,
            text.tau = text.tau, text.tau2 = text.tau2,
            #
            print.tau2 = print.tau2,
            print.tau2.ci = print.tau2.ci,
            print.tau = print.tau,
            print.tau.ci = print.tau.ci,
            #
            print.I2 = print.I2, text.I2 = text.I2,
            #
            print.df = TRUE, prediction.subgroup = FALSE)
  
  
  invisible(NULL)
}
