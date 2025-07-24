#' Print results of a cumulative meta-analysis
#' 
#' @description
#' Print results of a cumulative meta-analysis
#' 
#' @aliases print.metacum
#' 
#' @param x An object of class \code{\link{metacum}}.
#' @param prediction A logical indicating whether prediction intervals should
#'   be printed.
#' @param overall A logical indicating whether overall results should be
#'   printed.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf=TRUE}, results
#'   for \code{sm="OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param lab.NA A character string to label missing values.
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
#' @param digits.cid Minimal number of significant digits for
#'   CID / decision thresholds, see \code{print.default}.
#' @param digits.percent Minimal number of significant digits for
#'   probabilities, printed as percentages, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param print.stat A logical value indicating whether z- or t-value
#'   for test of treatment effect should be printed.
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
#' @param print.prob A logical specifying whether to print probabilities
#'   of clinically important benefit or harm.
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

print.metacum <- function(x,
                          #
                          prediction = x$prediction,
                          overall = x$overall,
                          backtransf = x$backtransf,
                          header = TRUE,
                          #
                          lab.NA = ".",
                          #
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.tau = gs("digits.tau"),
                          digits.I2 = gs("digits.I2"),
                          digits.cid = gs("digits.cid"),
                          digits.percent = 1,
                          #
                          big.mark = gs("big.mark"),
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
                          #
                          print.stat = FALSE,
                          print.tau2 = TRUE,
                          print.tau2.ci = FALSE,
                          print.tau = TRUE,
                          print.tau.ci = FALSE,
                          print.I2 = TRUE,
                          print.I2.ci = FALSE,
                          print.prob = TRUE,
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
  
  chkclass(x, c("metacum", "metainf"))
  x <- updateversion(x)
  #
  chklogical(prediction)
  chklogical(overall)
  chklogical(backtransf)
  chklogical(header)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.cid, min = 0, length = 1)
  chknumeric(digits.percent, min = 0, length = 1)
  #
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(print.stat)
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
  chklogical(print.prob)
  #
  chklogical(details.methods)
  #
  print.tau2.ci <- print.tau2 & print.tau2.ci
  print.tau.ci <- print.tau & print.tau.ci
  print.I2.ci <- print.I2 & print.I2.ci
  #
  sm <- x$sm
  #
  if (is.function(x$func.backtransf))
    fbt <- deparse(substitute(x$func.backtransf))
  else
    fbt <- x$func.backtransf
  ##
  abt <- x$args.backtransf
  
  
  #
  #
  # (3) Create data set
  #
  #
  
  # Get rid of warning 'no visible binding for global variable'
  #
  TE <- lower <- upper <- statistic <- pval <-
    tau2 <- lower.tau2 <- upper.tau2 <-
    tau <- lower.tau <- upper.tau <-
    I2 <- lower.I2 <- upper.I2 <-
    lower.predict <- upper.predict <-
    prop.cid.below.null <- prop.cid.above.null <-
    harmonic.mean <- n.harmonic.mean <- t.harmonic.mean <- NULL
  #
  dat <-
    with(x,
         data.frame(TE, lower, upper,
                    statistic, pval,
                    tau2, lower.tau2, upper.tau2,
                    tau, lower.tau, upper.tau,
                    I2, lower.I2, upper.I2,
                    lower.predict, upper.predict,
                    n.harmonic.mean,
                    t.harmonic.mean,
                    row.names = studlab))
  #
  avail.prop.cid.below.null <-
    !is.null(x$prop.cid.below.null) && !all(is.na(x$prop.cid.below.null))
  #
  avail.prop.cid.above.null <-
    !is.null(x$prop.cid.above.null) && !all(is.na(x$prop.cid.above.null))
  #
  avail.prop.cid <- avail.prop.cid.below.null | avail.prop.cid.above.null
  #
  if (avail.prop.cid.below.null)
    dat$prop.cid.below.null <- x$prop.cid.below.null
  else
    dat$prop.cid.below.null <- NA
  #
  if (avail.prop.cid.above.null)
    dat$prop.cid.above.null <- x$prop.cid.above.null
  else
    dat$prop.cid.above.null <- NA
  #
  if (overall) {
    space <- data.frame(TE = NA, lower = NA, upper = NA,
                        statistic = NA, pval = NA,
                        tau2 = NA, lower.tau2 = NA, upper.tau2 = NA,
                        tau = NA, lower.tau = NA, upper.tau = NA, 
                        I2 = NA, lower.I2 = NA, upper.I2 = NA,
                        lower.predict = NA, upper.predict = NA,
                        prop.cid.below.null = NA, prop.cid.above.null = NA,
                        n.harmonic.mean = NA, t.harmonic.mean = NA,
                        row.names = "")
    #
    dat.pooled <-
      with(x,
           data.frame(TE.pooled, lower.pooled, upper.pooled,
                      statistic.pooled, pval.pooled,
                      tau2.pooled, lower.tau2.pooled, upper.tau2.pooled,
                      tau.pooled, lower.tau.pooled, upper.tau.pooled,
                      I2.pooled, lower.I2.pooled, upper.I2.pooled,
                      lower.predict.pooled, upper.predict.pooled,
                      n.harmonic.mean.pooled, t.harmonic.mean.pooled)) %>%
      rename_with(~ gsub(".pooled", "", .x, fixed = TRUE))
    #
    rownames(dat.pooled) <- x$text.pooled
    #
    if (avail.prop.cid.below.null)
      dat.pooled$prop.cid.below.null <- x$prop.cid.below.null.pooled
    else
      dat.pooled$prop.cid.below.null <- NA
    #
    if (avail.prop.cid.above.null)
      dat.pooled$prop.cid.above.null <- x$prop.cid.above.null.pooled
    else
      dat.pooled$prop.cid.above.null <- NA
    #
    dat <- rbind(dat, space, dat.pooled)
  }
  #
  if (sm == "IRFT")
    dat %<>% rename(harmonic.mean = t.harmonic.mean) %>%
      select(-n.harmonic.mean)
  else
    dat %<>% rename(harmonic.mean = n.harmonic.mean) %>%
      select(-t.harmonic.mean)
  #
  if (backtransf) {
    dat %<>%
      mutate(
        TE = backtransf(TE, sm, harmonic.mean, harmonic.mean, fbt, abt),
        lower = backtransf(lower, sm, harmonic.mean, harmonic.mean, fbt, abt),
        upper = backtransf(upper, sm, harmonic.mean, harmonic.mean, fbt, abt),
        lower.predict =
          backtransf(lower.predict, sm, harmonic.mean, harmonic.mean, fbt, abt),
        upper.predict =
          backtransf(upper.predict, sm, harmonic.mean, harmonic.mean, fbt, abt),
      )    
  }
  #
  dat %<>% select(-harmonic.mean)
  #
  dat %<>%
    mutate(TE = formatN(TE, digits = digits, text.NA = lab.NA, big.mark = big.mark),
           lower = 
             if_else(is.na(lower) & is.na(upper), lab.NA,
                     formatCI(formatN(lower, digits = digits, text.NA = lab.NA,
                                      big.mark = big.mark),
                              formatN(upper, digits = digits, text.NA = lab.NA,
                                      big.mark = big.mark))),
           #
           statistic = formatN(statistic, digits = digits.stat, text.NA = lab.NA),
           #
           pval = formatPT(pval, digits = digits.pval, lab.NA = lab.NA,
                           scientific = scientific.pval,
                           zero = zero.pval, JAMA = JAMA.pval),
           #
           tau2 = formatPT(tau2, digits = digits.tau2, lab.NA = lab.NA,
                           big.mark = big.mark),
           lower.tau2 = 
             if_else(is.na(lower.tau2) & is.na(upper.tau2), lab.NA,
                     formatCI(formatPT(lower.tau2, digits = digits,
                                       lab.NA = lab.NA, big.mark = big.mark),
                              formatPT(upper.tau2, digits = digits,
                                       lab.NA = lab.NA, big.mark = big.mark))),
           #
           tau = formatPT(tau, digits = digits.tau, lab.NA = lab.NA,
                          big.mark = big.mark),
           lower.tau = 
             if_else(is.na(lower.tau) & is.na(upper.tau), lab.NA,
                     formatCI(formatPT(lower.tau, digits = digits,
                                       lab.NA = lab.NA, big.mark = big.mark),
                              formatPT(upper.tau, digits = digits,
                                       lab.NA = lab.NA, big.mark = big.mark))),
           #
           I2 = if_else(is.na(I2), lab.NA,
                        paste0(formatPT(100 * I2, digits = digits.I2,
                                        lab.NA = lab.NA), "%")),
           #
           lower.I2 =
             if_else(is.na(lower.I2), lab.NA,
                     paste0(formatPT(100 * lower.I2, digits = digits.I2,
                                     lab.NA = lab.NA), "%")),
           #
           upper.I2 =
             if_else(is.na(upper.I2), lab.NA,
                     paste0(formatPT(100 * upper.I2, digits = digits.I2,
                                     lab.NA = lab.NA), "%")),
           #
           lower.I2 =
             if_else(lower.I2 == lab.NA & upper.I2 == lab.NA, lab.NA,
                     formatCI(lower.I2, upper.I2)),
           #
           lower.predict = 
             if_else(is.na(lower.predict) & is.na(upper.predict), lab.NA,
                     formatCI(formatN(lower.predict, digits = digits,
                                      text.NA = lab.NA, big.mark = big.mark),
                              formatN(upper.predict, digits = digits,
                                      text.NA = lab.NA, big.mark = big.mark))),
           #
           prop.cid.below.null =
             if_else(is.na(prop.cid.below.null), lab.NA,
                     paste0(formatPT(100 * prop.cid.below.null,
                                     digits = digits.percent,
                                     lab.NA = lab.NA), "%")),
           #
           prop.cid.above.null =
             if_else(is.na(prop.cid.above.null), lab.NA,
                     paste0(formatPT(100 * prop.cid.above.null,
                                     digits = digits.percent,
                                     lab.NA = lab.NA), "%"))
    ) %>%
    select(-upper, -upper.predict, -upper.tau2, -upper.tau, -upper.I2)
  #
  # Use "" instead of lab.NA in empty row
  #
  sel.empty <- with(dat, TE == lab.NA & lower == lab.NA & I2 == lab.NA)
  #
  if (sum(sel.empty) > 0) {
    for (i in names(dat))
      dat[sel.empty, i] <- ""
  }
  #
  names(dat)[names(dat) == "TE"] <- smlab(sm, backtransf, x$pscale, x$irscale)
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
  if (all(is.na(c(x$pval, x$pval.pooled))))
    dat$pval <- NULL
  else
    names(dat)[names(dat) == "pval"] <- "p-value"
  #
  if (print.stat) {
    if (x$pooled == "random" & x$method.random.ci %in% c("HK", "KR"))
      names(dat)[names(dat) == "statistic"] <- "t"
    else
      names(dat)[names(dat) == "statistic"] <- "z"
  }
  else
    dat$statistic <- NULL
  #
  if (print.tau2)
    names(dat)[names(dat) == "tau2"] <- text.tau2
  else {
    dat$tau2 <- NULL
    dat$lower.tau2 <- NULL
  }
  #
  if (print.tau2.ci)
    names(dat)[names(dat) == "lower.tau2"] <-
    paste0(round(100 * x$level.hetstat, 1), "%-CI")
  else
    dat$lower.tau2 <- NULL
  #
  if (print.tau)
    names(dat)[names(dat) == "tau"] <- text.tau
  else
    dat$tau <- NULL
  #
  if (print.tau.ci)
    names(dat)[names(dat) == "lower.tau"] <-
    paste0(round(100 * x$level.hetstat, 1), "%-CI")
  else
    dat$lower.tau <- NULL
  #
  if (print.I2)
    names(dat)[names(dat) == "I2"] <- text.I2
  else
    dat$I2 <- NULL
  #
  if (print.I2.ci)
    names(dat)[names(dat) == "lower.I2"] <-
    paste0(round(100 * x$level.hetstat, 1), "%-CI")
  else
    dat$lower.I2 <- NULL
  #
  if (!avail.prop.cid.below.null)
    dat$prop.cid.below.null <- NULL
  #
  if (!avail.prop.cid.above.null)
    dat$prop.cid.above.null <- NULL
  #
  if (print.prob & avail.prop.cid.below.null) {
    names(dat)[names(dat) == "prop.cid.below.null"] <-
      paste0("P(",
             if (x$small.values == "desirable") "benefit" else "harm",
             ")")
  }
  else
    dat$prop.cid.below.null <- NULL
  #
  if (print.prob & avail.prop.cid.above.null) {
    names(dat)[names(dat) == "prop.cid.above.null"] <-
      paste0("P(",
             if (x$small.values == "desirable") "harm" else "benefit",
             ")")
  }
  else
    dat$prop.cid.above.null <- NULL
  #
  dat <- replaceNA(dat, "")
  
  
  #
  #
  # (4) Print results
  #
  #
  
  if (header)
    crtitle(x)
  #
  cat(paste0(if (inherits(x, "metacum")) "Cumulative " else "Leave-one-out ",
             "meta-analysis\n\n"))
  #
  prmatrix(dat, quote = FALSE, right = TRUE, ...)
  #
  if (details.methods) {
    details <-
      catmeth(x,
              x$pooled == "common", x$pooled == "random", prediction,
              overall, TRUE,
              #
              func.transf = x$func.transf,
              backtransf = backtransf,
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
    #
    if (avail.prop.cid)
      svd <- x$small.values == "desirable"
    #
    if (avail.prop.cid.below.null)
      cat(paste0("- Lower decision threshold (",
                 if (svd) "beneficial " else "harmful ",
                 "effects): ",
                 formatN(x$cid.below.null, digits = digits.cid, big.mark = big.mark),
                 "\n"))
    #
    if (avail.prop.cid.above.null)
      cat(paste0("- Upper decision threshold (",
                 if (svd) "harmful " else "beneficial ",
                 "effects): ",
                 formatN(x$cid.above.null, digits = digits.cid, big.mark = big.mark),
                 "\n"))
  }
  #
  invisible(NULL)
}
