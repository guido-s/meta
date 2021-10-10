#' Print meta-analysis results
#' 
#' @description
#' Print method for objects of class \code{meta}.
#'
#' R function cilayout can be utilised to change the layout to print
#' confidence intervals (both in printout from print.meta and
#' print.summary.meta function as well as in forest plots). The
#' default layout is "[lower; upper]". Another popular layout is
#' "(lower - upper)" which is used throughout an R session by using R
#' command \code{cilayout("(", " - ")}.
#' 
#' Argument \code{pscale} can be used to rescale single proportions or
#' risk differences, e.g. \code{pscale = 1000} means that proportions
#' are expressed as events per 1000 observations. This is useful in
#' situations with (very) low event probabilities.
#' 
#' Argument \code{irscale} can be used to rescale single rates or rate
#' differences, e.g. \code{irscale = 1000} means that rates are
#' expressed as events per 1000 time units, e.g. person-years. This is
#' useful in situations with (very) low rates. Argument \code{irunit}
#' can be used to specify the time unit used in individual studies
#' (default: "person-years"). This information is printed in summaries
#' and forest plots if argument \code{irscale} is not equal to 1.
#'
#' @aliases print.meta cilayout
#' 
#' @param x An object of class \code{meta}.
#' @param fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param overall A logical indicating whether overall summaries
#'   should be reported. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be reported.
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param test.subgroup.fixed A logical value indicating whether to
#'   print results of test for subgroup differences (based on fixed
#'   effect / common effect model).
#' @param test.subgroup.random A logical value indicating whether to
#'   print results of test for subgroup differences (based on random
#'   effects model).
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf=TRUE}, results
#'   for \code{sm="OR"} are printed as odds ratios rather than log
#'   odds ratios and results for \code{sm="ZCOR"} are printed as
#'   correlations rather than Fisher's z transformed correlations, for
#'   example.
#' @param pscale A numeric giving scaling factor for printing of
#'   single event probabilities or risk differences, i.e. if argument
#'   \code{sm} is equal to \code{"PLOGIT"}, \code{"PLN"},
#'   \code{"PRAW"}, \code{"PAS"}, \code{"PFT"}, or \code{"RD"}.
#' @param irscale A numeric defining a scaling factor for printing of
#'   single incidence rates or incidence rate differences, i.e. if
#'   argument \code{sm} is equal to \code{"IR"}, \code{"IRLN"},
#'   \code{"IRS"}, \code{"IRFT"}, or \code{"IRD"}.
#' @param irunit A character specifying the time unit used to
#'   calculate rates, e.g. person-years.
#' @param subgroup.name A character string with a name for the
#'   grouping variable.
#' @param print.subgroup.name A logical indicating whether the name of
#'   the grouping variable should be printed in front of the group
#'   labels.
#' @param sep.subgroup A character string defining the separator
#'   between label and levels of grouping variable.
#' @param nchar.subgroup A numeric specifying the number of characters
#'   to print from subgroup labels.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value of test for overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param digits.H Minimal number of significant digits for H
#'   statistic, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param print.I2 A logical specifying whether heterogeneity
#'   statistic I\eqn{^2} should be printed.
#' @param print.H A logical specifying whether heterogeneity statistic
#'   H should be printed.
#' @param print.Rb A logical specifying whether heterogeneity
#'   statistic R\eqn{_b} should be printed.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param text.Rb Text printed to identify heterogeneity statistic
#'   R\eqn{_b}.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param warn.backtransf A logical indicating whether a warning
#'   should be printed if backtransformed proportions and rates are
#'   below 0 and backtransformed proportions are above 1.
#' @param bracket A character with bracket symbol to print lower
#'   confidence interval: "[", "(", "\{", "".
#' @param separator A character string with information on separator
#'   between lower and upper confidence interval.
#' @param \dots Additional arguments (passed on to \code{prmatrix}).
#' 
#' @rdname print.meta
#' @method print meta
#' @export


print.meta <- function(x,
                       fixed = x$fixed,
                       random = x$random,
                       prediction = x$prediction,
                       overall = x$overall,
                       overall.hetstat = x$overall.hetstat,
                       ##
                       test.subgroup = x$test.subgroup,
                       test.subgroup.fixed = test.subgroup & fixed,
                       test.subgroup.random = test.subgroup & random,
                       ##
                       backtransf = x$backtransf,
                       pscale = x$pscale,
                       irscale = x$irscale,
                       irunit = x$irunit,
                       ##
                       subgroup.name = x$subgroup.name,
                       print.subgroup.name = x$print.subgroup.name,
                       sep.subgroup = x$sep.subgroup,
                       nchar.subgroup = 35,
                       ##
                       header = TRUE,
                       print.CMH = x$print.CMH,
                       ##
                       digits = gs("digits"),
                       digits.stat = gs("digits.stat"),
                       digits.pval = max(gs("digits.pval"), 2),
                       digits.pval.Q = max(gs("digits.pval.Q"), 2),
                       digits.Q = gs("digits.Q"),
                       digits.tau2 = gs("digits.tau2"),
                       digits.tau = gs("digits.tau"),
                       digits.H = gs("digits.H"),
                       digits.I2 = gs("digits.I2"),
                       ##
                       scientific.pval = gs("scientific.pval"),
                       big.mark = gs("big.mark"),
                       zero.pval = gs("zero.pval"),
                       JAMA.pval = gs("JAMA.pval"),
                       print.tau2 = TRUE,
                       print.tau = TRUE,
                       print.I2 = gs("print.I2"),
                       print.H = gs("print.H"),
                       print.Rb = gs("print.Rb"),
                       ##
                       text.tau2 = gs("text.tau2"),
                       text.tau = gs("text.tau"),
                       text.I2 = gs("text.I2"),
                       text.Rb = gs("text.Rb"),
                       ##
                       details.methods = TRUE,
                       ##
                       warn.backtransf = FALSE,
                       ...) {
  
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  if (inherits(x, "metacum") | inherits(x, "metainf"))
    return(invisible(NULL))
  ##
  is.metabind <- inherits(x, "metabind")
  is.netpairwise <- inherits(x, "netpairwise")
  
  
  ##
  ##
  ## (2) Check and set other arguments
  ##
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.H, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  chklogical(print.tau2)
  chklogical(print.tau)
  chklogical(print.I2)
  chklogical(print.H)
  chklogical(print.Rb)
  chkchar(text.tau2, length = 1)
  chkchar(text.tau, length = 1)
  chkchar(text.I2, length = 1)
  chkchar(text.Rb, length = 1)
  chklogical(warn.backtransf)
  is.prop <- is.prop(x$sm)
  is.rate <- is.rate(x$sm)
  ##
  if (!is.prop & x$sm != "RD")
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  if (!is.rate & x$sm != "IRD")
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  ##
  chklogical(fixed)
  chklogical(random)
  chklogical(prediction)
  overall <- replaceNULL(overall, TRUE)
  chklogical(overall)
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  test.subgroup <- replaceNULL(test.subgroup, TRUE)
  if (is.na(test.subgroup))
    test.subgroup <- FALSE
  chklogical(test.subgroup)
  test.subgroup.fixed <- replaceNULL(test.subgroup.fixed, test.subgroup)
  chklogical(test.subgroup.fixed)
  test.subgroup.random <- replaceNULL(test.subgroup.random, test.subgroup)
  chklogical(test.subgroup.random)
  ##
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  chklogical(header)
  ##
  chklogical(details.methods)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ##
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed", FALSE)
  chklogical(fixed)
  ##
  random <- deprecated(random, missing(random), args, "comb.random", FALSE)
  chklogical(random)
  ##
  subgroup.name <-
    deprecated(subgroup.name, missing(subgroup.name), args, "bylab", FALSE)
  ##
  print.subgroup.name <-
    deprecated(print.subgroup.name, missing(print.subgroup.name),
               args, "print.byvar", FALSE)
  ##
  sep.subgroup <-
    deprecated(sep.subgroup, missing(sep.subgroup), args, "byseparator",
               FALSE)
  ##
  nchar.subgroup <-
    deprecated(nchar.subgroup, missing(nchar.subgroup), args, "nchar.subgroup")
  chknumeric(nchar.subgroup, min = 1, length = 1)
  ##
  backtransf <-
    deprecated(backtransf, missing(backtransf), args, "logscale")
  if (is.untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  ##
  ## Additional settings
  ##
  if (!backtransf & pscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!backtransf & irscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  by <- !is.null(subgroup.name)
  if (by) {
    chklogical(print.subgroup.name)
    chkchar(sep.subgroup)
  }
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  k.all <- replaceNULL(x$k.all, length(x$TE))
  k <- x$k
  k.study <- replaceNULL(x$k.study, k)
  sm <- x$sm
  ##
  bip <- inherits(x, c("metabin", "metainc", "metaprop", "metarate"))
  metabin <- inherits(x, "metabin")
  metaprop <- inherits(x, "metaprop")
  ##
  print.ci <- attr(x, ".print.study.results.")
  if (!is.null(print.ci) && print.ci) {
    method.ci <- x$method.ci
    if (metaprop & !backtransf)
      method.ci <- "NAsm"
  }
  else
    method.ci <- NULL
  ##
  null.effect <- x$null.effect
  null.given <- !is.null(null.effect) && !is.na(null.effect)
  ##
  if (null.given & !backtransf) {
    ##
    if (sm %in% c("PFT", "PAS"))
      null.effect <- asin(sqrt(null.effect))
    else if (is.log.effect(sm))
      null.effect <- log(null.effect)
    else if (sm == c("PLOGIT"))
      null.effect <- log(null.effect / (1 - null.effect))
    else if (sm %in% c("IRS", "IRFT"))
      null.effect <- sqrt(null.effect)
    else if (sm == "ZCOR")
      null.effect <- 0.5 * log((1 + null.effect) / (1 - null.effect))
  }
  ##
  if (by) {
    k.w <- x$k.w
    ##
    prediction.w <- prediction & k.w >= 3
    prediction.w[is.na(prediction.w)] <- FALSE
    prediction.w <- any(prediction.w)
  }
  ##
  prediction <- prediction & k >= 3
  if (is.na(prediction))
    prediction <- FALSE
  ##
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is.mean(sm))
      sm.lab <- "mean"
    else if (is.prop) {
      if (pscale == 1)
        sm.lab <- "proportion"
      else
        sm.lab <- "events"
    }
    else if (is.rate) {
      if (irscale == 1)
        sm.lab <- "rate"
      else
        sm.lab <- "events"
    }
  }
  else
    if (is.relative.effect(sm))
      sm.lab <- paste0("log", sm)
  ##
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  ##
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  ##
  if (by)
    bylevs <- ifelse(nchar(x$bylevs) > nchar.subgroup,
                     paste0(substring(x$bylevs, 1, nchar.subgroup - 4), " ..."),
                     x$bylevs)
  ##
  if (is.null(x$text.fixed))
    text.fixed <- gs("text.fixed")
  else
    text.fixed <- x$text.fixed
  ##
  if (is.null(x$text.random))
    text.random <- gs("text.random")
  else
    text.random <- x$text.random
  ##
  if (is.null(x$text.predict))
    text.predict <- gs("text.predict")
  else
    text.predict <- x$text.predict
  ##
  if (substring(text.fixed, 1, 5) %in% c("Fixed", "Commo", "Rando")) {
    text.fixed.br <- tolower(text.fixed)
    text.random.br <- tolower(text.random)
  }
  else {
    text.fixed.br <- text.fixed
    text.random.br <- text.random
  }
  ##
  ci.lab <- paste0(round(100 * x$level.ma, 1), "%-CI")
  
  
  ##
  ##
  ## (4) Set and backtransform results of meta-analysis
  ##
  ##
  TE.fixed    <- x$TE.fixed
  lowTE.fixed <- x$lower.fixed
  uppTE.fixed <- x$upper.fixed
  ##
  TE.random    <- x$TE.random
  lowTE.random <- x$lower.random
  uppTE.random <- x$upper.random
  ##
  lowTE.predict <- x$lower.predict
  uppTE.predict <- x$upper.predict
  ##
  Q <- x$Q
  df.Q <- replaceNULL(x$df.Q, k - 1)
  pval.Q <- replaceNULL(x$pval.Q, pvalQ(Q, df.Q))
  ##
  if (!is.null(x$Q.CMH)) {
    Q.CMH <- x$Q.CMH
    df.Q.CMH <- replaceNULL(x$df.Q.CMH, 1)
    pval.Q.CMH <- replaceNULL(x$pval.Q.CMH, pvalQ(Q.CMH, df.Q.CMH))
  }
  ##
  if (x$method == "GLMM") {
    Q.LRT <- x$Q.LRT
    df.Q.LRT <- replaceNULL(x$df.Q.LRT, df.Q)
    pval.Q.LRT <- replaceNULL(x$pval.Q.LRT, pvalQ(Q.LRT, df.Q.LRT))
  }
  ##
  if (by) {
    TE.fixed.w      <- x$TE.fixed.w
    lowTE.fixed.w   <- x$lower.fixed.w
    uppTE.fixed.w   <- x$upper.fixed.w
    pval.fixed.w    <- x$pval.fixed.w
    harmonic.mean.w <- x$n.harmonic.mean.w
    ##
    TE.random.w     <- x$TE.random.w
    lowTE.random.w  <- x$lower.random.w
    uppTE.random.w  <- x$upper.random.w
    pval.random.w   <- x$pval.random.w
    ##
    lowTE.predict.w <- x$lower.predict.w
    uppTE.predict.w <- x$upper.predict.w
    ##
    Q.b.fixed <- x$Q.b.fixed
    Q.w.fixed <- x$Q.w.fixed
    Q.b.random <- x$Q.b.random
    Q.w.random <- x$Q.w.random
    ##
    Q.w <- x$Q.w
    ##
    df.Q.w <- replaceNULL(x$df.Q.w, sum((k.w - 1)[!is.na(x$Q.w)]))
    df.Q.b <- replaceNULL(x$df.Q.b, (k - 1) - sum((k.w - 1)[!is.na(x$Q.w)]))
    ##
    pval.Q.b.fixed  <- replaceNULL(x$pval.Q.b.fixed, pvalQ(Q.b.fixed, df.Q.b))
    pval.Q.w.fixed  <- replaceNULL(x$pval.Q.w.fixed, pvalQ(Q.w.fixed, df.Q.w))
    pval.Q.b.random <- replaceNULL(x$pval.Q.b.random, pvalQ(Q.b.random, df.Q.b))
    pval.Q.w.random <- replaceNULL(x$pval.Q.w.random, pvalQ(Q.w.random, df.Q.w))
  }
  ##
  if (backtransf) {
    if (sm == "IRFT") {
      if (is.metabind)
        harmonic.mean <- x$t.harmonic.mean.ma
      else
        harmonic.mean <- 1 / mean(1 / x$time)
    }
    else {
      if (is.metabind)
        harmonic.mean <- x$n.harmonic.mean.ma
      else
        harmonic.mean <- 1 / mean(1 / x$n)
    }
    ##
    TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                              harmonic.mean,
                              warn = overall & fixed & warn.backtransf)
    lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                              harmonic.mean,
                              warn = overall & fixed & warn.backtransf)
    uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                              harmonic.mean,
                              warn = overall & fixed & warn.backtransf)
    ##
    TE.random <- backtransf(TE.random, sm, "mean",
                            harmonic.mean,
                            warn = overall & random & warn.backtransf)
    lowTE.random <- backtransf(lowTE.random, sm, "lower",
                               harmonic.mean,
                               warn = overall & random & warn.backtransf)
    uppTE.random <- backtransf(uppTE.random, sm, "upper",
                               harmonic.mean,
                               warn = overall & random & warn.backtransf)
    ##
    lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                harmonic.mean,
                                warn = overall & prediction & warn.backtransf)
    uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                harmonic.mean,
                                warn = overall & prediction & warn.backtransf)
    ##
    if (by) {
      TE.fixed.w     <- backtransf(TE.fixed.w, sm, "mean",
                                   harmonic.mean.w,
                                   warn = overall & fixed &
                                     warn.backtransf)
      lowTE.fixed.w  <- backtransf(lowTE.fixed.w, sm, "lower",
                                   harmonic.mean.w,
                                   warn = overall & fixed &
                                     warn.backtransf)
      uppTE.fixed.w  <- backtransf(uppTE.fixed.w, sm, "upper",
                                   harmonic.mean.w,
                                   warn = overall & fixed &
                                     warn.backtransf)
      ##
      TE.random.w    <- backtransf(TE.random.w, sm, "mean",
                                   harmonic.mean.w,
                                   warn = overall & random &
                                     warn.backtransf)
      lowTE.random.w <- backtransf(lowTE.random.w, sm, "lower",
                                   harmonic.mean.w,
                                   warn = overall & random &
                                     warn.backtransf)
      uppTE.random.w <- backtransf(uppTE.random.w, sm, "upper",
                                   harmonic.mean.w,
                                   warn = overall & random &
                                     warn.backtransf)
      ##
      if (prediction.w) {
        lowTE.predict.w <- backtransf(lowTE.predict.w, sm, "lower",
                                      harmonic.mean.w,
                                      warn = warn.backtransf)
        uppTE.predict.w <- backtransf(uppTE.predict.w, sm, "upper",
                                      harmonic.mean.w,
                                      warn = warn.backtransf)
      }
    }
  }
  ##
  ## Apply argument 'pscale' to proportions / risk differences and
  ## 'irscale' to rates / incidence rate differences
  ##
  if (is.prop | sm == "RD" | is.rate | sm == "IRD") {
    if (is.prop | sm == "RD")
      scale <- pscale
    else if (is.rate | sm == "IRD")
      scale <- irscale
    ##
    TE.fixed    <- scale * TE.fixed
    lowTE.fixed <- scale * lowTE.fixed
    uppTE.fixed <- scale * uppTE.fixed
    ##
    TE.random    <- scale * TE.random
    lowTE.random <- scale * lowTE.random
    uppTE.random <- scale * uppTE.random
    ##
    lowTE.predict <- scale * lowTE.predict
    uppTE.predict <- scale * uppTE.predict
    ##
    if (by) {
      TE.fixed.w    <- scale * TE.fixed.w
      lowTE.fixed.w <- scale * lowTE.fixed.w
      uppTE.fixed.w <- scale * uppTE.fixed.w
      ##   
      TE.random.w    <- scale * TE.random.w
      lowTE.random.w <- scale * lowTE.random.w
      uppTE.random.w <- scale * uppTE.random.w
      ##   
      if (prediction.w) {
        lowTE.predict.w <- scale * lowTE.predict.w
        uppTE.predict.w <- scale * uppTE.predict.w
      }
    }
  }
  ##
  ## Round and round ...
  ##
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed <- x$pval.fixed
  sTE.fixed <- round(x$statistic.fixed, digits.stat)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$pval.random
  sTE.random <- round(x$statistic.random, digits.stat)
  ##
  lowTE.predict <- round(lowTE.predict, digits)
  uppTE.predict <- round(uppTE.predict, digits)
  ##
  if (by) {
    TE.fixed.w     <- round(TE.fixed.w, digits)
    lowTE.fixed.w  <- round(lowTE.fixed.w, digits)
    uppTE.fixed.w  <- round(uppTE.fixed.w, digits)
    ##
    TE.random.w    <- round(TE.random.w, digits)
    lowTE.random.w <- round(lowTE.random.w, digits)
    uppTE.random.w <- round(uppTE.random.w, digits)
    ##
    if (prediction.w) {
      lowTE.predict.w <- round(lowTE.predict.w, digits)
      uppTE.predict.w <- round(uppTE.predict.w, digits)
    }
    ##
    if (print.I2)
      I2.w <- round(100 * x$I2.w, digits.I2)
    ##
    if (print.Rb)
      Rb.w <- round(100 * x$Rb.w, digits.I2)
  }
  ##
  if (print.H) {
    H <- round(x$H, digits.H)
    lowH <- round(x$lower.H, digits.H)
    uppH <- round(x$upper.H, digits.H)
    if (all(!is.na(lowH) & !is.na(uppH)) && all(lowH == uppH)) {
      lowH <- NA
      uppH <- NA
    }
  }
  ##
  if (print.I2) {
    I2 <- round(100 * x$I2, digits.I2)
    lowI2 <- round(100 * x$lower.I2, digits.I2)
    uppI2 <- round(100 * x$upper.I2, digits.I2)
    print.I2.ci <- ((Q > k & k >= 2) | (Q <= k & k > 2)) &
      !(is.na(lowI2) | is.na(uppI2))
    if (is.na(print.I2.ci))
      print.I2.ci <- FALSE
    if (print.I2.ci && all(lowI2 == uppI2))
      print.I2.ci <- FALSE
  }
  else
    print.I2.ci <- FALSE
  ##
  if (print.Rb) {
    Rb <- round(100 * x$Rb, digits.I2)
    lowRb <- round(100 * x$lower.Rb, digits.I2)
    uppRb <- round(100 * x$upper.Rb, digits.I2)
    if (all(!is.na(lowRb) & !is.na(uppRb)) && all(lowRb == uppRb)) {
      lowRb <- NA
      uppRb <- NA
    }
  }
  ##
  three.level <- if (is.null(x$three.level)) FALSE else x$three.level
  is.glmm <- x$method == "GLMM"
  ##
  catobsev <- function(var1, var2 = NULL, type = "n", addrow = FALSE,
                       big.mark = gs("big.mark")) {
    if (type == "n") {
      txt <- "observations"
      idx <- "o"
    }
    else if (type == "e") {
      txt <- "events"
      idx <- "e"
    }
    ##
    if (!is.null(var1) & !is.null(var2)) {
      sum1 <- sum(var1, na.rm = TRUE)
      sum2 <- sum(var2, na.rm = TRUE)
      ##
      cat(paste0("Number of ", txt, ": ", idx, " = ",
                 format(sum1 + sum2, big.mark = big.mark),
                 ##" (", idx, ".e = ",
                 ##format(sum1, big.mark = big.mark),
                 ##", ", idx, ".c = ",
                 ##format(sum2, big.mark = big.mark),
                 ##")",
                 "\n"))
    }
    else if (!is.null(var1)) {
      cat(paste0("Number of ", txt, ": ", idx, " = ",
                 format(sum(var1, na.rm = TRUE),
                          big.mark = big.mark),
                 "\n"))
    }
    else if (!is.null(var2)) {
      cat(paste0("Number of ", txt, ": ", idx, " = ",
                 format(sum(var2, na.rm = TRUE), big.mark = big.mark),
                 "\n"))
    }
    ##
    if (addrow)
      cat("\n")
    ##
    invisible(NULL)
  }
  ##
  sel.n <- inherits(x, c("metacor", "metaprop", "metamean", "metarate"))
  ##
  sel.ev <- inherits(x, "metaprop")
  
  
  ##
  ##
  ## (5) Print result for meta-analysis
  ##
  ##
  if (header) {
    if (inherits(x, "metamiss"))
      cat("Sensitivity analysis for missing binary data\n\n")
    ##
    crtitle(x)
  }
  ##
  if (is.na(k.all)) {
    ## Do nothing
    return(invisible(NULL))
  }
  else if (k.all == 1) {
    ##
    ## Print results for a single study
    ##
    if (!is.metabind) {
      if (sel.n)
        catobsev(x$n, type = "n")
      else
        catobsev(x$n.e, x$n.c, type = "n")
      ##
      if (sel.ev)
        catobsev(x$event, type = "e", addrow = TRUE)
      else
        catobsev(x$event.e, x$event.c, type = "e", addrow = TRUE)
    }
    ##
    res <- cbind(formatN(TE.fixed, digits, "NA",
                         big.mark = big.mark),
                 formatCI(formatN(lowTE.fixed, digits, "NA",
                                  big.mark = big.mark),
                          formatN(uppTE.fixed, digits, "NA",
                                  big.mark = big.mark)),
                 if (null.given)
                   formatN(sTE.fixed, digits.stat, big.mark = big.mark),
                 if (null.given)
                   formatPT(pTE.fixed, digits = digits.pval,
                            scientific = scientific.pval,
                            zero = zero.pval, JAMA = JAMA.pval))
    dimnames(res) <- list(x$studlab,
                          c(sm.lab, ci.lab,
                            if (null.given) "z",
                            if (null.given) "p-value"))
    prmatrix(res, quote = FALSE, right = TRUE, ...)
    ## Print information on meta-analysis method:
    if (details.methods)
      catmeth(class = class(x),
              method =
                if (!metaprop | (overall & (fixed | random)) |
                    overall.hetstat | by)
                  x$method else "NoMA",
              sm = sm,
              k.all = k.all,
              sparse = ifelse(bip, x$sparse, FALSE),
              incr = if (bip) x$incr else FALSE,
              allincr = ifelse(bip, x$allincr, FALSE),
              addincr = ifelse(bip, x$addincr, FALSE),
              allstudies = x$allstudies,
              doublezeros = x$doublezeros,
              MH.exact = ifelse(metabin, x$MH.exact, FALSE),
              method.ci = method.ci,
              pooledvar = x$pooledvar,
              method.smd = x$method.smd,
              sd.glass = x$sd.glass,
              exact.smd = x$exact.smd,
              model.glmm = x$model.glmm,
              pscale = pscale,
              irscale = irscale,
              irunit = irunit,
              null.effect = if (null.given) null.effect else 0,
              big.mark = big.mark,
              digits = digits, digits.tau = digits.tau,
              text.tau = text.tau, text.tau2 = text.tau2,
              method.miss = x$method.miss,
              IMOR.e = x$IMOR.e, IMOR.c = x$IMOR.c,
              three.level = three.level)
  }
  else if (is.na(k)) {
    ## Do nothing
    return(invisible(NULL))
  }
  else {
    ##
    ##
    ## Print results for meta-analysis with more than one study
    ##
    ##
    if (overall & (fixed | random | prediction)) {
      if (!inherits(x, "trimfill")) {
        if (x$method == "MH" &&
            (inherits(x, c("metabin", "metainc")) &
             fixed & sm %in% c("RD", "IRD") &
             (!is.null(x$k.MH) == 1 && k != x$k.MH)))
          cat(paste0("Number of studies combined:   k.MH = ", x$k.MH,
                     " (", text.fixed.br, "), k = ",
                     format(k, big.mark = big.mark),
                     " (", text.random.br, ")\n"))
        else {
          if (k.study != k) {
            cat(paste0("Number of studies combined: n = ",
                       format(x$k.study, big.mark = big.mark), "\n"))
            cat(paste0("Number of estimates combined: k = ",
                       format(k, big.mark = big.mark), "\n"))
          }
          else
            cat(paste0("Number of studies combined: k = ",
                       format(k, big.mark = big.mark), "\n"))
        }
      }
      else
        cat(paste0("Number of studies combined: k = ",
                   format(k, big.mark = big.mark),
                   " (with ",
                   format(x$k0, big.mark = big.mark),
                   " added studies)\n"))
      ##
      if (!is.metabind) {
        if (sel.n)
          catobsev(x$n, type = "n")
        else
          catobsev(x$n.e, x$n.c, type = "n")
        ##
        if (sel.ev)
          catobsev(x$event, type = "e", addrow = TRUE)
        else
          catobsev(x$event.e, x$event.c, type = "e", addrow = TRUE)
      }
      ##
      res <- cbind(formatN(c(if (fixed) TE.fixed,
                             if (random) TE.random,
                             if (prediction) NA),
                           digits, "NA",
                           big.mark = big.mark),
                   formatCI(formatN(c(if (fixed) lowTE.fixed,
                                      if (random) lowTE.random,
                                      if (prediction) lowTE.predict),
                                    digits, "NA", big.mark = big.mark),
                            formatN(c(if (fixed) uppTE.fixed,
                                      if (random) uppTE.random,
                                      if (prediction) uppTE.predict),
                                    digits, "NA", big.mark = big.mark)),
                   if (null.given)
                     formatN(c(if (fixed) sTE.fixed,
                               if (random) sTE.random,
                               if (prediction) NA),
                             digits = digits.stat, big.mark = big.mark),
                   if (null.given)
                     formatPT(c(if (fixed) pTE.fixed,
                                if (random) pTE.random,
                                if (prediction) NA),
                              digits = digits.pval,
                              scientific = scientific.pval,
                              zero = zero.pval, JAMA = JAMA.pval))
      if (prediction) {
        res[dim(res)[1], 1] <- ""
        if (null.given)
          res[dim(res)[1], 3:4] <- ""
      }
      if (!is.null(x$hakn) && x$hakn) {
        if (fixed & random)
          zlab <- "z|t"
        else if (fixed & !random)
          zlab <- "z"
        else if (!fixed & random)
          zlab <- "t"
      }
      else
        zlab <- "z"
      ##
      if (prediction)
        if (x$level.ma == x$level.predict)
          lab.predict <- text.predict
        else
          lab.predict <- paste(text.predict,
                               paste0("(",
                                      round(100 * x$level.predict, 1),
                                      "%-PI)"))
      ##
      dimnames(res) <- list(c(if (fixed) text.fixed,
                              if (random) text.random,
                              if (prediction) lab.predict),
                            c(sm.lab, ci.lab,
                              if (null.given) zlab,
                              if (null.given) "p-value"))
      prmatrix(res, quote = FALSE, right = TRUE, ...)
      ##
      if (metabin && print.CMH) {
        Qdata <- cbind(formatN(round(Q.CMH, digits.Q), digits.Q, "NA",
                               big.mark = big.mark),
                       df.Q.CMH,
                       formatPT(pval.Q.CMH,
                                digits = digits.pval.Q,
                                scientific = scientific.pval,
                                zero = zero.pval, JAMA = JAMA.pval))
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        ##
        cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
    }
    else {
      if (k.study != k) {
        cat(paste0("Number of studies: ",
                   if (is.netpairwise) "k" else "n", " = ",
                   format(x$k.study, big.mark = big.mark), "\n"))
        cat(paste0("Number of ",
                   if (is.netpairwise)
                     "pairwise comparisons: m = " else "estimates: k = ",
                   format(k, big.mark = big.mark), "\n"))
      }
      else
        cat(paste0("Number of studies: k = ",
                   format(k, big.mark = big.mark), "\n"))
      ##
      if (!is.metabind) {
        if (sel.n)
          catobsev(x$n, type = "n")
        else
          catobsev(x$n.e, x$n.c, type = "n")
        ##
        if (sel.ev)
          catobsev(x$event, type = "e")
        else
          catobsev(x$event.e, x$event.c, type = "e")
      }
    }
    ##
    ## Print information on heterogeneity
    ##
    if (overall.hetstat) {
      cat("\nQuantifying heterogeneity:\n")
      ##
      print.tau2.ci <-
        print.tau2 & all(!(is.na(x$lower.tau2) | is.na(x$upper.tau2)))
      if (print.tau2.ci &&
          (all(x$lower.tau2 == 0) & all(x$upper.tau2 == 0)))
        print.tau2.ci <- FALSE
      ##
      print.tau.ci <-
        print.tau & all(!(is.na(x$lower.tau) | is.na(x$upper.tau)))
      if (print.tau.ci &&
          (all(x$lower.tau == 0) & all(x$upper.tau == 0)))
        print.tau.ci <- FALSE
      ##
      ##
      cathet(k,
             x$tau2, x$lower.tau2, x$upper.tau2,
             print.tau2, print.tau2.ci, text.tau2, digits.tau2,
             x$tau, x$lower.tau, x$upper.tau,
             print.tau, print.tau.ci, text.tau, digits.tau,
             x$sign.lower.tau, x$sign.upper.tau,
             I2, lowI2, uppI2,
             print.I2, print.I2.ci, text.I2, digits.I2,
             H, lowH, uppH,
             print.H, digits.H,
             Rb, lowRb, uppRb,
             print.Rb, text.Rb,
             big.mark,
             if (is.null(x$detail.tau)) "" else x$detail.tau)
      ##
      ## Print information on residual heterogeneity
      ##
      if (by & !is.metabind) {
        ##
        Q.resid <- x$Q.w.fixed
        k.resid <- x$df.Q.w + 1
        ##
        if (print.H) {
          H.resid <- round(x$H.resid, digits.H)
          lowH.resid <- round(x$lower.H.resid, digits.H)
          uppH.resid <- round(x$upper.H.resid, digits.H)
        }
        if (print.I2) {
          I2.resid <- round(100 * x$I2.resid, digits.I2)
          lowI2.resid <- round(100 * x$lower.I2.resid, digits.I2)
          uppI2.resid <- round(100 * x$upper.I2.resid, digits.I2)
          print.I2.ci <-
            ((Q.resid  > k.resid & k.resid >= 2) |
             (Q.resid <= k.resid & k.resid > 2)) &
            !(is.na(lowI2.resid) | is.na(uppI2.resid))
          ##
          if (is.na(print.I2.ci))
            print.I2.ci <- FALSE
        }
        ##
        if (!is.na(I2.resid)) {
          cat("\nQuantifying residual heterogeneity:\n")
          ##
          cathet(k.resid, 
                 x$tau2.resid, x$lower.tau2.resid, x$upper.tau2.resid,
                 x$tau.common, print.tau2.ci, text.tau2, digits.tau2,
                 x$tau.resid, x$lower.tau.resid, x$upper.tau.resid,
                 x$tau.common, print.tau.ci, text.tau, digits.tau,
                 "", "",
                 I2.resid, lowI2.resid, uppI2.resid,
                 print.I2, print.I2.ci, text.I2, digits.I2,
                 H.resid, lowH.resid, uppH.resid,
                 print.H, digits.H,
                 NA, NA, NA,
                 FALSE, text.Rb,
                 big.mark,
                 if (is.null(x$detail.tau)) "" else x$detail.tau)
        }
      }
      ##
      ## Test of heterogeneity
      ##
      if (fixed | random) {
        if (k > 1) {
          if (!is.glmm) {
            Qdata <- cbind(formatN(round(Q, digits.Q), digits.Q, "NA",
                                   big.mark = big.mark),
                           format(df.Q, big.mark = big.mark),
                           formatPT(pval.Q,
                                    digits = digits.pval.Q,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval))
            dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
          }
          else {
            Qdata <- cbind(formatN(round(c(Q, Q.LRT), digits.Q), digits.Q, "NA",
                                   big.mark = big.mark),
                           format(c(df.Q, df.Q.LRT), big.mark = big.mark),
                           formatPT(c(pval.Q, pval.Q.LRT),
                                    digits = digits.pval.Q,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval),
                           c("Wald-type", "Likelihood-Ratio"))
            dimnames(Qdata) <- list(rep("", 2),
                                    c("Q", "d.f.", "p-value", "Test"))
          }
          ##
          cat("\nTest of heterogeneity:\n")
          prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
        }
      }
    }
    else {
      print.tau2 <- FALSE
      print.tau2.ci <- FALSE
      print.tau <- FALSE
      print.tau.ci <- FALSE
    }
    ##
    ## Print information for subgroup analysis
    ##
    if (by) {
      if (is.metabind)
        anaunit <- "meta-analyses"
      else if (is.netpairwise)
        anaunit <- "pairwise comparisons"
      else
        anaunit <- "subgroups"
      ##
      if (fixed) {
        ##
        ## Subgroup analysis based on fixed effect model
        ##
        Tdata <- cbind(format(k.w, big.mark = big.mark),
                       formatN(TE.fixed.w, digits, "NA",
                               big.mark = big.mark),
                       formatCI(formatN(lowTE.fixed.w, digits, "NA",
                                        big.mark = big.mark),
                                formatN(uppTE.fixed.w, digits, "NA",
                                        big.mark = big.mark)),
                       formatN(round(Q.w, digits.Q), digits.Q,
                               big.mark = big.mark),
                       if (print.I2)
                         ifelse(is.na(I2.w),
                                "--",
                                paste0(formatN(I2.w, digits.I2), "%")),
                       if (print.Rb)
                         ifelse(is.na(Rb.w),
                                "--",
                                paste0(formatN(Rb.w, digits.I2), "%")),
                       if (!random)
                         ifelse(k.w == 1 & !x$tau.common, "--",
                                formatPT(x$tau.w^2,
                                         digits = digits.tau2,
                                         big.mark = big.mark,
                                         noblanks = TRUE)),
                       if (!random)
                         ifelse(k.w == 1 & !x$tau.common, "--",
                                formatPT(x$tau.w,
                                         digits = digits.tau,
                                         big.mark = big.mark,
                                         noblanks = TRUE))
                       )
        ##
        bylab.txt <- bylabel(subgroup.name, bylevs,
                             print.subgroup.name, sep.subgroup,
                             big.mark = big.mark)
        ##
        dimnames(Tdata) <- list(bylab.txt,
                                c(if (is.netpairwise)
                                    "  m" else "  k",
                                  sm.lab, ci.lab,
                                  "Q",
                                  if (print.I2) text.I2,
                                  if (print.Rb) text.Rb,
                                  if (!random) text.tau2,
                                  if (!random) text.tau)
                                )
        ##
        cat(paste0("\nResults for ", anaunit, " (",
                   text.fixed.br, "):\n"))
        ##
        prmatrix(Tdata, quote = FALSE, right = TRUE, ...)
        ##
        if (is.glmm & length(df.Q.b) > 1) {
          dfs.b <-
            rmSpace(paste(formatN(df.Q.b, digits = 0, big.mark = big.mark),
                          collapse = ", "), end = TRUE)
          Q.lab <- "F"
        }
        else {
          dfs.b <- formatN(df.Q.b, digits = 0, big.mark = big.mark)
          Q.lab <- "Q"
        }
        ##
        if (test.subgroup.fixed & !is.metabind) {
          cat(paste0("\nTest for subgroup differences (",
                     text.fixed.br, "):\n"))
          if (x$method == "MH") {
            Qdata <- cbind(formatN(round(Q.b.fixed, digits.Q), digits.Q, "NA",
                                   big.mark = big.mark),
                           formatN(dfs.b, digits = 0, big.mark = big.mark),
                           formatPT(pval.Q.b.fixed,
                                    digits = digits.pval.Q,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval))
            dimnames(Qdata) <- list("Between groups  ",
                                    c(Q.lab, "d.f.", "p-value"))
            prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
          }
          else {
            Qs  <- c(Q.b.fixed, Q.w.fixed)
            dfs <- c(dfs.b, formatN(df.Q.w, digits = 0, big.mark = big.mark))
            Q.lab <-
              ifelse(is.glmm && Q.lab == "F", "F/Q", Q.lab)
            ##
            pvals <- c(pval.Q.b.fixed, pval.Q.w.fixed)
            Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                   big.mark = big.mark),
                           dfs,
                           formatPT(pvals,
                                    digits = digits.pval.Q,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval))
            dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                    c(Q.lab, "d.f.", "p-value"))
            prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
          }
        }
      }
      ##
      if (random) {
        ##
        ## Subgroup analysis based on random effects model
        ##
        Tdata <- cbind(format(k.w, big.mark = big.mark),
                       formatN(TE.random.w, digits, "NA",
                               big.mark = big.mark),
                       formatCI(formatN(lowTE.random.w, digits, "NA",
                                        big.mark = big.mark),
                                formatN(uppTE.random.w, digits, "NA",
                                        big.mark = big.mark)),
                       ifelse(k.w == 1 & !x$tau.common, "--",
                              formatPT(x$tau.w^2,
                                       digits = digits.tau2,
                                       big.mark = big.mark,
                                       noblanks = TRUE)),
                       ifelse(k.w == 1 & !x$tau.common, "--",
                              formatPT(x$tau.w,
                                       digits = digits.tau,
                                       big.mark = big.mark,
                                       noblanks = TRUE)),
                       if (!fixed)
                         formatN(round(Q.w, digits.Q), digits.Q,
                                 big.mark = big.mark),
                       if (!fixed & print.I2)
                         ifelse(is.na(I2.w),
                                "--",
                                paste0(formatN(I2.w, digits.I2), "%")),
                       if (!fixed & print.Rb)
                         ifelse(is.na(Rb.w),
                                "--",
                                paste0(formatN(Rb.w, digits.I2,
                                               big.mark = big.mark), "%"))
                       )
        ##
        bylab.txt <- bylabel(subgroup.name, bylevs,
                             print.subgroup.name, sep.subgroup,
                             big.mark = big.mark)
        ##
        dimnames(Tdata) <- list(bylab.txt,
                                c(if (is.netpairwise)
                                    "  m" else "  k",
                                  sm.lab, ci.lab,
                                  text.tau2, text.tau,
                                  if (!fixed) "Q",
                                  if (!fixed & print.I2) text.I2,
                                  if (!fixed & print.Rb) text.Rb)
                                )
        ##
        cat(paste0("\nResults for ", anaunit, " (", text.random.br, "):\n"))
        ##
        prmatrix(Tdata, quote = FALSE, right = TRUE, ...)
        ##
        if ((three.level | is.glmm) & length(df.Q.b) > 1) {
          dfs.b <-
            rmSpace(paste(formatN(df.Q.b, digits = 0, big.mark = big.mark),
                          collapse = ", "), end = TRUE)
          Q.lab <- "F"
        }
        else {
          dfs.b <- formatN(df.Q.b, digits = 0, big.mark = big.mark)
          Q.lab <- "Q"
        }
        ##
        if (test.subgroup.random & !is.metabind & !is.na(Q.b.random)) {
          cat(paste0("\nTest for subgroup differences (",
                     text.random.br, "):\n"))
          if (is.na(Q.w.random)) {
            Qdata <- cbind(formatN(round(Q.b.random, digits.Q), digits.Q,
                                   "NA", big.mark = big.mark),
                           formatN(dfs.b, digits = 0, big.mark = big.mark),
                           formatPT(pval.Q.b.random,
                                    digits = digits.pval.Q,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval))
            dimnames(Qdata) <- list("Between groups  ",
                                    c(Q.lab, "d.f.", "p-value"))
          }
          else {
            Qs  <- c(Q.b.random, Q.w.random)
            dfs <- c(dfs.b, formatN(df.Q.w, digits = 0, big.mark = big.mark))
            Q.lab <-
              ifelse((three.level | is.glmm) && Q.lab == "F", "F/Q", Q.lab)
            ##
            pvals <- c(pval.Q.b.random, pval.Q.w.random)
            Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                   big.mark = big.mark),
                           dfs,
                           formatPT(pvals,
                                    digits = digits.pval.Q,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval))
            dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                    c(Q.lab, "d.f.", "p-value"))
          }
          prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
        }
      }
      ##
      if (prediction.w & !is.metabind) {
        ##
        ## Prediction intervals for subgroups
        ##
        Pdata <- cbind(formatCI(formatN(lowTE.predict.w, digits, "NA",
                                        big.mark = big.mark),
                                formatN(uppTE.predict.w, digits, "NA",
                                        big.mark = big.mark)))
        ##
        bylab.txt <-
          bylabel(subgroup.name, bylevs,
                  print.subgroup.name, sep.subgroup,
                  big.mark = big.mark)
        lab.predict <- paste0(round(100 * x$level.predict, 1), "%-PI")
        dimnames(Pdata) <- list(bylab.txt, lab.predict)
        ##
        cat("\nPrediction intervals for subgroups:\n")
        prmatrix(Pdata, quote = FALSE, right = TRUE, ...)
      }      
    }
    ##
    ## Print information on meta-analysis method:
    ##
    if (details.methods & (fixed | random | prediction))
      catmeth(class = class(x),
              method =
                if ((overall & (fixed | random)) |
                    overall.hetstat | by)
                  x$method else "NoMA",
              method.tau =
                if ((overall & random) | overall.hetstat | by)
                  x$method.tau else NULL,
              sm =
                if ((overall & (fixed | random)) |
                    overall.hetstat | by)
                  sm else "",
              k.all = k.all,
              hakn = !is.null(x$hakn) && (x$hakn & random),
              adhoc.hakn = !is.null(x$adhoc.hakn) &&
                (!is.null(x$seTE.random.hakn.orig) & x$adhoc.hakn != ""),
              tau.common = by & x$tau.common,
              tau.preset = x$tau.preset,
              sparse = ifelse(bip, x$sparse, FALSE),
              incr = if (bip) x$incr else FALSE,
              allincr = ifelse(bip, x$allincr, FALSE),
              addincr = ifelse(bip, x$addincr, FALSE),
              allstudies = x$allstudies,
              doublezeros = x$doublezeros,
              MH.exact = ifelse(metabin, x$MH.exact, FALSE),
              RR.Cochrane = ifelse(metabin, x$RR.Cochrane, FALSE),
              Q.Cochrane = ifelse(metabin, x$Q.Cochrane, TRUE),
              method.ci = method.ci,
              print.tau.ci = print.tau2.ci | print.tau.ci,
              method.tau.ci = x$method.tau.ci,
              pooledvar = x$pooledvar,
              method.smd = x$method.smd,
              sd.glass = x$sd.glass,
              exact.smd = x$exact.smd,
              model.glmm = x$model.glmm,
              pscale = pscale,
              irscale = irscale,
              irunit = irunit,
              null.effect = if (null.given) null.effect else 0,
              big.mark = big.mark,
              digits = digits, digits.tau = digits.tau,
              text.tau = text.tau, text.tau2 = text.tau2,
              method.miss = x$method.miss,
              IMOR.e = x$IMOR.e, IMOR.c = x$IMOR.c,
              three.level =
                if (is.null(x$three.level)) FALSE else x$three.level)
  }
  
  
  invisible(NULL)
}





#' @rdname print.meta
#' @export cilayout


cilayout <- function(bracket = "[", separator = "; ") {
  
  ibracket <- charmatch(bracket,
                        c("[", "(", "{", ""),
                        nomatch = NA)
  ##
  if (is.na(ibracket) | ibracket == 0)
    stop("No valid bracket type specified. ",
         "Admissible values: '[', '(', '{', '\"\"'")
  
  bracket <- c("[", "(", "{", "")[ibracket]
  
  setOption("CIbracket", bracket)
  setOption("CIseparator", separator)
  
  invisible(NULL)
}
