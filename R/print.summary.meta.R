#' Print detailed meta-analysis results
#' 
#' @description
#' Print method for objects of class \code{summary.meta}.
#' 
#' @aliases print.summary.meta
#' 
#' @param x An object of class \code{summary.meta}
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param common A logical indicating whether results of common effect
#'   meta-analysis should be printed.
#' @param random A logical indicating whether results of random
#'   effects meta-analysis should be printed.
#' @param details A logical indicating whether further details of
#'   individual studies should be printed.
#' @param ma A logical indicating whether the summary results of the
#'   meta-analysis should be printed.
#' @param overall A logical indicating whether overall summaries
#'   should be reported. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be reported.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf = TRUE}, results
#'   for \code{sm = "OR"} are printed as odds ratios rather than log
#'   odds ratios and results for \code{sm = "ZCOR"} are printed as
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
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   deviations and standard errors, see \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value of test for effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of test of treatment effect, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.df Minimal number of significant digits for degrees
#'   of freedom.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
#' @param digits.H Minimal number of significant digits for H
#'   statistic, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param digits.weight Minimal number of significant digits for
#'   weights, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param big.mark A character used as thousands separator.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau2.ci A logical value indicating whether to print
#'   the confidence interval of \eqn{\tau^2}.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param print.tau.ci A logical value indicating whether to print the
#'   confidence interval of \eqn{\tau}.
#' @param print.Q A logical value indicating whether to print the
#'   results of the test of heterogeneity.
#' @param print.I2 A logical specifying whether heterogeneity
#'   statistic I\eqn{^2} should be printed.
#' @param print.I2.ci A logical specifying whether confidence interval for
#'   heterogeneity statistic I\eqn{^2} should be printed.
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
#' @param truncate An optional vector used to truncate the printout of
#'   results for individual studies (must be a logical vector of same
#'   length as \code{x$TE} or contain numerical values).
#' @param text.truncate A character string printed if study results
#'   were truncated from the printout.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param warn.backtransf Deprecated argument (ignored).
#' @param \dots Additional arguments (passed on to
#'   \code{\link{print.meta}} called internally).
#' 
#' @details
#' Print method for objects of class \code{summary.meta} giving
#' detailed information on the meta-analysis.
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
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{summary.meta}}, \code{\link{update.meta}},
#'   \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}
#' 
#' @references Cooper H & Hedges LV (1994), \emph{The Handbook of
#'   Research Synthesis}.  Newbury Park, CA: Russell Sage Foundation.
#' 
#' Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D (2016), A new measure
#' of between-studies heterogeneity in meta-analysis.  \emph{Statistics in
#' Medicine}, \bold{35}, 3661--75.
#' 
#' Higgins JPT & Thompson SG (2002), Quantifying heterogeneity in a
#' meta-analysis.  \emph{Statistics in Medicine}, \bold{21}, 1539--58.
#' 
#' @keywords print
#' 
#' @examples
#' data(Fleiss1993cont)
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year))
#' sm1 <- summary(m1)
#' sm1
#' 
#' print(sm1, digits = 2)
#' 
#' \dontrun{
#' # Use unicode characters to print tau^2, tau, and I^2 
#' print(sm1,
#'   text.tau2 = "\u03c4\u00b2",
#'   text.tau = "\u03c4", text.I2 = "I\u00b2")
#' }
#' 
#' @method print summary.meta
#' @export


print.summary.meta <- function(x,
                               sortvar,
                               common = x$x$common,
                               random = x$x$random,
                               details = FALSE, ma = TRUE,
                               overall = x$overall & ma,
                               ##
                               backtransf = x$backtransf,
                               pscale = x$pscale,
                               irscale = x$irscale,
                               irunit = x$irunit,
                               #
                               digits = gs("digits"),
                               digits.se = gs("digits.se"),
                               digits.stat = gs("digits.stat"),
                               digits.pval = max(gs("digits.pval"), 2),
                               #
                               digits.tau2 = gs("digits.tau2"),
                               digits.tau = gs("digits.tau"),
                               #
                               digits.Q = gs("digits.Q"),
                               digits.df = gs("digits.df"),
                               digits.pval.Q = max(gs("digits.pval.Q"), 2),
                               #
                               digits.H = gs("digits.H"),
                               digits.I2 = gs("digits.I2"),
                               #
                               digits.prop = gs("digits.prop"),
                               digits.weight = gs("digits.weight"),
                               #
                               scientific.pval = gs("scientific.pval"),
                               zero.pval = gs("zero.pval"),
                               JAMA.pval = gs("JAMA.pval"),
                               #
                               big.mark = gs("big.mark"),
                               #
                               print.tau2 = gs("print.tau2"),
                               print.tau2.ci = gs("print.tau2.ci"),
                               print.tau = gs("print.tau"),
                               print.tau.ci = gs("print.tau.ci"),
                               #
                               print.Q = gs("print.Q"),
                               print.I2 = gs("print.I2"),
                               print.I2.ci = gs("print.I2.ci"),
                               print.H = gs("print.H"),
                               print.Rb = gs("print.Rb"),
                               #
                               text.tau2 = gs("text.tau2"),
                               text.tau = gs("text.tau"),
                               text.I2 = gs("text.I2"),
                               text.Rb = gs("text.Rb"),
                               #
                               truncate,
                               text.truncate = "*** Output truncated ***",
                               ##
                               details.methods = TRUE,
                               ##
                               warn.backtransf = FALSE,
                               ...
                               ) {
  
  
  ##
  ##
  ## (1) Check for summary.meta object
  ##
  ##
  chkclass(x, "summary.meta")
  ##
  k.all <- length(x$TE)
  ##
  x.meta <- updateversion(x$x)
  ##
  fbt <- x$x$func.backtransf
  abt <- x$x$args.backtransf
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  error <-
    try(sortvar <-
          catch("sortvar", mc, x.meta, sfsp),
        silent = TRUE)
  if (inherits(error, "try-error")) {
    sortvar <- catch("sortvar", mc, x$data,  NULL)
    if (isCol(x$data, ".subset"))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != k.all))
    stop("Number of studies in object 'x' and ",
         "argument 'sortvar' have different length.")
  if (!sort)
    sortvar <- 1:k.all
  ##
  chklogical(details)
  chklogical(ma)
  overall <- replaceNULL(overall, TRUE)
  chklogical(overall)
  ##
  if (is_untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  ##
  chklogical(details.methods)
  ##
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  if (!backtransf & pscale != 1 & !is_untransformed(x$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  if (!backtransf & irscale != 1 & !is_untransformed(x$sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  scale <- 1
  if (pscale != 1 || irscale != 1) {
    if (pscale != 1 && irscale != 1)
      stop("Provide either arguments 'pscale' or 'irscale' but not ",
           "both arguments.",
           call. = FALSE)
    if (pscale != 1)
      scale <- pscale
    else
      scale <- irscale
  }
  ##
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  chknumeric(digits.weight, min = 0, length = 1)
  ##
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  ##
  chklogical(print.tau2)
  chklogical(print.tau2.ci)
  chklogical(print.tau)
  chklogical(print.tau.ci)
  chklogical(print.I2)
  chklogical(print.I2.ci)
  chklogical(print.H)
  chklogical(print.Rb)
  chklogical(print.Q)
  ##
  chkchar(text.tau2, length = 1)
  chkchar(text.tau, length = 1)
  chkchar(text.I2, length = 1)
  chkchar(text.Rb, length = 1)
  ##
  ## Catch 'truncate' from meta-analysis object:
  ##
  missing.truncate <- missing(truncate)
  if (!missing.truncate) {
    truncate <- catch("truncate", mc, x.meta, sfsp)
    ##
    if (is.null(truncate))
      truncate <- catch("truncate", mc, x$data, sfsp)
    ##
    if (length(truncate) > k.all)
      stop("Length of argument 'truncate' is too long.",
           call. = FALSE)
    else if (length(truncate) < k.all) {
      if (is.numeric(truncate)) {
        if (any(is.na(truncate)) | max(truncate) > k.all | min(truncate) < 0)
          stop("Numeric values in argument 'truncate' must be between 1 and ",
               k.all, ".",
               call. = FALSE)
        truncate2 <- rep(FALSE, k.all)
        truncate2[truncate] <- TRUE
        truncate <- truncate2
      }
      else if (is.character(truncate)) {
        if (any(!(truncate %in% x$studlab)))
          stop("At least one value of argument 'truncate' does not ",
               "match a study label.",
               call. = FALSE)
        truncate2 <- rep(FALSE, k.all)
        truncate2[x$studlab %in% truncate] <- TRUE
        truncate <- truncate2
      }
      else
        stop("Argument 'truncate' must contain integers or study labels if ",
             "length differs from number of treatment effects.",
             call. = FALSE)
    }
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ##
  missing.common <- missing(common)
  common <- replaceNULL(common, x$comb.common)
  common <- deprecated(common, missing.common, args, "comb.fixed", FALSE)
  common <- deprecated(common, missing.common, args, "fixed", FALSE)
  chklogical(common)
  ##
  random <- replaceNULL(random, x$comb.random)
  random <- deprecated(random, missing(random), args, "comb.random", FALSE)
  chklogical(random)
  ##
  ## More checks ...
  ##
  cl <- paste0("update.meta() or ", class(x)[1], "()")
  addargs <- names(list(...))
  ##
  level <- x$level
  level.ma <- replaceNULL(x$level.ma, x$level.comb)
  level.predict <- x$level.predict
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  ci.lab <- paste0(round(100 * level, 1), "%-CI")
  ##
  sm <- x$sm
  #
  sm.lab <- smlab(sm, backtransf, pscale, irscale)
  #
  if (is.null(x$text.w.common))
    text.w.common <- paste0("%W(", gs("text.w.common"), ")")
  else
    text.w.common <- paste0("%W(", x$text.w.common, ")")
  ##
  if (is.null(x$text.w.random))
    text.w.random <- paste0("%W(", gs("text.w.random"), ")")
  else
    text.w.random <- paste0("%W(", x$text.w.random, ")")
  ##
  subgroup <- replaceNULL(x$subgroup, x$byvar)
  subgroup.name <- replaceNULL(x$subgroup.name, x$bylab)
  ##
  by <- !is.null(subgroup)
  three.level <- !is.null(x$three.level) && any(x$three.level)
  n_of_1 <- !is.null(x$cycles)
  
  
  ##
  ##
  ## (4) Print title and details
  ##
  ##
  is.metamiss <- inherits(x, "metamiss")
  show.imor <- is.metamiss &
    !is.null(x$IMOR.e) & !is.null(x$IMOR.c) &&
    (length(unique(x$IMOR.e)) != 1 | length(unique(x$IMOR.c)) != 1)
  ##
  if (is.metamiss)
    cat("Sensitivity analysis for missing binary data\n\n")
  ##
  crtitle(x)
  ##
  if (details) {
    if (is.metamiss) {
      res <- cbind(event.e = formatN(x$event.e, digits = 0,
                                     "NA", big.mark = big.mark),
                   noevent.e = formatN(x$n.e - x$event.e - x$miss.e,
                                       digits = 0,
                                       "NA", big.mark = big.mark),
                   miss.e = formatN(x$miss.e, digits = 0,
                                    "NA", big.mark = big.mark),
                   event.c = formatN(x$event.c, digits = 0,
                                     "NA", big.mark = big.mark),
                   noevent.c = formatN(x$n.c - x$event.c - x$miss.c,
                                       digits = 0,
                                       "NA", big.mark = big.mark),
                   miss.c = formatN(x$miss.c, digits = 0,
                                    "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metabin")) {
      res <- cbind(event.e = formatN(x$event.e, digits = 0,
                                     "NA", big.mark = big.mark),
                   n.e = formatN(x$n.e, digits = 0,
                                 "NA", big.mark = big.mark),
                   event.c = formatN(x$event.c, digits = 0,
                                     "NA", big.mark = big.mark),
                   n.c = formatN(x$n.c, digits = 0,
                                 "NA", big.mark = big.mark))
      ##
      if (pscale == 1) {
        res <- cbind(res,
                     p.e = formatN(round(x$event.e / x$n.e, digits.prop),
                                   digits.prop, big.mark = big.mark))
        res <- cbind(res,
                     p.c = formatN(round(x$event.c / x$n.c, digits.prop),
                                   digits.prop, big.mark = big.mark))
      }
      else {
        res <- cbind(res,
                     events.e = formatN(round(pscale * x$event.e / x$n.e,
                                              digits),
                                        digits, "NA", big.mark = big.mark))
        res <- cbind(res,
                     events.c = formatN(round(pscale * x$event.c / x$n.c,
                                              digits),
                                        digits, "NA", big.mark = big.mark))
      }
    }
    else if (inherits(x, "metacont")) {
      res <- cbind(n.e = formatN(x$n.e, digits = 0,
                                 "NA", big.mark = big.mark),
                   mean.e = formatN(round(x$mean.e, digits), digits,
                                    "NA", big.mark = big.mark),
                   sd.e = formatN(round(x$sd.e, digits.se), digits.se,
                                  "NA", big.mark = big.mark),
                   n.c = formatN(x$n.c, digits = 0,
                                 "NA", big.mark = big.mark),
                   mean.c = formatN(round(x$mean.c, digits), digits,
                                    "NA", big.mark = big.mark),
                   sd.c = formatN(round(x$sd.c, digits.se), digits.se,
                                  "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metacor")) {
      res <- cbind(cor = x$cor,
                   n = formatN(x$n, digits = 0,
                               "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metagen")) {
      res <- cbind(TE = formatN(round(x$TE, digits), digits,
                                "NA", big.mark = big.mark),
                   seTE = formatN(round(x$seTE, digits.se), digits.se,
                                  "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metainc")) {
      res <- cbind(event.e = formatN(x$event.e, digits = 0,
                                     "NA", big.mark = big.mark),
                   time.e = formatN(round(x$time.e, digits), digits,
                                    "NA", big.mark = big.mark),
                   event.c = formatN(x$event.c, digits = 0,
                                     "NA", big.mark = big.mark),
                   time.c = formatN(round(x$time.c, digits), digits,
                                    "NA", big.mark = big.mark))
      ##
      if (irscale == 1) {
        res <- cbind(res,
                     rate.e = formatN(round(x$event.e / x$time.e,
                                            digits.prop),
                                      big.mark = big.mark))
        res <- cbind(res,
                     rate.c = formatN(round(x$event.c / x$time.c,
                                            digits.prop),
                                      big.mark = big.mark))
      }
      else {
        res <- cbind(res,
                     events.e = formatN(round(irscale * x$event.e / x$n.e,
                                              digits),
                                        digits, "NA", big.mark = big.mark))
        res <- cbind(res,
                     events.c = formatN(round(irscale * x$event.c / x$n.c,
                                              digits),
                                        digits, "NA", big.mark = big.mark))
      }
    }
    else if (inherits(x, "metaprop")) {
      res <- cbind(event = formatN(x$event, digits = 0,
                                   "NA", big.mark = big.mark),
                   n = formatN(x$n, digits = 0,
                               "NA", big.mark = big.mark))
      if (pscale == 1)
        res <- cbind(res,
                     p = formatN(round(x$event / x$n, digits.prop),
                                 digits.prop, "NA", big.mark = big.mark))
      else
        res <- cbind(res,
                     events = formatN(round(pscale * x$event / x$n, digits),
                                      digits, "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metarate")) {
      res <- cbind(event = formatN(x$event, digits = 0,
                                   "NA", big.mark = big.mark),
                   time = formatN(x$time, digits = digits,
                                  "NA", big.mark = big.mark))
      if (irscale == 1)
        res <- cbind(res,
                     rate = formatN(round(x$event / x$time, digits.prop),
                                    digits.prop, "NA", big.mark = big.mark))
      else
        res <- cbind(res,
                     events = formatN(round(irscale * x$event / x$time,
                                            digits),
                                      digits, "NA", big.mark = big.mark))
      if (!is.null(x$n))
        res <- cbind(res,
                     n = formatN(x$n, digits = 0,
                                 "NA", big.mark = big.mark))
    }
    else {
      res <- cbind(TE = formatN(round(x$TE, digits), digits,
                                "NA", big.mark = big.mark),
                   seTE = formatN(round(x$seTE, digits), digits,
                                  "NA", big.mark = big.mark))
    }
    ##
    if (three.level)
      res <- cbind(res, cluster = as.character(x$cluster))
    ##
    if (n_of_1)
      res <- cbind(res, cycles = x$cycles)
    ##
    if (by)
      res <- cbind(res, subgroup = as.character(subgroup))
    #
    if (common & !is.null(x$weights.common))
      res <- cbind(res, weights.common = x$weights.common)
    #
    if (random & !is.null(x$weights.random))
      res <- cbind(res, weights.random = x$weights.random)
    #
    if ("weights.common" %in% colnames(res)) {
      if (!random)
        colnames(res)[colnames(res) == "weights.common"] <- "W"
      else
        colnames(res)[colnames(res) == "weights.common"] <- "W(common)"
    }
    #
    if ("weights.random" %in% colnames(res)) {
      if (!common)
        colnames(res)[colnames(res) == "weights.random"] <- "W"
      else
        colnames(res)[colnames(res) == "weights.random"] <- "W(random)"
    }
    #
    dimnames(res)[[1]] <- x$studlab
    ##
    if (!missing.truncate) {
      sortvar <- sortvar[truncate]
      res <- res[truncate, , drop = FALSE]
    }
    ##
    prmatrix(res[order(sortvar), , drop = FALSE],
             quote = FALSE, right = TRUE)
    if (!missing.truncate)
      cat(text.truncate, "\n")
    cat("\n")
  }
  
  
  ##
  ##
  ## (5) Print results for individual studies
  ##
  ##
  if (k.all == 1 &&
      !(inherits(x, c("metaprop", "metarate")) |
        (inherits(x, "metabin") && x$sm == "RR" && !x$RR.Cochrane &&
         !is_zero(x$TE - x$TE.common)) |
        (inherits(x, c("metacont", "metamean")) && x$method.ci == "t"))) {
    print.meta(x.meta,
               header = FALSE,
               digits = digits,
               backtransf = backtransf, pscale = pscale,
               irscale = irscale, irunit = irunit,
               #
               text.tau2 = text.tau2, text.tau = text.tau, text.I2 = text.I2,
               #
               scientific.pval = scientific.pval, big.mark = big.mark,
               zero.pval = zero.pval, JAMA.pval = JAMA.pval,
               #
               details.methods = details.methods,
               ...)
  }
  else {
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    if (k.all == 1 &&
        inherits(x, "metabin") && x$sm == "RR" && !x$RR.Cochrane &&
        !is_zero(x$TE - x$TE.common))
      x$method.ci <- "!RR.Cochrane"
    ##
    if (backtransf) {
      ## Freeman-Tukey Arcsin transformation
      if (metainf.metacum | inherits(x, "metabind")) {
        if (sm == "IRFT")
          harmonic.mean <- x$t.harmonic.mean
        else
          harmonic.mean <- x$n.harmonic.mean
      }
      else {
        if (sm == "IRFT")
          harmonic.mean <- x$time
        else
          harmonic.mean <- x$n
      }
      ##
      if (inherits(x, "metaprop"))
        TE <- x$event / x$n
      ##
      else if (inherits(x, "metarate"))
        TE <- x$event / x$time
      else
        TE <- backtransf(TE, sm, harmonic.mean, harmonic.mean, fbt, abt)
      ##
      lowTE <- backtransf(lowTE, sm, harmonic.mean, harmonic.mean, fbt, abt)
      uppTE <- backtransf(uppTE, sm, harmonic.mean, harmonic.mean, fbt, abt)
      ##
      TE <- scale * TE
      lowTE <- scale * lowTE
      uppTE <- scale * uppTE
      ##
      if (sm == "VE") {
        tmp.l <- lowTE
        lowTE <- uppTE
        uppTE <- tmp.l
      }
    }
    ##
    TE <- round(TE, digits)
    lowTE <- round(lowTE, digits)
    uppTE <- round(uppTE, digits)
    ##
    if (!metainf.metacum) {
      if (common) {
        w.common <- x$w.common
        ##
        if (is.matrix(w.common)) {
          keep <- !apply(w.common, 2, allNA)
          w.common <- w.common[, keep]
          text.w.common <- text.w.common[keep]
        }
        ##
        if (!is.null(w.common) & !all(is.na(w.common)) &&
            sum(w.common) > 0) {
          if (is.matrix(w.common))
            w.common.p <-
              round(apply(w.common, 2, calcPercent), digits.weight)
          else
            w.common.p <-
              round(calcPercent(w.common), digits.weight)
        }
        else
          w.common.p <- w.common
        ##
        if (!allNA(w.common.p) && all(w.common.p == 0, na.rm = TRUE))
          w.common.p[!is.na(w.common.p)] <- NA
      }
      ##
      if (random) {
        w.random <- x$w.random
        ##
        if (is.matrix(w.random)) {
          keep <- !apply(w.random, 2, allNA)
          w.random <- w.random[, keep]
          text.w.random <- text.w.random[keep]
        }
        ##
        if (!is.null(w.random) && !all(is.na(w.random)) &&
            sum(w.random) > 0) {
          if (is.matrix(w.random))
            w.random.p <-
              round(apply(w.random, 2, calcPercent), digits.weight)
          else
            w.random.p <-
              round(calcPercent(w.random), digits.weight)
        }
        else
          w.random.p <- w.random
        ##
        if (!allNA(w.random.p) && all(w.random.p == 0, na.rm = TRUE))
          w.random.p[!is.na(w.random.p)] <- NA
      }
    }
    ##
    if (metainf.metacum) {
      if (is.null(x$text.common))
        text.common <- gs("text.common")
      else
        text.common <- x$text.common
      ##
      if (is.null(x$text.random))
        text.random <- gs("text.random")
      else
        text.random <- x$text.random
      ##
      if (any(substring(text.common, 1, 5) %in% c("Fixed", "Commo"))) {
        text.common <- gsub("Fixed", "fixed", text.common)
        text.common <- gsub("Common", "common", text.common)
        text.common <- gsub("Effect", "effect", text.common)
      }
      ##
      if (any(substring(text.random, 1, 5) %in% c("Rando"))) {
        text.random <- gsub("Random", "random", text.random)
        text.random <- gsub("Effect", "effect", text.random)
      }
      ##
      is.random <- x$pooled == "random"
      ##
      I2 <- formatN(round(100 * x$I2, digits.I2), digits.I2, "")
      ##
      pval <- formatPT(x$pval, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval,
                       lab.NA = "")
      ##
      tau2 <- formatN(round(x$tau2, digits.tau2), digits.tau2, "",
                      big.mark = big.mark)
      tau <- formatN(round(x$tau, digits.tau), digits.tau, "",
                     big.mark = big.mark)
      ##
      res <- cbind(formatN(round(TE, digits), digits, "",
                           big.mark = big.mark),
                   formatCI(formatN(round(lowTE, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE, digits), digits, "NA",
                                    big.mark = big.mark)),
                   pval,
                   if (print.tau2) paste0(" ", tau2),
                   if (print.tau) paste0(" ", tau),
                   if (print.I2) paste0(" ", I2, ifelse(I2 == "", "", "%")))
      dimnames(res) <- list(x$studlab,
                            c(sm.lab, ci.lab, "p-value",
                              if (print.tau2) text.tau2,
                              if (print.tau) text.tau,
                              if (print.I2) text.I2))
      ##
      if (inherits(x, "metainf")) {
        if (!is.random)
          cat(paste0("Influential analysis (", text.common, ")\n"))
        else
          cat(paste0("Influential analysis (", text.random, ")\n"))
      }
      else if (inherits(x, "metacum")) {
        if (!is.random)
          cat(paste0("Cumulative meta-analysis (", text.common, ")\n"))
        else
          cat(paste0("Cumulative meta-analysis (", text.random, ")\n"))
      }
      cat("\n")
      prmatrix(res, quote = FALSE, right = TRUE, na.print = "--")
      #
      # Print information on summary method:
      #
      if (details.methods)
        catmeth(x,
                common, random, x$prediction, overall, x$overall.hetstat,
                #
                func.transf = x$func.transf,
                backtransf = backtransf, func.backtransf = x$func.backtransf,
                #
                big.mark = big.mark, digits = digits,
                digits.tau = digits.tau,
                text.tau = text.tau, text.tau2 = text.tau2, text.I2 = text.I2,
                #
                print.tau2 = print.tau2, print.tau = print.tau,
                print.I2 = print.I2)
    }
    else if (!(inherits(x, "metabind") && !x$show.studies)) {
      show.w.common  <-
        (overall | by) &
        (common && !all(is.na(w.common.p)))
      show.w.random <-
        (overall | by) &
        (random && !all(is.na(w.random.p)))
      ##
      res <- cbind(formatN(round(TE, digits), digits, "NA",
                           big.mark = big.mark),
                   formatCI(formatN(round(lowTE, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE, digits), digits, "NA",
                                    big.mark = big.mark)),
                   if (show.w.common)
                     formatN(w.common.p, digits.weight,
                             big.mark = big.mark),
                   if (show.w.random)
                     formatN(w.random.p, digits.weight,
                             big.mark = big.mark),
                   if (three.level) as.character(x$cluster),
                   if (n_of_1) x$cycles,
                   if (by) as.character(subgroup),
                   if (show.imor) round(x$IMOR.e, 4),
                   if (show.imor) round(x$IMOR.c, 4),
                   if (!is.null(x$exclude))
                     ifelse(is.na(x$exclude), "",
                     ifelse(x$exclude, "*", "")))
      ## Printout for a single proportion, mean difference or mean:
      if (k.all == 1) {
        ##
        print.stat <- FALSE
        print.pval <- FALSE
        ##
        if (!is.null(x$method.ci)) {
          if (x$method.ci == "t") {
            details.ci <-
              "Confidence interval based on t-distribution:\n\n"
            ##
            ## Add test statistic and p-value
            ##
            if (inherits(x, c("metacont", "metamean"))) {
              if (any(!is.na(x$statistic))) {
                res <- cbind(res,
                             formatN(x$statistic,
                                     digits = digits.stat,
                                     big.mark = big.mark))
                ##
                print.stat <- TRUE
              }
              ##
              if (any(!is.na(x$pval))) {
                res <- cbind(res,
                             formatPT(x$pval, digits = digits.pval,
                                      scientific = scientific.pval,
                                      zero = zero.pval, JAMA = JAMA.pval,
                                      lab.NA = ""))
                ##
                print.pval <- TRUE
              }
            }
          }            
          else if (x$method.ci == "CP") {
            details.ci <-
              "Clopper-Pearson confidence interval:\n\n"
            ##
            ## Add p-value of binomial test
            ##
            if (any(!is.na(x$pval))) {
              res <- cbind(res,
                           formatPT(x$pval, digits = digits.pval,
                                    scientific = scientific.pval,
                                    zero = zero.pval, JAMA = JAMA.pval,
                                    lab.NA = ""))
              ##
              print.pval <- TRUE
            }
          }
          else if (x$method.ci == "WS")
            details.ci <-
              "Wilson Score confidence interval:\n\n"
          else if (x$method.ci == "WSCC")
            details.ci <-
              "Wilson Score confidence interval with continuity correction:\n\n"
          else if (x$method.ci == "AC")
            details.ci <-
              "Agresti-Coull confidence interval:\n\n"
          else if (x$method.ci == "SA")
            details.ci <-
              "Simple approximation confidence interval:\n\n"
          else if (x$method.ci == "SACC")
            details.ci <-
              paste0("Simple approximation confidence interval with ",
                     "continuity correction:\n\n")
          else if (x$method.ci == "Poisson")
            details.ci <-
              "Exact Poisson confidence interval for individual studies:\n\n"
          else if (x$method.ci == "!RR.Cochrane")
            details.ci <-
              paste0("Continuity correction of 1*incr for sample sizes\n",
                     "(Hartung & Knapp, 2001, Stat Med, equation (18)):\n\n")
          ##
          if (x$method.ci != "NAsm") {
            if (inherits(x, "metacont")) {
              catobsev(x$n.e + x$n.c, type = "n", addrow = TRUE)
              x.meta$n.e <- x.meta$n.c <- NA
            }
            else if (inherits(x, "metamean")) {
              catobsev(x$n, type = "n", addrow = TRUE)
              x.meta$n <- NA
            }
            else if (x$method.ci == "!RR.Cochrane") {
              catobsev(x$n.e + x$n.c, type = "n")
              catobsev(x$event.e + x$event.c, type = "e", addrow = TRUE)
              x.meta$n.e <- x.meta$event.e <-
                x.meta$n.c <- x.meta$event.c <- NA
            }
            else {
              catobsev(x$n, type = "n")
              catobsev(x$event, type = "e", addrow = TRUE)
              x.meta$n <- x.meta$event <- NA
            }
            ##
            dimnames(res) <-
              list(x$studlab,
                   c(sm.lab, ci.lab,
                     if (show.w.common) text.w.common,
                     if (show.w.random) text.w.random,
                     if (three.level) "cluster",
                     if (n_of_1) "cycles",
                     if (by) subgroup.name,
                     if (!is.null(x$exclude)) "exclude",
                     if (print.stat) "t",
                     if (print.pval) "p-value")
                   )
            prmatrix(res, quote = FALSE, right = TRUE)
            cat("\n")
          }
        }
        if (ma) {
          if (inherits(x, c("metaprop", "metarate", "metacont", "metamean")))
            cat("Normal approximation confidence interval:")
          else if (!is.null(x$method.ci) && x$method.ci == "!RR.Cochrane")
            cat("Mantel-Haenszel method:")
        }
        else {
          if (!(x$method.ci %in% c("t", "NAsm"))) {
            if (pscale != 1)
              sm.details <- paste0("\n- Events per ", pscale, " observations")
            else
              sm.details <- ""
            ##
            if (x$method.ci == "CP" & any(!is.na(x$pval))) {
              if (pscale != 1)
                sm.details <-
                  paste0(sm.details,
                         "\n- Null hypothesis: effect is equal to ",
                         format(round(x$null.effect * pscale, digits),
                                scientific = FALSE, big.mark = big.mark),
                         " events per ",
                         format(pscale, scientific = FALSE,
                                big.mark = big.mark),
                         " observations")
              else
                sm.details <-
                  paste0(sm.details,
                         "\n- Null hypothesis: effect is equal to ",
                         format(x$null.effect, scientific = FALSE,
                                big.mark = big.mark))
            }
            ##
            if (sm.details != "")
              cat(paste0("Details:", sm.details, "\n"))
          }
        }
      }
      else {
        dimnames(res) <-
          list(x$studlab,
               c(sm.lab, ci.lab,
                 if (show.w.common) text.w.common,
                 if (show.w.random) text.w.random,
                 if (three.level) "cluster",
                 if (n_of_1) "cycles",
                 if (by) subgroup.name,
                 if (show.imor) "IMOR.e",
                 if (show.imor) "IMOR.c",
                 if (!is.null(x$exclude)) "exclude"))
        ##
        if (!missing.truncate) {
          sortvar <- sortvar[truncate]
          res <- res[truncate, , drop = FALSE]
        }
        ##
        prmatrix(res[order(sortvar), , drop = FALSE],
                 quote = FALSE, right = TRUE)
        if (!missing.truncate)
          cat(text.truncate, "\n")
      }
    }
    
    
    ##
    ##
    ## (6) Print result for meta-analysis
    ##
    ##
    if (ma & !metainf.metacum) {
      if (!all(is.na(x$k)))
        cat("\n")
      ##
      attr(x.meta, ".print.study.results.") <- max(k.all, na.rm = TRUE) > 1
      print.meta(x.meta,
                 header = FALSE,
                 digits = digits,
                 common = common, random = random,
                 overall = overall,
                 backtransf = backtransf, pscale = pscale,
                 irscale = irscale, irunit = irunit,
                 #
                 print.tau2 = print.tau2, print.tau2.ci = print.tau2.ci,
                 digits.tau2 = digits.tau2,
                 #
                 print.tau = print.tau, print.tau.ci = print.tau.ci,
                 digits.tau = digits.tau,
                 #
                 print.Q = print.Q, digits.Q = digits.Q,
                 print.I2 = print.I2, digits.I2 = digits.I2,
                 print.I2.ci = print.I2.ci,
                 print.H = print.H, digits.H = digits.H,
                 print.Rb = print.Rb,
                 #
                 scientific.pval = scientific.pval, big.mark = big.mark,
                 zero.pval = zero.pval, JAMA.pval = JAMA.pval,
                 #
                 text.tau2 = text.tau2, text.tau = text.tau,
                 text.I2 = text.I2, text.Rb = text.Rb,
                 #
                 details.methods = details.methods,
                 warn.deprecated = FALSE,
                 ...)
    }
  }
  
  
  invisible(NULL)
}
