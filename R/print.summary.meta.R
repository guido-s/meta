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
#' @param fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param details A logical indicating whether further details of
#'   individual studies should be printed.
#' @param ma A logical indicating whether the summary results of the
#'   meta-analysis should be printed.
#' @param overall A logical indicating whether overall summaries
#'   should be reported. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be reported.
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
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
#' @param digits.pval Minimal number of significant digits for p-value
#'   of test of treatment effect, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
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
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param truncate An optional vector used to truncate the printout of
#'   results for individual studies (must be a logical vector of same
#'   length as \code{x$TE} or contain numerical values).
#' @param text.truncate A character string printed if study results
#'   were truncated from the printout.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param warn.backtransf A logical indicating whether a warning
#'   should be printed if backtransformed proportions and rates are
#'   below 0 and backtransformed proportions are above 1.
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
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
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
#'                data = Fleiss1993cont, sm = "SMD",
#'                studlab = paste(study, year))
#' sm1 <- summary(m1)
#' sm1
#' 
#' print(sm1, digits = 2)
#' 
#' \dontrun{
#' # Use unicode characters to print tau^2, tau, and I^2 
#' print(sm1,
#'       text.tau2 = "\u03c4\u00b2",
#'       text.tau = "\u03c4", text.I2 = "I\u00b2")
#' }
#' 
#' @method print summary.meta
#' @export
#' @export print.summary.meta


print.summary.meta <- function(x,
                               sortvar,
                               fixed = x$x$fixed,
                               random = x$x$random,
                               prediction = x$prediction,
                               details = FALSE, ma = TRUE,
                               overall = x$overall,
                               overall.hetstat = x$overall.hetstat,
                               ##
                               backtransf = x$backtransf,
                               pscale = x$pscale,
                               irscale = x$irscale,
                               irunit = x$irunit,
                               ##
                               digits = gs("digits"),
                               digits.se = gs("digits.se"),
                               digits.pval = max(gs("digits.pval"), 2),
                               digits.tau2 = gs("digits.tau2"),
                               digits.tau = gs("digits.tau"),
                               digits.I2 = gs("digits.I2"),
                               digits.prop = gs("digits.prop"),
                               digits.weight = gs("digits.weight"),
                               ##
                               scientific.pval = gs("scientific.pval"),
                               zero.pval = gs("zero.pval"),
                               JAMA.pval = gs("JAMA.pval"),
                               ##
                               big.mark = gs("big.mark"),
                               ##
                               text.tau2 = gs("text.tau2"),
                               text.tau = gs("text.tau"),
                               text.I2 = gs("text.I2"),
                               ##
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
  x.meta <- x$x
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  mf <- match.call()
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x.meta, stringsAsFactors = FALSE),
                               enclos = sys.frame(sys.parent())),
               silent = TRUE)
  if (class(error) == "try-error") {
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
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
  chklogical(prediction)
  chklogical(details)
  chklogical(ma)
  overall <- replaceNULL(overall, TRUE)
  chklogical(overall)
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  ##
  if (is.untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  ##
  chklogical(details.methods)
  ##
  chklogical(warn.backtransf)
  ##
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  if (!backtransf & pscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate(x$sm) & x$sm != "IRD")
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  if (!backtransf & irscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
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
  chkchar(text.tau2, length = 1)
  chkchar(text.tau, length = 1)
  chkchar(text.I2, length = 1)
  tt2 <- text.tau2
  tt <- text.tau
  ti <- text.I2
  ##
  ## Catch 'truncate' from meta-analysis object:
  ##
  missing.truncate <- missing(truncate)
  if (!missing.truncate) {
    truncate <- eval(mf[[match("truncate", names(mf))]],
                     x.meta, enclos = sys.frame(sys.parent()))
    ##
    if (is.null(truncate))
      truncate <- eval(mf[[match("truncate", names(mf))]],
                       x$data, enclos = sys.frame(sys.parent()))
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
  fixed <- replaceNULL(fixed, x$comb.fixed)
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed", FALSE)
  chklogical(fixed)
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
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  mb.glmm <- inherits(x, "metabind") | x$method == "GLMM"
  ##
  prediction <- prediction & x$k >= 3
  if (is.na(prediction))
    prediction <- FALSE
  ##
  ci.lab <- paste0(round(100 * level, 1), "%-CI")
  ##
  sm <- x$sm
  ##
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is.mean(sm))
      sm.lab <- "mean"
    else if (is.prop(sm)) {
      if (pscale == 1)
        sm.lab <- "proportion"
      else
        sm.lab <- "events"
    }
    else if (is.rate(sm)) {
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
  if (is.null(x$text.w.fixed))
    text.w.fixed <- "%W(fixed)"
  else
    text.w.fixed <- paste0("%W(", x$text.w.fixed, ")")
  ##
  if (is.null(x$text.w.random))
    text.w.random <- "%W(random)"
  else
    text.w.random <- paste0("%W(", x$text.w.random, ")")
  ##
  subgroup <- replaceNULL(x$subgroup, x$byvar)
  subgroup.name <- replaceNULL(x$subgroup.name, x$bylab)
  ##
  by <- !is.null(subgroup.name)
  id <- !is.null(x$three.level) && x$three.level
  
  
  ##
  ##
  ## (4) Print title and details
  ##
  ##
  if (inherits(x, "metamiss"))
    cat("Sensitivity analysis for missing binary data\n\n")
  ##
  crtitle(x)
  ##
  if (details) {
    if (inherits(x, "metamiss")) {
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
    if (id)
      res <- cbind(res, id = as.character(x$id))
    ##
    if (by)
      res <- cbind(res, subgroup = as.character(subgroup))
    ##
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
  if (k.all == 1 & !inherits(x, "metaprop")) {
    print.meta(x.meta,
               header = FALSE,
               digits = digits,
               backtransf = backtransf, pscale = pscale,
               irscale = irscale, irunit = irunit, big.mark = big.mark,
               text.tau2 = text.tau2, text.tau = text.tau, text.I2 = text.I2,
               details.methods = details.methods,
               warn.backtransf = warn.backtransf,
               ...)
  }
  else {
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    method.ci <- x$method.ci
    ##
    if (inherits(x, "metaprop") & !backtransf) {
      ciTE <- ci(TE, seTE, level = level)
      lowTE <- ciTE$lower
      uppTE <- ciTE$upper
      ##
      method.ci <- "NAsm"
    }
    ##
    if (backtransf) {
      ## Freeman-Tukey Arcsin transformation
      if (metainf.metacum) {
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
      else {
        TE    <- backtransf(   TE, sm, "mean",  harmonic.mean, warn.backtransf)
        lowTE <- backtransf(lowTE, sm, "lower", harmonic.mean, warn.backtransf)
        uppTE <- backtransf(uppTE, sm, "upper", harmonic.mean, warn.backtransf)
      }
      ##
      if (is.prop(sm) | sm == "RD") {
        TE <- pscale * TE
        lowTE <- pscale * lowTE
        uppTE <- pscale * uppTE
      }
      ##
      if (is.rate(sm) | sm == "IRD") {
        TE <- irscale * TE
        lowTE <- irscale * lowTE
        uppTE <- irscale * uppTE
      }
    }
    ##
    TE <- round(TE, digits)
    lowTE <- round(lowTE, digits)
    uppTE <- round(uppTE, digits)
    ##
    if (!metainf.metacum) {
      if (fixed)
        if (!all(is.na(x$w.fixed)) && sum(x$w.fixed) > 0)
          w.fixed.p <- round(100 * x$w.fixed / sum(x$w.fixed, na.rm = TRUE),
                             digits.weight)
        else w.fixed.p <- x$w.fixed
      ##
      if (random)
        if (!is.null(x$w.random) & !all(is.na(x$w.random)) &&
            sum(x$w.random) > 0)
          w.random.p <- round(100 * x$w.random / sum(x$w.random, na.rm = TRUE),
                              digits.weight)
        else w.random.p <- x$w.random
    }
    ##
    if (metainf.metacum) {
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
                   paste0(" ", tau2),
                   paste0(" ", tau),
                   paste0(" ", I2, ifelse(I2 == "", "", "%")))
      dimnames(res) <- list(paste0(x$studlab, "  "),
                            c(sm.lab, ci.lab, "p-value",
                              text.tau2, text.tau, text.I2))
      ##
      if (inherits(x, "metainf")) {
        if (!is.random)
          cat("Influential analysis (Fixed effect model)\n")
        else
          cat("Influential analysis (Random effects model)\n")
      }
      else if (inherits(x, "metacum")) {
        if (!is.random)
          cat("Cumulative meta-analysis (Fixed effect model)\n")
        else
          cat("Cumulative meta-analysis (Random effects model)\n")
      }
      cat("\n")
      prmatrix(res, quote = FALSE, right = TRUE, na.print = "--")
      ## Print information on summary method:
      if (details.methods)
        catmeth(class = class(x),
                method = x$method,
                method.tau = x$method.tau,
                sm = sm,
                k.all = k.all,
                hakn = is.random & x$hakn,
                tau.preset = x$tau.preset,
                method.smd = x$method.smd,
                sd.glass = x$sd.glass,
                exact.smd = x$exact.smd,
                model.glmm = x$model.glmm,
                big.mark = big.mark,
                digits = digits, digits.tau = digits.tau,
                text.tau = text.tau, text.tau2 = text.tau2,
                method.miss = x$method.miss,
                IMOR.e = x$IMOR.e, IMOR.c = x$IMOR.c)
    }
    else if (!(inherits(x, "metabind") && !x$show.studies)) {
      show.w.fixed  <-
        (overall | by) & !mb.glmm &
        (fixed && !all(is.na(w.fixed.p)))
      show.w.random <-
        (overall | by) & !mb.glmm &
        (random && !all(is.na(w.random.p)))
      ##
      res <- cbind(formatN(round(TE, digits), digits, "NA",
                           big.mark = big.mark),
                   formatCI(formatN(round(lowTE, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE, digits), digits, "NA",
                                    big.mark = big.mark)),
                   if (show.w.fixed)
                     formatN(w.fixed.p, digits.weight,
                             big.mark = big.mark),
                   if (show.w.random)
                     formatN(w.random.p, digits.weight,
                             big.mark = big.mark),
                   if (id) as.character(x$id),
                   if (by) as.character(subgroup),
                   if (!is.null(x$exclude))
                     ifelse(is.na(x$exclude), "",
                     ifelse(x$exclude, "*", "")))
      ## Printout for a single proportion:
      if (k.all == 1) {
        ##
        if (!is.null(method.ci)) {
          if (method.ci == "CP")
            method.ci.details <-
              "Clopper-Pearson confidence interval:\n\n"
          else if (method.ci == "WS")
            method.ci.details <-
              "Wilson Score confidence interval:\n\n"
          else if (method.ci == "WSCC")
            method.ci.details <-
              "Wilson Score confidence interval with continuity correction:\n\n"
          else if (method.ci == "AC")
            method.ci.details <-
              "Agresti-Coull confidence interval:\n\n"
          else if (method.ci == "SA")
            method.ci.details <-
              "Simple approximation confidence interval:\n\n"
          else if (method.ci == "SACC")
            method.ci.details <-
              paste0("Simple approximation confidence interval with ",
                     "continuity correction:\n\n")
          else if (method.ci == "t")
            method.ci.details <-
              "Confidence interval based on t-distribution:\n\n"
          if (method.ci != "NAsm") {
            cat(method.ci.details)
            dimnames(res) <-
              list("",
                   c(sm.lab, ci.lab,
                     if (show.w.fixed) text.w.fixed,
                     if (show.w.random) text.w.random,
                     if (id) "id",
                     if (by) subgroup.name,
                     if (!is.null(x$exclude)) "exclude"))
            prmatrix(res, quote = FALSE, right = TRUE)
            cat("\n")
          }
        }
        cat("Normal approximation confidence interval:\n")
      }
      else {
        dimnames(res) <-
          list(x$studlab,
               c(sm.lab, ci.lab,
                 if (show.w.fixed) text.w.fixed,
                 if (show.w.random) text.w.random,
                 if (id) "id",
                 if (by) subgroup.name,
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
      if (!is.na(x$k))
        cat("\n")
      ##
      attr(x.meta, ".print.study.results.") <- k.all > 1
      print.meta(x.meta,
                 header = FALSE,
                 digits = digits,
                 fixed = fixed, random = random,
                 prediction = prediction,
                 overall = overall, overall.hetstat = overall.hetstat,
                 backtransf = backtransf, pscale = pscale,
                 irscale = irscale, irunit = irunit,
                 digits.tau2 = digits.tau2, digits.tau = digits.tau,
                 digits.I2 = digits.I2, big.mark = big.mark,
                 text.tau2 = text.tau2, text.tau = text.tau, text.I2 = text.I2,
                 details.methods = details.methods,
                 warn.backtransf = warn.backtransf,
                 warn.deprecated = FALSE,
                 ...)
    }
  }
  
  
  invisible(NULL)
}
