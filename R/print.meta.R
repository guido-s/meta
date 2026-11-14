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
#' @param common A logical indicating whether results for common
#'   effect meta-analysis should be printed.
#' @param random A logical indicating whether results for random
#'   effects meta-analysis should be printed.
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
#' @param test A named list controlling tests for subgroup differences; see
#'   \code{\link{meta-formatting}}. Defaults follow the setting in \code{x}
#'   and the models selected for printing.
#' @param prediction.subgroup A single logical or logical vector
#'   indicating whether / which prediction intervals should be printed
#'   for subgroups.
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
#' @param sep.subgroup A character string defining the separator
#'   between label and levels of grouping variable.
#' @param nchar.subgroup A numeric specifying the number of characters
#'   to print from subgroup labels.
#' @param sort.overall An optional vector used to sort meta-analysis
#'   results.
#' @param sort.tau An optional vector used to sort estimators of the
#'   between-study heterogeneity variance.
#' @param sort.het An optional vector used to sort heterogeneity
#'   statistics.
#' @param sort.Q An optional vector used to sort test of
#'   heterogeneity.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param print A named list controlling displayed information; see
#'   \code{\link{meta-formatting}}. Heterogeneity settings default to
#'   the corresponding \code{gs("print.*")} settings; \code{CMH} and
#'   \code{subgroup.name} default
#'   to the corresponding settings in \code{x}.
#' @param digits A named list controlling rounding; see
#'   \code{\link{meta-formatting}}. A scalar value is also accepted.
#'   The main default is \code{gs("digits")}; P-value defaults are
#'   \code{max(gs("digits.pval"), 2)} and \code{max(gs("digits.pval.Q"), 2)}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
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
#' @param warn.backtransf Deprecated argument (ignored).
#' @param func.backtransf A function used to back-transform results.
#' @param bracket A character with bracket symbol to print lower
#'   confidence interval: "[", "(", "\{", "".
#' @param separator A character string with information on separator
#'   between lower and upper confidence interval.
#' @param lower.blank A logical indicating whether blanks between left
#'   bracket and lower confidence limit should be printed.
#' @param upper.blank A logical indicating whether blanks between
#'   separator and upper confidence limit should be printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments, including individual \code{test.*},
#'   \code{print.*}, and \code{digits.*} settings. Remaining arguments are passed to
#'   \code{prmatrix}.
#' 
#' @rdname print.meta
#' @method print meta
#' @export

print.meta <- function(x,
                       common = x$common,
                       random = x$random,
                       prediction = x$prediction,
                       overall = x$overall,
                       overall.hetstat = x$overall.hetstat,
                       #
                       prediction.subgroup = x$prediction.subgroup,
                       #
                       backtransf = x$backtransf,
                       pscale = x$pscale,
                       irscale = x$irscale,
                       irunit = x$irunit,
                       #
                       subgroup.name = x$subgroup.name,
                       sep.subgroup = x$sep.subgroup,
                       nchar.subgroup = 35,
                       #
                       sort.overall = NULL,
                       sort.tau = NULL,
                       sort.het = NULL,
                       sort.Q = NULL,
                       #
                       header = TRUE,
                       #
                       test = NULL,
                       print = NULL,
                       digits = NULL,
                       #
                       big.mark = gs("big.mark"),
                       scientific.pval = gs("scientific.pval"),
                       zero.pval = gs("zero.pval"),
                       JAMA.pval = gs("JAMA.pval"),
                       #
                       text.tau2 = gs("text.tau2"),
                       text.tau = gs("text.tau"),
                       text.I2 = gs("text.I2"),
                       text.Rb = gs("text.Rb"),
                       #
                       details.methods = gs("details"),
                       #
                       warn.backtransf = FALSE,
                       func.backtransf = x$func.backtransf,
                       warn.deprecated = gs("warn.deprecated"),
                       ...) {
  
  
  #
  #
  # (1) Check for meta object
  #
  #
  chkclass(x, "meta")
  x <- updateversion(x)
  #
  if (inherits(x, "metainf"))
    return(invisible(NULL))
  #
  args <- list(...)
  #
  if (anyDuplicated(names(args)) || anyDuplicated(names(test)) ||
      anyDuplicated(names(print)) ||
      anyDuplicated(names(digits)))
    stop("Formatting arguments must not be duplicated.", call. = FALSE)
  #
  # Tests for subgroup differences
  #
  test.list <- merge_legacy(test, args, test_defaults(x, func = "print"))
  #
  test <- test.list$res
  args <- test.list$args
  #
  test.subgroup <- test[["subgroup"]]
  test.subgroup.common <- test[["subgroup.common"]]
  test.subgroup.random <- test[["subgroup.random"]]
  #
  # Columns shown when printing
  #
  print.list <-
    merge_legacy(print, args, print_defaults(x, func = "print"))
  #
  print.settings <- print.list$res
  args <- print.list$args
  #
  print.subgroup.name <- print.settings[["subgroup.name"]]
  print.CMH <- print.settings[["CMH"]]
  print.tau2 <- print.settings[["tau2"]]
  print.tau2.ci <- print.settings[["tau2.ci"]]
  print.tau <- print.settings[["tau"]]
  print.tau.ci <- print.settings[["tau.ci"]]
  print.Q <- print.settings[["Q"]]
  print.I2 <- print.settings[["I2"]]
  print.I2.ci <- print.settings[["I2.ci"]]
  print.H <- print.settings[["H"]]
  print.Rb <- print.settings[["Rb"]]
  #
  # Digits used to round results
  #
  digits <- as_named_list(digits)
  digits.list <-
    merge_legacy(digits, args, digits_defaults(func = "print"))
  #
  digits.settings <- digits.list$res
  args <- digits.list$args
  #
  digits <- digits.settings[["digits"]]
  digits.stat <- digits.settings[["stat"]]
  digits.pval <- digits.settings[["pval"]]
  digits.tau2 <- digits.settings[["tau2"]]
  digits.tau <- digits.settings[["tau"]]
  digits.Q <- digits.settings[["Q"]]
  digits.df <- digits.settings[["df"]]
  digits.pval.Q <- digits.settings[["pval.Q"]]
  digits.H <- digits.settings[["H"]]
  digits.I2 <- digits.settings[["I2"]]
  #
  is.netpairwise <- inherits(x, "netpairwise")
  #
  method <- x$method
  method.random <- x$method.random
  #
  by <- !is.null(x$subgroup)
  #
  if (by) {
    n.by <- length(x$subgroup.levels)
    #
    Q.resid <- replaceNULL(x$Q.resid)
    df.Q.resid <- replaceNULL(x$df.Q.resid)
    pval.Q.resid <- replaceNULL(x$pval.Q.resid)
  }
  #
  method.random.ci <- replaceNULL(x$method.random.ci, "")
  method.predict <- replaceNULL(x$method.predict, "")
  #
  if (is.function(func.backtransf))
    fbt <- deparse(substitute(func.backtransf))
  else
    fbt <- func.backtransf
  #
  abt <- x$args.backtransf
  #
  print.I2 <- print.I2 & ((overall.hetstat & !by) | by)
  #
  pairwise <- replaceNULL(x$pairwise, FALSE) & any(table(x$studlab) > 1)
  
  
  #
  #
  # (2) Check and set other arguments
  #
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.H, min = 0, length = 1)
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
  if (all(method == "LRP") & all(method.random == "LRP")) {
    print.tau2 <- FALSE
    print.tau <- FALSE
  }
  #
  chklogical(print.I2)
  chklogical(print.I2.ci)
  chklogical(print.H)
  chklogical(print.Rb)
  #
  chkchar(text.tau2, length = 1)
  chkchar(text.tau, length = 1)
  chkchar(text.I2, length = 1)
  chkchar(text.Rb, length = 1)
  chklogical(print.Q)
  #
  chklogical(warn.deprecated)
  #
  missing.pscale <- missing(pscale)
  pscale <- replaceNULL(pscale, 1)
  chknumeric_strict(pscale, min = 1, length = 1)
  #
  missing.irscale <- missing(irscale)
  irscale <- replaceNULL(irscale, 1)
  chknumeric_strict(irscale, min = 1, length = 1)
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  #
  chklogical(common)
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
  test.subgroup.common <- replaceNULL(test.subgroup.common,
                                      if (test.list$miss[["subgroup.common"]])
                                        test.subgroup & common else test.subgroup)
  chklogical(test.subgroup.common)
  test.subgroup.random <- replaceNULL(test.subgroup.random,
                                      if (test.list$miss[["subgroup.random"]])
                                        test.subgroup & random else test.subgroup)
  chklogical(test.subgroup.random)
  #
  missing.prediction.subgroup <- missing(prediction.subgroup)
  #
  if (by) {
    x$lower.predict.w <- list2mat(x$lower.predict.w)
    x$upper.predict.w <- list2mat(x$upper.predict.w)
    x$df.predict.w <- list2mat(x$df.predict.w)
    #
    if (!missing.prediction.subgroup)
      prediction.subgroup <- catch(enquo(prediction.subgroup), x)
    #
    prediction.subgroup <- replaceNULL(prediction.subgroup, FALSE)
    #
    if (length(prediction.subgroup) == 1) {
      if (is.matrix(x$lower.predict.w)) {
        prediction.subgroup.logical <-
          prediction.subgroup &
          apply(x$lower.predict.w, 1, notallNA) &
          apply(x$upper.predict.w, 1, notallNA)
      }
      else {
        prediction.subgroup.logical <-
          prediction.subgroup &
          notallNA(x$lower.predict.w) &
          notallNA(x$upper.predict.w)
        #
        prediction.subgroup.logical <-
          rep(prediction.subgroup.logical, n.by)
      }
    }
    else {
      chklength(prediction.subgroup, n.by,
                text = paste("Length of argument 'prediction.subgroup' must be",
                             "equal to 1 or number of subgroups."))
      prediction.subgroup.logical <- prediction.subgroup
    }
    #
    if (length(prediction.subgroup.logical) == 0)
      prediction.subgroup.logical <- FALSE
    #
    chklogical(prediction.subgroup.logical[1])
    #
    if (missing.prediction.subgroup) {
      if (any(method.predict %in% c("", "S")))
        prediction.w <- prediction.subgroup.logical
      else
        prediction.w <- prediction.subgroup.logical & !is.na(x$df.predict.w)
    }
    else
      prediction.w <- prediction.subgroup.logical
    #
    prediction.w[is.na(prediction.w)] <- FALSE
    prediction.w <- any(prediction.w)
  }
  else
    prediction.w <- FALSE
  #
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  #
  sort.tau <- catch(enquo(sort.tau), x)
  if (is.list(sort.tau))
    sort.tau <- unlist(sort.tau)
  #
  sort.het <- catch(enquo(sort.het), x)
  if (is.list(sort.het))
    sort.het <- unlist(sort.het)
  #
  sort.Q <- catch(enquo(sort.Q), x)
  if (is.list(sort.Q))
    sort.Q <- unlist(sort.Q)
  #
  chklogical(header)
  #
  chklogical(details.methods)
  #
  # Check for deprecated arguments in '...'
  #
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  #
  random <- deprecated(random, missing(random), args, "comb.random", FALSE)
  chklogical(random)
  #
  subgroup.name <-
    deprecated(subgroup.name, missing(subgroup.name), args, "bylab", FALSE)
  #
  print.subgroup.name <-
    deprecated(print.subgroup.name, print.list$miss[["subgroup.name"]],
               args, "print.byvar", FALSE)
  #
  sep.subgroup <-
    deprecated(sep.subgroup, missing(sep.subgroup), args, "byseparator",
               FALSE)
  #
  nchar.subgroup <-
    deprecated(nchar.subgroup, missing(nchar.subgroup), args, "nchar.subgroup",
               warn.deprecated)
  chknumeric(nchar.subgroup, min = 1, length = 1)
  #
  backtransf <-
    deprecated(backtransf, missing(backtransf), args, "logscale",
               warn.deprecated)
  if (is_untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  #
  # Additional settings
  #
  if (!backtransf & pscale != 1 & !is_untransformed(x$sm)) {
    if (!missing.pscale)
      warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!backtransf & irscale != 1 & !is_untransformed(x$sm)) {
    if (!missing.irscale)
      warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  #
  if (by) {
    chklogical(print.subgroup.name)
    chkchar(sep.subgroup)
  }
  else {
    if (!missing.prediction.subgroup)
      warning("Argument 'prediction.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
  }
  
  
  #
  #
  # (3) Some additional settings
  #
  #
  k.all <- replaceNULL(x$k.all, length(x$TE))
  k <- x$k
  k.study <- replaceNULL(x$k.study, k)
  sm <- x$sm
  #
  bip <- inherits(x, c("metabin", "metainc", "metaprop", "metarate"))
  metabin <- inherits(x, "metabin")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  #
  if (inherits(x, "trimfill")) {
    harmonic.mean <- NULL
  }
  else if (sm == "IRFT")
    harmonic.mean <- x$t.harmonic.mean
  else
    harmonic.mean <- replaceNULL(x$n.harmonic.mean)
  #
  # Null-effect display conversion is handled by catmeth()
  #
  null.given <- !is.null(x$null.effect) && !is.na(x$null.effect)
  #
  if (is.null(prediction) || is.na(prediction))
    prediction <- FALSE
  #
  sm.lab <- smlab(sm, backtransf, pscale, irscale)
  #
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  #
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  #
  if (by)
    subgroup.levels <-
      ifelse(nchar(x$subgroup.levels) > nchar.subgroup,
             paste0(substring(x$subgroup.levels, 1, nchar.subgroup - 4),
                    " ..."),
             x$subgroup.levels)
  #
  if (is.null(x$text.common))
    text.common <- gs("text.common")
  else
    text.common <- x$text.common
  #
  if (is.null(x$text.random))
    text.random <- gs("text.random")
  else
    text.random <- x$text.random
  #
  if (is.null(x$text.predict))
    text.predict <- gs("text.predict")
  else
    text.predict <- x$text.predict
  #
  if (any(substring(text.common, 1, 5) %in% c("Fixed", "Commo"))) {
    text.common.br <- gsub("Fixed", "fixed", text.common)
    text.common.br <- gsub("Common", "common", text.common.br)
    text.common.br <- gsub("Effect", "effect", text.common.br)
  }
  else
    text.common.br <- text.common
  #
  if (any(substring(text.random, 1, 5) %in% c("Rando"))) {
    text.random.br <- gsub("Random", "random", text.random)
    text.random.br <- gsub("Effect", "effect", text.random.br)
  }
  else
    text.random.br <- text.random
  #
  if (common | random)
    ci.lab <- percent_label(x$level.ma)
  else if (prediction)
    ci.lab <- percent_label(x$level.predict, "PI")
  else if (all(k.all == 1))
    ci.lab <- percent_label(x$level.ma)
  else
    ci.lab <- ""
  #
  details <- NULL
  
  
  #
  #
  # (4) Set and backtransform results of meta-analysis
  #
  #
  TE.common <- x$TE.common
  lowTE.common <- x$lower.common
  uppTE.common <- x$upper.common
  pTE.common <- x$pval.common
  sTE.common <- x$statistic.common
  #
  TE.random <- x$TE.random
  lowTE.random <- x$lower.random
  uppTE.random <- x$upper.random
  if (length(TE.random) == 1 &&
      length(TE.random) != length(lowTE.random))
    TE.random <- rep_len(TE.random, length(lowTE.random))
  #
  lowTE.predict <- x$lower.predict
  uppTE.predict <- x$upper.predict
  #
  Q <- unlist(x$Q)
  df.Q <- replaceNULL(unlist(x$df.Q), k - 1)
  pval.Q <- replaceNULL(unlist(x$pval.Q), pvalQ(Q, df.Q))
  #
  if (!is.null(x$Q.CMH)) {
    Q.CMH <- unlist(x$Q.CMH)
    df.Q.CMH <- replaceNULL(unlist(x$df.Q.CMH), 1)
    pval.Q.CMH <- replaceNULL(unlist(x$pval.Q.CMH), pvalQ(Q.CMH, df.Q.CMH))
  }
  #
  if (by) {
    k.w <- x$k.w
    #
    TE.common.w    <- list2mat(x$TE.common.w)
    lowTE.common.w <- list2mat(x$lower.common.w)
    uppTE.common.w <- list2mat(x$upper.common.w)
    pval.common.w  <- list2mat(x$pval.common.w)
    #
    TE.random.w    <- list2mat(x$TE.random.w, x$lower.random.w)
    lowTE.random.w <- list2mat(x$lower.random.w)
    uppTE.random.w <- list2mat(x$upper.random.w)
    pval.random.w  <- list2mat(x$pval.random.w)
    #
    if (is.matrix(lowTE.random.w) &&
        (!is.matrix(TE.random.w) ||
         (ncol(TE.random.w) == 1 && ncol(lowTE.random.w) > 1))) {
      TE.random.w <-
        matrix(TE.random.w,
               nrow = nrow(lowTE.random.w),
               ncol = ncol(lowTE.random.w))
      rownames(TE.random.w) <- rownames(lowTE.random.w)
      colnames(TE.random.w) <- colnames(lowTE.random.w)
    }
    #
    lowTE.predict.w <- x$lower.predict.w
    uppTE.predict.w <- x$upper.predict.w
    #
    harmonic.mean.w <- if (sm == "IRFT")
                         x$t.harmonic.mean.w
                       else
                         x$n.harmonic.mean.w
    #
    Q.b.common <- unlist(x$Q.b.common)
    names(Q.b.common) <- colnames(lowTE.common.w)
    Q.w.common <- x$Q.w.common
    #
    Q.b.random <- unlist(x$Q.b.random)
    #
    Q.w <- x$Q.w
    #
    df.Q.w <- replaceNULL(x$df.Q.w, sum((k.w - 1)[!is.na(x$Q.w)]))
    df.Q.b <- replaceNULL(x$df.Q.b, (k - 1) - sum((k.w - 1)[!is.na(x$Q.w)]))
    df.Q.b.common <- replaceNULL(x$df.Q.b.common, df.Q.b)
    df.Q.b.random <- replaceNULL(x$df.Q.b.random, df.Q.b)    
    #
    if (common && test.subgroup.common && length(Q.b.common) > 0L) {
      df.Q.b.common <- subgroup_test_df(df.Q.b.common, length(Q.b.common))
      pval.Q.b.common <-
        replaceNULL(x$pval.Q.b.common, pval_subgroup(Q.b.common, df.Q.b.common))
    }
    pval.Q.w.common <-
      replaceNULL(x$pval.Q.w.common, pvalQ(Q.w.common, df.Q.w))
    #
    if (random && test.subgroup.random && length(Q.b.random) > 0L) {
      df.Q.b.random <- subgroup_test_df(df.Q.b.random, length(Q.b.random))
      pval.Q.b.random <-
        unlist(replaceNULL(x$pval.Q.b.random, pval_subgroup(Q.b.random, df.Q.b.random)))
    }
  }
  #
  if (backtransf) {
    TE.common <-
      backtransf(TE.common, sm, harmonic.mean, harmonic.mean, fbt, abt)
    lowTE.common <-
      backtransf(lowTE.common, sm, harmonic.mean, harmonic.mean, fbt, abt)
    uppTE.common <-
      backtransf(uppTE.common, sm, harmonic.mean, harmonic.mean, fbt, abt)
    #
    TE.random <-
      backtransf(TE.random, sm, harmonic.mean, harmonic.mean, fbt, abt)
    lowTE.random <-
      backtransf(lowTE.random, sm, harmonic.mean, harmonic.mean, fbt, abt)
    uppTE.random <-
      backtransf(uppTE.random, sm, harmonic.mean, harmonic.mean, fbt, abt)
    #
    lowTE.predict <-
      backtransf(lowTE.predict, sm, harmonic.mean, harmonic.mean, fbt, abt)
    uppTE.predict <-
      backtransf(uppTE.predict, sm, harmonic.mean, harmonic.mean, fbt, abt)
    #
    if (by) {
      TE.common.w <-
        backtransf(TE.common.w, sm, harmonic.mean.w, harmonic.mean.w,
                   fbt, abt)
      lowTE.common.w <-
        backtransf(lowTE.common.w, sm, harmonic.mean.w, harmonic.mean.w,
                   fbt, abt)
      uppTE.common.w  <-
        backtransf(uppTE.common.w, sm, harmonic.mean.w, harmonic.mean.w)
      #
      TE.random.w <-
        backtransf(TE.random.w, sm, harmonic.mean.w, harmonic.mean.w, fbt, abt)
      lowTE.random.w <-
        backtransf(lowTE.random.w, sm, harmonic.mean.w, harmonic.mean.w,
                   fbt, abt)
      uppTE.random.w <-
        backtransf(uppTE.random.w, sm, harmonic.mean.w, harmonic.mean.w,
                   fbt, abt)
      #
      if (prediction.w) {
        lowTE.predict.w <-
          backtransf(lowTE.predict.w, sm, harmonic.mean.w, harmonic.mean.w,
                     fbt, abt)
        uppTE.predict.w <-
          backtransf(uppTE.predict.w, sm, harmonic.mean.w, harmonic.mean.w,
                     fbt, abt)
      }
    }
  }
  #
  # Apply argument 'pscale' to proportions / risk differences and
  # 'irscale' to rates / incidence rate differences
  #
  scale <- 1
  #
  if (pscale != 1 && irscale != 1)
    stop("Provide either arguments 'pscale' or 'irscale' but not both arguments.",
         call. = FALSE)
  #
  if ((is_prop(sm) || sm == "RD") && pscale != 1)
    scale <- pscale
  else if ((is_rate(sm) || sm == "IRD") && irscale != 1)
    scale <- irscale
  #
  if (scale != 1) {
    #
    TE.common    <- scale * TE.common
    lowTE.common <- scale * lowTE.common
    uppTE.common <- scale * uppTE.common
    #
    TE.random    <- scale * TE.random
    lowTE.random <- scale * lowTE.random
    uppTE.random <- scale * uppTE.random
    #
    lowTE.predict <- scale * lowTE.predict
    uppTE.predict <- scale * uppTE.predict
    #
    if (by) {
      TE.common.w    <- scale * TE.common.w
      lowTE.common.w <- scale * lowTE.common.w
      uppTE.common.w <- scale * uppTE.common.w
      #   
      TE.random.w    <- scale * TE.random.w
      lowTE.random.w <- scale * lowTE.random.w
      uppTE.random.w <- scale * uppTE.random.w
      #   
      if (prediction.w) {
        lowTE.predict.w <- scale * lowTE.predict.w
        uppTE.predict.w <- scale * uppTE.predict.w
      }
    }
  }
  #
  # Switch reversed lower and upper limits after back-transformation
  #
  if (backtransf) {
    wo <- !is.na(lowTE.common) & !is.na(uppTE.common) &
      lowTE.common > uppTE.common
    #
    tmp.l <- lowTE.common[wo]
    lowTE.common[wo] <- uppTE.common[wo]
    uppTE.common[wo] <- tmp.l
    #
    wo <- !is.na(lowTE.random) & !is.na(uppTE.random) &
      lowTE.random > uppTE.random
    #
    tmp.l <- lowTE.random[wo]
    lowTE.random[wo] <- uppTE.random[wo]
    uppTE.random[wo] <- tmp.l
    #
    wo <- !is.na(lowTE.predict) & !is.na(uppTE.predict) &
      lowTE.predict > uppTE.predict
    #
    tmp.l <- lowTE.predict[wo]
    lowTE.predict[wo] <- uppTE.predict[wo]
    uppTE.predict[wo] <- tmp.l
    #
    if (by) {
      wo <- !is.na(lowTE.common.w) & !is.na(uppTE.common.w) &
        lowTE.common.w > uppTE.common.w
      #
      tmp.l <- lowTE.common.w[wo]
      lowTE.common.w[wo] <- uppTE.common.w[wo]
      uppTE.common.w[wo] <- tmp.l
      #
      wo <- !is.na(lowTE.random.w) & !is.na(uppTE.random.w) &
        lowTE.random.w > uppTE.random.w
      #
      tmp.l <- lowTE.random.w[wo]
      lowTE.random.w[wo] <- uppTE.random.w[wo]
      uppTE.random.w[wo] <- tmp.l
      #
      if (prediction.w) {
        wo <- !is.na(lowTE.predict.w) & !is.na(uppTE.predict.w) &
          lowTE.predict.w > uppTE.predict.w
        #
        tmp.l <- lowTE.predict.w[wo]
        lowTE.predict.w[wo] <- uppTE.predict.w[wo]
        uppTE.predict.w[wo] <- tmp.l
      }
    }
  }
  #
  # Round and round ...
  #
  TE.common    <- round(TE.common, digits)
  lowTE.common <- round(lowTE.common, digits)
  uppTE.common <- round(uppTE.common, digits)
  sTE.common <- round(sTE.common, digits.stat)
  #
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$pval.random
  sTE.random <- round(x$statistic.random, digits.stat)
  #
  lowTE.predict <- round(lowTE.predict, digits)
  uppTE.predict <- round(uppTE.predict, digits)
  #
  if (by) {
    TE.common.w     <- round(TE.common.w, digits)
    lowTE.common.w  <- round(lowTE.common.w, digits)
    uppTE.common.w  <- round(uppTE.common.w, digits)
    #
    TE.random.w    <- round(TE.random.w, digits)
    lowTE.random.w <- round(lowTE.random.w, digits)
    uppTE.random.w <- round(uppTE.random.w, digits)
    #
    if (prediction.w) {
      lowTE.predict.w <- round(lowTE.predict.w, digits)
      uppTE.predict.w <- round(uppTE.predict.w, digits)
    }
    #
    if (print.I2)
      I2.w <- round(100 * x$I2.w, digits.I2)
    else
      I2.w <- NA
    #
    if (print.Rb)
      Rb.w <- round(100 * x$Rb.w, digits.I2)
    else
      Rb.w <- NA
  }
  #
  if (print.H) {
    H <- round(unlist(x$H), digits.H)
    lowH <- round(unlist(x$lower.H), digits.H)
    uppH <- round(unlist(x$upper.H), digits.H)
    if (all(!is.na(lowH) & !is.na(uppH)) && all(lowH == uppH)) {
      lowH <- NA
      uppH <- NA
    }
  }
  else {
    H <- lowH <- uppH <- NULL
  }
  #
  if (print.I2) {
    I2 <- round(100 * unlist(x$I2), digits.I2)
    lowI2 <- round(100 * unlist(x$lower.I2), digits.I2)
    uppI2 <- round(100 * unlist(x$upper.I2), digits.I2)
    #
    if (is.null(names(Q)))
      Q4I2 <- Q
    else
      Q4I2 <- Q[!grepl("LRT", names(Q))]
    #
    print.I2.ci <- print.I2.ci &
      ifelse(((Q4I2 > max(k, na.rm = TRUE) & max(k, na.rm = TRUE) >= 2) |
              (Q4I2 <= max(k, na.rm = TRUE) & max(k, na.rm = TRUE) > 2)) &
             !(is.na(lowI2) | is.na(uppI2)), TRUE, FALSE)
    #
    print.I2.ci[is.na(print.I2.ci)] <- FALSE
    print.I2.ci <- print.I2.ci & !(lowI2 == uppI2)
    print.I2.ci <- any(print.I2.ci)
  }
  else {
    I2 <- lowI2 <- uppI2 <- NULL
    print.I2.ci <- FALSE
  }
  #
  if (print.Rb) {
    Rb <- round(100 * unlist(x$Rb), digits.I2)
    lowRb <- round(100 * unlist(x$lower.Rb), digits.I2)
    uppRb <- round(100 * unlist(x$upper.Rb), digits.I2)
    #
    if (all(!is.na(lowRb) & !is.na(uppRb)) && all(lowRb == uppRb)) {
      lowRb <- rep(NA, length(Rb))
      uppRb <- rep(NA, length(Rb))
    }
  }
  else {
    Rb <- lowRb <- uppRb <- NULL
    print.Rb.ci <- FALSE    
  }
  #
  three.level <- if (is.null(x$three.level)) FALSE else any(x$three.level)
  is.glmm <- any(method == "GLMM") | any(method.random == "GLMM")
  #
  sel.n <- inherits(x, c("metacor", "metaprop", "metamean", "metarate"))
  #
  sel.ev <- inherits(x, c("metaprop", "metarate"))
  
  
  #
  #
  # (5) Print result for meta-analysis
  #
  #
  is.metamiss <- inherits(x, "metamiss")
  if (header) {
    if (is.metamiss)
      cat("Sensitivity analysis for missing binary data\n")
    #
    crtitle(x)
  }
  #
  if (all(is.na(k.all))) {
    # Do nothing
    return(invisible(NULL))
  }
  else if (all(k.all == 1)) {
    #
    # Print results for a single study
    #
    if (!(metaprop | metarate) | (overall & (common | random)) |
        overall.hetstat | by)
      method <- ifelse(method == "Peto", method, "")
    #
    if (is.metamiss)
      cat("\n")
    else if (sel.n)
      catobsev(if (is.null(x$exclude)) x$n else x$n[!x$exclude],
               type = "n")
    else
      catobsev(if (is.null(x$exclude)) x$n.e else x$n.e[!x$exclude],
               if (is.null(x$exclude)) x$n.c else x$n.c[!x$exclude],
               type = "n")
    #
    if (sel.ev)
      catobsev(if (is.null(x$exclude)) x$event else x$event[!x$exclude],
               type = "e", addrow = TRUE)
    else if (!is.metamiss)
      catobsev(if (is.null(x$exclude)) x$event.e else x$event.e[!x$exclude],
               if (is.null(x$exclude)) x$event.c else x$event.c[!x$exclude],
               type = "e", addrow = TRUE)
    #
    res <- cbind(formatN(TE.common, digits, "NA",
                         big.mark = big.mark),
                 formatCI(formatN(lowTE.common, digits, "NA",
                                  big.mark = big.mark),
                          formatN(uppTE.common, digits, "NA",
                                  big.mark = big.mark)),
                 if (null.given)
                   formatN(sTE.common, digits.stat, big.mark = big.mark),
                 if (null.given)
                   formatPT(pTE.common, digits = digits.pval,
                            scientific = scientific.pval,
                            zero = zero.pval, JAMA = JAMA.pval))
    dimnames(res) <- list(x$studlab,
                          c(sm.lab, ci.lab,
                            if (null.given) "z",
                            if (null.given) "p-value"))
    do.call(prmatrix, c(list(res, quote = FALSE, right = TRUE), args))
    #
    # Print information on meta-analysis method:
    #
    x$pscale <- pscale
    x$irscale <- irscale
    x$irunit <- irunit
    #
    if (!is.null(x$cycles))
      x$sd.n_of_1 <-
        formatPT(x$sd.n_of_1, digits = gs("digits.sd"), big.mark = big.mark)
    #
    if (details.methods)
      details <-
        catmeth(x,
                common, random, prediction, overall, overall.hetstat,
                #
                func.transf = x$func.transf,
                backtransf = backtransf, func.backtransf = fbt,
                #
                big.mark = big.mark, digits = digits,
                digits.tau = digits.tau,
                text.tau = text.tau, text.tau2 = text.tau2,
                #
                print.tau2 = FALSE)
  }
  else {
    #
    #
    # Print results for meta-analysis with more than one study
    #
    #
    zlab <- ""
    #
    if (header & is.metamiss)
      cat("\n")
    #
    if (all(!is.na(k))) {
      if (overall & (common | random | prediction)) {
        if (!inherits(x, "trimfill")) {
          if (any(method == "MH") &&
              (inherits(x, c("metabin", "metainc")) &
               common & sm %in% c("RD", "IRD") &
               (!is.null(x$k.MH) && any(k != x$k.MH, na.rm = TRUE))))
            cat(paste0("Number of studies:   k.MH = ",
                       cond(x$k.MH, digits = 0),
                       " (", text.common.br[1], "), k = ",
                       format(cond(k, digits = 0),
                              big.mark = big.mark),
                       " (", text.random.br[1], ")\n",
                       collapse = ""))
          else {
            if (any(k.study != k, na.rm = TRUE)) {
              cat(paste0("Number of ",
                         if (three.level) "clusters" else "studies",
                         ": n = ",
                         format(cond(x$k.study, digits = 0),
                                big.mark = big.mark), "\n",
                         collapse = ""))
              cat(paste0("Number of estimates: k = ",
                         format(cond(k, digits = 0),
                                big.mark = big.mark), "\n",
                         collapse = ""))
            }
            else
              cat(paste0("Number of studies: k = ",
                         format(cond(k, digits = 0),
                                big.mark = big.mark), "\n",
                         collapse = ""))
          }
        }
        else
          cat(paste0("Number of studies: k = ",
                     format(cond(k, digits = 0),
                            big.mark = big.mark),
                     " (with ",
                     format(cond(x$k0, digits = 0),
                            big.mark = big.mark),
                     " added studies)\n",
                     collapse = ""))
        #
        if (is.metamiss)
          cat("\n")
        else if (sel.n)
          catobsev(if (is.null(x$exclude)) x$n else x$n[!x$exclude],
                   type = "n")
        else
          catobsev(if (is.null(x$exclude)) x$n.e else x$n.e[!x$exclude],
                   if (is.null(x$exclude)) x$n.c else x$n.c[!x$exclude],
                   type = "n")
        #
        if (sel.ev)
          catobsev(if (is.null(x$exclude)) x$event else x$event[!x$exclude],
                   type = "e", addrow = TRUE)
        else if (!is.metamiss)
          catobsev(if (is.null(x$exclude)) x$event.e
                   else x$event.e[!x$exclude],
                   if (is.null(x$exclude)) x$event.c
                   else x$event.c[!x$exclude],
                   type = "e", addrow = TRUE)
        #
        res <- cbind(formatN(c(if (common) TE.common,
                               if (random) TE.random,
                               if (prediction)
                                 rep(NA, length(lowTE.predict))),
                             digits, "NA",
                             big.mark = big.mark),
                     formatCI(formatN(c(if (common) lowTE.common,
                                        if (random) lowTE.random,
                                        if (prediction) lowTE.predict),
                                      digits, "NA", big.mark = big.mark),
                              formatN(c(if (common) uppTE.common,
                                        if (random) uppTE.random,
                                        if (prediction) uppTE.predict),
                                      digits, "NA", big.mark = big.mark)),
                     if (null.given)
                       formatN(c(if (common) sTE.common,
                                 if (random) sTE.random,
                                 if (prediction)
                                   rep(NA, length(lowTE.predict))),
                               digits = digits.stat, big.mark = big.mark),
                     if (null.given)
                       formatPT(c(if (common) pTE.common,
                                  if (random) pTE.random,
                                  if (prediction)
                                    rep(NA, length(lowTE.predict))),
                                digits = digits.pval,
                                scientific = scientific.pval,
                                zero = zero.pval, JAMA = JAMA.pval))
        if (prediction) {
          seq <- (nrow(res) - length(lowTE.predict) + 1):nrow(res)
          res[seq, 1] <- ""
          if (null.given)
            res[seq, 3:4] <- ""
        }
        #
        methci.t <- method.random.ci %in% c("HK", "KR", "CR0", "CR1", "CR2")
        #
        if (any(methci.t)) {
          if ((common & random) | any(!methci.t))
            zlab <- "z|t"
          else if (common & !random)
            zlab <- "z"
          else if (!common & random)
            zlab <- "t"
        }
        else
          zlab <- "z"
        #
        if (prediction) {
          if (x$level.ma == x$level.predict)
            lab.predict <- text.predict
          else
            lab.predict <- paste0(text.predict, " (",
                                  percent_label(x$level.predict, "PI"), ")")
        }
        #
        dimnames(res) <- list(c(if (common) text.common,
                                if (random) text.random,
                                if (prediction) lab.predict),
                              c(sm.lab, ci.lab,
                                if (null.given) zlab,
                                if (null.given) "p-value"))
        #
        sort.overall <- setsort(sort.overall, nrow(res), "overall results")
        #
        if (prediction & !(common | random))
          res <- res[, ci.lab, drop = FALSE]
        #
        do.call(prmatrix,
                c(list(res[sort.overall, , drop = FALSE],
                       quote = FALSE, right = TRUE), args))
        #
        if (metabin && print.CMH) {
          Qdata <- cbind(formatN(round(Q.CMH, digits.Q), digits.Q, "NA",
                                 big.mark = big.mark),
                         df.Q.CMH,
                         formatPT(pval.Q.CMH,
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval,
                                  zero = zero.pval, JAMA = JAMA.pval))
          dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
          #
          cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
          do.call(prmatrix, c(list(Qdata, quote = FALSE, right = TRUE), args))
        }
      }
      else {
        if (any(k.study != k)) {
          cat(paste0("Number of studies: ",
                     if (is.netpairwise) "k" else "n", " = ",
                     format(cond(x$k.study, digits = 0),
                            big.mark = big.mark), "\n"))
          cat(paste0("Number of ",
                     if (is.netpairwise | pairwise)
                       "pairwise comparisons: m"
                     else
                       "estimates: k",
                     " = ",
                     format(cond(k, digits = 0),
                            big.mark = big.mark), "\n"))
        }
        else
          cat(paste0("Number of studies: k",
                     " = ",
                     format(cond(k, digits = 0),
                            big.mark = big.mark), "\n"))
        #
        if (is.metamiss)
          cat("\n")
        else if (sel.n)
          catobsev(if (is.null(x$exclude)) x$n else x$n[!x$exclude],
                   type = "n")
        else
          catobsev(if (is.null(x$exclude)) x$n.e else x$n.e[!x$exclude],
                   if (is.null(x$exclude)) x$n.c else x$n.c[!x$exclude],
                   type = "n")
        #
        if (sel.ev)
          catobsev(if (is.null(x$exclude)) x$event else x$event[!x$exclude],
                   type = "e")
        else if (!is.metamiss)
          catobsev(if (is.null(x$exclude)) x$event.e else x$event.e[!x$exclude],
                   if (is.null(x$exclude)) x$event.c else x$event.c[!x$exclude],
                   type = "e")
      }
    }
    #
    # Print information on heterogeneity
    #
    if (overall.hetstat) {
      print.tau2.ci <-
        print.tau2.ci & !all(is.na(x$lower.tau2) & is.na(x$upper.tau2))
      if (print.tau2.ci &&
          (all(x$lower.tau2 == 0, na.rm = TRUE) &
           all(x$upper.tau2 == 0, na.rm = TRUE)))
        print.tau2.ci <- FALSE
      #
      print.tau.ci <-
        print.tau.ci & !all(is.na(x$lower.tau) & is.na(x$upper.tau))
      if (print.tau.ci &&
          (all(x$lower.tau == 0, na.rm = TRUE) &
           all(x$upper.tau == 0, na.rm = TRUE)))
        print.tau.ci <- FALSE
      #
      cat(paste0("\nQuantifying heterogeneity",
                 if (!(is.null(x$level.hetstat) || is.na(x$level.hetstat)) &
                     (print.tau2.ci | print.tau.ci | print.I2.ci))
                   paste0(" (with ", percent_label(x$level.hetstat),
                          if (print.tau2.ci + print.tau.ci +
                              print.I2.ci + print.H * print.I2.ci > 1) "s",
                          ")"),
                 ":"))
      #
      if (sum(c(print.tau2, print.tau, print.I2, print.H, print.Rb)) > 1)
        cat("\n")
      else
        cat("")
      #
      cathet(k,
             x$method.tau, x$detail.tau,
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
             sort.tau, sort.het)
      #
      # Print information on residual heterogeneity
      #
      if (by) {
        #
        k.resid <- df.Q.resid + 1
        #
        if (print.H) {
          H.resid <- round(replaceNULL(x$H.resid), digits.H)
          lowH.resid <- round(replaceNULL(x$lower.H.resid), digits.H)
          uppH.resid <- round(replaceNULL(x$upper.H.resid), digits.H)
        }
        #
        I2.resid <- round(100 * replaceNULL(x$I2.resid), digits.I2)
        #
        if (print.I2) {
          lowI2.resid <- round(100 * replaceNULL(x$lower.I2.resid), digits.I2)
          uppI2.resid <- round(100 * replaceNULL(x$upper.I2.resid), digits.I2)
          print.I2.ci <-
            ((replaceNULL(Q.resid) > replaceNULL(k.resid) &
              replaceNULL(k.resid) >= 2) |
             (replaceNULL(Q.resid) <= replaceNULL(k.resid) &
              replaceNULL(k.resid) > 2)) &
            !(is.na(replaceNULL(lowI2.resid)) |
              is.na(replaceNULL(uppI2.resid)))
          #
          print.I2.ci <- unique(print.I2.ci)          
          if (is.na(print.I2.ci))
            print.I2.ci <- FALSE
        }
        #
        if (!is.na(replaceNULL(I2.resid))) {
          cat(paste0("\nQuantifying residual heterogeneity",
                     if (!(is.null(x$level.hetstat) || is.na(x$level.hetstat)) &
                         (print.tau2.ci | print.tau.ci | print.I2.ci))
                       paste0(" (with ", percent_label(x$level.hetstat),
                              if (print.tau2.ci + print.tau.ci +
                                  print.I2.ci + print.H * print.I2.ci > 1) "s",
                              ")"),
                     ":"))
          #
          if (sum(c(print.tau2, print.tau, print.I2, print.H, print.Rb)) > 1)
            cat("\n")
          else
            cat("")
          #
          cathet(k.resid, 
                 x$method.tau, x$detail.tau,
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
                 FALSE, "",
                 big.mark)
        }
      }
      #
      # Test of heterogeneity
      #
      if (print.Q & (common | random)) {
        if (any(k > 1, na.rm = TRUE)) {
          #
          Qdat <- qdat(Q, df.Q, pval.Q, x$hetlabel, text.common)
          #
          sort.Q <- setsort(sort.Q, nrow(Qdat), "tests of heterogeneity")
          #
          Qdata <- cbind(formatN(round(Qdat$Q, digits.Q), digits.Q, "NA",
                                 big.mark = big.mark),
                         format(Qdat$df.Q, big.mark = big.mark),
                         formatPT(Qdat$pval.Q,
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval,
                                  zero = zero.pval, JAMA = JAMA.pval))
          dimnames(Qdata) <- list(Qdat$names, c("Q", "d.f.", "p-value"))
          #
          Qdata <- Qdata[sort.Q, , drop = FALSE]
          #
          cat("\nTest of heterogeneity:\n")
          do.call(prmatrix, c(list(Qdata, quote = FALSE, right = TRUE), args))
        }
      }
      # Test of residual heterogeneity
      #
      if (by && print.Q && (common || random)) {
        if (length(Q.resid) == 1L &&
            length(df.Q.resid) == 1L &&
            !is.na(Q.resid) && !is.na(df.Q.resid) &&
            df.Q.resid > 0) {
          Qdata <- cbind(
            formatN(Q.resid, digits.Q, "NA", big.mark = big.mark),
            formatN(df.Q.resid, digits = 0, big.mark = big.mark),
            formatPT(pval.Q.resid,
                     digits = digits.pval.Q,
                     scientific = scientific.pval,
                     zero = zero.pval, JAMA = JAMA.pval))
          #
          dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
          #
          cat("\nTest of residual heterogeneity:\n")
          do.call(prmatrix, c(list(Qdata, quote = FALSE, right = TRUE), args))
        }
      }
    }
    else if (!by) {
      print.tau2 <- FALSE
      print.tau2.ci <- FALSE
      print.tau <- FALSE
      print.tau.ci <- FALSE
    }
    else {
      print.tau2.ci <- FALSE
      print.tau.ci <- FALSE
    }
    #
    # Print information for subgroup analysis
    #
    if (by) {
      if (is.netpairwise)
        anaunit <- "pairwise comparisons"
      else
        anaunit <- "subgroups"
      #
      if (common) {
        #
        # Subgroup analysis based on common effect model
        #
        if (is.vector(TE.common.w)) {
          nam <- names(TE.common.w)
          TE.common.w <- matrix(TE.common.w, ncol = 1)
          lowTE.common.w <- matrix(lowTE.common.w, ncol = 1)
          uppTE.common.w <- matrix(uppTE.common.w, ncol = 1)
          #
          if (!is.null(nam))
            rownames(TE.common.w) <- rownames(lowTE.common.w) <-
              rownames(uppTE.common.w) <- nam
        }
        #
        n.common <- ncol(TE.common.w)
        #
        for (i in seq_len(n.common)) {
          Tdata <- cbind(format(k.w, big.mark = big.mark),
                         formatN(TE.common.w[, i], digits, "NA",
                                 big.mark = big.mark),
                         formatCI(formatN(lowTE.common.w[, i], digits, "NA",
                                          big.mark = big.mark),
                                  formatN(uppTE.common.w[, i], digits, "NA",
                                          big.mark = big.mark)),
                         if (print.Q)
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
                         if (print.tau2 & !random & !is.list(x$tau.w))
                           ifelse(k.w == 1 & !x$tau.common, "--",
                                  formatPT(x$tau.w^2,
                                           digits = digits.tau2,
                                           big.mark = big.mark,
                                           noblanks = TRUE)),
                         if (print.tau & !random & !is.list(x$tau.w))
                           ifelse(k.w == 1 & !x$tau.common, "--",
                                  formatPT(x$tau.w,
                                           digits = digits.tau,
                                           big.mark = big.mark,
                                           noblanks = TRUE))
                         )
          #
          bylab.txt <- bylabel(subgroup.name, subgroup.levels,
                               print.subgroup.name, sep.subgroup,
                               big.mark = big.mark)
          #
          dimnames(Tdata) <- list(
            bylab.txt,
            c(if (is.netpairwise | pairwise) "  m" else "  k",
              sm.lab, ci.lab,
              if (print.Q) "Q",
              if (print.I2) text.I2,
              if (print.Rb) text.Rb,
              if (print.tau2 & !random & !is.list(x$tau.w)) text.tau2,
              if (print.tau & !random & !is.list(x$tau.w)) text.tau)
          )
          #
          cat(paste0("\nResults for ", anaunit,
                     paste0(" (", text.common.br[i], "):\n")))
          #
          do.call(prmatrix, c(list(Tdata, quote = FALSE, right = TRUE), args))
          #
          if (test.subgroup.common && i <= length(Q.b.common)) {
            df.i <- df.Q.b.common[[i]]
            dfs.b.i <- format_df(df.i, lab.NA = "NA", big.mark = big.mark)
            Q.lab <- if (length(df.i) == 2L) "F" else "Q"
            #
            cat(paste0("\nTest for subgroup differences (",
                       text.common.br[i], "):\n"))
            #
            if (any(method == "MH") | all(is.na(Q.w.common))) {
              Qdata <- cbind(formatN(round(Q.b.common[i], digits.Q),
                                     digits.Q, "NA", big.mark = big.mark),
                             formatN(dfs.b.i, digits = 0,
                                     big.mark = big.mark),
                             formatPT(pval.Q.b.common[i],
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval,
                                      zero = zero.pval, JAMA = JAMA.pval))
              dimnames(Qdata) <- list("Between groups",
                                      c(Q.lab, "d.f.", "p-value"))
              do.call(prmatrix, c(list(Qdata, quote = FALSE, right = TRUE), args))
            }
            else {
              Qs  <- c(Q.b.common[i], Q.w.common[i])
              dfs <- c(dfs.b.i,
                       formatN(df.Q.w, digits = 0, big.mark = big.mark))
              if (is.glmm & Q.lab == "F")
                Q.lab <- "F/Q"
              #
              pvals <- c(pval.Q.b.common[i], pval.Q.w.common[i])
              #
              Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             dfs,
                             formatPT(pvals,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval,
                                      zero = zero.pval, JAMA = JAMA.pval))
              dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                      c(Q.lab, "d.f.", "p-value"))
              do.call(prmatrix, c(list(Qdata, quote = FALSE, right = TRUE), args))
            }
          }
        }
      }
      #
      if (random) {
        #
        # Subgroup analysis based on random effects model
        #
        if (is.vector(TE.random.w)) {
          nam <- names(TE.random.w)
          TE.random.w <- matrix(TE.random.w, ncol = 1)
          lowTE.random.w <- matrix(lowTE.random.w, ncol = 1)
          uppTE.random.w <- matrix(uppTE.random.w, ncol = 1)
          #
          if (!is.null(nam))
            rownames(TE.random.w) <- rownames(lowTE.random.w) <-
              rownames(uppTE.random.w) <- nam
        }
        #
        n.random <- ncol(TE.random.w)
        #
        if (!is.list(x$tau.w))
          x$tau.w <- list(x$tau.w)
        #
        for (i in seq_len(n.random)) {
          tau.i <- length(x$tau.w) >= i
          #
          Tdata <- cbind(format(k.w, big.mark = big.mark),
                         formatN(TE.random.w[, i], digits, "NA",
                                 big.mark = big.mark),
                         formatCI(formatN(lowTE.random.w[, i], digits, "NA",
                                          big.mark = big.mark),
                                  formatN(uppTE.random.w[, i], digits, "NA",
                                          big.mark = big.mark)),
                         if (print.tau2 & (i == 1 | tau.i))
                           ifelse(k.w == 1 & !x$tau.common, "--",
                                  formatPT(x$tau.w[[i]]^2,
                                           digits = digits.tau2,
                                           big.mark = big.mark,
                                           noblanks = TRUE)),
                         if (print.tau & (i == 1 | tau.i))
                           ifelse(k.w == 1 & !x$tau.common, "--",
                                  formatPT(x$tau.w[[i]],
                                           digits = digits.tau,
                                           big.mark = big.mark,
                                           noblanks = TRUE)),
                         if (i == 1 & !common & print.Q)
                           formatN(round(Q.w, digits.Q), digits.Q,
                                   big.mark = big.mark),
                         if (i == 1 & !common & print.I2)
                           ifelse(is.na(I2.w),
                                  "--",
                                  paste0(formatN(I2.w, digits.I2), "%")),
                         if (i == 1 & !common & print.Rb)
                           ifelse(is.na(Rb.w),
                                  "--",
                                  paste0(formatN(Rb.w, digits.I2,
                                                 big.mark = big.mark), "%"))
                         )
          #
          bylab.txt <- bylabel(subgroup.name, subgroup.levels,
                               print.subgroup.name, sep.subgroup,
                               big.mark = big.mark)
          #
          dimnames(Tdata) <- list(bylab.txt,
                                  c(if (is.netpairwise | pairwise)
                                      "  m" else "  k",
                                    sm.lab, ci.lab,
                                    if (print.tau2 & (i == 1 | tau.i))
                                      text.tau2,
                                    if (print.tau & (i == 1 | tau.i))
                                      text.tau,
                                    if (i == 1 & !common & print.Q) "Q",
                                    if (i == 1 & !common & print.I2) text.I2,
                                    if (i == 1 & !common & print.Rb) text.Rb)
                                  )
          #
          cat(paste0("\nResults for ", anaunit,
                     paste0(" (", text.random.br[i], "):\n")))
          #
          do.call(prmatrix, c(list(Tdata, quote = FALSE, right = TRUE), args))
          #
          if (test.subgroup.random && i <= length(Q.b.random) && !is.na(Q.b.random[i])) {
            df.i <- df.Q.b.random[[i]]
            dfs.b.i <- format_df(df.i, lab.NA = "NA", big.mark = big.mark)
            Q.lab <- if (length(df.i) == 2L) "F" else "Q"
            #
            cat(paste0("\nTest for subgroup differences (",
                       text.random.br[i],
                       "):\n"))
            #
            Qdata <- cbind(
              formatN(round(Q.b.random[i], digits.Q),
                      digits.Q, "NA", big.mark = big.mark),
              dfs.b.i,
              formatPT(pval.Q.b.random[i],
                       digits = digits.pval.Q,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval))
            #
            dimnames(Qdata) <- list(
              "Between groups",
              c(Q.lab, "d.f.", "p-value"))
            #
            do.call(prmatrix, c(list(Qdata, quote = FALSE, right = TRUE), args))
          }
        }
      }
      #
      if (prediction.w) {
        #
        # Prediction intervals for subgroups
        #
        if (is.vector(lowTE.predict.w)) {
          nam <- names(lowTE.predict.w)
          lowTE.predict.w <- matrix(lowTE.predict.w, ncol = 1)
          uppTE.predict.w <- matrix(uppTE.predict.w, ncol = 1)
          #
          rownames(lowTE.predict.w) <- rownames(uppTE.predict.w) <- nam
        }
        #
        for (i in seq_len(ncol(lowTE.predict.w))) {
          Pdata <- cbind(formatCI(formatN(lowTE.predict.w[, i], digits, "NA",
                                          big.mark = big.mark),
                                  formatN(uppTE.predict.w[, i], digits, "NA",
                                          big.mark = big.mark)))
          #
          bylab.txt <-
            bylabel(subgroup.name, subgroup.levels,
                    print.subgroup.name, sep.subgroup,
                    big.mark = big.mark)
          lab.predict <- percent_label(x$level.predict, "PI")
          dimnames(Pdata) <- list(bylab.txt, lab.predict)
          #
          Pdata <- Pdata[prediction.subgroup.logical, , drop = FALSE]
          #
          cat(paste0("\n", text.predict[i], " for subgroups:\n"))
          do.call(prmatrix, c(list(Pdata, quote = FALSE, right = TRUE), args))
        }
      }
    }
    #
    # Print information on meta-analysis method:
    #
    x$pscale <- pscale
    x$irscale <- irscale
    x$irunit <- irunit
    #
    if (!is.null(x$cycles))
      x$sd.n_of_1 <-
      formatPT(x$sd.n_of_1, digits = gs("digits.sd"), big.mark = big.mark)
    #
    if (details.methods &
        (common | random | prediction | overall.hetstat | by))
      details <-
        catmeth(x,
                common, random, prediction, overall, overall.hetstat,
                #
                func.transf = x$func.transf,
                backtransf = backtransf, func.backtransf = fbt,
                #
                big.mark = big.mark, digits = digits,
                digits.tau = digits.tau,
                text.tau = text.tau, text.tau2 = text.tau2,
                #
                print.tau2 = print.tau2,
                print.tau2.ci = print.tau2 & print.tau2.ci,
                print.tau = print.tau,
                print.tau.ci = print.tau & print.tau.ci,
                #
                print.I2 = print.I2, text.I2 = text.I2,
                #
                print.df = overall, prediction.subgroup = prediction.subgroup)
  }
  
  invisible(details)
}





#' @rdname print.meta
#' @export cilayout

cilayout <- function(bracket = gs("CIbracket"),
                     separator = gs("CIseparator"),
                     lower.blank = gs("CIlower.blank"),
                     upper.blank = gs("CIupper.blank")) {

  chkchar(bracket, length = 1)
  chkchar(separator, length = 1)
  #
  bracket <- setchar(bracket, c("[", "(", "{", ""))
  #
  chklogical(lower.blank)
  chklogical(upper.blank)
  
  setOption("CIbracket", bracket)
  setOption("CIseparator", separator)
  setOption("CIlower.blank", lower.blank)
  setOption("CIupper.blank", upper.blank)
  
  invisible(NULL)
}
