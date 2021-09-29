#' Summary of meta-analysis results
#' 
#' @description
#' Summary method for objects of class \code{meta}.
#' 
#' @param object An object of class \code{meta}.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
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
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param bystud A logical indicating whether results of individual
#'   studies should be printed by grouping variable.
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Note, in R package \bold{meta}, version 3.0-0 some arguments have
#' been removed from R functions \code{\link{summary.meta}}
#' (arguments: byvar, level, level.comb, level.prediction) and
#' print.summary.meta (arguments: level, level.comb,
#' level.prediction). This functionality is now provided by R function
#' \code{\link{update.meta}} (or directly in meta-analysis functions,
#' e.g., \code{\link{metabin}}, \code{\link{metacont}},
#' \code{\link{metagen}}, \code{\link{metacor}}, and
#' \code{\link{metaprop}}).
#' 
#' Review Manager 5 (RevMan 5) is the current software used for
#' preparing and maintaining Cochrane Reviews
#' (\url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman}).
#' In RevMan 5, subgroup analyses can be defined and data from a
#' Cochrane review can be imported to Rusing the function \code{read.rm5}. If a
#' meta-analysis is then conducted using function \code{metacr}, information on
#' subgroups is available in R (components \code{byvar}, \code{bylab}, and
#' \code{print.byvar}, \code{byvar} in an object of class \code{"meta"}).
#' Accordingly, by using function \code{metacr} there is no need to define
#' subgroups in order to redo the statistical analysis conducted in the
#' Cochrane review.
#' 
#' Note, for an object of type \code{metaprop}, starting with version
#' 3.7-0 of meta, list elements \code{TE}, \code{lower} and
#' \code{upper} in element \code{study} correspond to transformed
#' proportions and confidence limits (regardless whether exact
#' confidence limits are calculated; argument \code{ciexact=TRUE} in
#' metaprop function). Accordingly, the following results are based on
#' the same transformation defined by argument \code{sm}: list
#' elements \code{TE}, \code{lower} and \code{upper} in elements
#' \code{study}, \code{fixed}, \code{random}, \code{within.fixed} and
#' \code{within.random}.
#' 
#' R function cilayout can be utilised to change the layout to print
#' confidence intervals (both in printout from print.meta and
#' print.summary.meta function as well as in forest plots). The
#' default layout is "[lower; upper]". Another popular layout is
#' "(lower - upper)" which is used throughout an R session by using R
#' command \code{cilayout("(", " - ")}.
#' 
#' Argument \code{pscale} can be used to rescale single proportions or
#' risk differences, e.g. \code{pscale=1000} means that proportions
#' are expressed as events per 1000 observations. This is useful in
#' situations with (very) low event probabilities.
#' 
#' Argument \code{irscale} can be used to rescale single rates or rate
#' differences, e.g. \code{irscale=1000} means that rates are
#' expressed as events per 1000 time units, e.g. person-years. This is
#' useful in situations with (very) low rates. Argument \code{irunit}
#' can be used to specify the time unit used in individual studies
#' (default: "person-years"). This information is printed in summaries
#' and forest plots if argument \code{irscale} is not equal to 1.
#' 
#' @return
#' An object of classes \code{summary.meta} and \code{meta}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{update.meta}}, \code{\link{metabin}},
#'   \code{\link{metacont}}, \code{\link{metagen}}
#' 
#' @references
#' Cooper H & Hedges LV (1994):
#' \emph{The Handbook of Research Synthesis}.
#' Newbury Park, CA: Russell Sage Foundation
#' 
#' Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D (2016):
#' A new measure of between-studies heterogeneity in meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 3661--75
#' 
#' Higgins JPT & Thompson SG (2002):
#' Quantifying heterogeneity in a meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{21}, 1539--58
#' 
#' @examples
#' data(Fleiss1993cont)
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                data = Fleiss1993cont, sm = "SMD",
#'                studlab = paste(study, year))
#' summary(m1)
#' 
#' summary(update(m1, byvar = c(1, 2, 1, 1, 2), bylab = "group"))
#' forest(update(m1, byvar = c(1, 2, 1, 1, 2), bylab = "group"))
#' 
#' \dontrun{
#' # Use unicode characters to print tau^2, tau, and I^2
#' print(summary(m1),
#'       text.tau2 = "\u03c4\u00b2", text.tau = "\u03c4", text.I2 = "I\u00b2")
#' }
#' 
#' @rdname summary.meta
#' @export
#' @export summary.meta


summary.meta <- function(object,
                         comb.fixed = object$comb.fixed,
                         comb.random = object$comb.random,
                         prediction = object$prediction,
                         overall = object$overall,
                         overall.hetstat = object$overall.hetstat,
                         test.subgroup = object$test.subgroup,
                         ##
                         backtransf = object$backtransf,
                         pscale = object$pscale,
                         irscale = object$irscale,
                         irunit = object$irunit,
                         ##
                         bylab = object$bylab,
                         print.byvar = object$print.byvar,
                         byseparator = object$byseparator,
                         bystud = FALSE,
                         ##
                         print.CMH = object$print.CMH,
                         ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(object, "meta")
  ##
  if (inherits(object, "metacum")) {
    warning("Summary method not defined for objects of class \"metacum\".")
    return(object)
  }
  ##
  if (inherits(object, "metainf")) {
    warning("Summary method not defined for objects of class \"metainf\".")
    return(object)
  }
  ##
  object <- updateversion(object)
  ##
  metaprop <- inherits(object, "metaprop")
  metarate <- inherits(object, "metarate")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  overall <- replaceNULL(overall, comb.fixed | comb.random)
  chklogical(overall)
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  test.subgroup <- replaceNULL(test.subgroup, TRUE)
  if (is.na(test.subgroup))
    test.subgroup <- FALSE
  chklogical(test.subgroup)
  ##
  chklogical(backtransf)
  ##
  chknumeric(pscale, length = 1)
  chknumeric(irscale, length = 1)
  ##
  if (!backtransf & pscale != 1 & !is.untransformed(object$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  if (!backtransf & irscale != 1 & !is.untransformed(object$sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  if (!is.null(print.byvar))
    chklogical(print.byvar)
  if (!is.null(byseparator))
    chkchar(byseparator)
  chklogical(bystud)
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  ##
  cl <- paste0("update.meta() or ", class(object)[1], "()")
  addargs <- names(list(...))
  ##
  fun <- "summary.meta"
  ##
  warnarg("byvar", addargs, fun, cl)
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  
  
  ##
  ##
  ## (3) Generate R object
  ##
  ##
  res <- object
  ##
  comb.fixed <- comb.fixed
  res$comb.random <- comb.random
  res$prediction <- prediction
  res$overall <- overall
  res$overall.hetstat <- overall.hetstat
  res$test.subgroup <- test.subgroup
  ##
  res$backtransf <- backtransf
  res$pscale <- pscale
  res$irscale <- irscale
  res$irunit <- irunit
  ##
  res$bylab <- bylab
  res$print.byvar <- print.byvar
  res$byseparator <- byseparator
  res$bystud <- bystud
  ##
  res$print.CMH <- print.CMH
  ##
  res$call.summary <- match.call()
  ##
  class(res) <- c("summary.meta", class(object))
  
  
  ##if (metaprop)
  ##  ci.f$harmonic.mean <- 1 / mean(1 / object$n)
  ##else if (metarate)
  ##  ci.f$harmonic.mean <- 1 / mean(1 / object$time)
  
  
  ##if (inherits(object, "metamiss")) {
  ##  res$n.e <- object$n.e + object$miss.e
  ##  res$n.c <- object$n.c + object$miss.c
  ##}
  ##
  ## Function netpairwise() from R package netmeta
  ##
  if (inherits(object, "netpairwise"))
    class(res) <- c(class(res), "is.netpairwise")
  ##
  attr(res, "class.orig") <- class(object)
  
  
  res
}
