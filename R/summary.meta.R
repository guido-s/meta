#' Summary of meta-analysis results
#' 
#' @description
#' Summary method for objects of class \code{meta}.
#' 
#' @param object An object of class \code{meta}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Review Manager 5 (RevMan 5) is the current software used for
#' preparing and maintaining Cochrane Reviews
#' (\url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman}).
#' In RevMan 5, subgroup analyses can be defined and data from a
#' Cochrane review can be imported to Rusing the function
#' \code{read.rm5}. If a meta-analysis is then conducted using
#' function \code{metacr}, information on subgroups is available in R
#' (components \code{subgroup}, \code{subgroup.name}, and
#' \code{print.subgroup.name}, \code{subgroup} in an object of class
#' \code{"meta"}).  Accordingly, by using function \code{metacr} there
#' is no need to define subgroups in order to redo the statistical
#' analysis conducted in the Cochrane review.
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
#' @seealso \code{\link{print.summary.meta}}, \code{\link{metabin}},
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
#' summary(update(m1, subgroup = c(1, 2, 1, 1, 2), subgroup.name = "group"))
#' forest(update(m1, subgroup = c(1, 2, 1, 1, 2), subgroup.name = "group"))
#' 
#' \dontrun{
#' # Use unicode characters to print tau^2, tau, and I^2
#' print(summary(m1),
#'       text.tau2 = "\u03c4\u00b2", text.tau = "\u03c4", text.I2 = "I\u00b2")
#' }
#' 
#' @method summary meta
#' @export
#' @export summary.meta


summary.meta <- function(object, ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(object, "meta")
  ##
  if (inherits(object, c("metacum", "metainf")))
    return(object)
  ##
  metaprop <- inherits(class(object), "metaprop")
  metarate <- inherits(class(object), "metarate")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  addargs <- names(list(...))
  ##
  if (length(addargs) > 0)
    warning("Additional arguments provided in '...' are ignored.",
            call. = FALSE)
  
  
  ##
  ##
  ## (3) Results for individual studies
  ##
  ##
  object$df <- replaceNULL(object$df, NA)
  method.ci <- replaceNULL(object$method.ci, "")
  object$statistic <- replaceNULL(object$statistic, object$zval)
  ##
  ci.study <- list(TE = object$TE,
                   seTE = object$seTE,
                   lower = object$lower,
                   upper = object$upper,
                   statistic = object$statistic,
                   p = object$pval,
                   level = object$level,
                   df = object$df)
  ##
  if (metaprop) {
    ci.study$event <- object$event
    ci.study$n <- object$n
  }
  
  
  ##
  ##
  ## (4) Results for meta-analysis
  ##
  ##
  ci.f <- list(TE = object$TE.fixed,
               seTE = object$seTE.fixed,
               lower = object$lower.fixed,
               upper = object$upper.fixed,
               statistic = object$statistic.fixed,
               p = object$pval.fixed,
               level = object$level.ma)
  if (metaprop)
    ci.f$harmonic.mean <- 1 / mean(1 / object$n)
  else if (metarate)
    ci.f$harmonic.mean <- 1 / mean(1 / object$time)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               statistic = object$statistic.random,
               p = object$pval.random,
               level = object$level.ma,
               df = if (!is.null(object$df.hakn)) object$df.hakn else NA)
  if (metaprop)
    ci.r$harmonic.mean <- 1 / mean(1 / object$n)
  else if (metarate)
    ci.r$harmonic.mean <- 1 / mean(1 / object$time)
  ##
  ci.p <- list(TE = NA,
               seTE = object$seTE.predict,
               lower = object$lower.predict,
               upper = object$upper.predict,
               statistic = NA,
               p = NA,
               level = object$level.predict,
               df = object$k - 2)
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- object
  ##
  res$fixed <- ci.f
  res$random <- ci.r
  res$predict <- ci.p
  ##
  ## Add results from subgroup analysis
  ##
  if (length(object$subgroup) > 0) {
    ##
    n.by <- length(object$bylevs)
    ##
    ci.fixed.w <- list(TE = object$TE.fixed.w,
                       seTE = object$seTE.fixed.w,
                       lower = object$lower.fixed.w,
                       upper = object$upper.fixed.w,
                       statistic = object$statistic.fixed.w,
                       p = object$pval.fixed.w,
                       level = object$level.ma,
                       harmonic.mean = object$n.harmonic.mean.w)
    ##
    if (metarate)
      ci.fixed.w$harmonic.mean <- object$t.harmonic.mean.w
    ##
    ci.random.w <- list(TE = object$TE.random.w,
                        seTE = object$seTE.random.w,
                        lower = object$lower.random.w,
                        upper = object$upper.random.w,
                        statistic = object$statistic.random.w,
                        p = object$pval.random.w,
                        level = object$level.ma,
                        df = object$df.hakn.w,
                        harmonic.mean = object$n.harmonic.mean.w)
    ##
    ci.predict.w <- list(TE = rep(NA, n.by),
                         seTE = object$seTE.predict.w,
                         lower = object$lower.predict.w,
                         upper = object$upper.predict.w,
                         statistic = rep(NA, n.by),
                         p = rep(NA, n.by),
                         level = object$level.predict,
                         df = object$k.w - 2,
                         harmonic.mean = object$n.harmonic.mean.w)
    ##
    if (metarate)
      ci.random.w$harmonic.mean <- object$t.harmonic.mean.w
    ## 
    res$within.fixed    <- ci.fixed.w
    res$within.random   <- ci.random.w
    res$within.predict  <- ci.predict.w
    ##
    if (is.null(res$test.subgroup))
      res$test.subgroup <- TRUE
  }
  ##
  res$x <- object
  ##
  class(res) <- c("summary.meta", class(object))
  ##
  attr(res, "class.orig") <- class(object)
  
  
  res
}
