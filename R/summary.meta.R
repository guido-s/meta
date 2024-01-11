#' Summary of meta-analysis results
#' 
#' @description
#' Summary method for objects of class \code{meta}.
#' 
#' @param object An object of class \code{meta}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Summary method for objects of class \code{meta}.
#' 
#' @return
#' An object of classes \code{summary.meta} and \code{meta} (see
#' \code{\link{meta-object}}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#'   data = Fleiss1993cont, studlab = paste(study, year), sm = "SMD")
#' summary(m1)
#' 
#' summary(update(m1, subgroup = c(1, 2, 1, 1, 2), subgroup.name = "group"))
#' forest(update(m1, subgroup = c(1, 2, 1, 1, 2), subgroup.name = "group"))
#' 
#' \dontrun{
#' # Use unicode characters to print tau^2, tau, and I^2
#' print(summary(m1),
#'   text.tau2 = "\u03c4\u00b2", text.tau = "\u03c4", text.I2 = "I\u00b2")
#' }
#' 
#' @method summary meta
#' @export


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
  object <- updateversion(object)
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
  object$df <- replaceNULL(object$df, Inf)
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
  ci.c <- list(TE = object$TE.common,
               seTE = object$seTE.common,
               lower = object$lower.common,
               upper = object$upper.common,
               statistic = object$statistic.common,
               p = object$pval.common,
               level = object$level.ma)
  if (metaprop)
    ci.c$harmonic.mean <- 1 / mean(1 / object$n)
  else if (metarate)
    ci.c$harmonic.mean <- 1 / mean(1 / object$time)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               statistic = object$statistic.random,
               p = object$pval.random,
               level = object$level.ma,
               df = object$df.random)
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
               df = object$df.predict)
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- object
  ##
  res$common <- ci.c
  res$random <- ci.r
  res$predict <- ci.p
  ##
  ## Backward compatibility
  ##
  res$fixed <- ci.c
  ##
  ## Add results from subgroup analysis
  ##
  if (length(object$subgroup) > 0) {
    ##
    n.subgroups <- length(object$subgroup.levels)
    ##
    ci.common.w <- list(TE = object$TE.common.w,
                        seTE = object$seTE.common.w,
                        lower = object$lower.common.w,
                        upper = object$upper.common.w,
                        statistic = object$statistic.common.w,
                        p = object$pval.common.w,
                        level = object$level.ma,
                        harmonic.mean = object$n.harmonic.mean.w)
    ##
    if (metarate)
      ci.common.w$harmonic.mean <- object$t.harmonic.mean.w
    ##
    ci.random.w <- list(TE = object$TE.random.w,
                        seTE = object$seTE.random.w,
                        lower = object$lower.random.w,
                        upper = object$upper.random.w,
                        statistic = object$statistic.random.w,
                        p = object$pval.random.w,
                        level = object$level.ma,
                        df = object$df.random.w,
                        harmonic.mean = object$n.harmonic.mean.w)
    ##
    ci.predict.w <- list(TE = rep(NA, n.subgroups),
                         seTE = object$seTE.predict.w,
                         lower = object$lower.predict.w,
                         upper = object$upper.predict.w,
                         statistic = rep(NA, n.subgroups),
                         p = rep(NA, n.subgroups),
                         level = object$level.predict,
                         df = object$df.predict.w,
                         harmonic.mean = object$n.harmonic.mean.w)
    ##
    if (metarate)
      ci.random.w$harmonic.mean <- object$t.harmonic.mean.w
    ## 
    res$within.common  <- ci.common.w
    res$within.random  <- ci.random.w
    res$within.predict <- ci.predict.w
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
