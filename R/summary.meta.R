#' Summary of meta-analysis results
#' 
#' @description
#' Summary method for objects of class \code{meta}.
#' 
#' @param x An object of class \code{summary.meta}.
#' @param object An object of class \code{meta}.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
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
#' @param bylab.nchar A numeric specifying the number of characters to
#'   print from label for the grouping variable.
#' @param bystud A logical indicating whether results of individual
#'   studies should be printed by grouping variable.
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.zval Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.H Minimal number of significant digits for H
#'   statistic, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param print.I2 A logical specifying whether heterogeneity
#'   statistic I^2 should be printed.
#' @param warn A logical indicating whether the use of
#'   \code{summary.meta} in connection with \code{metacum} or
#'   \code{metainf} should result in a warning.
#' @param warn.backtransf A logical indicating whether a warning
#'   should be printed if backtransformed proportions and rates are
#'   below 0 and backtransformed proportions are above 1.
#' @param print.H A logical specifying whether heterogeneity statistic
#'   H should be printed.
#' @param print.Rb A logical specifying whether heterogeneity
#'   statistic Rb should be printed.
#' @param text.tau2 Text printed to identify between-study variance
#'   tau^2.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I^2.
#' @param text.Rb Text printed to identify heterogeneity statistic Rb.
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
#' (\url{http://community.cochrane.org/tools/review-production-tools/revman-5}).
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
#' A list is returned by the function \code{summary.meta} with the
#' following elements:
#' \item{study}{Results for individual studies (a list with elements
#'   TE, seTE, lower, upper, z, p, level, df).}
#' \item{fixed}{Results for fixed effect model (a list with elements
#'   TE, seTE, lower, upper, z, p, level, df).}  #
#' \item{random}{Results for random effects model (a list with
#'   elements TE, seTE, lower, upper, z, p, level, df).}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic Q.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{se.tau2}{Standard error of between-study variance.}
#' \item{C}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups.}
#' \item{H}{Heterogeneity statistic H (a list with elements TE, lower,
#'   upper).}
#' \item{I2}{Heterogeneity statistic I2 (a list with elements TE,
#'   lower, upper), see Higgins & Thompson (2002).}
#' \item{Rb}{Heterogeneity statistic Rb (a list with elements TE,
#'   lower, upper), see Crippa et al. (2016).}  # \item{k.all}{Total
#'   number of studies.}
#' \item{Q.CMH}{Cochran-Mantel-Haenszel test statistic for overall
#'   effect.}
#' \item{sm}{A character string indicating underlying summary
#'   measure.}
#' \item{method}{A character string with the pooling method.}
#' \item{call}{Function call.}
#' \item{ci.lab}{Label for confidence interval.}
#' \item{hakn}{A logical indicating whether method by Hartung and
#'   Knapp was used.}
#' 
#' \item{method.tau}{A character string indicating which method is
#'   used to estimate the between-study variance tau-squared.}
#' 
#' \item{tau.common}{A logical indicating whether tau-squared is
#'   assumed to be the same across subgroups.}
#' 
#' \item{within.fixed}{Result for fixed effect model within groups (a
#'   list with elements TE, seTE, lower, upper, z, p, level, df,
#'   harmonic.mean) - if \code{byvar} is not missing.}
#' 
#' \item{within.random}{Result for random effects model within groups
#'   (a list with elements TE, seTE, lower, upper, z, p, level, df,
#'   harmonic.mean) - if \code{byvar} is not missing.}
#' 
#' \item{k.w}{Number of studies combined within groups - if
#'   \code{byvar} is not missing.}
#' \item{Q.w}{Heterogeneity statistic Q within groups - if
#'   \code{byvar} is not missing.}
#' \item{Q.b.fixed}{Heterogeneity statistic Q between groups (based on
#'   fixed effect model) - if \code{byvar} is not missing.}
#' \item{Q.b.random}{Heterogeneity statistic Q between groups (based
#'   on random effects model) - if \code{byvar} is not missing.}
#' \item{tau.w}{Square-root of between-study variance within subgroups
#'   - if \code{byvar} is not missing.}
#' \item{C.w}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups.}
#' \item{H.w}{Heterogeneity statistic H within subgroups (a list with
#'   elements TE, lower, upper) - if \code{byvar} is not missing.}
#' \item{I2.w}{Heterogeneity statistic I2 within subgroups (a list
#'   with elements TE, lower, upper) - if \code{byvar} is not
#'   missing.}
#' \item{Rb.w}{Heterogeneity statistic Rb within subgroups (a list
#'   with elements TE, lower, upper) - if \code{byvar} is not
#'   missing.}
#' \item{H.resid}{Statistic H for residual heterogeneity (a list with
#'   elements TE, lower, upper) - if \code{byvar} is not missing.}
#' \item{I2.resid}{Statistic I2 for residual heterogeneity (a list
#'   with elements TE, lower, upper) - if \code{byvar} is not
#'   missing.}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{complab}{Comparison label.}
#' \item{outclab}{Outcome label.}
#' \item{data}{Original data (set) used to create meta object.}
#' \item{subset}{Information on subset of original data used in
#'   meta-analysis.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{comb.fixed, comb.random, print.CMH}{As defined above.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
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
#' data(Fleiss93cont)
#' m1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
#'                data = Fleiss93cont, sm = "SMD",
#'                studlab = paste(study, year))
#' summary(m1)
#' 
#' summary(update(m1, byvar = c(1, 2, 1, 1, 2), bylab = "group"))
#' forest(update(m1, byvar = c(1, 2, 1, 1, 2), bylab = "group"))
#' 
#' \dontrun{
#' # Use unicode characters to print tau^2 and I^2
#' print(summary(m1), text.tau2 = "\u03c4\u00b2", text.I2 = "I\u00b2")
#' }
#' 
#' @rdname summary.meta
#' @export
#' @export summary.meta


summary.meta <- function(object,
                         comb.fixed = object$comb.fixed,
                         comb.random = object$comb.random,
                         prediction = object$prediction,
                         backtransf = object$backtransf,
                         pscale = object$pscale,
                         irscale = object$irscale,
                         irunit = object$irunit,
                         bylab = object$bylab,
                         print.byvar = object$print.byvar,
                         byseparator = object$byseparator,
                         bystud = FALSE,
                         print.CMH = object$print.CMH,
                         warn = object$warn,
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
  if (length(warn) == 0)
    warn <- gs("warn")
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
  ##
  chklogical(backtransf)
  ##
  chknumeric(pscale, single = TRUE)
  chknumeric(irscale, single = TRUE)
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
  chklogical(warn)
  ##
  cl <- paste("update.meta() or ", class(object)[1], "()", sep = "")
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
  ## (3) Results for individual studies
  ##
  ##
  ci.study <- list(TE = object$TE,
                   seTE = object$seTE,
                   lower = object$lower,
                   upper = object$upper,
                   z = object$zval,
                   p = object$pval,
                   level = object$level,
                   df = NA)
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
               z = object$zval.fixed,
               p = object$pval.fixed,
               level = object$level.comb)
  if (metaprop)
    ci.f$harmonic.mean <- mean(1 / object$n)
  else if (metarate)
    ci.f$harmonic.mean <- mean(1 / object$time)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               z = object$zval.random,
               p = object$pval.random,
               level = object$level.comb,
               df = if (!is.null(object$df.hakn)) object$df.hakn else NA)
  if (metaprop)
    ci.r$harmonic.mean <- mean(1 / object$n)
  else if (metarate)
    ci.r$harmonic.mean <- mean(1 / object$time)
  ##
  ci.H <- list(TE = object$H, lower = object$lower.H, upper = object$upper.H)
  ##
  ci.I2 <- list(TE = object$I2, lower = object$lower.I2, upper = object$upper.I2)
  ##
  ci.Rb <- list(TE = object$Rb, lower = object$lower.Rb, upper = object$upper.Rb)
  ##
  ci.H.resid <- list(TE = object$H.resid,
                     lower = object$lower.H.resid,
                     upper = object$upper.H.resid)
  ##
  ci.I2.resid <- list(TE = object$I2.resid,
                      lower = object$lower.I2.resid,
                      upper = object$upper.I2.resid)
  ##
  ci.p <- list(TE = NA,
               seTE = object$seTE.predict,
               lower = object$lower.predict,
               upper = object$upper.predict,
               z = NA,
               p = NA,
               level = object$level.predict,
               df = object$k - 2)
  ##  
  ci.lab <- paste(round(100 * object$level.comb, 1), "%-CI", sep = "")
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- list(study = ci.study,
              fixed = ci.f, random = ci.r,
              predict = ci.p,
              k = object$k, Q = object$Q, df.Q = object$df.Q,
              Q.LRT = object$Q.LRT,
              tau = object$tau, H = ci.H, I2 = ci.I2, Rb = ci.Rb,
              H.resid = ci.H.resid, I2.resid = ci.I2.resid,
              tau.preset = object$tau.preset,
              k.all = length(object$TE),
              Q.CMH = object$Q.CMH,
              k.MH = object$k.MH,
              sm = object$sm, method = object$method,
              call = match.call(),
              ci.lab = ci.lab,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              prediction = prediction)
  ##  
  res$se.tau2    <- object$se.tau2
  res$hakn       <- object$hakn
  res$df.hakn    <- object$df.hakn
  res$method.tau <- object$method.tau
  res$TE.tau     <- object$TE.tau
  res$C          <- object$C
  ##
  ## Add results from subgroup analysis
  ##
  if (length(object$byvar) > 0) {
    ##
    ci.fixed.w <- list(TE = object$TE.fixed.w,
                       seTE = object$seTE.fixed.w,
                       lower = object$lower.fixed.w,
                       upper = object$upper.fixed.w,
                       z = object$zval.fixed.w,
                       p = object$pval.fixed.w,
                       level = object$level.comb,
                       harmonic.mean = object$n.harmonic.mean.w)
    ##
    if (metarate)
      ci.fixed.w$harmonic.mean <- object$t.harmonic.mean.w
    ##
    ci.random.w <- list(TE = object$TE.random.w,
                        seTE = object$seTE.random.w,
                        lower = object$lower.random.w,
                        upper = object$upper.random.w,
                        z = object$zval.random.w,
                        p = object$pval.random.w,
                        level = object$level.comb,
                        df = object$df.hakn.w,
                        harmonic.mean = object$n.harmonic.mean.w)
    ##
    if (metarate)
      ci.random.w$harmonic.mean <- object$t.harmonic.mean.w
    ##
    ci.H <- list(TE = object$H.w, lower = object$lower.H.w, upper = object$upper.H.w)
    ci.I2 <- list(TE = object$I2.w, lower = object$lower.I2.w, upper = object$upper.I2.w)
    ci.Rb <- list(TE = object$Rb.w, lower = object$lower.Rb.w, upper = object$upper.Rb.w)
    ## 
    res$within.fixed    <- ci.fixed.w
    res$within.random   <- ci.random.w
    res$k.w             <- object$k.w
    res$Q.w             <- object$Q.w
    res$Q.w.fixed       <- object$Q.w.fixed
    res$Q.w.random      <- object$Q.w.random
    res$df.Q.w          <- object$df.Q.w
    res$pval.Q.w        <- object$pval.Q.w
    res$Q.b.fixed       <- object$Q.b.fixed
    res$Q.b.random      <- object$Q.b.random
    res$df.Q.b          <- object$df.Q.b
    res$pval.Q.b.fixed  <- object$pval.Q.b.fixed
    res$pval.Q.b.random <- object$pval.Q.b.random
    res$tau.w           <- object$tau.w
    res$C.w             <- object$C.w
    res$H.w             <- ci.H
    res$I2.w            <- ci.I2
    res$Rb.w            <- ci.Rb
    res$bylab           <- bylab
    res$tau.common      <- object$tau.common
    res$bylevs          <- object$bylevs
  }
  ##
  class(res) <- "summary.meta"
  ##
  if (inherits(object, "metabin")) {
    res$sparse      <- object$sparse
    res$incr        <- object$incr
    res$allincr     <- object$allincr
    res$addincr     <- object$addincr
    res$allstudies  <- object$allstudies
    res$doublezeros <- object$doublezeros
    res$MH.exact    <- object$MH.exact
    ##
    res$model.glmm   <- object$model.glmm
    res$.glmm.fixed  <- object$.glmm.fixed
    res$.glmm.random <- object$.glmm.random
    ##
    class(res) <- c(class(res), "metabin")
  }
  ##
  if (inherits(object, "metacont")) {
    res$pooledvar  <- object$pooledvar
    res$method.smd <- object$method.smd
    res$sd.glass   <- object$sd.glass
    res$exact.smd  <- object$exact.smd
    ##
    class(res) <- c(class(res), "metacont")
  }
  ##
  if (inherits(object, "metacor")) {
    res$cor <- object$cor
    res$n   <- object$n
    ##
    res$null.effect <- object$null.effect
    ##
    class(res) <- c(class(res), "metacor")
  }
  ##
  if (inherits(object, "metagen")) {
    res$n.e <- object$n.e
    res$n.c <- object$n.c
    ##
    res$null.effect <- object$null.effect
    ##
    class(res)  <- c(class(res), "metagen")
  }
  ##
  if (inherits(object, "metainc")) {
    class(res)  <- c(class(res), "metainc")
    res$sparse  <- object$sparse
    res$incr    <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
    ##
    res$model.glmm   <- object$model.glmm
    res$.glmm.fixed  <- object$.glmm.fixed
    res$.glmm.random <- object$.glmm.random
  }
  ##
  if (inherits(object, "metamean")) {
    res$n    <- object$n
    res$mean <- object$mean
    res$sd   <- object$sd
    ##
    res$null.effect <- object$null.effect
    ##
    class(res)  <- c(class(res), "metamean")
  }
  ##
  if (metaprop) {
    res$event <- object$event
    res$n     <- object$n
    ##
    res$sparse  <- object$sparse
    res$incr    <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
    ##
    res$null.effect <- object$null.effect
    ##
    res$method.ci <- object$method.ci
    ##
    res$model.glmm   <- object$model.glmm
    res$.glmm.fixed  <- object$.glmm.fixed
    res$.glmm.random <- object$.glmm.random
    ##
    class(res) <- c(class(res), "metaprop")
  }
  ##
  if (is.prop(object$sm)) {
    res$event <- object$event
    res$n     <- object$n
    ##
    res$null.effect <- object$null.effect    
  }
  ##
  if (metarate) {
    res$event <- object$event
    res$time  <- object$time
    ##
    res$sparse  <- object$sparse
    res$incr    <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
    ##
    res$null.effect <- object$null.effect
    ##
    res$model.glmm   <- object$model.glmm
    res$.glmm.fixed  <- object$.glmm.fixed
    res$.glmm.random <- object$.glmm.random
    ##
    class(res) <- c(class(res), "metarate")
  }
  ##
  if (is.rate(object$sm)) {
    res$event     <- object$event
    res$time      <- object$time
    ##
    res$null.effect <- object$null.effect    
  }
  ##
  if (inherits(object, "trimfill")) {
    res$object <- object
    res$k0     <- object$k0
    ##
    class(res) <- c(class(res), "trimfill")
  }
  ##
  if (inherits(object, "metabind")) {
    res$null.effect <- object$null.effect    
    ##
    class(res) <- c(class(res), "metabind")
  }
  ##
  ## Function metamiss() from R package metasens
  ##
  if (inherits(object, "metamiss")) {
    ##
    res$event.e <- object$event.e
    res$miss.e <- object$miss.e
    res$n.e <- object$n.e + object$miss.e
    ##
    res$event.c <- object$event.c
    res$miss.c <- object$miss.c
    res$n.c <- object$n.c + object$miss.c
    ##
    res$IMOR.e <- object$IMOR.e
    res$IMOR.c <- object$IMOR.c
    ##
    res$method.miss <- object$method.miss
    res$small.values <- object$small.values
    ##
    res$incr <- object$incr
    res$p.e <- object$p.e
    res$p.c <- object$p.c
    ##
    res$pmiss.e <- object$pmiss.e
    res$pmiss.c <- object$pmiss.c
    ##
    res$p.star.e <- object$p.star.e
    res$p.star.c <- object$p.star.c
    ##
    res$var.p.star.e <- object$var.p.star.e
    res$var.p.star.c <- object$var.p.star.c
    ##
    class(res) <- c(class(res), "metamiss")
  }
  ##
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  ##
  res$print.byvar <- print.byvar
  res$byseparator <- byseparator
  res$print.CMH   <- print.CMH
  ##
  res$data   <- object$data
  res$subset <- object$subset
  ##
  res$backtransf <- backtransf
  res$pscale <- pscale
  res$irscale <- irscale
  res$irunit  <- irunit
  ##
  res$version <- object$version
  if (is.null(res$version))
    res$version <- packageDescription("meta")$Version
  ##
  res$version.metafor <- object$version.metafor
  
  
  res
}





#' @rdname summary.meta
#' @method print summary.meta
#' @export
#' @export print.summary.meta


print.summary.meta <- function(x,
                               comb.fixed = x$comb.fixed,
                               comb.random = x$comb.random,
                               prediction = x$prediction,
                               print.byvar = x$print.byvar,
                               byseparator = x$byseparator,
                               print.CMH = x$print.CMH,
                               header = TRUE,
                               backtransf = x$backtransf,
                               pscale = x$pscale,
                               irscale = x$irscale,
                               irunit = x$irunit,
                               bylab.nchar = 35,
                               digits = gs("digits"),
                               digits.zval = gs("digits.zval"),
                               digits.pval = max(gs("digits.pval"), 2),
                               digits.pval.Q = max(gs("digits.pval.Q"), 2),
                               digits.Q = gs("digits.Q"),
                               digits.tau2 = gs("digits.tau2"),
                               digits.H = gs("digits.H"),
                               digits.I2 = gs("digits.I2"),
                               scientific.pval = gs("scientific.pval"),
                               big.mark = gs("big.mark"),
                               print.I2 = gs("print.I2"),
                               print.H = gs("print.H"),
                               print.Rb = gs("print.Rb"),
                               text.tau2 = gs("text.tau2"),
                               text.I2 = gs("text.I2"),
                               text.Rb = gs("text.Rb"),
                               warn.backtransf = FALSE,
                               ...) {
  
  
  ##
  ##
  ## (1) Check for summary.meta object
  ##
  ##
  chkclass(x, "summary.meta")
  ##
  if (inherits(x, "metacum") | inherits(x, "metainf"))
    return(invisible(NULL))
  ##
  by <- !is.null(x$bylab)
  
  
  ##
  ##
  ## (2) Check and set other arguments
  ##
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.H, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  chklogical(scientific.pval)
  ##
  if (is.untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  ##
  chklogical(print.I2)
  chklogical(print.H)
  chklogical(print.Rb)
  chkchar(text.tau2)
  chkchar(text.I2)
  chkchar(text.Rb)
  chklogical(warn.backtransf)
  is.prop <- is.prop(x$sm)
  is.rate <- is.rate(x$sm)
  ##
  if (!is.prop & x$sm != "RD")
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, single = TRUE)
  else
    pscale <- 1
  if (!backtransf & pscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate & x$sm != "IRD")
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, single = TRUE)
  else
    irscale <- 1
  if (!backtransf & irscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  if (by) {
    chklogical(print.byvar)
    chkchar(byseparator)
  }
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  chklogical(header)
  chknumeric(bylab.nchar)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  cl <- paste("update.meta() or ", class(x)[1], "()", sep = "")
  addargs <- names(list(...))
  ##
  fun <- "print.summary.meta"
  ##
  warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  k.all <- length(x$study$TE)
  k <- x$k
  sm <- x$sm
  ##
  bip <- inherits(x, c("metabin", "metainc", "metaprop", "metarate"))
  null.given <- (inherits(x, c("metacor", "metagen", "metamean",
                               "metaprop", "metarate")) |
                 is.prop(sm) | is.rate(sm) | is.cor(sm) | is.mean(sm))
  ##
  null.effect <- x$null.effect
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
      sm.lab <- paste("log", sm, sep = "")
  ##
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  ##
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  ##
  if (by)
    bylevs <- ifelse(nchar(x$bylevs) > bylab.nchar,
                     paste(substring(x$bylevs, 1, bylab.nchar - 4), " ...", sep = ""),
                     x$bylevs)
  
  
  ##
  ##
  ## (4) Set and backtransform results of meta-analysis
  ##
  ##
  TE.fixed    <- x$fixed$TE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  lowTE.predict <- x$predict$lower
  uppTE.predict <- x$predict$upper
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
    TE.fixed.w     <- x$within.fixed$TE
    lowTE.fixed.w  <- x$within.fixed$lower
    uppTE.fixed.w  <- x$within.fixed$upper
    pval.fixed.w   <- x$within.fixed$p
    harmonic.mean.w <- x$within.fixed$harmonic.mean
    TE.random.w    <- x$within.random$TE
    lowTE.random.w <- x$within.random$lower
    uppTE.random.w <- x$within.random$upper
    pval.random.w   <- x$within.random$p
    ##
    Q.b.fixed <- x$Q.b.fixed
    Q.w.fixed <- x$Q.w.fixed
    Q.b.random <- x$Q.b.random
    Q.w.random <- x$Q.w.random
    ##
    Q.w <- x$Q.w
    ##
    k.w <- x$k.w
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
    if (sm %in% c("IR", "IRLN", "IRS", "IRFT"))
      harmonic.mean <- 1 / mean(1 / x$time)
    else
      harmonic.mean <- 1 / mean(1 / x$n)
    ##
    TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                              harmonic.mean, warn = comb.fixed & warn.backtransf)
    lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                              harmonic.mean, warn = comb.fixed & warn.backtransf)
    uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                              harmonic.mean, warn = comb.fixed & warn.backtransf)
    ##
    TE.random <- backtransf(TE.random, sm, "mean",
                            harmonic.mean, warn = comb.random & warn.backtransf)
    lowTE.random <- backtransf(lowTE.random, sm, "lower",
                               harmonic.mean, warn = comb.random & warn.backtransf)
    uppTE.random <- backtransf(uppTE.random, sm, "upper",
                               harmonic.mean, warn = comb.random & warn.backtransf)
    ##
    lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                harmonic.mean, warn = prediction & warn.backtransf)
    uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                harmonic.mean, warn = prediction & warn.backtransf)
    ##
    if (by) {
      TE.fixed.w     <- backtransf(TE.fixed.w, sm, "mean",
                                   harmonic.mean.w, warn = comb.fixed & warn.backtransf)
      lowTE.fixed.w  <- backtransf(lowTE.fixed.w, sm, "lower",
                                   harmonic.mean.w, warn = comb.fixed & warn.backtransf)
      uppTE.fixed.w  <- backtransf(uppTE.fixed.w, sm, "upper",
                                   harmonic.mean.w, warn = comb.fixed & warn.backtransf)
      ##
      TE.random.w    <- backtransf(TE.random.w, sm, "mean",
                                   harmonic.mean.w, warn = comb.random & warn.backtransf)
      lowTE.random.w <- backtransf(lowTE.random.w, sm, "lower",
                                   harmonic.mean.w, warn = comb.random & warn.backtransf)
      uppTE.random.w <- backtransf(uppTE.random.w, sm, "upper",
                                   harmonic.mean.w, warn = comb.random & warn.backtransf)
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
    }
  }
  ##
  ## Round and round ...
  ##
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed <- x$fixed$p
  zTE.fixed <- round(x$fixed$z, digits.zval)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$random$p
  zTE.random <- round(x$random$z, digits.zval)
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
    if (print.I2)
      I2.w <- round(100 * x$I2.w$TE, digits.I2)
    ##
    if (print.Rb)
      Rb.w <- round(100 * x$Rb.w$TE, digits.I2)
  }
  ##
  if (print.H) {
    H <- round(x$H$TE, digits.H)
    lowH <- round(x$H$lower, digits.H)
    uppH <- round(x$H$upper, digits.H)
  }
  ##
  if (print.I2) {
    I2 <- round(100 * x$I2$TE, digits.I2)
    lowI2 <- round(100 * x$I2$lower, digits.I2)
    uppI2 <- round(100 * x$I2$upper, digits.I2)
    print.ci.I2 <- ((Q > k & k >= 2) | (Q <= k & k > 2)) &
      !(is.na(lowI2) | is.na(uppI2))
    if (is.na(print.ci.I2))
      print.ci.I2 <- FALSE
  }
  else
    print.ci.I2 <- FALSE
  ##
  if (print.Rb) {
    Rb <- round(100 * x$Rb$TE, digits.I2)
    lowRb <- round(100 * x$Rb$lower, digits.I2)
    uppRb <- round(100 * x$Rb$upper, digits.I2)
  }
  
  
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
    res <- cbind(formatN(TE.fixed, digits, "NA",
                         big.mark = big.mark),
                 formatCI(formatN(lowTE.fixed, digits, "NA",
                                  big.mark = big.mark),
                          formatN(uppTE.fixed, digits, "NA",
                                  big.mark = big.mark)),
                 formatN(zTE.fixed, digits.zval, big.mark = big.mark),
                 formatPT(pTE.fixed, digits = digits.pval,
                          scientific = scientific.pval))
    dimnames(res) <- list("", c(sm.lab, x$ci.lab, "z", "p-value"))
    prmatrix(res, quote = FALSE, right = TRUE, ...)
    ## Print information on summary method:
    catmeth(class = class(x),
            method = x$method,
            sm = sm,
            k.all = k.all,
            sparse = ifelse(bip, x$sparse, FALSE),
            incr = if (bip) x$incr else FALSE,
            allincr = ifelse(bip, x$allincr, FALSE),
            addincr = ifelse(bip, x$addincr, FALSE),
            allstudies = x$allstudies,
            doublezeros = x$doublezeros,
            method.ci = x$method.ci,
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
            digits = digits, digits.tau2 = digits.tau2,
            method.miss = x$method.miss, IMOR.e = x$IMOR.e, IMOR.c = x$IMOR.c)
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
    if (comb.fixed | comb.random | prediction) {
      if (!inherits(x, "trimfill")) {
        if (x$method == "MH" &&
            (inherits(x, c("metabin", "metainc")) &
             comb.fixed & sm %in% c("RD", "IRD") &
             (!is.null(x$k.MH) == 1 && k != x$k.MH)))
          cat(paste("Number of studies combined:   k.MH = ", x$k.MH,
                    " (fixed effect), k = ", format(k, big.mark = big.mark),
                    " (random effects)\n\n", sep = ""))
        else
          cat(paste("Number of studies combined: k = ",
                    format(k, big.mark = big.mark), "\n\n", sep = ""))
      }
      else
        cat(paste("Number of studies combined: k = ",
                  format(k, big.mark = big.mark),
                  " (with ",
                  format(x$k0, big.mark = big.mark),
                  " added studies)\n\n", sep = ""))
      ##
      res <- cbind(formatN(c(if (comb.fixed) TE.fixed,
                             if (comb.random) TE.random,
                             if (prediction) NA),
                           digits, "NA",
                           big.mark = big.mark),
                   formatCI(formatN(c(if (comb.fixed) lowTE.fixed,
                                      if (comb.random) lowTE.random,
                                      if (prediction) lowTE.predict),
                                    digits, "NA", big.mark = big.mark),
                            formatN(c(if (comb.fixed) uppTE.fixed,
                                      if (comb.random) uppTE.random,
                                      if (prediction) uppTE.predict),
                                    digits, "NA", big.mark = big.mark)),
                   formatN(c(if (comb.fixed) zTE.fixed,
                             if (comb.random) zTE.random,
                             if (prediction) NA),
                           digits = digits.zval, big.mark = big.mark),
                   formatPT(c(if (comb.fixed) pTE.fixed,
                              if (comb.random) pTE.random,
                              if (prediction) NA),
                            digits = digits.pval,
                            scientific = scientific.pval))
      if (prediction)
        res[dim(res)[1], c(1,3:4)] <- ""
      if (!is.null(x$hakn) && x$hakn) {
        if (comb.fixed & comb.random)
          zlab <- "z|t"
        else if (comb.fixed & !comb.random)
          zlab <- "z"
        else if (!comb.fixed & comb.random)
          zlab <- "t"
      }
      else
        zlab <- "z"
      ##
      dimnames(res) <- list(c(if (comb.fixed) "Fixed effect model",
                              if (comb.random) "Random effects model",
                              if (prediction) "Prediction interval"),  
                            c(sm.lab, x$ci.lab, zlab, "p-value"))
      prmatrix(res, quote = FALSE, right = TRUE, ...)
      ##
      if (inherits(x, "metabin") && print.CMH) {
        Qdata <- cbind(formatN(round(Q.CMH, digits.Q), digits.Q, "NA",
                               big.mark = big.mark),
                       df.Q.CMH,
                       formatPT(pval.Q.CMH,
                                digits = digits.pval.Q,
                                scientific = scientific.pval))
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        ##
        cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
    }
    else
      cat(paste("Number of studies: k = ", k, "\n", sep = ""))
    ##
    ## Print information on heterogeneity
    ##
    cat("\nQuantifying heterogeneity:\n")
    ##
    cathet(k, 
           TRUE, text.tau2, x$tau, digits.tau2, big.mark,
           print.H, H, lowH, uppH, digits.H,
           print.I2, print.ci.I2, text.I2,
           I2, lowI2, uppI2, digits.I2,
           print.Rb, text.Rb, Rb, lowRb, uppRb)
    ##
    ## Print information on residual heterogeneity
    ##
    if (by & !inherits(x, "metabind")) {
      ##
      Q.resid <- x$Q.w.fixed
      k.resid <- x$df.Q.w + 1
      ##
      if (print.H) {
        H.resid <- round(x$H.resid$TE, digits.H)
        lowH.resid <- round(x$H.resid$lower, digits.H)
        uppH.resid <- round(x$H.resid$upper, digits.H)
      }
      if (print.I2) {
        I2.resid <- round(100 * x$I2.resid$TE, digits.I2)
        lowI2.resid <- round(100 * x$I2.resid$lower, digits.I2)
        uppI2.resid <- round(100 * x$I2.resid$upper, digits.I2)
        print.ci.I2 <-
          ((Q.resid  > k.resid & k.resid >= 2) |
           (Q.resid <= k.resid & k.resid > 2)) &
          !(is.na(lowI2.resid) | is.na(uppI2.resid))
        ##
        if (is.na(print.ci.I2))
          print.ci.I2 <- FALSE
      }
      ##
      if (!is.na(I2.resid)) {
        cat("\nQuantifying residual heterogeneity:\n")
        ##
        cathet(k.resid, 
               x$tau.common, text.tau2, unique(x$tau.w),
               digits.tau2, big.mark,
               print.H, H.resid, lowH.resid, uppH.resid, digits.H,
               print.I2, print.ci.I2, text.I2,
               I2.resid, lowI2.resid, uppI2.resid, digits.I2,
               FALSE, text.Rb, NA, NA, NA)
      }
    }
    ##
    ## Test of heterogeneity
    ##
    if (comb.fixed | comb.random) {
      if (k > 1) {
        if (x$method != "GLMM") {
          Qdata <- cbind(formatN(round(Q, digits.Q), digits.Q, "NA",
                                 big.mark = big.mark),
                         format(df.Q, big.mark = big.mark),
                         formatPT(pval.Q,
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval))
          dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        }
        else {
          Qdata <- cbind(formatN(round(c(Q, Q.LRT), digits.Q), digits.Q, "NA",
                                 big.mark = big.mark),
                         format(c(df.Q, df.Q.LRT), big.mark = big.mark),
                         formatPT(c(pval.Q, pval.Q.LRT),
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval),
                         c("Wald-type", "Likelihood-Ratio"))
          dimnames(Qdata) <- list(rep("", 2),
                                  c("Q", "d.f.", "p-value", "Test"))
        }
        ##
        cat("\nTest of heterogeneity:\n")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
      ##
      if (by) {
        ##
        ## Print information for subgroup analysis
        ##
        if (comb.fixed) {
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
                         ifelse(k.w == 1, "--",
                                formatPT(x$tau.w^2,
                                         digits = digits.tau2,
                                         big.mark = big.mark,
                                         noblanks = TRUE)),
                         if (print.I2)
                           ifelse(is.na(I2.w),
                                  "--",
                                  paste(formatN(I2.w, digits.I2),
                                        "%", sep = "")),
                         if (print.Rb)
                           ifelse(is.na(Rb.w),
                                  "--",
                                  paste(formatN(Rb.w, digits.I2),
                                        "%", sep = ""))
                         )
          ##
          bylab <- bylabel(x$bylab, bylevs, print.byvar, byseparator,
                           big.mark = big.mark)
          ##
          dimnames(Tdata) <- list(bylab,
                                  c("  k", sm.lab, x$ci.lab,
                                    "Q", text.tau2,
                                    if (print.I2) text.I2,
                                    if (print.Rb) text.Rb)
                                  )
          if (inherits(x, "metabind"))
            cat("\nResults for meta-analyses (fixed effect model):\n")
          else
            cat("\nResults for subgroups (fixed effect model):\n")
          prmatrix(Tdata, quote = FALSE, right = TRUE, ...)
          ##
          if (!inherits(x, "metabind")) {
            cat("\nTest for subgroup differences (fixed effect model):\n")
            if (x$method == "MH") {
              Qdata <- cbind(formatN(round(Q.b.fixed, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             format(df.Q.b, big.mark = big.mark),
                             formatPT(pval.Q.b.fixed,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list("Between groups  ",
                                      c("Q", "d.f.", "p-value"))
              prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
            }
            else {
              Qs  <- c(Q.b.fixed, Q.w.fixed)
              dfs <- c(df.Q.b, df.Q.w)
              pvals <- c(pval.Q.b.fixed, pval.Q.w.fixed)
              Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             format(dfs, big.mark = big.mark),
                             formatPT(pvals,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                      c("Q", "d.f.", "p-value"))
              prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
            }
          }
        }
        ##
        if (comb.random) {
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
                         formatN(round(Q.w, digits.Q), digits.Q,
                                 big.mark = big.mark),
                         ifelse(k.w == 1, "--",
                                formatPT(x$tau.w^2,
                                         digits = digits.tau2,
                                         big.mark = big.mark,
                                         noblanks = TRUE)),
                         if (print.I2)
                           ifelse(is.na(I2.w),
                                  "--",
                                  paste(formatN(I2.w, digits.I2),
                                        "%", sep = "")),
                         if (print.Rb)
                           ifelse(is.na(Rb.w),
                                  "--",
                                  paste(formatN(Rb.w, digits.I2,
                                                big.mark = big.mark),
                                        "%", sep = ""))
                         )
          ##
          bylab <- bylabel(x$bylab, bylevs, print.byvar, byseparator,
                           big.mark = big.mark)
          ##
          dimnames(Tdata) <- list(bylab,
                                  c("  k", sm.lab, x$ci.lab,
                                    "Q", text.tau2,
                                    if (print.I2) text.I2,
                                    if (print.Rb) text.Rb)
                                  )
          if (inherits(x, "metabind"))
            cat("\nResults for meta-analyses (random effects model):\n")
          else
            cat("\nResults for subgroups (random effects model):\n")
          prmatrix(Tdata, quote = FALSE, right = TRUE, ...)
          ##
          if (!inherits(x, "metabind")) {
            cat("\nTest for subgroup differences (random effects model):\n")
            if (is.na(Q.w.random)) {
              Qdata <- cbind(formatN(round(Q.b.random, digits.Q), digits.Q,
                                     "NA", big.mark = big.mark),
                             format(df.Q.b, big.mark = big.mark),
                             formatPT(pval.Q.b.random,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list("Between groups  ",
                                      c("Q", "d.f.", "p-value"))
            }
            else {
              Qs  <- c(Q.b.random, Q.w.random)
              dfs <- c(df.Q.b, df.Q.w)
              pvals <- c(pval.Q.b.random, pval.Q.w.random)
              Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             format(dfs, big.mark = big.mark),
                             formatPT(pvals,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                      c("Q", "d.f.", "p-value"))
            }
            prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
          }
        }
      }
    }
    ##
    ## Print information on summary method:
    ##
    if (comb.fixed | comb.random | prediction)
      catmeth(class = class(x),
              method = x$method,
              method.tau = if (comb.random) x$method.tau else "",
              sm = sm,
              k.all = k.all,
              hakn = !is.null(x$hakn) && (x$hakn & comb.random),
              tau.common = by & x$tau.common,
              tau.preset = x$tau.preset,
              sparse = ifelse(bip, x$sparse, FALSE),
              incr = if (bip) x$incr else FALSE,
              allincr = ifelse(bip, x$allincr, FALSE),
              addincr = ifelse(bip, x$addincr, FALSE),
              allstudies = x$allstudies,
              doublezeros = x$doublezeros,
              MH.exact = ifelse(inherits(x, "metabin"), x$MH.exact, FALSE),
              method.ci = x$method.ci,
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
              digits = digits, digits.tau2 = digits.tau2,
              method.miss = x$method.miss, IMOR.e = x$IMOR.e, IMOR.c = x$IMOR.c)
  }
  
  
  invisible(NULL)
}
