#' Update a meta-analysis object
#' 
#' @description
#' Update an existing meta-analysis object.
#' 
#' @param object An object of class \code{meta}.
#' @param data Dataset.
#' @param subset Subset.
#' @param studlab Study label.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param cluster An optional vector specifying which estimates come
#'   from the same cluster resulting in the use of a three-level
#'   meta-analysis model.
#' @param rho Assumed correlation of estimates within a cluster.
#' @param cycles A numeric vector with the number of cycles per patient / study
#'   in n-of-1 trials (see \code{\link{metagen}}).
#' @param weights A single numeric or vector with user-specified weights.
#' @param weights.common User-specified weights (common effect model).
#' @param weights.random User-specified weights (random effects model).
#' @param method A character string indicating which method is to be
#'   used for pooling of studies (see \code{\link{metabin}},
#'   \code{\link{metainc}}, \code{\link{metaprop}} and
#'   \code{\link{metarate}}).
#' @param sm A character string indicating which summary measure is
#'   used for pooling.
#' @param incr Information on increment added to cell frequencies of
#'   studies with zero cell counts (see \code{\link{metabin}},
#'   \code{\link{metainc}}, \code{\link{metaprop}} and
#'   \code{\link{metarate}}).
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (see \code{\link{metabin}},
#'   \code{\link{metainc}}, \code{\link{metaprop}} and
#'   \code{\link{metarate}}).
#' @param allstudies A logical indicating if studies with zero or all
#'   events in both groups are to be included in the meta-analysis
#'   (applies only to \code{\link{metabin}} object with \code{sm}
#'   equal to \code{"RR"} or \code{"OR"}).
#' @param incr.e Continuity correction in experimental group (see
#'   \code{\link{metabin}} and \code{\link{metainc}}).
#' @param incr.c Continuity correction in control group (see
#'   \code{\link{metabin}} and \code{\link{metainc}}).
#' @param incr.event Continuity correction (see
#'   \code{\link{metaprop}} and \code{\link{metarate}}).
#' @param MH.exact A logical indicating if \code{incr} is not to be
#'   added to all cell frequencies for studies with a zero cell count
#'   to calculate the pooled estimate based on the Mantel-Haenszel
#'   method (applies only to \code{\link{metabin}} object).
#' @param RR.Cochrane A logical indicating which method to use as
#'   continuity correction for the risk ratio (see
#'   \code{\link{metabin}}).
#' @param Q.Cochrane A logical indicating which method to use to
#'   calculate the heterogeneity statistic Q (see
#'   \code{\link{metabin}}).
#' @param model.glmm A character string indicating which GLMM model
#'   should be used (see \code{\link{metabin}} and
#'   \code{\link{metainc}}).
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param common A logical indicating whether a common effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param overall A logical indicating whether overall summaries
#'   should be reported. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be reported.
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
#' @param method.common.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   common effect estimate (see \code{\link{meta-package}}).
#' @param method.random.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   random effects estimate (see \code{\link{meta-package}}).
#' @param adhoc.hakn.ci A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied in the case
#'   of an arbitrarily small Hartung-Knapp variance estimate (see
#'   \code{\link{meta-package}}).
#' @param method.predict A character string indicating which method is
#'   used to calculate a prediction interval (see
#'   \code{\link{meta-package}}).
#' @param adhoc.hakn.pi A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied for
#'   prediction interval (see \code{\link{meta-package}}).
#' @param seed.predict A numeric value used as seed to calculate
#'   bootstrap prediction interval (see \code{\link{meta-package}}).
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau} (see \code{\link{meta-package}}).
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau} (see \code{\link{meta-package}}).
#' @param level.hetstat The level used to calculate confidence intervals
#'   for heterogeneity statistics.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param detail.tau Detail on between-study variance estimate.
#' @param method.I2 A character string indicating which method is
#'   used to estimate the heterogeneity statistic I\eqn{^2}. Either
#'   \code{"Q"} or \code{"tau2"}, can be abbreviated
#'   (see \code{\link{meta-package}}).
#' @param method.bias A character string indicating which test for
#'   funnel plot asymmetry is to be used, can be abbreviated. See
#'   function \code{\link{metabias}.}
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE}, results for \code{sm = "OR"} are printed as odds ratios
#'   rather than log odds ratios and results for \code{sm = "ZCOR"}
#'   are printed as correlations rather than Fisher's z transformed
#'   correlations, for example.
#' @param func.backtransf A function used to back-transform results
#'   (see \code{\link{metagen}}).
#' @param args.backtransf An optional list to provide additional
#'   arguments to \code{func.backtransf} (see \code{\link{metagen}}).
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
#' @param text.common A character string used in printouts and forest
#'   plot to label the pooled common effect estimate.
#' @param text.random A character string used in printouts and forest
#'   plot to label the pooled random effects estimate.
#' @param text.predict A character string used in printouts and forest
#'   plot to label the prediction interval.
#' @param text.w.common A character string used to label weights of
#'   common effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
#' @param label.left Graph label on left side of null effect in forest plot.
#' @param label.right Graph label on right side of null effect in forest plot.
#' @param col.label.left The colour of the graph label on the left side of
#'   the null effect.
#' @param col.label.right The colour of the graph label on the right side of
#'   the null effect.
#' @param n.e Number of observations in experimental group (only for
#'   \code{\link{metagen}} or \code{\link{metainc}} object).
#' @param n.c Number of observations in control group (only for
#'   \code{\link{metagen}} or \code{\link{metainc}} object).
#' @param method.mean A character string indicating which method to
#'   use to approximate the mean from the median and other statistics
#'   (see \code{\link{metacont}} and \code{\link{metamean}}).
#' @param method.sd A character string indicating which method to use
#'   to approximate the standard deviation from sample size, median,
#'   interquartile range and range (see \code{\link{metacont}} and
#'   \code{\link{metamean}}).
#' @param approx.mean.e Approximation method to estimate means in
#'   experimental group (see \code{\link{metacont}}).
#' @param approx.mean.c Approximation method to estimate means in
#'   control group (see \code{\link{metacont}}).
#' @param approx.sd.e Approximation method to estimate standard
#'   deviations in experimental group (see \code{\link{metacont}}).
#' @param approx.sd.c Approximation method to estimate standard
#'   deviations in control group (see \code{\link{metacont}}).
#' @param approx.mean Approximation method to estimate means (see
#'   \code{\link{metamean}}).
#' @param approx.sd Approximation method to estimate standard
#'   deviations (see \code{\link{metamean}}).
#' @param approx.TE Approximation method to estimate treatment
#'   estimate (see \code{\link{metagen}}).
#' @param approx.seTE Approximation method to estimate standard error
#'   (see \code{\link{metagen}}).
#' @param pooledvar A logical indicating if a pooled variance should
#'   be used for the mean difference or ratio of means (see
#'   \code{\link{metacont}}).
#' @param method.smd A character string indicating which method is
#'   used to estimate the standardised mean difference (see
#'   \code{\link{metacont}}).
#' @param sd.glass A character string indicating which standard
#'   deviation is used in the denominator for Glass' method to
#'   estimate the standardised mean difference (only for metacont
#'   object with \code{sm = "SMD"}). Either \code{"control"} using the
#'   standard deviation in the control group (\code{sd.c}) or
#'   \code{"experimental"} using the standard deviation in the
#'   experimental group (\code{sd.e}), can be abbreviated.
#' @param exact.smd A logical indicating whether exact formulae should
#'   be used in estimation of the standardised mean difference and its
#'   standard error.
#' @param method.ci A character string indicating which method is used
#'   to calculate confidence intervals for individual studies. Either
#'   \code{"z"}, \code{"t"}, \code{"WS"}, \code{"WSCC"}, \code{"AC"},
#'   \code{"SA"}, \code{"SACC"}, \code{"NAsm"}, or \code{"Poisson"},
#'   can be abbreviated. See functions \code{\link{metacont}},
#'   \code{\link{metaprop}} and \code{\link{metarate}}.
#' @param subgroup An optional vector to conduct a meta-analysis with
#'   subgroups.
#' @param subgroup.name A character string with a name for the
#'   subgroup variable.
#' @param print.subgroup.name A logical indicating whether the name of
#'   the subgroup variable should be printed in front of the group
#'   labels.
#' @param sep.subgroup A character string defining the separator
#'   between name of subgroup variable and subgroup label.
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param prediction.subgroup A logical indicating whether prediction
#'   intervals should be printed for subgroups.
#' @param seed.predict.subgroup A numeric vector providing seeds to
#'   calculate bootstrap prediction intervals within subgroups. Must
#'   be of same length as the number of subgroups.
#' @param byvar Deprecated argument (replaced by 'subgroup').
#' @param id Deprecated argument (replaced by 'cluster').
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param left A logical indicating whether studies are supposed to be
#'   missing on the left or right side of the funnel plot. If NULL,
#'   the linear regression test for funnel plot symmetry (i.e.,
#'   function \code{metabias(..., method = "linreg")}) is used to
#'   determine whether studies are missing on the left or right side.
#' @param ma.common A logical indicating whether a common effect or
#'   random effects model is used to estimate the number of missing
#'   studies.
#' @param type A character indicating which method is used to estimate
#'   the number of missing studies. Either \code{"L"} or \code{"R"}.
#' @param n.iter.max Maximum number of iterations to estimate number
#'   of missing studies.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param verbose A logical indicating whether to print information on
#'   updates of older meta versions.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}, respectively.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @details
#' Wrapper function to update an existing meta-analysis object which
#' was created with R function \code{\link{metabin}},
#' \code{\link{metacont}}, \code{\link{metacor}},
#' \code{\link{metagen}}, \code{\link{metainc}},
#' \code{\link{metamean}}, \code{\link{metaprop}}, or
#' \code{\link{metarate}}. More details on function arguments are
#' available in help files of respective R functions.
#' 
#' This function can also be used for objects of class 'trimfill',
#' 'metacum', and 'metainf'.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metabin"},
#' \code{"metacont"}, \code{"metacor"}, \code{"metainc"},
#' \code{"metagen"}, \code{"metamean"}, \code{"metaprop"}, or
#' \code{"metarate"} (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metacor}}, \code{\link{metagen}},
#'   \code{\link{metainc}}, \code{\link{metamean}},
#'   \code{\link{metaprop}}, \code{\link{metarate}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, studlab = paste(study, year), sm = "SMD")
#' m1
#' 
#' # Change summary measure (from 'SMD' to 'MD')
#' #
#' update(m1, sm = "MD")
#' 
#' # Restrict analysis to subset of studies
#' #
#' update(m1, subset = 1:2)
#' 
#' # Use different levels for confidence intervals
#' #
#' m2 <- update(m1, level = 0.66, level.ma = 0.99)
#' print(m2, digits = 2)
#' forest(m2)
#' 
#' @method update meta
#' @export

update.meta <- function(object, 
                        data = object$data,
                        subset, studlab, exclude, cluster,
                        rho = object$rho,
                        cycles,
                        #
                        weights = NULL,
                        weights.common = object$weights.common,
                        weights.random = object$weights.random,
                        #
                        method,
                        sm = object$sm,
                        incr = object$incr,
                        method.incr = object$method.incr,
                        allstudies = object$allstudies,
                        #
                        incr.e, incr.c, incr.event,
                        #
                        MH.exact = object$MH.exact,
                        RR.Cochrane = object$RR.Cochrane,
                        Q.Cochrane = object$Q.Cochrane,
                        model.glmm = object$model.glmm,
                        #
                        level = object$level,
                        level.ma = object$level.ma,
                        #
                        common = object$common,
                        random = object$random,
                        overall = object$overall,
                        overall.hetstat = object$overall.hetstat,
                        method.common.ci = object$method.common.ci,
                        method.random.ci = object$method.random.ci,
                        adhoc.hakn.ci = object$adhoc.hakn.ci,
                        method.predict = object$method.predict,
                        adhoc.hakn.pi = object$adhoc.hakn.pi,
                        seed.predict = object$seed.predict,
                        method.tau = object$method.tau,
                        method.tau.ci = object$method.tau.ci,
                        level.hetstat = object$level.hetstat,
                        tau.preset = object$tau.preset,
                        TE.tau = object$TE.tau,
                        tau.common = object$tau.common,
                        detail.tau = object$detail.tau,
                        #
                        method.I2 = object$method.I2,
                        #
                        prediction = object$prediction,
                        level.predict = object$level.predict,
                        null.effect = object$null.effect,
                        method.bias = object$method.bias,
                        ##
                        backtransf = object$backtransf,
                        func.backtransf = object$func.backtransf,
                        args.backtransf = object$args.backtransf,    
                        #
                        pscale = object$pscale,
                        irscale = object$irscale,
                        irunit = object$irunit,
                        ##
                        text.common = object$text.common,
                        text.random = object$text.random,
                        text.predict = object$text.predict,
                        text.w.common = object$text.w.common,
                        text.w.random = object$text.w.random,
                        ##
                        title = object$title,
                        complab = object$complab,
                        outclab = object$outclab,
                        #
                        label.e = object$label.e,
                        label.c = object$label.c,
                        label.left = object$label.left,
                        label.right = object$label.right,
                        col.label.left = object$col.label.left,
                        col.label.right = object$col.label.right,
                        #
                        n.e = object$n.e,
                        n.c = object$n.c,
                        ##
                        method.mean = object$method.mean,
                        method.sd = object$method.sd,
                        ##
                        approx.mean.e = object$approx.mean.e,
                        approx.mean.c = object$approx.mean.c,
                        approx.sd.e = object$approx.sd.e,
                        approx.sd.c = object$approx.sd.c,
                        ##
                        approx.mean = object$approx.mean,
                        approx.sd = object$approx.sd,
                        ##
                        approx.TE = object$approx.TE,
                        approx.seTE = object$approx.seTE,
                        ##
                        pooledvar = object$pooledvar,
                        method.smd = object$method.smd,
                        sd.glass = object$sd.glass,
                        exact.smd = object$exact.smd,
                        method.ci = object$method.ci,
                        ##
                        subgroup,
                        subgroup.name = object$subgroup.name,
                        print.subgroup.name = object$print.subgroup.name,
                        sep.subgroup = object$sep.subgroup,
                        test.subgroup = object$test.subgroup,
                        prediction.subgroup = object$prediction.subgroup,
                        seed.predict.subgroup = object$seed.predict.subgroup,
                        ##
                        byvar, id,
                        ##
                        print.CMH = object$print.CMH,
                        keepdata = TRUE,
                        ##
                        left = object$left,
                        ma.common = object$ma.common,
                        type = object$type,
                        n.iter.max = object$n.iter.max,
                        ##
                        warn = FALSE, warn.deprecated = gs("warn.deprecated"),
                        verbose = FALSE,
                        ##
                        control = object$control,
                        ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and update older meta objects
  ##
  ##
  chkclass(object, "meta")
  suitable <-
    chksuitable(object, "Update",
                c("metabind", "metaadd", "metamerge", "metamiss"),
                check.mlm = FALSE,
                stop = FALSE, status = "possible")
  if (!suitable)
    return(object)
  ##
  metabin  <- inherits(object, "metabin")
  metacont <- inherits(object, "metacont")
  metacor  <- inherits(object, "metacor")
  metagen  <- inherits(object, "metagen")
  metainc  <- inherits(object, "metainc")
  metamean <- inherits(object, "metamean")
  metaprop <- inherits(object, "metaprop")
  metarate <- inherits(object, "metarate")
  ##
  missing.method.common.ci <- missing(method.common.ci)
  missing.method.random.ci <- missing(method.random.ci)
  missing.adhoc.hakn.ci <- missing(adhoc.hakn.ci)
  missing.text.random <- missing(text.random)
  ##
  missing.method.predict <- missing(method.predict)
  missing.adhoc.hakn.pi <- missing(adhoc.hakn.pi)
  missing.text.predict <- missing(text.predict)
  #
  missing.method.bias <- missing(method.bias)
  missing.method.incr <- missing(method.incr)
  #
  chklogical(verbose)
  ##
  if (update_needed(object$version, 3, 2, verbose)) {
    ##
    ## Changes for meta objects with version < 3.2
    ##
    object$subset <- NULL
    ##
    object$data <- data.frame(.studlab = object$studlab,
                              .exclude = rep_len(FALSE,
                                                 length(object$studlab)))
    ##
    if (!is.null(object$byvar))
      object$data$.byvar <- object$byvar
    ##
    if (metabin) {
      object$data$.event.e <- object$event.e
      object$data$.n.e <- object$n.e
      object$data$.event.c <- object$event.c
      object$data$.n.c <- object$n.c
    }
    ##
    if (metacont) {
      object$data$.n.e <- object$n.e
      object$data$.mean.e <- object$mean.e
      object$data$.sd.e <- object$sd.e
      object$data$.n.c <- object$n.c
      object$data$.mean.c <- object$mean.c
      object$data$.sd.c <- object$sd.c
    }
    ##
    if (metacor) {
      object$data$.cor <- object$cor
      object$data$.n <- object$n
    }
    ##
    if (metagen) {
      object$data$.TE <- object$TE
      object$data$.seTE <- object$seTE
    }
    ##
    if (metaprop) {
      object$data$.event <- object$event
      object$data$.n <- object$n
    }
  }
  ##
  if (update_needed(object$version, 4, 8, verbose)) {
    ##
    ## Changes for meta objects with version < 4.8
    ##
    if (metabin | metainc | metaprop | metarate)
      object$data$.incr <- object$incr
    ##
    if (metabin | metainc)
      if (object$method == "MH")
        object$k.MH <- sum(object$w.fixed > 0)
      else
        object$k.MH <- NA
  }
  ##
  if (update_needed(object$version, 5, 0, verbose)) {
    ##
    ## Changes for meta objects with version < 5.0
    ##
    object$fixed <- object$comb.fixed
    object$random <- object$comb.random
    object$level.ma <- object$level.comb
    ##
    object$detail.tau <- ""
    ##
    if (!is.null(object$byvar)) {
      object$data$.subgroup <- object$byvar
      object$subgroup.name <- object$bylab
      object$print.subgroup.name <- object$print.byvar
      object$sep.subgroup <- object$byseparator
    }
  }
  #
  if (update_needed(object$version, 5, 5, verbose)) {
    ##
    ## Changes for meta objects with version < 5.5
    ##
    object$common <- object$fixed
    object$w.common <- object$w.fixed
    ##
    object$TE.common <- object$TE.fixed
    object$seTE.common <- object$seTE.fixed
    object$lower.common <- object$lower.fixed
    object$upper.common <- object$upper.fixed
    object$statistic.common <- object$statistic.fixed
    object$pval.common <- object$pval.fixed
    object$zval.common <- object$zval.fixed
    ##
    object$text.common <- object$text.fixed
    object$text.w.common <- object$text.w.fixed
    ##
    if (!is.null(object$pooled) && object$pooled == "fixed")
      object$pooled <- "common"
    ##
    if (!is.null(object$byvar)) {
      object$TE.common.w <- object$TE.fixed.w
      object$seTE.common.w <- object$seTE.fixed.w
      object$lower.common.w <- object$lower.fixed.w
      object$upper.common.w <- object$upper.fixed.w
      object$statistic.common.w <- object$statistic.fixed.w
      object$pval.common.w <- object$pval.fixed.w
      object$zval.common.w <- object$zval.fixed.w
      object$w.common.w <- object$w.fixed.w
      ##
      object$Q.w.common <- object$Q.w.fixed
      object$pval.Q.w.common <- object$pval.Q.w.fixed
      object$Q.b.common <- object$Q.b.fixed
      object$pval.Q.b.common <- object$pval.Q.b.fixed
    }
    ##
    method.incr <- gs("method.incr")
    if (is.logical(object$addincr) && object$addincr)
      method.incr <- "all"
    else if (is.logical(object$allincr) && object$allincr)
      method.incr <- "if0all"
    ##
    if (!is.null(object$id))
      object$cluster <- object$id
      object$data$.cluster <- object$data$.id
  }
  #
  if (update_needed(object$version, 6, 0, verbose)) {
    ##
    ## Changes for meta objects with version < 6.0
    ##
    object$method.random.ci <- if (object$hakn) "HK" else "classic"
    object$adhoc.hakn.ci <- object$adhoc.hakn
    object$df.random <- object$df.hakn
    object$seTE.hakn.ci <- object$seTE.hakn
    object$seTE.hakn.adhoc.ci <- object$seTE.hakn.adhoc
    ##
    object$method.predict <- "HTS"
    object$adhoc.hakn.pi <- ""
    object$df.predict <- object$k - 2
    object$seTE.hakn.pi <- NA
    object$seTE.hakn.adhoc.pi <- NA
    ##
    object$seTE.kero <- NA
    object$df.kero <- NA
    ##
    if (!is.null(object$subgroup)) {
      object$bylevs <- object$subgroup.levels
      object$df.random.w <- object$df.hakn.w
      object$df.predict.w <- object$k.w - 2
      ##
      object$seTE.hakn.ci.w <-
        object$seTE.hakn.adhoc.ci.w <-
          object$seTE.hakn.pi.w <-
            object$seTE.hakn.adhoc.pi.w <- NA
      ##
      object$df.Q.b.random <- object$df.Q.b.common <- object$df.Q.b
    }
  }
  #
  if (update_needed(object$version, 6, 1, verbose)) {
    ##
    ## Changes for meta objects with version < 6.1
    ##
    object$transf <- TRUE
  }
  #
  if (update_needed(object$version, 6, 5, verbose)) {
    ##
    ## Changes for meta objects with version < 6.5
    ##
    if (length(object$TE.random) == 1 & length(object$lower.random) > 1) {
      object$TE.random <- rep(object$TE.random, length(object$lower.random))
      names(object$TE.random) <- names(object$lower.random)
    }
    ##
    if (length(object$seTE.random) == 1 & length(object$lower.random) > 1) {
      object$seTE.random <- rep(object$seTE.random, length(object$lower.random))
      names(object$seTE.random) <- names(object$lower.random)
    }
    ##
    object$method.random <- object$method
    object$method.random[object$method.random %in% c("MH", "Cochran")] <-
      "Inverse"
    ##
    if (!is.null(object$Q.LRT)) {
      object$Q <- c(object$Q, object$Q.LRT)
      object$df.Q <- c(object$df.Q, object$df.Q.LRT)
      object$pval.Q <- c(object$pval.Q, object$pval.Q.LRT)
      names(object$Q) <- c("Wald", "LRT")
    }
    ##
    if (metacont) {
      object$data$.approx.mean.e <-
        setVal(object$data, ".approx.mean.e", object$approx.mean.e)
      object$data$.approx.mean.c <-
        setVal(object$data, ".approx.mean.c", object$approx.mean.c)
      object$data$.approx.sd.e <-
        setVal(object$data, ".approx.sd.e", object$approx.sd.e)
      object$data$.approx.sd.c <-
        setVal(object$data, ".approx.sd.c", object$approx.sd.c)
    }
    ##
    if (metamean) {
      object$data$.approx.mean <-
        setVal(object$data, ".approx.mean", object$approx.mean)
      object$data$.approx.sd <-
        setVal(object$data, ".approx.sd", object$approx.sd)
    }
    ##
    object$hetlabel <- object$label
    ##
    if (length(object$tau) > 1)
      names(object$tau) <- object$detail.tau
    if (length(object$tau2) > 1)
      names(object$tau2) <- object$detail.tau
    if (length(object$I2) > 1)
      names(object$I2) <- object$detail.tau
    ##
    object$seed.predict <- NULL
    if (!is.null(object$byvar))
      object$seed.predict.subgroup <- NULL
  }
  #
  if (update_needed(object$version, 7, 0, verbose)) {
    ##
    ## Changes for meta objects with version < 7.0
    ##
    object$rho <- 0
    ##
    if (inherits(object, c("metacum", "metainf"))) {
      object$label.e <- replaceNULL(object$label.e, "")
      object$label.c <- replaceNULL(object$label.c, "")
    }
    ##
    if (inherits(object, "metaprop") && object$method.ci != "NAsm") {
      if (object$sm == "PLOGIT") {
        object$lower <- p2logit(object$lower)
        object$upper <- p2logit(object$upper)
      }
      ##
      else if (object$sm == "PAS") {
        object$lower <- p2asin(object$lower)
        object$upper <- p2asin(object$upper)
      }
      ##
      else if (object$sm == "PFT") {
        lower.ev <- object$n * object$lower 
        upper.ev <- object$n * object$upper 
        ##
        object$lower <-
          0.5 * (asin(sqrt(lower.ev / object$n)) +
                 asin(sqrt((lower.ev + 1) / object$n)))
        object$upper <-
          0.5 * (asin(sqrt(upper.ev / object$n)) +
                 asin(sqrt((upper.ev + 1) / object$n)))
      }
      ##
      else if (object$sm == "PLN") {
        object$lower <- log(object$lower)
        object$upper <- log(object$upper)
      }
    }
    ##
    if (inherits(object, "metarate") && object$method.ci != "NAsm") {
      if (object$sm == "IRLN") {
        object$lower <- log(object$lower)
        object$upper <- log(object$upper)
      }
      else if (object$sm == "IRS") {
        object$lower <- sqrt(object$lower)
        object$upper <- sqrt(object$upper)
      }
      ##
      else if (object$sm == "IRFT") {
        lower.ev <- object$time * object$lower 
        upper.ev <- object$time * object$upper 
        ##
        object$lower <-
          0.5 * (sqrt(lower.ev / object$time) +
                 sqrt((lower.ev + 1) / object$time))
        object$upper <-
          0.5 * (sqrt(upper.ev / object$time) +
                 sqrt((upper.ev + 1) / object$time))
      }
      ##
      if (inherits(object, "metabind")) {
        object$with.subgroups <- any(object$is.subgroup)
        ##
        if (object$with.subgroups) {
          object$data$Q.b.common <- object$data$Q.b
          object$data$Q.b.random <- object$data$Q.b
          ##
          object$data$pval.Q.b.common <- object$data$pval.Q.b
          object$data$pval.Q.b.random <- object$data$pval.Q.b
        }
      }
    }
    #
    object$seTE.kero <-
      replaceNA(object$seTE.kero, kenwardroger(object$w.random)$se)
    object$df.kero <-
      replaceNA(object$df.kero, kenwardroger(object$w.random)$df)
  }
  #
  if (update_needed(object$version, 8, 0, verbose)) {
    #
    # Changes for meta objects with version < 8.0
    #
    object$method.I2 <- "Q"
    #
    object$col.label.left <- replaceNULL(object$col.label.left, "black")
    object$col.label.right <- replaceNULL(object$col.label.left, "black")
    #
    object$level.hetstat <- object$level.ma
  }
  #
  if (update_needed(object$version, 8, 1, verbose)) {
    #
    # Changes for meta objects with version < 8.1
    #
    object$pairwise <- FALSE
    #
    if (object$keepdata) {
      if (inherits(object, c("metabin", "metainc"))) {
        if (isCol(object$data, ".subset")) {
          object$incr.e <- object$data$.incr[object$data$.subset]
          object$incr.c <- object$data$.incr[object$data$.subset]
        }
        else {
          object$incr.e <- object$data$.incr
          object$incr.c <- object$data$.incr
        }
        #
        object$data$.incr.e <- object$data$.incr
        object$data$.incr.c <- object$data$.incr
        object$data$.incr <- NULL
      }
      #
      if (inherits(object, "metainc"))
        object$incr.event <- NULL
      #
      if (inherits(object, c("metaprop", "metarate"))) {
        object$data$.incr.event <- object$data$.incr
        object$data$.incr <- NULL
      }
    }
  }
  #
  if (update_needed(object$version, 8, 2, verbose)) {
    ##
    ## Changes for meta objects with version < 8.2
    ##
    object$method.common.ci <- "classic"
  }
  
  
  ##
  ##
  ## (2) Check arguments
  ##
  ##
  if (missing(method))
    method <- object$method
  ##
  pscale <- replaceNULL(pscale, 1)
  irscale <- replaceNULL(irscale, 1)
  irunit <- replaceNULL(irunit, "")
  ##
  tau.common <- replaceNULL(tau.common, gs("tau.common"))
  method.I2 <- replaceNULL(method.I2, gs("method.I2"))
  sep.subgroup <- replaceNULL(sep.subgroup, gs("sep.subgroup"))
  ##
  if (!backtransf & pscale != 1 & !is_untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!backtransf & irscale != 1 & !is_untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  args <- list(...)
  nam.args <- names(args)
  #
  chklogical(warn.deprecated)
  ##
  level.ma <- deprecated(level.ma, missing(level.ma), args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                      warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                      warn.deprecated)
  chklogical(common)
  ##
  missing.random <- missing(random)
  random <- deprecated(random, missing.random, args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  #
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  method.random.ci <-
    deprecated(method.random.ci, missing.method.random.ci,
               args, "hakn", warn.deprecated)
  if (is.logical(method.random.ci))
    if (method.random.ci)
      method.random.ci <- "HK"
    else
      method.random.ci <- "classic"
  method.random.ci <- setchar(method.random.ci, gs("meth4random.ci"))
  ##
  missing.subgroup.name <- missing(subgroup.name)
  subgroup.name <-
    deprecated(subgroup.name, missing.subgroup.name, args, "bylab",
               warn.deprecated)
  ##
  print.subgroup.name <-
    deprecated(print.subgroup.name, missing(print.subgroup.name),
               args, "print.byvar", warn.deprecated)
  print.subgroup.name <-
    replaceNULL(print.subgroup.name, gs("print.subgroup.name"))
  chklogical(print.subgroup.name)
  ##
  sep.subgroup <-
    deprecated(sep.subgroup, missing(sep.subgroup), args, "byseparator",
               warn.deprecated)
  if (!is.null(sep.subgroup))
    chkchar(sep.subgroup, length = 1)
  ##
  test.subgroup <- replaceNULL(test.subgroup, gs("test.subgroup"))
  prediction.subgroup <-
    replaceNULL(prediction.subgroup, gs("prediction.subgroup"))
  ##
  missing.method.incr <- missing(method.incr)
  #
  # Ignore deprecated arguments 'addincr' and 'allincr'
  #
  txt.ignore <- "(deprecated); use argument 'method.incr'"
  #
  addincr <- allincr <- NULL
  if (!is.na(charmatch("addincr", nam.args)))
    warn_ignore_input(addincr, TRUE, txt.ignore)
  if (!is.na(charmatch("allincr", nam.args)))
    warn_ignore_input(allincr, TRUE, txt.ignore)
  ##
  ## Some more checks
  ##
  overall <- replaceNULL(overall, common | random)
  overall.hetstat <- replaceNULL(overall.hetstat, common | random)
  chklogical(overall)
  chklogical(overall.hetstat)
  
  
  ##
  ##
  ## (3) Update trim-and-fill object
  ##
  ##
  if (inherits(object, "trimfill")) {
    ##
    rmfilled <- function(x) {
      ##
      if (!is.null(object[[x]]))
        res <- object[[x]][!object$trimfill]
      else
        res <- NULL
      ##
      res
    }
    ##
    tfnames <- c("TE", "seTE",
                 "studlab",
                 "n.e", "n.c",
                 "event.e", "event.c",
                 "mean.e", "mean.c", "sd.e", "sd.c",
                 "n", "event", "cor")
    ##
    for (i in tfnames)
      object[[i]] <- rmfilled(i)
    ##
    oldclass <- object$class.x
    ##
    object$class <- object$class[object$class != "trimfill"]
    #
    res <- trimfill(object,
                    left = left, ma.common = ma.common,
                    type = type, n.iter.max = n.iter.max,
                    level = level, level.ma = level.ma,
                    common = common, random = random,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    prediction = prediction | !missing.method.predict,
                    level.predict = level.predict,
                    silent = TRUE,
                    ...)
    ##
    res$call.object <- object$call
    res$call <- match.call()
    res$class.x <- oldclass
    ##
    return(res)
  }
  
  
  ##
  ##
  ## (4) Update metacum or metainf object
  ##
  ##
  if (inherits(object, "metacum") | inherits(object, "metainf")) {
    ##
    res <- object
    ##
    res$common <- ifelse(res$pooled == "common", TRUE, FALSE)
    res$random <- ifelse(res$pooled == "random", TRUE, FALSE)
    ##
    res$call.object <- object$call
    res$call <- match.call()
    res$version <- packageDescription("meta")$Version
    ##
    return(res)
  }
  ##
  if (is.null(object$data)) {
    warning("Necessary data not available. Please, recreate ",
            "meta-analysis object without option 'keepdata = FALSE'.")
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (5) Catch variables
  ##
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  ## Catch argument 'subset'
  ##
  if (!missing(subset))
    subset <- catch("subset", mc, data, sfsp)
  else
    subset <- catch2(object, "subset", fromobject = TRUE)
  ##
  ## Catch argument 'studlab'
  ##
  if (!missing(studlab))
    studlab <- catch("studlab", mc, data, sfsp)
  else
    studlab <- catch2(object, "studlab")
  ##
  ## Catch argument 'exclude'
  ##
  if (!missing(exclude))
    exclude <- catch("exclude", mc, data, sfsp)
  else
    exclude <- catch2(object, "exclude")
  ##
  ## Catch argument 'cluster'
  ##
  missing.cluster <- missing(cluster)
  missing.id <- missing(id)
  ##
  if (!missing.cluster | !missing.id) {
    cluster <- catch("cluster", mc, data, sfsp)
    id <- catch("id", mc, data, sfsp)
    ...cluster <-
      deprecated2(cluster, missing.cluster, id, missing.id,
                  warn.deprecated)
    ##
    data$.cluster <- ...cluster
  }
  else
    ...cluster <- catch2(object, "cluster")
  #
  # Catch argument 'cycles'
  #
  if (!missing(cycles))
    cycles <- catch("cycles", mc, data, sfsp)
  else {
    if (isCol(object$data, ".cycles"))
      cycles <- catch(".cycles", object$data, object$data, sfsp)
    else
      cycles <- NULL
  }
  #
  avail.cycles <- !is.null(cycles)
  ##
  ## Catch argument 'incr'
  ##
  missing.incr <- missing(incr)
  #
  # Catch argument 'incr.e'
  #
  missing.incr.e <- missing(incr.e)
  #
  if (!missing.incr.e)
    incr.e <- catch("incr.e", mc, data, sfsp)
  else
    incr.e <- catch2(object, "incr.e", gs("incr"))
  #
  avail.incr.e <- !missing.incr.e & !is.null(incr.e)
  #
  # Catch argument 'incr.c'
  #
  missing.incr.c <- missing(incr.c)
  #
  if (!missing.incr.c)
    incr.c <- catch("incr.c", mc, data, sfsp)
  else
    incr.c <- catch2(object, "incr.c", gs("incr"))
  #
  avail.incr.c <- !missing.incr.c & !is.null(incr.c)
  #
  # Catch argument 'n.e'
  #
  if (!missing(n.e))
    n.e <- catch("n.e", mc, data, sfsp)
  else
    n.e <- catch2(object, "n.e")
  #
  # Catch argument 'n.c'
  #
  if (!missing(n.c))
    n.c <- catch("n.c", mc, data, sfsp)
  else
    n.c <- catch2(object, "n.c")
  ##
  ## Catch argument 'approx.mean.e'
  ##
  if (!missing(approx.mean.e))
    approx.mean.e <- catch("approx.mean.e", mc, data, sfsp)
  else
    approx.mean.e <- catch2(object, "approx.mean.e")
  ##
  ## Catch argument 'approx.mean.c'
  ##
  if (!missing(approx.mean.c))
    approx.mean.c <- catch("approx.mean.c", mc, data, sfsp)
  else
    approx.mean.c <- catch2(object, "approx.mean.c")
  ##
  ## Catch argument 'approx.sd.e'
  ##
  if (!missing(approx.sd.e))
    approx.sd.e <- catch("approx.sd.e", mc, data, sfsp)
  else
    approx.sd.e <- catch2(object, "approx.sd.e")
  ##
  ## Catch argument 'approx.sd.c'
  ##
  if (!missing(approx.sd.c))
    approx.sd.c <- catch("approx.sd.c", mc, data, sfsp)
  else
    approx.sd.c <- catch2(object, "approx.sd.c")
  ##
  ## Catch argument 'subgroup'
  ##
  missing.subgroup <- missing(subgroup)
  missing.byvar <- missing(byvar)
  ##
  if (!missing.subgroup | !missing.byvar) {
    subgroup <- catch("subgroup", mc, data, sfsp)
    byvar <- catch("byvar", mc, data, sfsp)
    subgroup <-
      deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                  warn.deprecated)
    ##
    data$.subgroup <- subgroup
  }
  else if (isCol(object$data, ".subgroup"))
    subgroup <- object$data$.subgroup
  else
    subgroup <- NULL
  ##
  if (missing.subgroup.name & is.null(subgroup.name)) {
    if (!missing.subgroup)
      subgroup.name <- byvarname("subgroup", mc)
    else if (!missing.byvar)
      subgroup.name <- byvarname("byvar", mc)
  }
  ##
  missing.sm <- missing(sm)
  ##
  if (!is.null(subgroup.name))
    chkchar(subgroup.name, length = 1)
  #
  # Check arguments for common effect model
  #
  if (!missing.method.common.ci & missing.common)
    common <- TRUE
  ##
  ## Check arguments for random effects model(s)
  ##
  if (!(missing.method.random.ci & missing.adhoc.hakn.ci) & missing.random)
    random <- TRUE
  ##
  if (!missing.method.random.ci) {
    if (missing.adhoc.hakn.ci)
      adhoc.hakn.ci <- rep("", length(method.random.ci))
    else if (length(method.random.ci) != length(adhoc.hakn.ci))
      stop("Arguments 'method.random.ci' and 'adhoc.hakn.ci' must be of ",
           "same length.",
           call. = FALSE)
  }
  if (!missing.adhoc.hakn.ci) {
    if (missing.method.random.ci)
      missing.method.random.ci <- rep("HK", length(adhoc.hakn.ci))
    else if (length(method.random.ci) != length(adhoc.hakn.ci))
      stop("Arguments 'method.random.ci' and 'adhoc.hakn.ci' must be of ",
           "same length.",
           call. = FALSE)
  }
  ##
  if (!missing.method.random.ci | !missing.text.random) {
    if (length(method.random.ci) != length(text.random)) {
      if (!missing.method.random.ci) {
        if (!(length(text.random) == 1 && text.random == gs("text.random")))
          warning("Setting argument 'text.random' to default as number of ",
                  "random effects \n   methods changed by ",
                  "argument 'method.random.ci'.",
                  call. = FALSE)
        text.random <- gs("text.random")
      }
      if (!missing.text.random)
        stop("Argument 'text.random' must be of same length as \n   ",
             "number of random effects methods specified by setting for ",
             "'method.random.ci'.",
             call. = FALSE)
    }
  }
  ##
  ## Check variables for prediction interval(s)
  ##
  if (!missing.method.predict) {
    if (missing.adhoc.hakn.pi)
      adhoc.hakn.pi <- rep("", length(method.predict))
    else if (length(method.predict) != length(adhoc.hakn.pi))
      stop("Arguments 'method.predict' and 'adhoc.hakn.pi' must be of ",
           "same length.",
           call. = FALSE)
  }
  if (!missing.adhoc.hakn.pi) {
    if (missing.method.predict)
      missing.method.predict <- rep("HK", length(adhoc.hakn.pi))
    else if (length(method.predict) != length(adhoc.hakn.pi))
      stop("Arguments 'method.predict' and 'adhoc.hakn.pi' must be of ",
           "same length.",
           call. = FALSE)
  }
  ##
  if (!missing.method.predict | !missing.text.predict) {
    if (length(method.predict) != length(text.predict)) {
      if (!missing.method.predict) {
        if (!(length(text.predict) == 1 && text.predict == gs("text.predict")))
          warning("Setting argument 'text.predict' to default as number of ",
                  "prediction intervals \n   changed by ",
                  "argument 'method.predict'.",
                  call. = FALSE)
        text.predict <- gs("text.predict")
      }
      if (!missing.text.predict)
        stop("Argument 'text.predict' must be of same length as \n   ",
             "number of prediction intervals specified by setting for ",
             "'method.predict'.",
             call. = FALSE)
    }
  }
  
  
  ##
  ##
  ## (6) Update meta object
  ##
  ##
  method.predict <- replaceVal(method.predict, "", gs("method.predict"))
  ##
  if (metabin) {
    sm <- setchar(sm, gs("sm4bin"))
    method <- setchar(method, gs("meth4bin"))
    ##
    if (!is.null(...cluster))
      method <- "Inverse"
    ##
    if (method %in% c("GLMM", "LRP") & !missing.sm & sm != "OR")
      warning("Summary measure 'sm = \"OR\" used as 'method = \"",
              method, "\".")
    #
    if (method == "GLMM") {
      sm <- "OR"
      method.tau <- "ML"
      model.glmm <- replaceNULL(model.glmm, gs("model.glmm"))
    }
    #
    if (method == "LRP") {
      sm <- "OR"
      method.tau <- "DL"
    }
    #
    if (!(method == "MH" & method.tau == "DL" &
          (sm %in% c("OR", "RR", "RD", "DOR"))))
      Q.Cochrane <- FALSE
    #
    if (sm == "DOR" & missing.method.bias)
      method.bias <- "Deeks"
    else if (sm == "OR" & missing.method.bias)
      method.bias <- "Harbord"
    #
    if (!missing.method.incr)
      method.incr <- setchar(method.incr, gs("meth4incr"))
    else {
      if (!missing.incr & missing.incr.e & missing.incr.c)
        method.incr <- gs("method.incr")
    }
    #
    if (method.incr == "user") {
      incr.orig <- incr
      incr <- NULL
    }
    else {
      incr.e <- NULL
      incr.c <- NULL
      #
      if (sm == "ASD") {
        if (!missing.incr && any(incr != 0))
          warning("Note, no continuity correction considered for ",
                  "arcsine difference (sm = \"ASD\").",
                  call. = FALSE)
        incr <- 0
      }
      #
      if (method == "Peto") {
        if (!missing.incr && any(incr != 0))
          warning("Note, no continuity correction considered for ",
                  "method = \"Peto\".",
                  call. = FALSE)
        incr <- 0
      }
      #
      RR.Cochrane <- replaceNULL(RR.Cochrane, gs("RR.cochrane"))
    }
    #
    m <- metabin(event.e = object$data$.event.e,
                 n.e = object$data$.n.e,
                 event.c = object$data$.event.c,
                 n.c = object$data$.n.c,
                 studlab = studlab,
                 ##
                 data = data, subset = subset, exclude = exclude,
                 cluster = ...cluster, rho = rho,
                 #
                 weights = weights,
                 weights.common = weights.common,
                 weights.random = weights.random,
                 #
                 method = method,
                 sm = sm,
                 #
                 incr = incr, method.incr = method.incr,
                 allstudies = allstudies,
                 incr.e = incr.e, incr.c = incr.c,
                 #
                 MH.exact = MH.exact, RR.Cochrane = RR.Cochrane,
                 Q.Cochrane = Q.Cochrane, model.glmm = model.glmm,
                 ##
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 method.common.ci = method.common.ci,
                 method.random.ci = method.random.ci,
                 adhoc.hakn.ci = adhoc.hakn.ci,
                 method.predict = method.predict,
                 adhoc.hakn.pi = adhoc.hakn.pi,
                 seed.predict = seed.predict,
                 ##
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 level.hetstat = level.hetstat,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 detail.tau = detail.tau,
                 #
                 method.I2 = method.I2,
                 #
                 prediction = prediction | !missing.method.predict,
                 level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 backtransf = backtransf, pscale = pscale,
                 ##
                 text.common = text.common, text.random = text.random,
                 text.predict = text.predict,
                 text.w.common = text.w.common, text.w.random = text.w.random,
                 ##
                 title = title, complab = complab, outclab = outclab,
                 #
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 col.label.right = col.label.right,
                 col.label.left = col.label.left,
                 #
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 seed.predict.subgroup = seed.predict.subgroup,
                 print.CMH = print.CMH,
                 ##
                 warn = warn, warn.deprecated = FALSE,
                 ##
                 control = control,
                 ...)
    #
    if (method.incr == "user")
      m$incr <- incr.orig
  }
  ##
  if (metacont) {
    if (!isCol(object$data, ".approx.mean.e"))
      mean.e <- object$data$.mean.e
    else
      mean.e <-
        setNA_ifnot(object$data$.mean.e, object$data$.approx.mean.e, "")
    ##
    if (!isCol(object$data, ".approx.mean.c"))
      mean.c <- object$data$.mean.c
    else
      mean.c <-
        setNA_ifnot(object$data$.mean.c, object$data$.approx.mean.c, "")
    ##
    if (!isCol(object$data, ".approx.sd.e"))
      sd.e <- object$data$.sd.e
    else
      sd.e <-
        setNA_ifnot(object$data$.sd.e, object$data$.approx.sd.e, "")
    ##
    if (!isCol(object$data, ".approx.sd.c"))
      sd.c <- object$data$.sd.c
    else
      sd.c <-
        setNA_ifnot(object$data$.sd.c, object$data$.approx.sd.c, "")
    ##
    m <- metacont(n.e = object$data$.n.e,
                  mean.e = object$data$.mean.e,
                  sd.e = object$data$.sd.e,
                  n.c = object$data$.n.c,
                  mean.c = object$data$.mean.c,
                  sd.c = object$data$.sd.c,
                  studlab = studlab,
                  ##
                  data = data, subset = subset, exclude = exclude,
                  cluster = ...cluster, rho = rho,
                  #
                  weights = weights,
                  weights.common = weights.common,
                  weights.random = weights.random,
                  #
                  median.e = setVal(object$data, ".median.e"),
                  q1.e = setVal(object$data, ".q1.e"),
                  q3.e = setVal(object$data, ".q3.e"),
                  min.e = setVal(object$data, ".min.e"),
                  max.e = setVal(object$data, ".max.e"),
                  ##
                  median.c = setVal(object$data, ".median.c"),
                  q1.c = setVal(object$data, ".q1.c"),
                  q3.c = setVal(object$data, ".q3.c"),
                  min.c = setVal(object$data, ".min.c"),
                  max.c = setVal(object$data, ".max.c"),
                  ##
                  method.mean =
                    replaceVal(replaceNULL(method.mean, "Luo"), "", "Luo"),
                  method.sd =
                    replaceVal(replaceNULL(method.sd, "Shi"), "", "Shi"),
                  #
                  approx.mean.e = approx.mean.e,
                  approx.mean.c = approx.mean.c,
                  approx.sd.e = approx.sd.e,
                  approx.sd.c = approx.sd.c,
                  ##
                  sm = sm,
                  pooledvar = replaceNA(pooledvar, gs("pooledvar")),
                  method.smd = replaceVal(method.smd, "", gs("method.smd")),
                  sd.glass = replaceVal(sd.glass, "", gs("sd.glass")),
                  exact.smd = replaceNA(exact.smd, gs("exact.smd")),
                  ##
                  method.ci = replaceNULL(method.ci, gs("method.ci.cont")),
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  method.common.ci = method.common.ci,
                  method.random.ci = method.random.ci,
                  adhoc.hakn.ci = adhoc.hakn.ci,
                  method.predict = method.predict,
                  adhoc.hakn.pi = adhoc.hakn.pi,
                  seed.predict = seed.predict,
                  ##
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  level.hetstat = level.hetstat,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  detail.tau = detail.tau,
                  #
                  method.I2 = method.I2,
                  #
                  prediction = prediction | !missing.method.predict,
                  level.predict = level.predict,
                  ##
                  method.bias = method.bias,
                  ##
                  text.common = text.common, text.random = text.random,
                  text.predict = text.predict,
                  text.w.common = text.w.common, text.w.random = text.w.random,
                  ##
                  title = title, complab = complab, outclab = outclab,
                  #
                  label.e = label.e, label.c = label.c,
                  label.right = label.right, label.left = label.left,
                  col.label.right = col.label.right,
                  col.label.left = col.label.left,
                  #
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  seed.predict.subgroup = seed.predict.subgroup,
                  ##
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control)
  }
  ##
  if (metacor)
    m <- metacor(cor = object$data$.cor,
                 n = object$data$.n,
                 studlab = studlab,
                 ##
                 data = data, subset = subset, exclude = exclude,
                 cluster = ...cluster, rho = rho,
                 #
                 weights = weights,
                 weights.common = weights.common,
                 weights.random = weights.random,
                 #
                 sm = sm,
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 method.common.ci = method.common.ci,
                 method.random.ci = method.random.ci,
                 adhoc.hakn.ci = adhoc.hakn.ci,
                 method.predict = method.predict,
                 adhoc.hakn.pi = adhoc.hakn.pi,
                 seed.predict = seed.predict,
                 ##
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 level.hetstat = level.hetstat,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 detail.tau = detail.tau,
                 #
                 method.I2 = method.I2,
                 #
                 prediction = prediction | !missing.method.predict,
                 level.predict = level.predict,
                 ##
                 null.effect = null.effect,
                 ##
                 method.bias = method.bias,
                 ##
                 backtransf = backtransf,
                 ##
                 text.common = text.common, text.random = text.random,
                 text.predict = text.predict,
                 text.w.common = text.w.common, text.w.random = text.w.random,
                 ##
                 title = title, complab = complab, outclab = outclab,
                 #
                 label.right = label.right, label.left = label.left,
                 col.label.right = col.label.right,
                 col.label.left = col.label.left,
                 #
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 seed.predict.subgroup = seed.predict.subgroup,
                 ##
                 warn.deprecated = FALSE,
                 ##
                 control = control)
  ##
  if (metagen) {
    ...n.e <- n.e
    ...n.c <- n.c
    #
    rm(n.e)
    rm(n.c)
    #
    if (missing(approx.TE)) {
      if (isCol(object$data, ".approx.TE"))
        approx.TE <- object$data$.approx.TE
      else
        approx.TE <- NULL
    }
    if (missing(approx.seTE)) {
      if (isCol(object$data, ".approx.seTE"))
        approx.seTE <- object$data$.approx.seTE
      else
        approx.seTE <- NULL
    }
    #
    m <- metagen(TE = object$data$.TE,
                 seTE =
                   if (!avail.cycles) object$data$.seTE
                   else replaceNULL(object$data$.seTE.orig, object$data$.seTE),
                 studlab = studlab,
                 ##
                 data = data, subset = subset, exclude = exclude,
                 cluster = ...cluster, rho = rho,
                 #
                 cycles = cycles,
                 #
                 weights = weights,
                 weights.common = weights.common,
                 weights.random = weights.random,
                 #
                 sm = sm,
                 method.ci = method.ci,
                 level = level,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 level.ma = level.ma,
                 method.common.ci = method.common.ci,
                 method.random.ci = method.random.ci,
                 adhoc.hakn.ci = adhoc.hakn.ci,
                 method.predict = method.predict,
                 adhoc.hakn.pi = adhoc.hakn.pi,
                 seed.predict = seed.predict,
                 ##
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 level.hetstat = level.hetstat,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 detail.tau = detail.tau,
                 #
                 method.I2 = method.I2,
                 #
                 prediction = prediction | !missing.method.predict,
                 level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 n.e = ...n.e, n.c = ...n.c,
                 ##
                 pval = setVal(object$data, ".pval"),
                 df = setVal(object$data, ".df"),
                 lower = setVal(object$data, ".lower"),
                 upper = setVal(object$data, ".upper"),
                 level.ci = setVal(object$data, ".level.ci"),
                 ##
                 median = setVal(object$data, ".median"),
                 q1 = setVal(object$data, ".q1"),
                 q3 = setVal(object$data, ".q3"),
                 min = setVal(object$data, ".min"),
                 max = setVal(object$data, ".max"),
                 ##
                 method.mean =
                   replaceVal(replaceNULL(method.mean, "Luo"), "", "Luo"),
                 method.sd =
                   replaceVal(replaceNULL(method.sd, "Shi"), "", "Shi"),
                 ##
                 approx.TE = approx.TE,
                 approx.seTE = approx.seTE,
                 ##
                 transf = TRUE,
                 backtransf = backtransf,
                 func.transf = object$func.transf,
                 func.backtransf = func.backtransf,
                 args.transf = object$args.transf,
                 args.backtransf = args.backtransf,         
                 #
                 pscale = pscale,
                 irscale = irscale, irunit = irunit,
                 ##
                 text.common = text.common, text.random = text.random,
                 text.predict = text.predict,
                 text.w.common = text.w.common, text.w.random = text.w.random,
                 ##
                 title = title, complab = complab, outclab = outclab,
                 #
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 col.label.right = col.label.right,
                 col.label.left = col.label.left,
                 #
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 seed.predict.subgroup = seed.predict.subgroup,
                 ##
                 warn = warn, warn.deprecated = FALSE,
                 ##
                 control = control)
  }
  ##
  if (metainc) {
    ...n.e <- n.e
    ...n.c <- n.c
    #
    rm(n.e)
    rm(n.c)
    #
    sm <- setchar(sm, gs("sm4inc"))
    method <- setchar(method, gs("meth4inc"))
    ##
    if (!is.null(...cluster))
      method <- "Inverse"
    ##
    if (method == "GLMM" & !missing.sm & !(sm %in% c("IRR", "VE")))
      warning("Summary measure 'sm = \"IRR\" used as 'method = \"GLMM\".")
    ##
    if (method == "GLMM") {
      if (sm != "VE")
        sm <- "IRR"
      method.tau <- "ML"
      model.glmm <- replaceNULL(model.glmm, "UM.FS")
    }
    #
    if (!missing.method.incr)
      method.incr <- setchar(method.incr, gs("meth4incr"))
    else {
      if (!missing.incr & missing.incr.e & missing.incr.c)
        method.incr <- gs("method.incr")
    }
    #
    if (method.incr == "user") {
      incr.orig <- incr
      incr <- NULL
    }
    else {
      incr.e <- NULL
      incr.c <- NULL
    }
    #
    m <- metainc(event.e = object$data$.event.e,
                 time.e = object$data$.time.e,
                 event.c = object$data$.event.c,
                 time.c = object$data$.time.c,
                 studlab = studlab,
                 ##
                 data = data, subset = subset, exclude = exclude,
                 cluster = ...cluster, rho = rho,
                 #
                 weights = weights,
                 weights.common = weights.common,
                 weights.random = weights.random,
                 #
                 method = method,
                 sm = sm,
                 #
                 incr = incr, method.incr = method.incr,
                 incr.e = incr.e, incr.c = incr.c,
                 #
                 model.glmm = model.glmm,
                 ##
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 method.common.ci = method.common.ci,
                 method.random.ci = method.random.ci,
                 adhoc.hakn.ci = adhoc.hakn.ci,
                 method.predict = method.predict,
                 adhoc.hakn.pi = adhoc.hakn.pi,
                 seed.predict = seed.predict,
                 ##
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 level.hetstat = level.hetstat,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 detail.tau = detail.tau,
                 #
                 method.I2 = method.I2,
                 #
                 prediction = prediction | !missing.method.predict,
                 level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 n.e = ...n.e, n.c = ...n.c,
                 ##
                 backtransf = backtransf, irscale = irscale, irunit = irunit,
                 ##
                 text.common = text.common, text.random = text.random,
                 text.predict = text.predict,
                 text.w.common = text.w.common, text.w.random = text.w.random,
                 ##
                 title = title, complab = complab, outclab = outclab,
                 #
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 col.label.right = col.label.right,
                 col.label.left = col.label.left,
                 #
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 seed.predict.subgroup = seed.predict.subgroup,
                 ##
                 warn = warn, warn.deprecated = FALSE,
                 ##
                 control = control,
                 ...)
    #
    if (method.incr == "user")
      m$incr <- incr.orig
  }
  ##
  if (metamean)
    m <- metamean(n = object$data$.n,
                  mean = object$data$.mean,
                  sd = object$data$.sd,
                  studlab = studlab,
                  ##
                  data = data, subset = subset, exclude = exclude,
                  cluster = ...cluster, rho = rho,
                  #
                  weights = weights,
                  weights.common = weights.common,
                  weights.random = weights.random,
                  #
                  median = setVal(object$data, ".median"),
                  q1 = setVal(object$data, ".q1"),
                  q3 = setVal(object$data, ".q3"),
                  min = setVal(object$data, ".min"),
                  max = setVal(object$data, ".max"),
                  ##
                  method.mean = replaceVal(method.mean, "", "Luo"),
                  method.sd = replaceVal(method.sd, "", "Shi"),
                  ##
                  approx.mean = approx.mean,
                  approx.sd = approx.sd,
                  ##
                  sm = sm,
                  method.ci = replaceNULL(method.ci, gs("method.ci.cont")),
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  method.common.ci = method.common.ci,
                  method.random.ci = method.random.ci,
                  adhoc.hakn.ci = adhoc.hakn.ci,
                  method.predict = method.predict,
                  adhoc.hakn.pi = adhoc.hakn.pi,
                  seed.predict = seed.predict,
                  ##
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  level.hetstat = level.hetstat,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  detail.tau = detail.tau,
                  #
                  method.I2 = method.I2,
                  #
                  prediction = prediction | !missing.method.predict,
                  level.predict = level.predict,
                  ##
                  null.effect = null.effect,
                  ##
                  method.bias = method.bias,
                  ##
                  backtransf = backtransf,
                  ##
                  text.common = text.common, text.random = text.random,
                  text.predict = text.predict,
                  text.w.common = text.w.common, text.w.random = text.w.random,
                  ##
                  title = title, complab = complab, outclab = outclab,
                  #
                  label.right = label.right, label.left = label.left,
                  col.label.right = col.label.right,
                  col.label.left = col.label.left,
                  #
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  seed.predict.subgroup = seed.predict.subgroup,
                  ##
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control)
  ##
  if (metaprop) {
    sm <- setchar(sm, gs("sm4prop"))
    method <- setchar(method, gs("meth4prop"))
    ##
    if (!is.null(...cluster))
      method <- "Inverse"
    if (method == "GLMM" & !missing.sm & sm != "PLOGIT")
      warning("Summary measure 'sm = \"PLOGIT\" used as 'method = \"GLMM\".")
    ##
    if (method == "GLMM") {
      sm <- "PLOGIT"
      method.tau <- "ML"
    }
    ##
    method.ci <- replaceNULL(method.ci, gs("method.ci.prop"))
    ##
    m <- metaprop(event = object$data$.event,
                  n = object$data$.n,
                  studlab = studlab,
                  ##
                  data = data, subset = subset, exclude = exclude,
                  cluster = ...cluster, rho = rho,
                  #
                  weights = weights,
                  weights.common = weights.common,
                  weights.random = weights.random,
                  #
                  method = method,
                  sm = sm,
                  incr = incr,
                  method.incr = method.incr,
                  ##
                  method.ci = method.ci,
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  method.common.ci = method.common.ci,
                  method.random.ci = method.random.ci,
                  adhoc.hakn.ci = adhoc.hakn.ci,
                  method.predict = method.predict,
                  adhoc.hakn.pi = adhoc.hakn.pi,
                  seed.predict = seed.predict,
                  ##
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  level.hetstat = level.hetstat,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  detail.tau = detail.tau,
                  #
                  method.I2 = method.I2,
                  #
                  prediction = prediction | !missing.method.predict,
                  level.predict = level.predict,
                  ##
                  null.effect = null.effect,
                  ##
                  method.bias = method.bias,
                  ##
                  backtransf = backtransf, pscale = pscale,
                  ##
                  text.common = text.common, text.random = text.random,
                  text.predict = text.predict,
                  text.w.common = text.w.common, text.w.random = text.w.random,
                  ##
                  title = title, complab = complab, outclab = outclab,
                  #
                  label.right = label.right, label.left = label.left,
                  col.label.right = col.label.right,
                  col.label.left = col.label.left,
                  #
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  seed.predict.subgroup = seed.predict.subgroup,
                  ##
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control,
                  ...)
  }
  ##
  if (metarate) {
    sm <- setchar(sm, gs("sm4rate"))
    method <- setchar(method, gs("meth4rate"))
    ##
    if (!is.null(...cluster))
      method <- "Inverse"
    ##
    if (method == "GLMM" & !missing.sm & sm != "IRLN")
      warning("Summary measure 'sm = \"IRLN\" used as 'method = \"GLMM\".")
    ##
    if (method == "GLMM") {
      sm <- "IRLN"
      method.tau <- "ML"
    }
    ##
    method.ci <- replaceNULL(method.ci, gs("method.ci.rate"))
    ##
    m <- metarate(event = object$data$.event,
                  time = object$data$.time,
                  studlab = studlab,
                  ##
                  data = data, subset = subset, exclude = exclude,
                  cluster = ...cluster, rho = rho,
                  #
                  weights = weights,
                  weights.common = weights.common,
                  weights.random = weights.random,
                  #
                  method = method,
                  sm = sm,
                  incr = incr,
                  method.incr = method.incr,
                  ##
                  method.ci = method.ci,
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  method.common.ci = method.common.ci,
                  method.random.ci = method.random.ci,
                  adhoc.hakn.ci = adhoc.hakn.ci,
                  method.predict = method.predict,
                  adhoc.hakn.pi = adhoc.hakn.pi,
                  seed.predict = seed.predict,
                  ##
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  level.hetstat = level.hetstat,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  detail.tau = detail.tau,
                  #
                  method.I2 = method.I2,
                  #
                  prediction = prediction | !missing.method.predict,
                  level.predict = level.predict,
                  ##
                  null.effect = null.effect,
                  ##
                  method.bias = method.bias,
                  ##
                  backtransf = backtransf, irscale = irscale, irunit = irunit,
                  ##
                  text.common = text.common, text.random = text.random,
                  text.predict = text.predict,
                  text.w.common = text.w.common, text.w.random = text.w.random,
                  ##
                  title = title, complab = complab, outclab = outclab,
                  #
                  label.right = label.right, label.left = label.left,
                  col.label.right = col.label.right,
                  col.label.left = col.label.left,
                  #
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  seed.predict.subgroup = seed.predict.subgroup,
                  ##
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control,
                  ...)
  }
  ##
  m$call.object <- object$call
  m$call <- match.call()
  ##
  if (!is.null(object$rob)) {
    if (!is.null(m$subset))
      m$rob <- object$rob[m$subset, ]
    else
      m$rob <- object$rob
  }
  ##
  if (!keepdata)
    m$data <- NULL
  
  m
}
