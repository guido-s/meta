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
#' @param method A character string indicating which method is to be
#'   used for pooling of studies; see \code{\link{metabin}} and
#'   \code{\link{metainc}} function for admissible values.
#' @param sm A character string indicating which summary measure is
#'   used for pooling.
#' @param incr Either a numerical value or vector which can be added
#'   to each cell frequency for studies with a zero cell count or the
#'   character string \code{"TA"} which stands for treatment arm
#'   continuity correction.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, or \code{"all"}).
#' @param allstudies A logical indicating if studies with zero or all
#'   events in both groups are to be included in the meta-analysis
#'   (applies only if \code{sm} is equal to \code{"RR"} or
#'   \code{"OR"}).
#' @param MH.exact A logical indicating if \code{incr} is not to be
#'   added to all cell frequencies for studies with a zero cell count
#'   to calculate the pooled estimate based on the Mantel-Haenszel
#'   method.
#' @param RR.Cochrane A logical indicating if 2*\code{incr} instead of
#'   1*\code{incr} is to be added to \code{n.e} and \code{n.c} in the
#'   calculation of the risk ratio (i.e., \code{sm="RR"}) for studies
#'   with a zero cell. This is used in RevMan 5, the program for
#'   preparing and maintaining Cochrane reviews.
#' @param Q.Cochrane A logical indicating if the Mantel-Haenszel
#'   estimate is used in the calculation of the heterogeneity
#'   statistic Q which is implemented in RevMan 5, the program for
#'   preparing and maintaining Cochrane reviews.
#' @param model.glmm A character string indicating which GLMM model
#'   should be used.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.ma The level used to calculate confidence intervals
#'   for pooled estimates.
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
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param adhoc.hakn A character string indicating whether an \emph{ad
#'   hoc} variance correction should be applied in the case of an
#'   arbitrarily small Hartung-Knapp variance estimate.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"PM"},
#'   \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
#'   \code{"HE"}, or \code{"EB"}, can be abbreviated. See function
#'   \code{\link{metagen}}.
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau}. Either \code{"QP"}, \code{"BJ"}, or \code{"J"}, can
#'   be abbreviated.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance \eqn{\tau^2}.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param method.bias A character string indicating which test for
#'   funnel plot asymmetry is to be used, can be abbreviated. See
#'   function \code{\link{metabias}.}
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE}, results for \code{sm = "OR"} are printed as odds ratios
#'   rather than log odds ratios and results for \code{sm = "ZCOR"}
#'   are printed as correlations rather than Fisher's z transformed
#'   correlations, for example.
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
#' @param label.left Graph label on left side of forest plot.
#' @param label.right Graph label on right side of forest plot.
#' @param n.e Number of observations in experimental group. (only for
#'   metagen object)
#' @param n.c Number of observations in control group. (only for
#'   metagen object)
#' @param pooledvar A logical indicating if a pooled variance should
#'   be used for the mean difference (only for metacont object with
#'   \code{sm = "MD"}).
#' @param method.smd A character string indicating which method is
#'   used to estimate the standardised mean difference (only for
#'   metacont object with \code{sm = "SMD"}). Either \code{"Hedges"}
#'   for Hedges' g (default), \code{"Cohen"} for Cohen's d, or
#'   \code{"Glass"} for Glass' delta, can be abbreviated.
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
#' \code{"metarate"}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
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
#' @export update.meta


update.meta <- function(object, 
                        data = object$data,
                        subset, studlab, exclude, cluster,
                        ##
                        method = object$method,
                        sm = object$sm,
                        incr,
                        method.incr = object$method.incr,
                        allstudies = object$allstudies,
                        MH.exact = object$MH.exact,
                        RR.Cochrane = object$RR.Cochrane,
                        Q.Cochrane = object$Q.Cochrane,
                        model.glmm = object$model.glmm,
                        level = object$level,
                        level.ma = object$level.ma,
                        common = object$common,
                        random = object$random,
                        overall = object$overall,
                        overall.hetstat = object$overall.hetstat,
                        hakn = object$hakn,
                        adhoc.hakn = object$adhoc.hakn,
                        method.tau = object$method.tau,
                        method.tau.ci = object$method.tau.ci,
                        tau.preset = object$tau.preset,
                        TE.tau = object$TE.tau,
                        tau.common = object$tau.common,
                        prediction = object$prediction,
                        level.predict = object$level.predict,
                        null.effect = object$null.effect,
                        method.bias = object$method.bias,
                        ##
                        backtransf = object$backtransf,
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
                        label.e = object$label.e,
                        label.c = object$label.c,
                        label.left = object$label.left,
                        label.right = object$label.right,
                        n.e = object$n.e,
                        n.c = object$n.c,
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
    if (!is.null(object$byvar)) {
      object$data$.subgroup <- object$byvar
      object$subgroup.name <- object$bylab
      object$print.subgroup.name <- object$print.byvar
      object$sep.subgroup <- object$byseparator
    }
  }
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
  
  
  ##
  ##
  ## (2) Check arguments
  ##
  ##
  if (!backtransf & pscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!backtransf & irscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
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
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  missing.subgroup.name <- missing(subgroup.name)
  subgroup.name <-
    deprecated(subgroup.name, missing.subgroup.name, args, "bylab",
               warn.deprecated)
  ##
  print.subgroup.name <-
    deprecated(print.subgroup.name, missing(print.subgroup.name),
               args, "print.byvar", warn.deprecated)
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
  addincr <-
    deprecated(method.incr, missing.method.incr, args, "addincr",
               warn.deprecated)
  allincr <-
    deprecated(method.incr, missing.method.incr, args, "allincr",
               warn.deprecated)
  if (missing.method.incr) {
    if (is.logical(addincr) && addincr)
      method.incr <- "all"
    else if (is.logical(allincr) && allincr)
      method.incr <- "if0all"
  }
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
    res <- trimfill(object,
                    left = left, ma.common = ma.common,
                    type = type, n.iter.max = n.iter.max,
                    level = level, level.ma = level.ma,
                    common = common, random = random,
                    hakn = hakn, adhoc.hakn = adhoc.hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    prediction = prediction, level.predict = level.predict,
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
  missing.subset  <- missing(subset)
  ##
  if (!missing.subset)
    subset <- catch("subset", mc, data, sfsp)
  else {
    if (!is.null(object$subset))
      subset <- object$subset
    else if (isCol(object$data, ".subset"))
      subset <- object$data$.subset
    else
      subset <- NULL
  }
  ##
  ## Catch argument 'studlab'
  ##
  missing.studlab <- missing(studlab)
  ##
  if (!missing.studlab)
    studlab <- catch("studlab", mc, data, sfsp)
  else if (isCol(object$data, ".studlab"))
    studlab <- object$data$.studlab
  else
    studlab <- NULL
  ##
  ## Catch argument 'exclude'
  ##
  missing.exclude <- missing(exclude)
  ##
  if (!missing.exclude)
    exclude <- catch("exclude", mc, data, sfsp)
  else if (isCol(object$data, ".exclude"))
    exclude <- object$data$.exclude
  else
    exclude <- NULL
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
  else if (isCol(object$data, ".cluster"))
    ...cluster <- object$data$.cluster
  else
    ...cluster <- NULL
  ##
  ## Catch argument 'incr'
  ##
  missing.incr <- missing(incr)
  ##
  if (!missing.incr)
    incr <- catch("incr", mc, data, sfsp)
  else {
    if (isCol(object$data, ".incr"))
      incr <- object$data$.incr
    else
      incr <- gs("incr")
  }
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
  
  
  ##
  ##
  ## (6) Update meta object
  ##
  ##
  if (metabin) {
    sm <- setchar(sm, gs("sm4bin"))
    method <- setchar(method, gs("meth4bin"))
    ##
    if (!is.null(...cluster))
      method <- "Inverse"
    ##
    if (method == "GLMM" & !missing.sm & sm != "OR")
      warning("Summary measure 'sm = \"OR\" used as 'method = \"GLMM\".")
    ##
    if (sm == "ASD") {
      if (!missing.incr && any(incr != 0))
        warning("Note, no continuity correction considered for ",
                "arcsine difference (sm = \"ASD\").",
                call. = FALSE)
      incr <- 0
      object$data$.incr <- 0
    }
    ##
    if (method == "Peto") {
      if (!missing.incr && any(incr != 0))
        warning("Note, no continuity correction considered for ",
                "method = \"Peto\".",
                call. = FALSE)
      incr <- 0
      object$data$.incr <- 0
    }
    ##
    if (method == "GLMM") {
      sm <- "OR"
      method.tau <- "ML"
      model.glmm <- replaceNULL(model.glmm, gs("model.glmm"))
    }
    ##
    RR.Cochrane <- replaceNULL(RR.Cochrane, gs("RR.cochrane"))
    ##
    if (method != "MH" |
        method.tau != "DL" |
        !(sm %in% c("OR", "RR", "RD", "DOR")))
      Q.Cochrane <- FALSE
    ##
    m <- metabin(event.e = object$data$.event.e,
                 n.e = object$data$.n.e,
                 event.c = object$data$.event.c,
                 n.c = object$data$.n.c,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 cluster = ...cluster,
                 ##
                 data = data, subset = subset,
                 ##
                 method = method,
                 sm = sm,
                 incr = incr,
                 method.incr = method.incr,
                 allstudies = allstudies,
                 MH.exact = MH.exact, RR.Cochrane = RR.Cochrane,
                 Q.Cochrane = Q.Cochrane, model.glmm = model.glmm,
                 ##
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn, adhoc.hakn = adhoc.hakn,
                 method.tau = method.tau,
                 method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
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
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 ##
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 print.CMH = print.CMH,
                 ##
                 keepdata = keepdata,
                 warn = warn, warn.deprecated = FALSE,
                 ##
                 control = control,
                 ...)
  }
  ##
  if (metacont) {
    method.ci <- replaceNULL(method.ci, gs("method.ci.cont"))
    ##
    m <- metacont(n.e = object$data$.n.e,
                  mean.e = object$data$.mean.e,
                  sd.e = object$data$.sd.e,
                  n.c = object$data$.n.c,
                  mean.c = object$data$.mean.c,
                  sd.c = object$data$.sd.c,
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  cluster = ...cluster,
                  ##
                  data = data, subset = subset,
                  ##
                  sm = sm, pooledvar = pooledvar,
                  method.smd = method.smd, sd.glass = sd.glass,
                  exact.smd = exact.smd,
                  ##
                  method.ci = method.ci,
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn, adhoc.hakn = adhoc.hakn,
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
                  ##
                  method.bias = method.bias,
                  ##
                  text.common = text.common, text.random = text.random,
                  text.predict = text.predict,
                  text.w.common = text.w.common, text.w.random = text.w.random,
                  ##
                  title = title, complab = complab, outclab = outclab,
                  label.e = label.e, label.c = label.c,
                  label.right = label.right, label.left = label.left,
                  ##
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  ##
                  keepdata = keepdata,
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control)
  }
  ##
  if (metacor)
    m <- metacor(cor = object$data$.cor,
                 n = object$data$.n,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 cluster = ...cluster,
                 ##
                 data = data, subset = subset,
                 ##
                 sm = sm,
                 ##
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn, adhoc.hakn = adhoc.hakn,
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
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
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 ##
                 keepdata = keepdata,
                 warn.deprecated = FALSE,
                 ##
                 control = control)
  ##
  if (metagen) {
    data.m <- data
    add.e <- FALSE
    add.c <- FALSE
    ##
    if ("n.e" %in% names(data)) {
      add.e <- TRUE
      data.m <- data.m[, names(data.m) != "n.e"]
    }
    if ("n.c" %in% names(data)) {
      add.c <- TRUE
      data.m <- data.m[, names(data.m) != "n.c"]
    }
    ##
    m <- metagen(TE = object$data$.TE,
                 seTE = object$data$.seTE,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 cluster = ...cluster,
                 ##
                 data = data.m, subset = subset,
                 ##
                 sm = sm,
                 ##
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn, adhoc.hakn = adhoc.hakn,
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 n.e = n.e, n.c = n.c,
                 ##
                 backtransf = backtransf, pscale = pscale,
                 irscale = irscale, irunit = irunit,
                 ##
                 text.common = text.common, text.random = text.random,
                 text.predict = text.predict,
                 text.w.common = text.w.common, text.w.random = text.w.random,
                 ##
                 title = title, complab = complab, outclab = outclab,
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 ##
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 ##
                 keepdata = keepdata,
                 warn = warn, warn.deprecated = FALSE,
                 ##
                 control = control)
    if (add.e)
      m$data$n.e <- data$n.e
    if (add.c)
      m$data$n.c <- data$n.c
    if (add.e | add.c)
      m$data <- m$data[, names(data)]
  }
  ##
  if (metainc) {
    sm <- setchar(sm, gs("sm4inc"))
    method <- setchar(method, gs("meth4inc"))
    ##
    if (!is.null(...cluster))
      method <- "Inverse"
    ##
    if (method == "GLMM" & !missing.sm & sm != "IRR")
      warning("Summary measure 'sm = \"IRR\" used as 'method = \"GLMM\".")
    ##
    data.m <- data
    add.e <- FALSE
    add.c <- FALSE
    ##
    if ("n.e" %in% names(data)) {
      add.e <- TRUE
      data.m <- data.m[, names(data.m) != "n.e"]
    }
    if ("n.c" %in% names(data)) {
      add.c <- TRUE
      data.m <- data.m[, names(data.m) != "n.c"]
    }
    ##
    if (method == "GLMM") {
      sm <- "IRR"
      method.tau <- "ML"
      model.glmm <- replaceNULL(model.glmm, "UM.FS")
    }
    ##
    m <- metainc(event.e = object$data$.event.e,
                 time.e = object$data$.time.e,
                 event.c = object$data$.event.c,
                 time.c = object$data$.time.c,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 cluster = ...cluster,
                 ##
                 data = data, subset = subset,
                 ##
                 method = method,
                 sm = sm,
                 incr = incr,
                 method.incr = method.incr,
                 model.glmm = model.glmm,
                 ##
                 level = level, level.ma = level.ma,
                 common = common, random = random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn, adhoc.hakn = adhoc.hakn,
                 method.tau = method.tau,
                 method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 n.e = n.e, n.c = n.c,
                 ##
                 backtransf = backtransf, irscale = irscale, irunit = irunit,
                 ##
                 text.common = text.common, text.random = text.random,
                 text.predict = text.predict,
                 text.w.common = text.w.common, text.w.random = text.w.random,
                 ##
                 title = title, complab = complab, outclab = outclab,
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 ##
                 subgroup = subgroup, subgroup.name = subgroup.name,
                 print.subgroup.name = print.subgroup.name,
                 sep.subgroup = sep.subgroup,
                 test.subgroup = test.subgroup,
                 prediction.subgroup = prediction.subgroup,
                 ##
                 keepdata = keepdata,
                 warn = warn, warn.deprecated = FALSE,
                 ##
                 control = control,
                 ...)
    if (add.e)
      m$data$n.e <- data$n.e
    if (add.c)
      m$data$n.c <- data$n.c
    if (add.e | add.c)
      m$data <- m$data[, names(data)]
  }
  ##
  if (metamean) {
    method.ci <- replaceNULL(method.ci, gs("method.ci.cont"))
    ##
    m <- metamean(n = object$data$.n,
                  mean = object$data$.mean,
                  sd = object$data$.sd,
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  cluster = ...cluster,
                  ##
                  data = data, subset = subset,
                  ##
                  sm = sm,
                  ##
                  method.ci = method.ci,
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn, adhoc.hakn = adhoc.hakn,
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
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
                  ##
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  ##
                  keepdata = keepdata,
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control)
  }
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
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  cluster = ...cluster,
                  ##
                  data = data, subset = subset,
                  ##
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
                  hakn = hakn, adhoc.hakn = adhoc.hakn,
                  method.tau = method.tau,
                  method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
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
                  ##
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  ##
                  keepdata = keepdata,
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
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  cluster = ...cluster,
                  ##
                  data = data, subset = subset, method = method,
                  ##
                  sm = sm,
                  incr = incr,
                  method.incr = method.incr,
                  ##
                  method.ci = method.ci,
                  level = level, level.ma = level.ma,
                  common = common, random = random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn, adhoc.hakn = adhoc.hakn,
                  method.tau = method.tau,
                  method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
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
                  subgroup = subgroup, subgroup.name = subgroup.name,
                  print.subgroup.name = print.subgroup.name,
                  sep.subgroup = sep.subgroup,
                  test.subgroup = test.subgroup,
                  prediction.subgroup = prediction.subgroup,
                  ##
                  keepdata = keepdata,
                  warn = warn, warn.deprecated = FALSE,
                  ##
                  control = control,
                  ...)
  }
  ##  
  m$call.object <- object$call
  m$call <- match.call()
  
  
  m
}
