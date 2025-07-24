#' Meta-analysis of single incidence rates
#' 
#' @description
#' Calculation of an overall incidence rate from studies reporting a
#' single incidence rate. Inverse variance method and generalised
#' linear mixed model (GLMM) are available for pooling. For GLMMs, the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} (Viechtbauer 2010) is called internally.
#' 
#' @param event Number of events.
#' @param time Person time at risk.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., event and time.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param cluster An optional vector specifying which estimates come
#'   from the same cluster resulting in the use of a three-level
#'   meta-analysis model.
#' @param rho Assumed correlation of estimates within a cluster.
#' @param weights A single numeric or vector with user-specified weights.
#' @param weights.common User-specified weights (common effect model).
#' @param weights.random User-specified weights (random effects model).
#' @param n Number of observations.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"} and
#'   \code{"GLMM"}, can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"IR"}, \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}) is
#'   to be used for pooling of studies, see Details.
#' @param incr A numeric which is added to the event number of studies
#'   with zero events, i.e., studies with an incidence rate of 0. Or a
#'   numeric vector with the continuity correction for each study.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, or \code{"all"}), see Details.
#' @param method.ci A character string indicating whether to use
#'   approximate normal ("NAsm") or exact Poisson ("Poisson")
#'   confidence limits.
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
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param method.predict A character string indicating which method is
#'   used to calculate a prediction interval (see
#'   \code{\link{meta-package}}).
#' @param adhoc.hakn.pi A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied for
#'   prediction interval (see \code{\link{meta-package}}).
#' @param seed.predict A numeric value used as seed to calculate
#'   bootstrap prediction interval (see \code{\link{meta-package}}).
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param method.bias A character string indicating which test is to
#'   be used. Either \code{"Begg"}, \code{"Egger"}, or
#'   \code{"Thompson"}, can be abbreviated. See function
#'   \code{\link{metabias}}.
#' @param backtransf A logical indicating whether results for
#'   transformed rates (argument \code{sm != "IR"}) should be back
#'   transformed in printouts and plots. If TRUE (default), results
#'   will be presented as incidence rates; otherwise transformed rates
#'   will be shown.
#' @param irscale A numeric defining a scaling factor for printing of
#'   rates.
#' @param irunit A character string specifying the time unit used to
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
#' @param label.left Graph label on left side of null effect in forest plot.
#' @param label.right Graph label on right side of null effect in forest plot.
#' @param col.label.left The colour of the graph label on the left side of
#'   the null effect.
#' @param col.label.right The colour of the graph label on the right side of
#'   the null effect.
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
#' @param hakn Deprecated argument (replaced by 'method.random.ci').
#' @param adhoc.hakn Deprecated argument (replaced by 'adhoc.hakn.ci').
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies or if estimation problems exist in fitting a GLMM).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}, respectively.
#' @param \dots Additional arguments passed on to
#'   \code{\link[metafor]{rma.glmm}} function and to catch deprecated
#'   arguments.
#' 
#' @details
#' This function provides methods for common effect and random effects
#' meta-analysis of single incidence rates to calculate an overall
#' rate. Note, you should use R function \code{\link{metainc}} to
#' compare incidence rates of pairwise comparisons instead of using
#' \code{metarate} for each treatment arm separately which will break
#' randomisation in randomised controlled trials.
#'
#' The following transformations of incidence rates are implemented to
#' calculate an overall rate:
#' 
#' \itemize{
#' \item Log transformation (\code{sm = "IRLN"}, default)
#' \item Square root transformation (\code{sm = "IRS"})
#' \item Freeman-Tukey Double arcsine transformation (\code{sm =
#'   "IRFT"})
#' \item No transformation (\code{sm = "IR"})
#' }
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain the transformed incidence rates. In printouts and plots
#' these values are back transformed if argument \code{backtransf =
#' TRUE} (default).
#'
#' By default (argument \code{method = "Inverse"}), the inverse
#' variance method (Borenstein et al., 2010) is used for pooling by
#' calling \code{\link{metagen}} internally. A random intercept
#' Poisson regression model (Stijnen et al., 2010) can be utilised
#' instead with argument \code{method = "GLMM"} which calls the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor}.
#' 
#' A three-level random effects meta-analysis model (Van den Noortgate
#' et al., 2013) is utilised if argument \code{cluster} is used and at
#' least one cluster provides more than one estimate. Internally,
#' \code{\link[metafor]{rma.mv}} is called to conduct the analysis and
#' \code{\link[metafor]{weights.rma.mv}} with argument \code{type =
#' "rowsum"} is used to calculate random effects weights.
#' 
#' Default settings are utilised for several arguments (assignments
#' using \code{\link{gs}} function). These defaults can be changed for
#' the current R session using the \code{\link{settings.meta}}
#' function.
#' 
#' Furthermore, R function \code{\link{update.meta}} can be used to
#' rerun a meta-analysis with different settings.
#' 
#' \subsection{Continuity correction}{
#'
#' Three approaches are available to apply a continuity correction:
#' \itemize{
#' \item Only studies with a zero cell count (\code{method.incr =
#'   "only0"})
#' \item All studies if at least one study has a zero cell count
#'   (\code{method.incr = "if0all"})
#' \item All studies irrespective of zero cell counts
#'   (\code{method.incr = "all"})
#' }
#' 
#' If the summary measure (argument \code{sm}) is equal to "IR" or
#' "IRLN", the continuity correction is applied if a study has zero
#' events, i.e., an incidence rate of 0.
#'
#' By default, 0.5 is used as continuity correction (argument
#' \code{incr}). This continuity correction is used both to calculate
#' individual study results with confidence limits and to conduct
#' meta-analysis based on the inverse variance method.
#'
#' For the Freeman-Tukey (Freeman & Tukey, 1950) and square root
#' transformation as well as GLMMs no continuity correction is used.
#' Furthermore, the value of \code{incr} is not considered for Poisson
#' confidence intervals for individual studies (\code{method.ci = "Poisson"}).
#' }
#'
#' \subsection{Subgroup analysis}{
#' 
#' Argument \code{subgroup} can be used to conduct subgroup analysis for
#' a categorical covariate. The \code{\link{metareg}} function can be
#' used instead for more than one categorical covariate or continuous
#' covariates.
#' }
#' 
#' \subsection{Specify the null hypothesis of test for an overall effect}{
#'
#' Argument \code{null.effect} can be used to specify the rate used
#' under the null hypothesis in a test for an overall effect.
#'
#' By default (\code{null.effect = NA}), no hypothesis test is
#' conducted as it is unclear which value is a sensible choice for the
#' data at hand.  An overall rate of 2, for example, could be tested
#' by setting argument \code{null.effect = 2}.
#'
#' Note, all tests for an overall effect are two-sided with the
#' alternative hypothesis that the effect is unequal to
#' \code{null.effect}.
#' }
#' 
#' \subsection{Exclusion of studies from meta-analysis}{
#'
#' Arguments \code{subset} and \code{exclude} can be used to exclude
#' studies from the meta-analysis. Studies are removed completely from
#' the meta-analysis using argument \code{subset}, while excluded
#' studies are shown in printouts and forest plots using argument
#' \code{exclude} (see Examples in \code{\link{metagen}}).
#' Meta-analysis results are the same for both arguments.
#' }
#' 
#' \subsection{Presentation of meta-analysis results}{
#' 
#' Internally, both common effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{common} and \code{random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{common} and \code{random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{random = FALSE}.
#' 
#' Argument \code{irscale} can be used to rescale rates, e.g.
#' \code{irscale = 1000} means that rates are expressed as events per
#' 1000 time units, e.g. person-years. This is useful in situations
#' with (very) low rates. Argument \code{irunit} can be used to
#' specify the time unit used in individual studies (default:
#' "person-years"). This information is printed in summaries and
#' forest plots if argument \code{irscale} is not equal to 1.
#'
#' A prediction interval will only be shown if \code{prediction =
#' TRUE}.
#' }
#' 
#' @return
#' An object of class \code{c("metarate", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{metacont}}, \code{\link{metagen}},
#'   \code{\link{print.meta}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
#' A basic introduction to fixed-effect and random-effects models for
#' meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{1}, 97--111
#' 
#' Freeman MF & Tukey JW (1950):
#' Transformations related to the angular and the square root.
#' \emph{Annals of Mathematical Statistics},
#' \bold{21}, 607--11
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67
#'
#' Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca J (2013):
#' Three-level meta-analysis of dependent effect sizes.
#' \emph{Behavior Research Methods},
#' \bold{45}, 576--94
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' @examples
#' # Apply various meta-analysis methods to estimate incidence rates
#' #
#' m1 <- metarate(4:1, c(10, 20, 30, 40))
#' m2 <- update(m1, sm = "IR")
#' m3 <- update(m1, sm = "IRS")
#' m4 <- update(m1, sm = "IRFT")
#' #
#' m1
#' m2
#' m3
#' m4
#' #
#' forest(m1)
#' forest(m1, irscale = 100)
#' forest(m1, irscale = 100, irunit = "person-days")
#' forest(m1, backtransf = FALSE)
#' \dontrun{
#' forest(m2)
#' forest(m3)
#' forest(m4)
#' }
#' 
#' m5 <- metarate(40:37, c(100, 200, 300, 400), sm = "IRFT")
#' m5
#' 
#' @export metarate


metarate <- function(event, time, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     cluster = NULL, rho = 0,
                     #
                     weights = NULL,
                     weights.common = weights, weights.random = weights,
                     #
                     n = NULL,
                     ##
                     method = "Inverse",
                     sm = gs("smrate"),
                     ##
                     incr = gs("incr"), method.incr = gs("method.incr"),
                     ##
                     method.ci = gs("method.ci.rate"), level = gs("level"),
                     ##
                     common = gs("common"),
                     random = gs("random") | !is.null(tau.preset),
                     overall = common | random,
                     overall.hetstat =
                       if (is.null(gs("overall.hetstat")))
                         common | random
                       else
                         gs("overall.hetstat"),   
                     prediction = gs("prediction") | !missing(method.predict),
                     ##
                     method.tau,
                     method.tau.ci = gs("method.tau.ci"),
                     level.hetstat = gs("level.hetstat"),
                     tau.preset = NULL, TE.tau = NULL,
                     tau.common = gs("tau.common"),
                     detail.tau = NULL,
                     #
                     method.I2 = gs("method.I2"),
                     #
                     level.ma = gs("level.ma"),
                     method.common.ci = gs("method.common.ci"),
                     method.random.ci = gs("method.random.ci"),
                     adhoc.hakn.ci = gs("adhoc.hakn.ci"),
                     ##
                     level.predict = gs("level.predict"),
                     method.predict = gs("method.predict"),
                     adhoc.hakn.pi = gs("adhoc.hakn.pi"),
                     seed.predict = NULL,
                     ##
                     null.effect = NA,
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     irscale = 1, irunit = "person-years",
                     ##
                     text.common = gs("text.common"),
                     text.random = gs("text.random"),
                     text.predict = gs("text.predict"),
                     text.w.common = gs("text.w.common"),
                     text.w.random = gs("text.w.random"),
                     ##
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
                     #
                     label.left = gs("label.left"),
                     label.right = gs("label.right"),
                     col.label.left = gs("col.label.left"),
                     col.label.right = gs("col.label.right"),
                     #
                     subgroup, subgroup.name = NULL,
                     print.subgroup.name = gs("print.subgroup.name"),
                     sep.subgroup = gs("sep.subgroup"),
                     test.subgroup = gs("test.subgroup"),
                     prediction.subgroup = gs("prediction.subgroup"),
                     seed.predict.subgroup = NULL,
                     ##
                     byvar, hakn, adhoc.hakn,
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"), warn.deprecated = gs("warn.deprecated"),
                     ##
                     control = NULL,
                     ...
                     ) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  
  chknumeric(rho, min = -1, max = 1)
  ##
  missing.method <- missing(method)
  chknull(sm)
  sm <- setchar(sm, gs("sm4rate"))
  ##
  chklevel(level)
  ##
  missing.method.tau <- missing(method.tau)
  if (missing.method.tau)
    method.tau <- if (method == "GLMM") "ML" else gs("method.tau")
  method.tau <- setchar(method.tau, gs("meth4tau"))
  ##
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
  method.tau.ci <- setchar(method.tau.ci, gs("meth4tau.ci"))
  missing.tau.common <- missing(tau.common)
  tau.common <- replaceNULL(tau.common, FALSE)
  chklogical(tau.common)
  #
  method.I2 <- setchar(method.I2, gs("meth4i2"))
  #
  chklogical(prediction)
  chklevel(level.predict)
  ##
  missing.method.predict <- missing(method.predict)
  ##
  method.tau <-
    set_method_tau(method.tau, missing.method.tau,
                 method.predict, missing.method.predict)
  method.predict <-
    set_method_predict(method.predict, missing.method.predict,
                     method.tau, missing.method.tau)
  ##
  if (any(method.predict == "NNF"))
    is_installed_package("pimeta", argument = "method.predict", value = "NNF")
  ##
  adhoc.hakn.pi <- setchar(replaceNA(adhoc.hakn.pi, ""), gs("adhoc4hakn.pi"))
  #
  chknumeric(null.effect, length = 1)
  ##
  method.bias <- setmethodbias(method.bias)
  ##
  chklogical(backtransf)
  ##
  chknumeric(irscale, length = 1)
  if (!backtransf & irscale != 1 & !is_untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.",
            call. = FALSE)
    irscale <- 1
  }
  chkchar(irunit)
  ##
  if (!is.null(text.common))
    chkchar(text.common, length = 1)
  if (!is.null(text.random))
    chkchar(text.random)
  if (!is.null(text.predict))
    chkchar(text.predict)
  if (!is.null(text.w.common))
    chkchar(text.w.common, length = 1)
  if (!is.null(text.w.random))
    chkchar(text.w.random, length = 1)
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metarate"
  ##
  method <- setchar(method, gs("meth4rate"))
  is.glmm <- method == "GLMM"
  #
  missing.method.common.ci <- missing(method.common.ci)
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  if (method != "Inverse" & method.common.ci == "IVhet") {
    if (!missing.method.common.ci)
      warning("Argument 'method.common.ci = \"IVhet\"' only available ",
              "if 'method = \"Inverse\".",
              call. = FALSE)
    method.common.ci <- "classic"
  }
  #
  missing.method.incr <- missing(method.incr)
  method.incr <- setchar(method.incr, gs("meth4incr"))
  ##
  method.ci <- setchar(method.ci, gs("ci4rate"))
  chklogical(warn)
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
  method.random.ci <-
    deprecated2(method.random.ci, missing(method.random.ci),
                hakn, missing(hakn),
                warn.deprecated)
  if (is.logical(method.random.ci))
    if (method.random.ci)
      method.random.ci <- "HK"
    else
      method.random.ci <- "classic"
  method.random.ci <- setchar(method.random.ci, gs("meth4random.ci"))
  ##
  adhoc.hakn.ci <-
    deprecated2(adhoc.hakn.ci, missing(adhoc.hakn.ci),
                adhoc.hakn, missing(adhoc.hakn), warn.deprecated)
  adhoc.hakn.ci <- setchar(replaceNA(adhoc.hakn.ci, ""), gs("adhoc4hakn.ci"))
  #
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
  addincr <-
    deprecated(method.incr, missing.method.incr, args, "addincr",
               warn.deprecated)
  allincr <-
    deprecated(method.incr, missing.method.incr, args, "allincr",
               warn.deprecated)
  if (missing.method.incr) {
    method.incr <- gs("method.incr")
    ##
    if (is.logical(addincr) && addincr)
      method.incr <- "all"
    else if (is.logical(allincr) && allincr)
      method.incr <- "if0all"
  }
  ##
  addincr <- allincr <- FALSE
  if (method.incr == "all")
    addincr <- TRUE
  else if (method.incr == "if0all")
    allincr <- TRUE
  ##
  ## Some more checks
  ##
  chklogical(overall)
  chklogical(overall.hetstat)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch 'event', 'time' and 'n' from data:
  ##
  event <- catch("event", mc, data, sfsp)
  chknull(event)
  k.All <- length(event)
  ##
  time <- catch("time", mc, data, sfsp)
  chknull(time)
  ##
  n <- catch("n", mc, data, sfsp)
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- catch("incr", mc, data, sfsp)
  chknumeric(incr, min = 0)
  ##
  ## Catch 'studlab', 'subgroup', 'subset', 'exclude' and 'cluster'
  ## from data:
  ##
  studlab <- catch("studlab", mc, data, sfsp)
  studlab <- setstudlab(studlab, k.All)
  ##
  missing.subgroup <- missing(subgroup)
  subgroup <- catch("subgroup", mc, data, sfsp)
  missing.byvar <- missing(byvar)
  byvar <- catch("byvar", mc, data, sfsp)
  ##
  subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                          warn.deprecated)
  by <- !is.null(subgroup)
  ##
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  ##
  exclude <- catch("exclude", mc, data, sfsp)
  missing.exclude <- is.null(exclude)
  ##
  cluster <- catch("cluster", mc, data, sfsp)
  with.cluster <- !is.null(cluster)
  #
  # Catch 'weights', 'weights.common', and 'weights.random' from data:
  #
  if (!missing(weights))
    weights <- catch("weights", mc, data, sfsp)
  if (!missing(weights.common))
    weights.common <- catch("weights.common", mc, data, sfsp)
  if (!missing(weights.random))
    weights.random <- catch("weights.random", mc, data, sfsp)
  #
  if (!is.null(weights) & is.null(weights.common))
    weights.common <- weights
  #
  if (!is.null(weights) & is.null(weights.random))
    weights.random <- weights
  #
  usw.common <- !is.null(weights.common)
  usw.random <- !is.null(weights.random)
  #
  if (usw.common)
    chknumeric(weights.common, min = 0)
  #
  if (usw.random)
    chknumeric(weights.random, min = 0)
  #
  if (usw.common & method != "Inverse")
    stop("User-specified weights for the common effect model only implemented ",
         "for the inverse variance method (method = \"Inverse\").",
         call. = FALSE)
  #
  if (usw.random & method == "GLMM")
    stop("User-specified weights for the random effects model not implemented ",
         "for generalized linear mixed models (method = \"GLMM\").",
         call. = FALSE)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  
  chklength(time, k.All, fun)
  if (!is.null(n))
    chklength(n, k.All, fun)
  chklength(studlab, k.All, fun)
  #
  if (with.cluster)
    chklength(cluster, k.All, fun)
  #
  if (usw.common) {
    if (length(weights.common) == 1)
      weights.common <- rep(weights.common, k.All)
    else
      chklength(weights.common, k.All, fun)
  }
  #
  if (usw.random) {
    if (length(weights.random) == 1)
      weights.random <- rep(weights.random, k.All)
    else
      chklength(weights.random, k.All, fun)
  }
  #
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (by) {
    chklength(subgroup, k.All, fun)
    chklogical(test.subgroup)
    chklogical(prediction.subgroup)
  }
  ##
  ## Additional checks
  ##
  if (!by & tau.common) {
    warning("Value for argument 'tau.common' set to FALSE as ",
            "argument 'subgroup' is missing.",
            call. = FALSE)
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as ",
            "argument tau.preset is not NULL.",
            call. = FALSE)
    tau.common <- TRUE
  }
  
  
  ##
  ##
  ## (4) Subset, exclude studies, and subgroups
  ##
  ##
  
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
  ##
  if (!missing.exclude) {
    if ((is.logical(exclude) & (sum(exclude) > k.All)) ||
        (length(exclude) > k.All))
      stop("Length of argument 'exclude' is larger than number of studies.")
    ##
    exclude2 <- rep(FALSE, k.All)
    exclude2[exclude] <- TRUE
    exclude <- exclude2
  }
  else
    exclude <- rep(FALSE, k.All)
  
  
  ##
  ##
  ## (5) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  
  if (keepdata) {
    if (nulldata)
      data <- data.frame(.event = event)
    else
      data$.event <- event
    ##
    data$.time <- time
    data$.studlab <- studlab
    if (!is.null(n))
      data$.n <- n
    ##
    data$.incr <- NA
    ##
    if (by)
      data$.subgroup <- subgroup
    ##
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
    ##
    if (!missing.exclude)
      data$.exclude <- exclude
    ##
    if (with.cluster)
      data$.id <- data$.cluster <- cluster
    #
    if (usw.common)
      data$.weights.common <- weights.common
    #
    if (usw.random)
      data$.weights.random <- weights.random
  }
  
    
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  
  if (!missing.subset) {
    event <- event[subset]
    time  <- time[subset]
    studlab <- studlab[subset]
    if (!is.null(n))
      n <- n[subset]
    ##
    cluster <- cluster[subset]
    exclude <- exclude[subset]
    #
    weights.common <- weights.common[subset]
    weights.random <- weights.random[subset]
    #
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (by)
      subgroup <- subgroup[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1) {
    common <- FALSE
    random <- FALSE
    prediction <- FALSE
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(event, 0)
  chknumeric(time, 0, zero = TRUE)
  if (!is.null(n))
    chknumeric(n, 0, zero = TRUE)
  ##
  ## Recode integer as numeric:
  ##
  event <- int2num(event)
  time  <- int2num(time)
  if (!is.null(n))
    n <- int2num(n)
  ##
  ## Check for whole numbers
  ##
  if (method.ci != "NAsm") {
    if (any(!is_wholenumber(event), na.rm = TRUE)) {
      warning("Normal approximation confidence interval ",
              "(argument method.ci = \"NAsm\") used as\n",
              "at least one number of events contains a non-integer value.",
              call. = FALSE)
      method.ci <- "NAsm"
    }
  }
  ##
  if (by) {
    chkmiss(subgroup)
    ##
    if (missing.subgroup.name & is.null(subgroup.name)) {
      if (!missing.subgroup)
        subgroup.name <- byvarname("subgroup", mc)
      else if (!missing.byvar)
        subgroup.name <- byvarname("byvar", mc)
    }
  }
  ##
  if (!is.null(subgroup.name))
    chkchar(subgroup.name, length = 1)
  
  
  #
  #
  # (7) Continuity correction
  #
  #
  
  sel <- switch(sm,
                IR   = event == 0,
                IRLN = event == 0,
                IRS  = rep(FALSE, length(event)),
                IRFT = rep(FALSE, length(event)))
  #
  sparse <- any(sel, na.rm = TRUE)
  #
  # No need to add anything to cell counts for arcsine transformation
  #
  if (addincr | method.incr == "user")
    incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
  else {
    if (sparse) {
      if (allincr)
        incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
      else
        incr.event <- incr * sel
    }
    else
      incr.event <- rep(0, k.all)
  }
  #
  if (keepdata) {
    if (missing.subset) {
      data$.incr <- incr.event
    }
    else {
      data$.incr <- NA
      #
      data$.incr[subset] <- incr.event
    }
  }
  
  
  ##
  ##
  ## (8) Calculate results for individual studies
  ##
  ##
  
  if (sm == "IR") {
    TE <- (event + incr.event) / time
    seTE <- sqrt(TE / time)
    transf.null.effect <- null.effect
  }
  else if (sm == "IRLN") {
    TE <- log((event + incr.event) / time)
    seTE <- sqrt(1 / (event + incr.event))
    transf.null.effect <- log(null.effect)
  }
  else if (sm == "IRS") {
    TE <- sqrt(event / time)
    seTE <- sqrt(1 / (4 * time))
    transf.null.effect <- sqrt(null.effect)
  }
  else if (sm == "IRFT") {
    TE <- 0.5 * (sqrt(event / time) + sqrt((event + 1) / time))
    seTE <- sqrt(1 / (4 * time))
    transf.null.effect <- sqrt(null.effect)
  }
  ##
  ## Calculate confidence intervals
  ##
  if (method.ci == "NAsm")
    ci.study <- ci(TE, seTE, level = level)
  else
    ci.study <- ciPoisson(event, time, level, null.effect)
  ##
  lower.study <- ci.study$lower
  upper.study <- ci.study$upper
  ##
  if (method.ci != "NAsm") {
    if (sm == "IRLN") {
      lower.study <- log(lower.study)
      upper.study <- log(upper.study)
    }
    else if (sm == "IRS") {
      lower.study <- sqrt(lower.study)
      upper.study <- sqrt(upper.study)
    }
    ##
    else if (sm == "IRFT") {
      lower.ev <- time * lower.study 
      upper.ev <- time * upper.study 
      ##
      lower.study <-
        0.5 * (sqrt(lower.ev / time) + sqrt((lower.ev + 1) / time))
      upper.study <-
        0.5 * (sqrt(upper.ev / time) + sqrt((upper.ev + 1) / time))
    }
  }
  
  
  ##
  ##
  ## (9) Additional checks for three-level model
  ##
  ##
  
  three.level <- FALSE
  sel.ni <- !is.infinite(TE) & !is.infinite(seTE)
  ##
  ## Only conduct three-level meta-analysis if variable 'cluster'
  ## contains duplicate values after removing inestimable study
  ## results standard errors
  ##
  if (with.cluster &&
      length(unique(cluster[sel.ni])) != length(cluster[sel.ni]))
    three.level <- TRUE
  ##
  if (three.level) {
    chkmlm(method.tau, missing.method.tau, method.predict,
           method, missing.method)
    ##
    common <- FALSE
    method <- "Inverse"
    is.glmm <- FALSE
    ##
    if (!(method.tau %in% c("REML", "ML")))
      method.tau <- "REML"
  }
  
  
  ##
  ##
  ## (10) Additional checks for GLMM
  ##
  ##
  
  if (is.glmm) {
    chkglmm(sm, method.tau, method.random.ci, method.predict,
            adhoc.hakn.ci, adhoc.hakn.pi,
            "IRLN")
    ##
    if (sparse & warn &
        ((!missing(incr) & any(incr != 0)) | allincr | addincr))
        warning("Note, for method = \"GLMM\", continuity correction only ",
                "used to calculate individual study results.",
                call. = FALSE)
    ##
    if (!is.null(TE.tau)) {
      if (warn)
        warning("Argument 'TE.tau' not considered for GLMM.",
                call. = FALSE)
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)) {
      if (warn)
        warning("Argument 'tau.preset' not considered for GLMM.",
                call. = FALSE)
      tau.preset <- NULL
    }
  }
  
  
  ##
  ##
  ## (11) Do meta-analysis
  ##
  ##
  
  k <- sum(!is.na(event[!exclude]) & !is.na(time[!exclude]))
  ##
  for (i in seq_along(method.random.ci))
    if (k == 1 & method.random.ci[i] == "HK")
      method.random.ci[i] <- "classic"
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               cluster = cluster, rho = rho,
               #
               weights.common = weights.common,
               weights.random = weights.random,
               #
               sm = sm,
               level = level,
               ##
               common = common,
               random = random,
               overall = overall,
               overall.hetstat = overall.hetstat,
               prediction = prediction,
               ##
               method.tau = if (is.glmm) "DL" else method.tau,
               method.tau.ci = if (is.glmm) "" else method.tau.ci,
               level.hetstat = level.hetstat,
               tau.preset = tau.preset,
               TE.tau = TE.tau,
               tau.common = FALSE,
               detail.tau = detail.tau,
               #
               method.I2 = method.I2,
               #
               level.ma = level.ma,
               method.common.ci = method.common.ci,
               method.random.ci = method.random.ci,
               adhoc.hakn.ci = adhoc.hakn.ci,
               ##
               level.predict = level.predict,
               method.predict = method.predict,
               adhoc.hakn.pi = adhoc.hakn.pi,
               seed.predict = seed.predict,
               ##
               null.effect = transf.null.effect,
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
               label.left = label.left, label.right = label.right,
               col.label.left = col.label.left,
               col.label.right = col.label.right,
               #
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  #
  # Estimate common tau-squared across subgroups
  #
  if (by & tau.common & !is.glmm)
    hcc <- hetcalc(TE, seTE, method.tau, "", TE.tau,
                   method.I2, level.hetstat, subgroup, control)
  
  
  ##
  ##
  ## (12) Generate R object
  ##
  ##
  
  res <- list(event = event, time = time,
              n = n,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              method.incr = method.incr,
              sparse = sparse,
              method.ci = method.ci,
              incr.event = incr.event)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$pscale <- NULL
  m$irscale <- NULL
  m$irunit <- NULL
  m$method.ci <- NULL
  m$method.mean <- NULL
  m$approx.TE <- NULL
  m$approx.seTE <- NULL
  ##
  m$label.e <- ""
  m$label.c <- ""
  m$label.left <- ""
  m$label.right <- ""
  ##
  if (method.ci == "exact") {
    m$statistic <- rep(NA, length(m$statistic))
    m$pval <- ci.study$p
  }    
  ##
  if (is.glmm) {
    m$method.tau <- method.tau
    #
    m$seTE.hakn.ci <- m$seTE.hakn.adhoc.ci <-
      m$seTE.hakn.pi <- m$seTE.hakn.adhoc.pi <-
        m$seTE.kero <- NA
    ##
    m$text.random <- gsub("(HK)", "(T)", m$text.random, fixed = TRUE)
  }
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Run GLMM and add data
  ##
  if (is.glmm & k > 0) {
    res$method <- "GLMM"
    res$method.random <- "GLMM"
    ##
    list.rate <-
      list(xi = event[!exclude], ti = time[!exclude], measure = "IRLN")
    ##
    use.random <-
      sum(!exclude) > 1 &
      sum(event[!exclude], na.rm = TRUE) > 0 &
      any(event[!exclude] != time[!exclude])
    ##
    res.glmm <-
      runGLMM(list.rate,
              method.tau = method.tau,
              method.random.ci = method.random.ci,
              level = level.ma,
              control = control, use.random = use.random,
              warn = warn)
    ##
    res <- addGLMM(res, res.glmm, method.I2, transf.null.effect)
    ##
    if (by) {
      n.subgroups <- length(unique(subgroup[!exclude]))
      if (n.subgroups > 1)
        subgroup.glmm <-
          factor(subgroup[!exclude], bylevs(subgroup[!exclude]))
      ##
      hcc <-
        hccGLMM(
          res,
          runGLMM(list.rate,
                  method.tau = method.tau,
                  method.random.ci = method.random.ci,
                  level = level.ma,
                  data =
                    if (n.subgroups > 1)
                      list(data = data.frame(subgroup.glmm))
                    else
                      NULL,
                  mods =
                    if (n.subgroups > 1)
                      as.call(~ subgroup.glmm)
                    else
                      NULL,
                  control = control, use.random = use.random,
                  warn = warn)$glmm.random[[1]],
          method.I2
        )
    }
  }
  ##
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$irscale <- irscale
  res$irunit  <- irunit
  #
  res$pairwise <- FALSE
  #
  res$call <- match.call()
  #
  if (keepdata) {
    res$data <- data
    if (!missing.subset)
      res$subset <- subset
  }
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (by) {
    res$subgroup <- subgroup
    res$subgroup.name <- subgroup.name
    res$print.subgroup.name <- print.subgroup.name
    res$sep.subgroup <- sep.subgroup
    res$test.subgroup <- test.subgroup
    res$prediction.subgroup <- prediction.subgroup
    res$tau.common <- tau.common
    ##
    if (!tau.common) {
      res <- c(res, subgroup(res, seed = seed.predict.subgroup))
      if (res$three.level)
        res <- setNA3(res)
    }
    else if (!is.null(tau.preset))
      res <-
        c(res, subgroup(res, tau.preset, seed = seed.predict.subgroup))
    else {
      if (is.glmm)
        res <- c(res,
                 subgroup(res, NULL,
                          factor(res$subgroup, bylevs(res$subgroup)), ...))
      else if (res$three.level)
        res <- c(res,
                 subgroup(res, NULL,
                          factor(res$subgroup, bylevs(res$subgroup))))
      else
        res <-
          c(res, subgroup(res, hcc$tau.resid, seed = seed.predict.subgroup))
    }
    ##
    if (tau.common && is.null(tau.preset))
      res <- addHet(res, hcc, !is.glmm)
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
    ##
    res <- setNAwithin(res, res$three.level | is.glmm)
  }
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")

  
  res
}
