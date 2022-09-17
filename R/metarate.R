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
#' @param n Number of observations.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"} and
#'   \code{"GLMM"}, can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"IR"}, \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}) is
#'   to be used for pooling of studies, see Details.
#' @param incr A numeric which is added to the event number of studies
#'   with zero events, i.e., studies with an incidence rate of 0.
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
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
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
#' @param hakn Deprecated argument (replaced by 'method.random.ci').
#' @param adhoc.hakn Deprecated argument (replaced by 'adhoc.hakn.ci').
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether the addition of
#'   \code{incr} to studies with zero events should result in a
#'   warning.
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
#' By default (argument \code{method = "Inverse"}), the inverse
#' variance method (Borenstein et al., 2010) is used for pooling by
#' calling \code{\link{metagen}} internally. A random intercept
#' Poisson regression model (Stijnen et al., 2010) can be utilised
#' instead with argument \code{method = "GLMM"} which calls the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor}.
#' 
#' A three-level random effects meta-analysis model (Van den Noortgate
#' et al., 2013) is utilized if argument \code{cluster} is used and at
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
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
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
                     cluster = NULL,
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
                     overall.hetstat = common | random,
                     prediction = gs("prediction") | !missing(method.predict),
                     ##
                     method.tau,
                     method.tau.ci = gs("method.tau.ci"),
                     tau.preset = NULL, TE.tau = NULL,
                     tau.common = gs("tau.common"),
                     ##
                     level.ma = gs("level.ma"),
                     method.random.ci = gs("method.random.ci"),
                     adhoc.hakn.ci = gs("adhoc.hakn.ci"),
                     ##
                     level.predict = gs("level.predict"),
                     method.predict = gs("method.predict"),
                     adhoc.hakn.pi = gs("adhoc.hakn.pi"),
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
                     ##
                     subgroup, subgroup.name = NULL,
                     print.subgroup.name = gs("print.subgroup.name"),
                     sep.subgroup = gs("sep.subgroup"),
                     test.subgroup = gs("test.subgroup"),
                     prediction.subgroup = gs("prediction.subgroup"),
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
  tau.common <- replaceNULL(tau.common, FALSE)
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  missing.method.predict <- missing(method.predict)
  ##
  method.tau <-
    setmethodtau(method.tau, missing.method.tau,
                 method.predict, missing.method.predict)
  method.predict <-
    setmethodpredict(method.predict, missing.method.predict,
                     method.tau, missing.method.tau)
  ##
  if (method.predict == "NNF")
    is.installed.package("pimeta", argument = "method.predict", value = "NNF")
  ##
  missing.adhoc.hakn.pi <- missing(adhoc.hakn.pi)
  adhoc.hakn.pi <- setchar(adhoc.hakn.pi, gs("adhoc4hakn.pi"))
  ##
  chknumeric(null.effect, length = 1)
  ##
  method.bias <- setmethodbias(method.bias)
  ##
  chklogical(backtransf)
  ##
  chknumeric(irscale, length = 1)
  if (!backtransf & irscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  chkchar(irunit)
  ##
  if (!is.null(text.common))
    chkchar(text.common, length = 1)
  if (!is.null(text.random))
    chkchar(text.random, length = 1)
  if (!is.null(text.predict))
    chkchar(text.predict, length = 1)
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
  ##
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
  missing.adhoc.hakn.ci <- missing(adhoc.hakn.ci)
  adhoc.hakn.ci <-
    deprecated2(adhoc.hakn.ci, missing(adhoc.hakn.ci),
                adhoc.hakn, missing(adhoc.hakn), warn.deprecated)
  adhoc.hakn.ci <- setchar(adhoc.hakn.ci, gs("adhoc4hakn.ci"))
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
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(time, k.All, fun)
  if (!is.null(n))
    chklength(n, k.All, fun)
  chklength(studlab, k.All, fun)
  if (with.cluster)
    chklength(cluster, k.All, fun)
  ##
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
    warning("Value for argument 'tau.common' set to FALSE as argument 'subgroup' is missing.")
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as argument tau.preset is not NULL.")
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
    data$.incr <- incr
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
    exclude <- exclude[subset]
    ##
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
    if (any(!is.wholenumber(event), na.rm = TRUE)) {
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
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                IR   = event == 0,
                IRLN = event == 0,
                IRS  = rep(FALSE, length(event)),
                IRFT = rep(FALSE, length(event)))
  ##
  sparse <- any(sel, na.rm = TRUE)
  ##
  ## No need to add anything to cell counts for arcsine transformation
  ##
  if (addincr)
    incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
  else
    if (sparse)
      if (allincr)
        incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
      else
        incr.event <- incr * sel
  else
    incr.event <- rep(0, k.all)
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
  if (method.ci == "NAsm") {
    if (sm == "IRLN") {
      lower.study <- exp(lower.study)
      upper.study <- exp(upper.study)
    }
    else if (sm == "IRS") {
      lower.study <- lower.study^2
      upper.study <- upper.study^2
    }
    ##
    else if (sm == "IRFT") {
      lower.study <-
        asin2ir(lower.study, time, value = "lower", warn = FALSE)
      upper.study <-
        asin2ir(upper.study, time, value = "upper", warn = FALSE)
    }
    ##
    lower.study[lower.study < 0] <- 0
  }
  
  
  ##
  ##
  ## (8) Additional checks for three-level model
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
    common <- FALSE
    ##
    if (method != "Inverse") {
      if (!missing.method)
        warning("Inverse variance method used in three-level model.",
                call. = FALSE)
      method <- "Inverse"
      is.glmm <- FALSE
    }
    ##
    if (!(method.tau %in% c("REML", "ML"))) {
      if (!missing(method.tau))
        warning("For three-level model, argument 'method.tau' set to \"REML\".",
                call. = FALSE)
      method.tau <- "REML"
    }
    ##
    if (by & !tau.common) {
      if (!missing(tau.common))
        warning("For three-level model, argument 'tau.common' set to ",
                "\"TRUE\".",
                call. = FALSE)
      tau.common <- TRUE
    }
  }
  
  
  ##
  ##
  ## (9) Additional checks for GLMM
  ##
  ##
  if (is.glmm & sm != "IRLN")
    stop("Generalised linear mixed models only possible with ",
         "argument 'sm = \"IRLN\"'.")
  ##
  if (is.glmm & method.tau != "ML")
    stop("Generalised linear mixed models only possible with ",
         "argument 'method.tau = \"ML\"'.")
  ##
  if (is.glmm & method.random.ci == "KR")
    stop("Kenward-Roger method for random effects meta-analysis not ",
         "available for GLMMs.",
         call. = FALSE)
  ##
  if (is.glmm & method.predict == "KR")
    stop("Kenward-Roger method for prediction interval not ",
         "available for GLMMs.",
         call. = FALSE)
  ##
  if (is.glmm & method.predict == "NNF")
    stop("Bootstrap method for prediction interval not ",
         "available for GLMMs.",
         call. = FALSE)
  ##
  if (is.glmm & method.random.ci == "HK" & adhoc.hakn.ci != "") {
    if (!missing.adhoc.hakn.ci)
      warning("Ad hoc correction for Hartung-Knapp method not ",
              "available for GLMMs.",
              call. = FALSE)
    adhoc.hakn.ci <- ""
  }
  ##
  if (is.glmm & method.predict == "HK" & adhoc.hakn.pi != "") {
    if (!missing.adhoc.hakn.pi)
      warning("Ad hoc Hartung-Knapp correction ffor prediction interval not ",
              "available for GLMMs.",
              call. = FALSE)
    adhoc.hakn.pi <- ""
  }
  ##
  if (is.glmm & sparse)
    if ((!missing(incr) & any(incr != 0)) |
        allincr | addincr)
      warning("Note, for method = \"GLMM\", continuity correction only ",
              "used to calculate individual study results.",
              call. = FALSE)
  ##
  if (is.glmm) {
    if (!is.null(TE.tau)) {
      if (warn)
        warning("Argument 'TE.tau' not considered for GLMM.",
                call. = FALSE)
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)) {
      if (warn)
        warning("Argument 'tau.preset' not considered for GLMM.")
      tau.preset <- NULL
    }
  }
  
  
  ##
  ##
  ## (10) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event[!exclude]) & !is.na(time[!exclude]))
  ##
  if (k == 1 & method.random.ci != "classic")
    method.random.ci <- "classic"
  ##
  if (is.glmm & k > 0) {
    glmm.common <-
      runNN(rma.glmm,
            list(xi = event[!exclude], ti = time[!exclude],
                 method = "FE",
                 test = ifelse(method.random.ci == "HK", "t", "z"),
                 level = 100 * level.ma,
                 measure = "IRLN", control = control,
                 ...))
    ##
    TE.common   <- as.numeric(glmm.common$b)
    seTE.common <- as.numeric(glmm.common$se)
    ##
    w.common <- rep(NA, length(event))
  }
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               cluster = cluster,
               ##
               sm = sm,
               level = level,
               ##
               common = common,
               random = random,
               overall = overall,
               overall.hetstat = overall.hetstat,
               prediction = prediction,
               ##
               method.tau = method.tau, method.tau.ci = method.tau.ci,
               tau.preset = tau.preset,
               TE.tau = TE.tau,
               tau.common = FALSE,
               ##
               level.ma = level.ma,
               method.random.ci = method.random.ci,
               adhoc.hakn.ci = adhoc.hakn.ci,
               ##
               level.predict = level.predict,
               method.predict = method.predict,
               adhoc.hakn.pi = adhoc.hakn.pi,
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
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (by & tau.common & !is.glmm) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   TE.tau, level.ma, subgroup, control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
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
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Add data
  ##
  if (is.glmm & k > 0) {
    ##
    ci.f <- ci(TE.common, seTE.common, level = level.ma,
               null.effect = transf.null.effect)
    ##
    res$method <- "GLMM"
    ##
    res$TE.common <- TE.common
    res$seTE.common <- seTE.common
    res$w.common <- w.common
    res$lower.common <- ci.f$lower
    res$upper.common <- ci.f$upper
    res$statistic.common <- ci.f$statistic
    res$pval.common <- ci.f$statistic
    res$zval.common <- ci.f$statistic
    ##
    if (sum(!exclude) > 1 &
        sum(event[!exclude], na.rm = TRUE) > 0 &
        any(event[!exclude] != time[!exclude]))
      glmm.random <-
        runNN(rma.glmm,
              list(xi = event[!exclude], ti = time[!exclude],
                   method = method.tau,
                   test = ifelse(method.random.ci == "HK", "t", "z"),
                   level = 100 * level.ma,
                   measure = "IRLN", control = control,
                   ...))
    else {
      ##
      ## Fallback to common effect model due to small number of studies
      ## or zero events (or number of events equal to time at risk)
      ##
      glmm.random <- glmm.common
    }
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.ma,
               null.effect = transf.null.effect,
               df = if (method.random.ci == "HK") k - 1)
    ##
    res$w.random <- rep(NA, length(event))
    ##
    res$TE.random <- TE.random
    res$seTE.random <- seTE.random
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p
    res$zval.random <- ci.r$statistic
    ##
    ## Prediction interval
    ##
    res$upper.predict <- res$lower.predict <- res$seTE.predict <- NA
    ##
    tau2.calc <- if (is.na(glmm.random$tau2)) 0 else glmm.random$tau2
    seTE.predict <- sqrt(seTE.random^2 + tau2.calc)
    ##
    if (method.predict == "HTS" && k >= 3)
      ci.p <- ci(TE.random, seTE.predict, level.predict, k - 2)
    else if (method.predict == "S")
      ci.p <- ci(TE.random, seTE.predict, level.predict)
    else
      ci.p <- list(lower = NA, upper = NA)
    ##
    res$seTE.predict <- seTE.predict
    res$lower.predict <- ci.p$lower
    res$upper.predict <- ci.p$upper
    ##
    res$Q <- if (glmm.random$k > 1) glmm.random$QE.Wld else 0
    res$df.Q <- glmm.random$QE.df
    res$pval.Q <- pvalQ(res$Q, res$df.Q)
    ##
    res$Q.LRT <- if (glmm.random$k > 1) glmm.random$QE.LRT else 0
    res$df.Q.LRT <- res$df.Q
    res$pval.Q.LRT <- pvalQ(res$Q.LRT, res$df.Q.LRT)
    ##
    if (k > 1) {
      res$tau <- sqrt(glmm.random$tau2)
      res$tau2 <- glmm.random$tau2
      res$se.tau2 <- glmm.random$se.tau2
    }
    else
      res$se.tau2 <- NA
    ##
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    ##
    res$lower.tau <- NA
    res$upper.tau <- NA
    ##
    res$method.tau.ci <- ""
    res$sign.lower.tau <- ""
    res$sign.upper.tau <- ""
    ##
    H <- calcH(res$Q, res$df.Q, level.ma)
    res$H <- H$TE
    res$lower.H <- H$lower
    res$upper.H <- H$upper
    ##
    I2 <- isquared(res$Q, res$df.Q, level.ma)
    res$I2 <- I2$TE
    res$lower.I2 <- I2$lower
    res$upper.I2 <- I2$upper
    ##
    res$Rb <- NA
    res$lower.Rb <- NA
    res$upper.Rb <- NA
    ##
    res$.glmm.common  <- glmm.common
    res$.glmm.random <- glmm.random
    res$version.metafor <- packageDescription("metafor")$Version
    ##
    if (by) {
      n.subgroups <- length(unique(subgroup[!exclude]))
      if (n.subgroups > 1)
        subgroup.glmm <- factor(subgroup[!exclude], bylevs(subgroup[!exclude]))
      else
        subgroup.glmm <- NA
      ##
      glmm.random.by <-
        try(suppressWarnings(
          runNN(rma.glmm,
                list(xi = event[!exclude],
                     ti = time[!exclude],
                     mods =
                       if (n.subgroups > 1)
                         as.call(~ subgroup.glmm) else NULL,
                     method = method.tau,
                     test = ifelse(method.random.ci == "HK", "t", "z"),
                     level = 100 * level.ma,
                     measure = "IRLN", control = control,
                     data = data.frame(subgroup.glmm),
                     ...))),
          silent = TRUE)
      ##
      if ("try-error" %in% class(glmm.random.by))
        if (grepl(paste0("Number of parameters to be estimated is ",
                         "larger than the number of observations"),
                  glmm.random.by)) {
          glmm.random.by <-
            suppressWarnings(
              runNN(rma.glmm,
                    list(xi = event[!exclude],
                         ti = time[!exclude],
                         mods =
                           if (n.subgroups > 1)
                             as.call(~ subgroup.glmm) else NULL,
                         method = "FE",
                         test = ifelse(method.random.ci == "HK", "t", "z"),
                         level = 100 * level.ma,
                         measure = "IRLN", control = control,
                         data = data.frame(subgroup.glmm),
                         ...)))
        }
        else
          stop(glmm.random.by)
      ##
      Q.r <- glmm.random.by$QE.Wld
      df.Q.r <- glmm.random.by$k - glmm.random.by$p
      ##
      H.r  <- calcH(Q.r, df.Q.r, level.ma)
      I2.r <- isquared(Q.r, df.Q.r, level.ma)
      ##
      hcc <- list(tau2.resid = glmm.random.by$tau2,
                  lower.tau2.resid = NA,
                  upper.tau2.resid = NA,
                  ##
                  tau.resid = sqrt(glmm.random.by$tau2),
                  lower.tau.resid = NA,
                  upper.tau.resid = NA,
                  sign.lower.tau.resid = "",
                  sign.upper.tau.resid = "",
                  ##
                  Q.resid = Q.r,
                  df.Q.resid = df.Q.r,
                  pval.Q.resid = pvalQ(Q.r, df.Q.r),
                  ##
                  H.resid = H.r$TE,
                  lower.H.resid = H.r$lower,
                  upper.H.resid = H.r$upper,
                  ##
                  I2.resid = I2.r$TE,
                  lower.I2.resid = I2.r$lower,
                  upper.I2.resid = I2.r$upper
                  )
    }
  }
  ##
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$irscale <- irscale
  res$irunit  <- irunit
  ##
  res$call <- match.call()
  res$allincr <- allincr
  res$addincr <- addincr
  ##
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
      res <- c(res, subgroup(res))
      if (res$three.level)
        res <- setNA3(res)
    }
    else if (!is.null(tau.preset))
      res <- c(res, subgroup(res, tau.preset))
    else {
      if (is.glmm)
        res <- c(res, subgroup(res, NULL,
                               factor(res$subgroup, bylevs(res$subgroup)), ...))
      else if (res$three.level)
        res <- c(res,
                 subgroup(res, NULL,
                          factor(res$subgroup, bylevs(res$subgroup))))
      else
        res <- c(res, subgroup(res, hcc$tau.resid))
    }
    ##
    if (tau.common && is.null(tau.preset)) {
      res$Q.w.random <- hcc$Q.resid
      res$df.Q.w.random <- hcc$df.Q.resid
      res$pval.Q.w.random <- hcc$pval.Q.resid
      ##
      res$tau2.resid <- hcc$tau2.resid
      res$lower.tau2.resid <- hcc$lower.tau2.resid
      res$upper.tau2.resid <- hcc$upper.tau2.resid
      ##
      res$tau.resid <- hcc$tau.resid
      res$lower.tau.resid <- hcc$lower.tau.resid
      res$upper.tau.resid <- hcc$upper.tau.resid
      res$sign.lower.tau.resid <- hcc$sign.lower.tau.resid
      res$sign.upper.tau.resid <- hcc$sign.upper.tau.resid
      ##
      res$Q.resid <- hcc$Q.resid
      res$df.Q.resid <- hcc$df.Q.resid
      res$pval.Q.resid <- hcc$pval.Q.resid
      ##
      res$H.resid <- hcc$H.resid
      res$lower.H.resid <- hcc$lower.H.resid
      res$upper.H.resid <- hcc$upper.H.resid
      ##
      res$I2.resid <- hcc$I2.resid
      res$lower.I2.resid <- hcc$lower.I2.resid
      res$upper.I2.resid <- hcc$upper.I2.resid
    }
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  ## Backward compatibility
  ##
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$w.fixed <- res$w.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  res$zval.fixed <- res$zval.common
  ##
  if (by) {
    res$byvar <- subgroup
    res$bylab <- subgroup.name
    res$print.byvar <- print.subgroup.name
    res$byseparator <- sep.subgroup
  }
  ##
  class(res) <- c(fun, "meta")

  
  res
}
