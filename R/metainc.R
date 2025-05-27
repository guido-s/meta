#' Meta-analysis of incidence rates
#' 
#' @description
#' Calculation of common effect and random effects estimates (incidence
#' rate ratio or incidence rate difference) for meta-analyses with
#' event counts.  Mantel-Haenszel, Cochran, inverse variance method,
#' and generalised linear mixed model (GLMM) are available for
#' pooling. For GLMMs, the \code{\link[metafor]{rma.glmm}} function
#' from R package \bold{metafor} (Viechtbauer 2010) is called
#' internally.
#' 
#' @param event.e Number of events in experimental group or an R object
#'   created with \code{\link{pairwise}}.
#' @param time.e Person time at risk in experimental group.
#' @param event.c Number of events in control group.
#' @param time.c Person time at risk in control group.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., event.e, time.e, event.c, and time.c.
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
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"MH"},
#'   \code{"Inverse"}, \code{"Cochran"}, or \code{"GLMM"} can be
#'   abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"IRR"}, \code{"IRD"}, \code{"IRSD"}, or \code{"VE"}) is to
#'   be used for pooling of studies, see Details.
#' @param incr A numerical value which is added to cell frequencies
#'   for studies with a zero cell count or a numeric vector with the continuity
#'   correction for each study, see Details.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, \code{"all"}, or \code{"user"}), see Details.
#' @param incr.e Continuity correction in experimental group, see Details.
#' @param incr.c Continuity correction in control group, see Details.
#' @param model.glmm A character string indicating which GLMM should
#'   be used. One of \code{"UM.FS"}, \code{"UM.RS"}, and
#'   \code{"CM.EL"}, see Details.
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
#' @param method.bias A character string indicating which test is to
#'   be used. Either \code{"Begg"}, \code{"Egger"}, or
#'   \code{"Thompson"}, can be abbreviated. See function
#'   \code{\link{metabias}}.
#' @param n.e Number of observations in experimental group (optional).
#' @param n.c Number of observations in control group (optional).
#' @param backtransf A logical indicating whether results for
#'   incidence rate ratio (\code{sm = "IRR"}) and vaccine efficacy or
#'   vaccine effectiveness (\code{sm = "VE"}) should be back
#'   transformed in printouts and plots. If TRUE (default), results
#'   will be presented as incidence rate ratios or vaccine efficacy /
#'   effectiveness; otherwise log incidence rate ratios or log vaccine
#'   rate ratios will be shown.
#' @param irscale A numeric defining a scaling factor for printing of
#'   incidence rate differences.
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
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
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
#' @param prediction.subgroup A logical indicating whether prediction
#'   intervals should be printed for subgroups.
#' @param seed.predict.subgroup A numeric vector providing seeds to
#'   calculate bootstrap prediction intervals within subgroups. Must
#'   be of same length as the number of subgroups.
#' @param byvar Deprecated argument (replaced by 'subgroup').
#' @param hakn Deprecated argument (replaced by 'method.random.ci').
#' @param adhoc.hakn Deprecated argument (replaced by
#'   'adhoc.hakn.ci').
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
#' Calculation of common and random effects estimates for meta-analyses
#' comparing two incidence rates.
#' 
#' The following measures of treatment effect are available:
#' \itemize{
#' \item Incidence Rate Ratio (\code{sm = "IRR"})
#' \item Incidence Rate Difference (\code{sm = "IRD"})
#' \item Square root transformed Incidence Rate Difference (\code{sm =
#'   "IRSD"})
#' \item Vaccine efficacy or vaccine effectiveness (\code{sm = "VE"})
#' }
#'
#' Note, log incidence rate ratio (logIRR) and log vaccine ratio
#' (logVR) are mathematical identical, however, back-transformed
#' results differ as vaccine efficacy or effectiveness is defined as
#' \code{VE = 100 * (1 - IRR)}.
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
#' \subsection{Meta-analysis method}{
#' 
#' By default, both common effect and random effects models are
#' considered (see arguments \code{common} and \code{random}). If
#' \code{method} is \code{"MH"} (default), the Mantel-Haenszel method
#' is used to calculate the common effect estimate (Greenland &
#' Robbins, 1985); if \code{method} is \code{"Inverse"}, inverse
#' variance weighting is used for pooling; if \code{method} is
#' \code{"Cochran"}, the Cochran method is used for pooling
#' (Bayne-Jones, 1964, Chapter 8). For these three methods, the random
#' effects estimate is always based on the inverse variance method.
#' 
#' A distinctive and frequently overlooked advantage of incidence
#' rates is that individual patient data (IPD) can be extracted from
#' count data. Accordingly, statistical methods for IPD, i.e.,
#' generalised linear mixed models, can be utilised in a meta-analysis
#' of incidence rate ratios (Stijnen et al., 2010). These methods are
#' available (argument \code{method = "GLMM"}) for the common effect
#' and random effects model by calling the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} internally.
#'
#' Three different GLMMs are available for meta-analysis of incidence
#' rate ratios using argument \code{model.glmm} (which corresponds to
#' argument \code{model} in the \code{\link[metafor]{rma.glmm}}
#' function):
#' \tabular{cl}{
#' 1. \tab Poisson regression model with fixed study effects (default)
#'  \cr
#'  \tab (\code{model.glmm = "UM.FS"}, i.e., \bold{U}nconditional
#'  \bold{M}odel - \bold{F}ixed \bold{S}tudy effects) \cr
#' 2. \tab Mixed-effects Poisson regression model with random study
#'  effects \cr
#'  \tab (\code{model.glmm = "UM.RS"}, i.e., \bold{U}nconditional
#'  \bold{M}odel - \bold{R}andom \bold{S}tudy effects) \cr
#' 3. \tab Generalised linear mixed model (conditional Poisson-Normal)
#'  \cr
#'  \tab (\code{model.glmm = "CM.EL"}, i.e., \bold{C}onditional
#'   \bold{M}odel - \bold{E}xact \bold{L}ikelihood)
#' }
#'
#' Details on these three GLMMs as well as additional arguments which
#' can be provided using argument '\code{\dots}' in \code{metainc}
#' are described in \code{\link[metafor]{rma.glmm}} where you can also
#' find information on the iterative algorithms used for estimation.
#' Note, regardless of which value is used for argument
#' \code{model.glmm}, results for two different GLMMs are calculated:
#' common effect model (with fixed treatment effect) and random effects
#' model (with random treatment effects).
#' }
#' 
#' \subsection{Continuity correction}{
#'
#' Four approaches are available to apply a continuity correction:
#' \itemize{
#' \item Only studies with a zero cell count (\code{method.incr =
#'   "only0", default})
#' \item All studies if at least one study has a zero cell count
#'   (\code{method.incr = "if0all"})
#' \item All studies irrespective of zero cell counts
#'   (\code{method.incr = "all"})
#' \item Use values provided in arguments \code{incr.e} and \code{incr.c}
#'   (\code{method.incr = "user"})
#' }
#'
#' For studies with a zero cell count, by default, 0.5 is added to all
#' cell frequencies of these studies (argument \code{incr}). This
#' continuity correction is used both to calculate individual study
#' results with confidence limits and to conduct meta-analysis based
#' on the inverse variance method. For Mantel-Haenszel method, Cochran
#' method, and GLMMs, nothing is added to zero cell counts.
#' Accordingly, estimates for these methods are not defined if the
#' number of events is zero in all studies either in the experimental
#' or control group.
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
#' A prediction interval will only be shown if \code{prediction =
#' TRUE}.
#' }
#' 
#' @return
#' An object of class \code{c("metainc", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{metabin}},
#'   \code{\link{update.meta}}, \code{\link{print.meta}}
#' 
#' @references
#' Bayne-Jones S et al. (1964):
#' Smoking and Health: Report of the Advisory Committee to the Surgeon
#' General of the United States.
#' U-23 Department of Health, Education, and Welfare.
#' Public Health Service Publication No. 1103.
#' 
#' Greenland S & Robins JM (1985):
#' Estimation of a common effect parameter from sparse follow-up data.
#' \emph{Biometrics},
#' \bold{41}, 55--68
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
#' data(smoking)
#' m1 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'   data = smoking, studlab = study)
#' print(m1, digits = 2)
#' 
#' m2 <- update(m1, method = "Cochran")
#' print(m2, digits = 2)
#' 
#' data(lungcancer)
#' m3 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'   data = lungcancer, studlab = study)
#' print(m3, digits = 2)
#' 
#' # Redo Cochran meta-analysis with inflated standard errors
#' #
#' # All cause mortality
#' #
#' TEa <- log((smoking$d.smokers/smoking$py.smokers) /
#'   (smoking$d.nonsmokers/smoking$py.nonsmokers))
#' seTEa <- sqrt(1 / smoking$d.smokers + 1 / smoking$d.nonsmokers +
#'   2.5 / smoking$d.nonsmokers)
#' metagen(TEa, seTEa, sm = "IRR", studlab = smoking$study)
#' 
#' # Lung cancer mortality
#' #
#' TEl <- log((lungcancer$d.smokers/lungcancer$py.smokers) /
#'   (lungcancer$d.nonsmokers/lungcancer$py.nonsmokers))
#' seTEl <- sqrt(1 / lungcancer$d.smokers + 1 / lungcancer$d.nonsmokers +
#'   2.25 / lungcancer$d.nonsmokers)
#' metagen(TEl, seTEl, sm = "IRR", studlab = lungcancer$study)
#'
#' \dontrun{
#' # Meta-analysis using generalised linear mixed models
#' # (only if R packages 'metafor' and 'lme4' are available)
#' 
#' # Poisson regression model (fixed study effects)
#' #
#' m4 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'   data = smoking, studlab = study, method = "GLMM")
#' m4
#' 
#' # Mixed-effects Poisson regression model (random study effects)
#' #
#' update(m4, model.glmm = "UM.RS", nAGQ = 1)
#' #
#' # Generalised linear mixed model (conditional Poisson-Normal)
#' #
#' update(m4, model.glmm = "CM.EL")
#' }
#' 
#' @export metainc


metainc <- function(event.e, time.e, event.c, time.c, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    cluster = NULL, rho = 0,
                    #
                    weights = NULL,
                    weights.common = weights, weights.random = weights,
                    #
                    method = if (sm == "IRSD") "Inverse" else "MH",
                    sm = gs("sminc"),
                    incr = gs("incr"), method.incr = gs("method.incr"),
                    incr.e = if (length(incr) > 1) incr else NULL,
                    incr.c = if (length(incr) > 1) incr else NULL,
                    model.glmm = "UM.FS",
                    ##
                    level = gs("level"),
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
                    method.tau =
                      ifelse(!is.na(charmatch(tolower(method), "glmm",
                                              nomatch = NA)),
                             "ML", gs("method.tau")),
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
                    method.bias = gs("method.bias"),
                    ##
                    n.e = NULL, n.c = NULL,
                    ##
                    backtransf = if (sm == "IRSD") FALSE else gs("backtransf"),
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
                    label.e = gs("label.e"), label.c = gs("label.c"),
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
  ## (1) Check arguments
  ##
  ##
  
  args <- list(...)
  nam.args <- names(args)
  #
  missing.sm <- missing(sm)
  missing.subgroup <- missing(subgroup)
  missing.byvar <- missing(byvar)
  missing.overall <- missing(overall)
  missing.overall.hetstat <- missing(overall.hetstat)
  missing.test.subgroup <- missing(test.subgroup)
  #
  missing.event.c <- missing(event.c)
  missing.time.e <- missing(time.e)
  missing.time.c <- missing(time.c)
  missing.n.e <- missing(n.e)
  missing.n.c <- missing(n.c)
  #
  missing.studlab <- missing(studlab)
  #
  missing.incr <- missing(incr)
  avail.incr <- !missing.incr && !is.null(incr)
  #
  missing.method.incr <- missing(method.incr)
  avail.method.incr <- !missing.method.incr && !is.null(method.incr)
  #
  missing.incr.e <- missing(incr.e)
  missing.incr.c <- missing(incr.c)
  #
  missing.method.tau <- missing(method.tau)
  missing.tau.common <- missing(tau.common)
  missing.method.predict <- missing(method.predict)
  missing.method <- missing(method)
  missing.level.ma <- missing(level.ma)
  missing.common <- missing(common)
  missing.random <- missing(random)
  missing.method.common.ci <- missing(method.common.ci)
  missing.method.random.ci <- missing(method.random.ci)
  #
  missing.hakn <- missing(hakn)
  missing.adhoc.hakn.ci <- missing(adhoc.hakn.ci)
  missing.adhoc.hakn <- missing(adhoc.hakn)
  #
  missing.subgroup.name <- missing(subgroup.name)
  missing.print.subgroup.name <- missing(print.subgroup.name)
  missing.sep.subgroup <- missing(sep.subgroup)
  #
  missing.label.e <- missing(label.e)
  missing.label.c <- missing(label.c)
  missing.complab <- missing(complab)
  #
  missing.cluster <- missing(cluster)
  #
  chknumeric(rho, min = -1, max = 1)
  ##
  chknull(sm)
  sm <- setchar(sm, gs("sm4inc"))
  ##
  chklevel(level)
  ##
  method.tau <- setchar(method.tau, gs("meth4tau"))
  ##
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
  method.tau.ci <- setchar(method.tau.ci, gs("meth4tau.ci"))
  ##
  tau.common <- replaceNULL(tau.common, FALSE)
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  method.predict <- setchar(method.predict, gs("meth4pi"))
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
  chklogical(backtransf)
  ##
  chknumeric(irscale, length = 1)
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
  fun <- "metainc"
  ##
  if (sm != "IRD" & irscale != 1) {
    warning("Argument 'irscale' only considered for ",
            "incidence rate differences.",
            call. = FALSE)
    irscale <- 1
  }
  ##
  method <- setchar(method, gs("meth4inc"))
  #
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
  method.incr <- setchar(method.incr, gs("meth4incr"))
  ##
  is.glmm <- method == "GLMM"
  ##
  model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL"))
  chklogical(warn)
  ##
  ## Check for deprecated arguments in '...'
  ##
  chklogical(warn.deprecated)
  ##
  level.ma <- deprecated(level.ma, missing.level.ma, args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
  #
  method.I2 <- setchar(method.I2, gs("meth4i2"))
  #
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing.random, args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  method.random.ci <-
    deprecated2(method.random.ci, missing.method.random.ci,
                hakn, missing.hakn,
                warn.deprecated)
  if (is.logical(method.random.ci))
    if (method.random.ci)
      method.random.ci <- "HK"
    else
      method.random.ci <- "classic"
  method.random.ci <- setchar(method.random.ci, gs("meth4random.ci"))
  ##
  adhoc.hakn.ci <-
    deprecated2(adhoc.hakn.ci, missing.adhoc.hakn.ci,
                adhoc.hakn, missing.adhoc.hakn, warn.deprecated)
  adhoc.hakn.ci <- setchar(replaceNA(adhoc.hakn.ci, ""), gs("adhoc4hakn.ci"))
  #
  subgroup.name <-
    deprecated(subgroup.name, missing.subgroup.name, args, "bylab",
               warn.deprecated)
  ##
  print.subgroup.name <-
    deprecated(print.subgroup.name, missing.print.subgroup.name,
               args, "print.byvar", warn.deprecated)
  print.subgroup.name <-
    replaceNULL(print.subgroup.name, gs("print.subgroup.name"))
  chklogical(print.subgroup.name)
  ##
  sep.subgroup <-
    deprecated(sep.subgroup, missing.sep.subgroup, args, "byseparator",
               warn.deprecated)
  if (!is.null(sep.subgroup))
    chkchar(sep.subgroup, length = 1)
  ##
  ## Some more checks
  ##
  chklogical(overall)
  chklogical(overall.hetstat)
  #
  # Ignore deprecated arguments 'addincr' and 'allincr'
  #
  txt.ignore <- "(deprecated); use argument 'method.incr'"
  #
  if (!is.na(charmatch("addincr", nam.args)))
    warn_ignore_input(addincr, TRUE, txt.ignore)
  if (!is.na(charmatch("allincr", nam.args)))
    warn_ignore_input(allincr, TRUE, txt.ignore)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (nulldata) {
    data <- sfsp
    data.pairwise <- FALSE
  }
  else
    data.pairwise <- inherits(data, "pairwise")
  #
  # Catch 'event.e', 'time.e', 'event.c', 'time.c', 'n.e', 'n.c',
  # 'incr.e', 'incr.c', studlab', and 'subgroup' from data:
  #
  event.e <- catch("event.e", mc, data, sfsp)
  chknull(event.e)
  #
  if (inherits(event.e, "pairwise")) {
    is.pairwise <- TRUE
    #
    type <- attr(event.e, "type")
    if (type != "count")
      stop("Wrong type for pairwise() object: '", type, "'.", call. = FALSE)
    #
    txt.ignore <- "as first argument is a pairwise object"
    #
    warn_ignore_input(event.c, !missing.event.c, txt.ignore)
    warn_ignore_input(time.e, !missing.time.e, txt.ignore)
    warn_ignore_input(time.c, !missing.time.c, txt.ignore)
    warn_ignore_input(n.e, !missing.n.e, txt.ignore)
    warn_ignore_input(n.c, !missing.n.c, txt.ignore)
    warn_ignore_input(incr.e, !missing.incr.e, txt.ignore)
    warn_ignore_input(incr.c, !missing.incr.c, txt.ignore)
    #
    warn_ignore_input(studlab, !missing.studlab, txt.ignore)
    #
    if (missing.sm)
      sm <- attr(event.e, "sm")
    #
    if (missing.method)
      method <- attr(event.e, "method")
    #
    if (!avail.method.incr & !avail.incr)
      method.incr <- "user"
    #
    missing.sm <- FALSE
    missing.method <- FALSE
    #
    avail.method.incr <- TRUE
    missing.incr <- FALSE
    missing.method.incr <- FALSE
    #
    reference.group <- attr(event.e, "reference.group")
    #
    studlab <- event.e$studlab
    #
    treat1 <- event.e$treat1
    treat2 <- event.e$treat2
    #
    event.c <- event.e$event2
    #
    time.e <- event.e$time1
    time.c <- event.e$time2
    #
    n.e <- event.e$n1
    n.c <- event.e$n2
    #
    avail.n.e <- !is.null(n.e)
    avail.n.c <- !is.null(n.c)
    #
    incr.e <- event.e$incr1
    incr.c <- event.e$incr2
    #
    if (avail.incr | method.incr != "user") {
      incr.e <- NULL
      incr.c <- NULL
    }
    #
    if (!avail.incr)
      incr <- attr(event.e, "incr")
    #
    pairdata <- event.e
    data <- event.e
    nulldata <- FALSE
    #
    event.e <- event.e$event1
    #
    wo <- treat1 == reference.group
    #
    if (any(wo)) {
      ttreat1 <- treat1
      treat1[wo] <- treat2[wo]
      treat2[wo] <- ttreat1[wo]
      #
      tevent.e <- event.e
      event.e[wo] <- event.c[wo]
      event.c[wo] <- tevent.e[wo]
      #
      ttime.e <- time.e
      time.e[wo] <- time.c[wo]
      time.c[wo] <- ttime.e[wo]
      #
      if (avail.n.e & avail.n.c) {
        tn.e <- n.e
        n.e[wo] <- n.c[wo]
        n.c[wo] <- tn.e[wo]
      }
      #
      if (!(is.null(incr.e) | is.null(incr.c))) {
        tincr.e <- incr.e
        incr.e[wo] <- incr.c[wo]
        incr.c[wo] <- tincr.e[wo]
      }
    }
    #
    if (missing.subgroup & missing.byvar) {
      subgroup <- paste(treat1, treat2, sep = " vs ")
      #
      if (length(unique(subgroup)) == 1) {
        if (missing.complab)
          complab <- unique(subgroup)
        #
        if (missing.label.e)
          label.e <- unique(treat1)
        if (missing.label.c)
          label.c <- unique(treat2)
        #
        subgroup <- NULL
      }
      else {
        if (missing.overall)
          overall <- FALSE
        if (missing.overall.hetstat)
          overall.hetstat <- FALSE
        if (missing.test.subgroup)
          test.subgroup <- FALSE
      }
    }
    else {
      subgroup <- catch("subgroup", mc, data, sfsp)
      byvar <- catch("byvar", mc, data, sfsp)
      #
      subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                              warn.deprecated)
    }
  }
  else {
    is.pairwise <- FALSE
    #
    if (missing.sm && !is.null(data) && !is.null(attr(data, "sm")))
      sm <- attr(data, "sm")
    #
    time.e <- catch("time.e", mc, data, sfsp)
    n.e <- catch("n.e", mc, data, sfsp)
    #
    event.c <- catch("event.c", mc, data, sfsp)
    time.c <- catch("time.c", mc, data, sfsp)
    n.c <- catch("n.c", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    subgroup <- catch("subgroup", mc, data, sfsp)
    byvar <- catch("byvar", mc, data, sfsp)
    #
    subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                            warn.deprecated)
    #
    if (!missing.incr)
      incr <- catch("incr", mc, data, sfsp)
    #
    if (!missing.incr.e)
      incr.e <- catch("incr.e", mc, data, sfsp)
    if (!missing.incr.c)
      incr.c <- catch("incr.c", mc, data, sfsp)
  }
  #
  method.bias <- setmethodbias(method.bias)
  #
  chknumeric(incr, min = 0)
  #
  avail.n.e <- !is.null(n.e)
  avail.n.c <- !is.null(n.c)
  #
  avail.incr.e <- !is.null(incr.e)
  avail.incr.c <- !is.null(incr.c)
  avail.incr.both <- avail.incr.e & avail.incr.c
  #
  if (avail.incr.e + avail.incr.c == 1)
    stop("Arguments 'incr.e' and 'incr.c' are required together.",
         call. = FALSE)
  #
  if (avail.incr.e)
    chknumeric(incr.e, min = 0, NA.ok = FALSE)
  #
  if (avail.incr.c)
    chknumeric(incr.c, min = 0, NA.ok = FALSE)
  #
  if (avail.incr.both) {
    if (!avail.method.incr)
      method.incr <- "user"
    #
    txt.ignore <- "as arguments 'incr.e' and 'incr.c' are provided"
    #
    warn_set_input(method.incr, method.incr != "user", '"user"', txt.ignore)
    #
    if (length(incr) != length(incr.e) | any(incr != incr.e))
      warn_ignore_input(incr, avail.incr, txt.ignore)
  }
  else {
    addincr <- allincr <- FALSE
    #
    if (method.incr == "all")
      addincr <- TRUE
    else if (method.incr == "if0all")
      allincr <- TRUE
  }
  #
  k.All <- length(event.e)
  #
  chknull(time.e)
  chknull(event.c)
  chknull(time.c)
  #
  studlab <- setstudlab(studlab, k.All)
  #
  by <- !is.null(subgroup)
  #
  # Catch 'subset', 'exclude' and 'cluster' from data:
  #
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
  #
  # Check variable values
  #
  chknumeric(event.e, 0)
  chknumeric(time.e, 0, zero = TRUE)
  chknumeric(event.c, 0)
  chknumeric(time.c, zero = TRUE)
  #
  # Recode integer as numeric:
  #
  event.e <- int2num(event.e)
  time.e  <- int2num(time.e)
  event.c <- int2num(event.c)
  time.c  <- int2num(time.c)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  
  chklength(time.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(time.c, k.All, fun)
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
  #
  if (avail.incr.e) {
    if (length(incr.e) == 1)  
      incr.e <- rep_len(incr.e, k.All)
    else
      chklength(incr.e, k.All, fun)
  }
  #
  if (avail.incr.c) {
    if (length(incr.c) == 1)  
      incr.c <- rep_len(incr.c, k.All)
    else
      chklength(incr.c, k.All, fun)
  }
  #
  if (avail.n.e)
    chklength(n.e, k.All, fun)
  #
  if (avail.n.c)
    chklength(n.c, k.All, fun)
  #
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
  #
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as ",
            "argument tau.preset is not NULL.",
            call. = FALSE)
    tau.common <- TRUE
  }
  #
  if (method.incr == "user" & !avail.incr.both)
    stop("Arguments 'incr.e' and 'incr.c' must be provided if ",
         "'method.incr = \"user\".",
         call. = FALSE)
  
  
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
      data <- data.frame(.studlab = studlab)
    else
      data$.studlab <- studlab
    ##
    data$.event.e <- event.e
    data$.time.e <- time.e
    data$.event.c <- event.c
    data$.time.c <- time.c
    #
    data$.incr.e <- NA
    data$.incr.c <- NA
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
    #
    if (with.cluster)
      data$.id <- data$.cluster <- cluster
    #
    if (usw.common)
      data$.weights.common <- weights.common
    #
    if (usw.random)
      data$.weights.random <- weights.random
    #
    if (avail.n.e)
      data$.n.e <- n.e
    if (avail.n.e)
      data$.n.c <- n.c
  }  
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  
  if (!missing.subset) {
    event.e <- event.e[subset]
    time.e <- time.e[subset]
    event.c <- event.c[subset]
    time.c <- time.c[subset]
    studlab <- studlab[subset]
    ##
    cluster <- cluster[subset]
    exclude <- exclude[subset]
    #
    weights.common <- weights.common[subset]
    weights.random <- weights.random[subset]
    #
    if (avail.n.e)
      n.e <- n.e[subset]
    if (avail.n.c)
      n.c <- n.c[subset]
    #
    incr.e <- incr.e[subset]
    incr.c <- incr.c[subset]
    #
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (by)
      subgroup <- subgroup[subset]
  }
  #
  if (missing.subgroup & is.pairwise & by) {
    if (length(unique(subgroup)) == 1) {
      by <- FALSE
      #
      if (missing.complab)
        complab <- unique(subgroup)
      #
      subgroup <- NULL
      #
      if (keepdata)
        data$.subgroup <- NULL
      #
      if (missing.overall)
        overall <- TRUE
      if (missing.overall.hetstat)
        overall.hetstat <- TRUE
    }
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event.e)
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
                IRD = event.e == 0 | event.c == 0,
                IRR = event.e == 0 | event.c == 0,
                VE = event.e == 0 | event.c == 0,
                IRSD = event.e == 0 | event.c == 0)
  #
  # Sparse computation
  #
  sparse <- any(sel, na.rm = TRUE)
  #
  if (avail.incr.both) {
    chknumeric(incr.e, min = 0, NA.ok = FALSE)
    chknumeric(incr.c, min = 0, NA.ok = FALSE)
  }
  else {
    if (addincr)
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
    incr.e <- incr.c <- incr.event
  }
  #
  if (keepdata) {
    if (missing.subset) {
      data$.incr.e <- incr.e
      data$.incr.c <- incr.c
    }
    else {
      data$.incr.e <- NA
      data$.incr.c <- NA
      #
      data$.incr.e[subset] <- incr.e
      data$.incr.c[subset] <- incr.c
    }
  }
  
  
  ##
  ##
  ## (8) Calculate results for individual studies
  ##
  ##
  
  if (sm %in% c("IRR", "VE")) {
    TE <- log(((event.e + incr.e) / time.e) / ((event.c + incr.c) / time.c))
    seTE <- sqrt(1 / (event.e + incr.e) + 1 / (event.c + incr.c))
  }
  else if (sm == "IRD") {
    TE <- event.e / time.e - event.c / time.c
    seTE <- sqrt((event.e + incr.e) / time.e^2 + (event.c + incr.c) / time.c^2)
  }
  else if (sm == "IRSD") {
    TE <- sqrt(event.e / time.e) - sqrt(event.c / time.c)
    seTE <- sqrt(0.25 / time.e + 0.25 / time.c)
  }
  #
  # Set NaN to NA
  #
  TE[is.nan(TE)] <- NA
  
  
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
  ## (9) Additional checks for GLMMs
  ##
  ##
  
  if (is.glmm) {
    chkglmm(sm, method.tau, method.random.ci, method.predict,
            adhoc.hakn.ci, adhoc.hakn.pi,
            c("IRR", "VE"))
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
  ## (10) Do meta-analysis
  ##
  ##
  
  k <- sum(!is.na(event.e[!exclude]) & !is.na(event.c[!exclude]) &
           !is.na(time.e[!exclude]) & !is.na(time.c[!exclude]))
  ##
  for (i in seq_along(method.random.ci))
    if (k == 1 & method.random.ci[i] == "HK")
      method.random.ci[i] <- "classic"
  ##
  if (method == "MH") {
    ##
    ## Greenland, Robins (1985)
    ## 
    x.k <- event.e
    y.k <- event.c
    n.k <- time.e
    m.k <- time.c
    ##
    N.k <- n.k + m.k
    t.k <- x.k + y.k
    ##
    if (sm %in% c("IRR", "VE")) {
      D <- n.k * m.k * t.k / N.k^2
      R <- x.k * m.k / N.k
      S <- y.k * n.k / N.k
      ##
      D[exclude] <- R[exclude] <- S[exclude] <- 0
      ##
      w.common <- S
      TE.common <- log(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE))
      seTE.common <- sqrt(sum(D, na.rm = TRUE) / (sum(R, na.rm = TRUE) *
                                                  sum(S, na.rm = TRUE)))
    }
    else if (sm == "IRD") {
      L <- (x.k * m.k^2 + y.k * n.k^2) / N.k^2
      S <- n.k * m.k / N.k
      ##
      L[exclude] <- S[exclude] <- 0
      ##
      w.common <- S
      TE.common <- weighted.mean(TE, w.common, na.rm = TRUE)
      seTE.common <- sqrt(sum(L, na.rm = TRUE) / sum(S, na.rm = TRUE)^2)
    }
  }
  ##
  else if (method == "Cochran") {
    ##
    ## Smoking and Health - Report of the Advisory Committee to the
    ## Surgeon General of the Public Health Service,
    ## Chapter 8
    ## 
    if (sm %in% c("IRR", "VE")) {
      w.common <- event.c * time.e / time.c
      w.common[exclude] <- 0
      TE.common <- weighted.mean(TE, w.common)
      seTE.common <- sqrt(1 / sum(event.e) + 1 / sum(event.c))
    }
    else if (sm == "IRD") {
      warning("Cochran method only available for ",
              "Incidence Rate Ratio (sm = \"IRR\") ",
              "and Vaccine Efficacy / Effectiveness (sm = \"VE\")",
              call. = FALSE)
      return(NULL)
    }
  }
  else if (is.glmm) {
    list.inc <- list(x1i = event.e[!exclude], t1i = time.e[!exclude],
                     x2i = event.c[!exclude], t2i = time.c[!exclude],
                     measure = "IRR", model = model.glmm)
    ##
    use.random <-
      sum(!exclude) > 1 &
      !((sum(event.e[!exclude], na.rm = TRUE) == 0 &
         sum(event.c[!exclude], na.rm = TRUE) == 0) |
        (!any(event.e[!exclude] != time.e[!exclude]) |
         !any(event.c[!exclude] != time.c[!exclude])))
    ##
    res.glmm <-
      runGLMM(list.inc,
              method.tau = method.tau,
              method.random.ci = method.random.ci,
              level = level.ma,
              control = control, use.random = use.random,
              warn = warn)
    ##
    TE.common   <- as.numeric(res.glmm$glmm.common$b)
    seTE.common <- as.numeric(res.glmm$glmm.common$se)
    ##
    w.common <- rep(NA, length(event.e))
  }
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
               label.e = label.e, label.c = label.c,
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
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   if (method == "Inverse") TE.tau else TE.common,
                   method.I2, level.hetstat, subgroup, control)
  
  
  ##
  ##
  ## (11) Generate R object
  ##
  ##
  
  res <- list(event.e = event.e, time.e = time.e,
              event.c = event.c, time.c = time.c,
              method = method, method.random = method,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              method.incr = method.incr,
              sparse = sparse,
              incr.e = incr.e,
              incr.c = incr.c,
              k.MH = if (method == "MH") sum(w.common > 0) else NA)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$method <- NULL
  m$method.random <- NULL
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
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$n.e <- n.e
  res$n.c <- n.c
  res$TE.tau <- TE.tau
  ##
  res$irscale <- irscale
  res$irunit  <- irunit
  #
  if (is.pairwise | data.pairwise) {
    res$pairwise <- TRUE
    res$k.study <- length(unique(res$studlab[!is.na(res$TE)]))
  }
  #
  res$call <- match.call()
  #
  if (method %in% c("MH", "Cochran", "GLMM")) {
    res <- ci2meta(res, ci.c = ci(TE.common, seTE.common, level = level.ma))
    res$w.common <- w.common
  }
  ##
  if (is.glmm) {
    res$method.tau <- method.tau
    #
    res <- addGLMM(res, res.glmm, method.I2)
    res$model.glmm <- model.glmm
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
          runGLMM(list.inc,
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
    res$n.w <- NULL
    res$event.w <- NULL
    ##
    if (!avail.n.e)
      res$n.e.w <- NULL
    if (!avail.n.c)
      res$n.c.w <- NULL
    ##
    res$n.harmonic.mean.w <- NULL
    res$t.harmonic.mean.w <- NULL
    ##
    res <- setNAwithin(res, res$three.level | is.glmm)
  }
  ##
  ## Mantel-Haenszel and Cochran method are common effect methods
  ##
  if (res$method.random %in% c("MH", "Cochran"))
    res$method.random <- "Inverse"
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
