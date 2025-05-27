#' Meta-analysis of correlations
#' 
#' @description
#' Calculation of common effect and random effects estimates for
#' meta-analyses with correlations; inverse variance weighting is used
#' for pooling.
#' 
#' @param cor Correlations.
#' @param n Number of observations.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., cor and n.
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
#' @param sm A character string indicating which summary measure
#'   (\code{"ZCOR"} or \code{"COR"}) is to be used for pooling of
#'   studies.
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
#' @param backtransf A logical indicating whether results for Fisher's
#'   z transformed correlations (\code{sm = "ZCOR"}) should be back
#'   transformed in printouts and plots. If TRUE (default), results
#'   will be presented as correlations; otherwise Fisher's z
#'   transformed correlations will be shown.
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
#' @param adhoc.hakn Deprecated argument (replaced by 'adhoc.hakn.ci').
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}}.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' This function conducts common effect and random effects meta-analysis of
#' correlations based either on Fisher's z transformation of correlations
#' (\code{sm = "ZCOR"}) or direct combination of (untransformed) correlations
#' (\code{sm = "COR"}) (see Cooper et al., 2009, p264-5 and p273-4). Note, the
#' input to argument \code{cor} is always correlations and not Fisher's z
#' transformed correlations if \code{sm = "ZCOR"}.
#' 
#' Only few statisticians would advocate the use of untransformed correlations
#' unless sample sizes are very large (see Cooper et al., 2009, p265). The
#' artificial example given below shows that the smallest study gets the largest
#' weight if correlations are combined directly because the correlation is
#' closest to 1.
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
#' \code{common} and \code{random}. E.g. functions
#' \code{\link{print.meta}} and \code{\link{forest.meta}} will not
#' print results for the random effects model if \code{random =
#' FALSE}.
#'
#' A prediction interval will only be shown if \code{prediction =
#' TRUE}.
#' }
#' 
#' @note
#' The function \code{\link{metagen}} is called internally to
#' calculate individual and overall treatment estimates and standard
#' errors.
#' 
#' @return
#' An object of class \code{c("metacor", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{metacont}}, \code{\link{metagen}},
#'   \code{\link{print.meta}}
#' 
#' @references
#' Cooper H, Hedges LV, Valentine JC (2009):
#' \emph{The Handbook of Research Synthesis and Meta-Analysis},
#' 2nd Edition.
#' New York: Russell Sage Foundation
#'
#' Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca J (2013):
#' Three-level meta-analysis of dependent effect sizes.
#' \emph{Behavior Research Methods},
#' \bold{45}, 576--94
#' 
#' @examples
#' m1 <- metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
#' 
#' # Print correlations (back transformed from Fisher's z
#' # transformation)
#' #
#' summary(m1)
#' 
#' # Print Fisher's z transformed correlations 
#' #
#' print(summary(m1), backtransf = FALSE)
#' 
#' # Forest plot with back transformed correlations
#' #
#' forest(m1)
#' 
#' # Forest plot with Fisher's z transformed correlations
#' #
#' forest(m1, backtransf = FALSE)
#' 
#' m2 <- update(m1, sm = "cor")
#' summary(m2)
#'
#' \dontrun{
#' # Identical forest plots (as back transformation is the identity
#' # transformation)
#' forest(m2)
#' forest(m2, backtransf = FALSE)
#' }
#' 
#' @export metacor


metacor <- function(cor, n, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    cluster = NULL, rho = 0,
                    #
                    weights = NULL,
                    weights.common = weights, weights.random = weights,
                    #
                    sm = gs("smcor"),
                    level = gs("level"),
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
                    method.tau = gs("method.tau"),
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
                    null.effect = 0,
                    ##
                    method.bias = gs("method.bias"),
                    ##
                    backtransf = gs("backtransf"),
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
                    byvar, adhoc.hakn,
                    ##
                    keepdata = gs("keepdata"),
                    warn.deprecated = gs("warn.deprecated"),
                    ##
                    control = NULL,
                    ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chknumeric(rho, min = -1, max = 1)
  ##
  chknull(sm)
  chklevel(level)
  #
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  missing.method.tau <- missing(method.tau)
  method.tau <- setchar(method.tau, gs("meth4tau"))
  ##
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
    deprecated(method.random.ci, missing(method.random.ci),
               args, "hakn", warn.deprecated)
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
  ## Some more checks
  ##
  chklogical(overall)
  chklogical(overall.hetstat)
  ##
  ## Additional arguments / checks for metacor objects
  ##
  fun <- "metacor"
  sm <- setchar(sm, gs("sm4cor"))
  
  
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
  ## Catch 'cor' and 'n' from data:
  ##
  cor <- catch("cor", mc, data, sfsp)
  chknull(cor)
  k.All <- length(cor)
  ##
  n <- catch("n", mc, data, sfsp)
  chknull(n)
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
  subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar)
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
  ## Additional checks
  ##
  if (!by & tau.common) {
    warning("Value for argument 'tau.common' set to FALSE as ",
            "argument 'subgroup' is missing.")
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as ",
            "argument tau.preset is not NULL.")
    tau.common <- TRUE
  }
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
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(n, k.All, fun)
  chklength(studlab, k.All, fun)
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
  if (by) {
    chklength(subgroup, k.All, fun)
    chklogical(test.subgroup)
    chklogical(prediction.subgroup)
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
      data <- data.frame(.cor = cor)
    else
      data$.cor <- cor
    ##
    data$.n <- n
    data$.studlab <- studlab
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
    cor <- cor[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    ##
    cluster <- cluster[subset]
    exclude <- exclude[subset]
    #
    weights.common <- weights.common[subset]
    weights.random <- weights.random[subset]
    #
    if (by)
      subgroup <- subgroup[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(cor)
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
  chknumeric(cor, -1, 1)
  chknumeric(n, 0, zero = TRUE)
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
  if (sm == "ZCOR") {
    TE   <- cor2z(cor)
    seTE <- sqrt(1 / (n - 3))
    transf.null.effect <- 0.5 * log((1 + null.effect) / (1 - null.effect))
  }
  if (sm == "COR") {
    TE <- cor
    seTE <- sqrt((1 - cor^2)^2 / (n - 1))
    transf.null.effect <- null.effect
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
    chkmlm(method.tau, missing.method.tau, method.predict)
    ##
    common <- FALSE
    ##
    if (!(method.tau %in% c("REML", "ML")))
      method.tau <- "REML"
  }
  
  
  ##
  ##
  ## (9) Do meta-analysis
  ##
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
               method.tau = method.tau, method.tau.ci = method.tau.ci,
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
               warn = FALSE,
               ##
               control = control)
  #
  # Estimate common tau-squared across subgroups
  #
  if (by & tau.common)
    hcc <- hetcalc(TE, seTE, method.tau, "", TE.tau,
                   method.I2, level.hetstat, subgroup, control)
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(cor = cor, n = n)
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
  m$warn <- NULL
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Add data
  ##
  res$pairwise <- FALSE
  #
  res$call <- match.call()
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
      if (res$three.level)
        res <- c(res,
                 subgroup(res, NULL,
                          factor(res$subgroup, bylevs(res$subgroup))))
      else
        res <-
          c(res, subgroup(res, hcc$tau.resid, seed = seed.predict.subgroup))
    }
    ##
    if (tau.common && is.null(tau.preset))
      res <- addHet(res, hcc)
    ##
    res$event.w <- NULL
    ##
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$n.harmonic.mean.w <- NULL
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    ##
    res$time.e.w <- NULL
    res$time.c.w <- NULL
    res$t.harmonic.mean.w <- NULL
    ##
    res <- setNAwithin(res, res$three.level)
  }
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
