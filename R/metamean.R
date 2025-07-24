#' Meta-analysis of single means
#' 
#' @description
#' Calculation of an overall mean from studies reporting a single mean
#' using the inverse variance method for pooling; inverse variance
#' weighting is used for pooling.
#' 
#' @param n Number of observations.
#' @param mean Estimated mean.
#' @param sd Standard deviation.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information.
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
#' @param median Median (used to estimate the mean and standard
#'   deviation).
#' @param q1 First quartile (used to estimate the mean and standard
#'   deviation).
#' @param q3 Third quartile (used to estimate the mean and standard
#'   deviation).
#' @param min Minimum (used to estimate the mean and standard
#'   deviation).
#' @param max Maximum (used to estimate the mean and standard
#'   deviation).
#' @param method.mean A character string indicating which method to
#'   use to approximate the mean from the median and other statistics
#'   (see Details).
#' @param method.sd A character string indicating which method to use
#'   to approximate the standard deviation from sample size, median,
#'   interquartile range and range (see Details).
#' @param approx.mean Approximation method to estimate means (see
#'   Details).
#' @param approx.sd Approximation method to estimate standard
#'   deviations (see Details).
#' @param method.ci A character string indicating which method is used
#'   to calculate confidence intervals for individual studies, see
#'   Details.
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
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots for \code{sm = "MLN"}. If
#'   TRUE (default), results will be presented as means; otherwise
#'   logarithm of means will be shown.
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
#' @param sm A character string indicating which summary measure
#'   (\code{"MRAW"} or \code{"MLN"}) is to be used for pooling of
#'   studies.
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
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard deviations).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}}.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' Common effect and random effects meta-analysis of single means to
#' calculate an overall mean; inverse variance weighting is used for
#' pooling. Note, you should use R function \code{\link{metacont}} to
#' compare means of pairwise comparisons instead of using
#' \code{metamean} for each treatment arm separately which will break
#' randomisation in randomised controlled trials.
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
#' The following transformations of means are implemented to
#' calculate an overall mean:
#' \itemize{
#' \item Raw, i.e. untransformed, means (\code{sm = "MRAW"}, default)
#' \item Log transformed means (\code{sm = "MLN"})
#' }
#' 
#' Calculations are conducted on the log scale if \code{sm =
#' "MLN"}. Accordingly, list elements \code{TE}, \code{TE.common}, and
#' \code{TE.random} contain the logarithm of means. In printouts and
#' plots these values are back transformed if argument
#' \code{backtransf = TRUE} (default).
#' 
#' \subsection{Approximate means from sample sizes, medians and other statistics}{
#' 
#' Missing means can be derived from
#' \enumerate{
#' \item sample size, median, interquartile range and range (arguments
#'   \code{n}, \code{median}, \code{q1}, \code{q3}, \code{min}, and
#'   \code{max}),
#' \item sample size, median and interquartile range (arguments
#'   \code{n}, \code{median}, \code{q1}, and \code{q3}), or
#' \item sample size, median and range (arguments \code{n},
#'   \code{median}, \code{min}, and \code{max}).
#' }
#' 
#' By default, methods described in Luo et al. (2018) are utilised
#' (argument \code{method.mean = "Luo"}):
#' \itemize{
#' \item equation (15) if sample size, median, interquartile range and 
#'   range are available,
#' \item equation (11) if sample size, median and interquartile range
#'   are available,
#' \item equation (7) if sample size, median and range are available.
#' }
#' 
#' Instead the methods described in Wan et al. (2014) are used if
#' argument \code{method.mean = "Wan"}:
#' \itemize{
#' \item equation (10) if sample size, median, interquartile range and 
#'   range are available,
#' \item equation (14) if sample size, median and interquartile range
#'   are available,
#' \item equation (2) if sample size, median and range are available.
#' }
#'
#' The following methods are also available to estimate means from
#' quantiles or ranges if R package \bold{estmeansd} is installed:
#' \itemize{
#' \item Method for Unknown Non-Normal Distributions (MLN) approach
#'   (Cai et al. (2021), argument \code{method.mean = "Cai"}),
#' \item Quantile Estimation (QE) method (McGrath et al. (2020),
#'   argument \code{method.mean = "QE-McGrath"})),
#' \item Box-Cox (BC) method (McGrath et al. (2020),
#'   argument \code{method.mean = "BC-McGrath"})).
#' }
#'
#' By default, missing means are replaced successively using
#' interquartile ranges and ranges (if available), interquartile
#' ranges (if available) and finally ranges. Argument
#' \code{approx.mean} can be used to overwrite this behaviour for each
#' individual study and treatment arm:
#' \itemize{
#' \item use means directly (entry \code{""} in argument
#'   \code{approx.mean});
#' \item median, interquartile range and range (\code{"iqr.range"});
#' \item median and interquartile range (\code{"iqr"});
#' \item median and range (\code{"range"}).
#' }
#' }
#'
#' \subsection{Approximate standard deviations from sample sizes, medians and other statistics}{
#' 
#' Missing standard deviations can be derived from
#' \enumerate{
#' \item sample size, median, interquartile range and range (arguments
#'   \code{n}, \code{median}, \code{q1}, \code{q3}, \code{min}, and
#'   \code{max}),
#' \item sample size, median and interquartile range (arguments
#'   \code{n}, \code{median}, \code{q1} and \code{q3}), or
#' \item sample size, median and range (arguments \code{n},
#'   \code{median}, \code{min} and \code{max}).
#' }
#' 
#' Wan et al. (2014) describe methods to estimate the standard
#' deviation from the sample size, median and additional
#' statistics. Shi et al. (2020) provide an improved estimate of the
#' standard deviation if the interquartile range and range are
#' available in addition to the sample size and median. Accordingly,
#' equation (11) in Shi et al. (2020) is the default (argument
#' \code{method.sd = "Shi"}), if the median, interquartile range and
#' range are provided. The method by Wan et al. (2014) is used if
#' argument \code{method.sd = "Wan"} and, depending on the sample
#' size, either equation (12) or (13) is used. If only the
#' interquartile range or range is available, equations (15) / (16)
#' and (7) / (9) in Wan et al. (2014) are used, respectively.
#'
#' The following methods are also available to estimate standard
#' deviations from quantiles or ranges if R package \bold{estmeansd}
#' is installed:
#' \itemize{
#' \item Method for Unknown Non-Normal Distributions (MLN) approach
#'   (Cai et al. (2021), argument \code{method.mean = "Cai"}),
#' \item Quantile Estimation (QE) method (McGrath et al. (2020),
#'   argument \code{method.mean = "QE-McGrath"})),
#' \item Box-Cox (BC) method (McGrath et al. (2020),
#'   argument \code{method.mean = "BC-McGrath"})).
#' }
#'
#' By default, missing standard deviations are replaced successively
#' using these method, i.e., interquartile ranges and ranges are used
#' before interquartile ranges before ranges. Argument
#' \code{approx.sd} can be used to overwrite this default for each
#' individual study and treatment arms:
#' \itemize{
#' \item sample size, median, interquartile range and range
#'   (\code{"iqr.range"});
#' \item sample size, median and interquartile range (\code{"iqr"});
#' \item sample size, median and range (\code{"range"}).
#' }
#' }
#' 
#' \subsection{Confidence intervals for individual studies}{
#' 
#' For untransformed means (argument \code{sm = "MRAW"}), the
#' confidence interval for individual studies can be based on the
#' \itemize{
#' \item standard normal distribution (\code{method.ci = "z"}, default), or
#' \item t-distribution (\code{method.ci = "t"}).
#' }
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
#' An object of class \code{c("metamean", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{metamean}}, \code{\link{metagen}}
#' 
#' @references
#' Cai S, Zhou J, Pan J (2021):
#' Estimating the sample mean and standard deviation from order
#' statistics and sample size in meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{30}, 2701--2719
#' 
#' Luo D, Wan X, Liu J, Tong T (2018):
#' Optimally estimating the sample mean from the sample size, median,
#' mid-range, and/or mid-quartile range.
#' \emph{Statistical Methods in Medical Research},
#' \bold{27}, 1785--805
#'
#' McGrath S, Zhao X, Steele R, et al. and the DEPRESsion Screening
#' Data (DEPRESSD) Collaboration (2020):
#' Estimating the sample mean and standard deviation from commonly
#' reported quantiles in meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{29}, 2520--2537
#' 
#' Shi J, Luo D, Weng H, Zeng X-T, Lin L, Chu H, et al. (2020):
#' Optimally estimating the sample standard deviation from the
#' five-number summary.
#' \emph{Research Synthesis Methods}.
#' 
#' Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca
#' J (2013):
#' Three-level meta-analysis of dependent effect sizes.
#' \emph{Behavior Research Methods},
#' \bold{45}, 576--94
#'
#' Wan X, Wang W, Liu J, Tong T (2014):
#' Estimating the sample mean and standard deviation from the sample
#' size, median, range and/or interquartile range.
#' \emph{BMC Medical Research Methodology},
#' \bold{14}, 135
#' 
#' @examples
#' m1 <- metamean(rep(100, 3), 1:3, rep(1, 3))
#' m1
#' 
#' m2 <- update(m1, sm = "MLN")
#' m2
#' 
#' # With test for overall mean equal to 2
#' #
#' update(m1, null.effect = 2)
#' update(m2, null.effect = 2)
#' 
#' # Print results without back-transformation
#' #
#' update(m1, backtransf = FALSE)
#' update(m2, backtransf = FALSE)
#' update(m1, null.effect = 2, backtransf = FALSE)
#' update(m2, null.effect = 2, backtransf = FALSE)
#' 
#' @export metamean


metamean <- function(n, mean, sd, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     cluster = NULL, rho = 0,
                     #
                     weights = NULL,
                     weights.common = weights, weights.random = weights,
                     #
                     median, q1, q3, min, max,
                     method.mean = "Luo", method.sd = "Shi",
                     approx.mean, approx.sd,
                     ##
                     sm = gs("smmean"),
                     method.ci = gs("method.ci.cont"),
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
                     null.effect = NA,
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
                     title = gs("title"), complab = gs("complab"), outclab = "",
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
                     warn = gs("warn"), warn.deprecated = gs("warn.deprecated"),
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
  ##
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
  adhoc.hakn.pi <- setchar(adhoc.hakn.pi, gs("adhoc4hakn.pi"))
  ##
  chknumeric(null.effect, length = 1)
  ##
  method.bias <- setmethodbias(method.bias)
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
  ## Additional arguments / checks for metamean objects
  ##
  fun <- "metamean"
  sm <- setchar(sm, gs("sm4mean"))
  if (sm != "MRAW")
    method.ci <- "z"
  method.ci <- setchar(method.ci, gs("ci4cont"))
  ##
  method.mean <-
    setchar(method.mean, c("Luo", "Wan", "Cai", "QE-McGrath", "BC-McGrath"))
  method.sd <-
    setchar(method.sd, c("Shi", "Wan", "Cai", "QE-McGrath", "BC-McGrath"))
  ##
  if (method.mean %in% c("Cai", "QE-McGrath", "BC-McGrath"))
    is_installed_package("estmeansd", argument = "method.mean",
                         value = method.mean)
  if (method.sd %in% c("Cai", "QE-McGrath", "BC-McGrath"))
    is_installed_package("estmeansd", argument = "method.sd",
                         value = method.sd)
  ##
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
  #
  method.I2 <- setchar(method.I2, gs("meth4i2"))
  #
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
  ## Catch 'n', 'mean', and 'sd' from data:
  ##
  missing.mean <- missing(mean)
  missing.sd <- missing(sd)
  missing.median <- missing(median)
  missing.q1 <- missing(q1)
  missing.q3 <- missing(q3)
  missing.min <- missing(min)
  missing.max <- missing(max)
  ##
  if (missing.mean & missing.median)
    stop("Provide either argument 'mean' or 'median'.",
         call. = FALSE)
  ##
  if (missing.sd &
      !((!missing.q1 & !missing.q3) |
        (!missing.min & !missing.max)))
    stop("Provide either argument 'sd' and ",
         "arguments 'q1' & 'q3' or 'min & 'max'.",
         call. = FALSE)
  ##
  n <- catch("n", mc, data, sfsp)
  chknull(n)
  k.All <- length(n)
  ##
  mean <- catch("mean", mc, data, sfsp)
  if (!missing.mean)
    chknull(mean)
  else
    mean <- rep(NA, k.All)
  ##
  sd <- catch("sd", mc, data, sfsp)
  if (!missing.sd)
    chknull(sd)
  else
    sd <- rep(NA, k.All)
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
  ##
  ## Catch 'median', 'q1', 'q3', 'min', 'max', 'approx.mean', and
  ## 'approx.sd', from data:
  ##
  median <- catch("median", mc, data, sfsp)
  ##
  q1 <- catch("q1", mc, data, sfsp)
  ##
  q3 <- catch("q3", mc, data, sfsp)
  ##
  min <- catch("min", mc, data, sfsp)
  ##
  max <- catch("max", mc, data, sfsp)
  ##
  avail.median <- !(missing.median || is.null(median))
  avail.q1 <- !(missing.q1 || is.null(q1))
  avail.q3 <- !(missing.q3 || is.null(q3))
  avail.min <- !(missing.min || is.null(min))
  avail.max <- !(missing.max || is.null(max))
  ##
  missing.approx.mean <- missing(approx.mean)
  approx.mean <- catch("approx.mean", mc, data, sfsp)
  avail.approx.mean <-
    !(missing.approx.mean || is.null(approx.mean)) && any(approx.mean != "")
  #
  missing.approx.sd <- missing(approx.sd)
  approx.sd <- catch("approx.sd", mc, data, sfsp)
  avail.approx.sd <-
    !(missing.approx.sd || is.null(approx.sd)) && any(approx.sd != "")
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  
  chklength(mean, k.All, fun)
  chklength(sd, k.All, fun)
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
  if (avail.median)
    chklength(median, k.All, fun)
  if (avail.q1)
    chklength(q1, k.All, fun)
  if (avail.q3)
    chklength(q3, k.All, fun)
  if (avail.min)
    chklength(min, k.All, fun)
  if (avail.max)
    chklength(max, k.All, fun)
  ##
  if (avail.approx.mean) {
    if (length(approx.mean) == 1)
      rep_len(approx.mean, k.All)
    else
      chklength(approx.mean, k.All, fun)
    ##
    approx.mean <- setchar(approx.mean, c("", "iqr.range", "iqr", "range"))
  }
  else
    approx.mean <- rep_len("", k.All)
  #
  if (avail.approx.sd) {
    if (length(approx.sd) == 1)
      rep_len(approx.sd, k.All)
    else
      chklength(approx.sd, k.All, fun)
    ##
    approx.sd <- setchar(approx.sd, c("", "iqr.range", "iqr", "range"))
  }
  else
    approx.sd <- rep_len("", k.All)
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
            "argument 'subgroup' is missing.")
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as ",
            "argument tau.preset is not NULL.")
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
      data <- data.frame(.n = n)
    else
      data$.n <- n
    ##
    data$.mean <- mean
    data$.sd <- sd
    data$.studlab <- studlab
    ##
    if (avail.median)
      data$.median <- median
    if (avail.q1)
      data$.q1 <- q1
    if (avail.q3)
      data$.q3 <- q3
    if (avail.min)
      data$.min <- min
    if (avail.max)
      data$.max <- max
    #
    data$.approx.mean <- approx.mean
    data$.approx.sd <- approx.sd
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
    n <- n[subset]
    mean <- mean[subset]
    sd <- sd[subset]
    studlab <- studlab[subset]
    ##
    cluster <- cluster[subset]
    exclude <- exclude[subset]
    #
    weights.common <- weights.common[subset]
    weights.random <- weights.random[subset]
    #
    if (avail.median)
      median <- median[subset]
    if (avail.q1)
      q1 <- q1[subset]
    if (avail.q3)
      q3 <- q3[subset]
    if (avail.min)
      min <- min[subset]
    if (avail.max)
      max <- max[subset]
    #
    approx.mean <- approx.mean[subset]
    approx.sd <- approx.sd[subset]
    ##
    if (by)
      subgroup <- subgroup[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(n)
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
  chknumeric(n)
  chknumeric(mean)
  chknumeric(sd)
  ##
  if (avail.median)
    chknumeric(median)
  if (avail.q1)
    chknumeric(q1)
  if (avail.q3)
    chknumeric(q3)
  if (avail.min)
    chknumeric(min)
  if (avail.max)
    chknumeric(max)
  ##
  ## Recode integer as numeric:
  ##
  n    <- int2num(n)
  mean <- int2num(mean)
  sd   <- int2num(sd)
  ##
  if (avail.median)
    median <- int2num(median)
  if (avail.q1)
    q1 <- int2num(q1)
  if (avail.q3)
    q3 <- int2num(q3)
  if (avail.min)
    min <- int2num(min)
  if (avail.max)
    max <- int2num(max)
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
  ## (7) Calculate means from other information
  ##
  ##
  
  if (!avail.approx.mean) {
    ##
    ## (a) Use IQR and range
    ##
    sel.NA <- is.na(mean)
    if (any(sel.NA) & avail.median &
        avail.q1 & avail.q3 &
        avail.min & avail.max) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3) &
        !is.na(min) & !is.na(max)
      approx.mean[j] <- "iqr.range"
      ##
      mean[j] <- mean_sd_iqr_range(n[j], median[j], q1[j], q3[j],
                                   min[j], max[j], method.mean)$mean
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA <- is.na(mean)
    if (any(sel.NA) & avail.median & avail.q1 & avail.q3) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3)
      approx.mean[j] <- "iqr"
      mean[j] <- mean_sd_iqr(n[j], median[j], q1[j], q3[j],
                             method.mean)$mean
    }
    ##
    ## (c) Use range
    ##
    sel.NA <- is.na(mean)
    if (any(sel.NA) & avail.median & avail.min & avail.max) {
      j <- sel.NA & !is.na(median) & !is.na(min) & !is.na(max)
      approx.mean[j] <- "range"
      mean[j] <- mean_sd_range(n[j], median[j], min[j], max[j],
                               method.mean)$mean
    }
  }
  else {
    j <- 0
    for (i in approx.mean) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        mean[j] <- mean_sd_iqr_range(n[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
      else if (i == "iqr")
        mean[j] <- mean_sd_iqr(n[j], median[j], q1[j], q3[j],
                               method.mean)$mean
      else if (i == "range")
        mean[j] <- mean_sd_range(n[j], median[j], min[j], max[j],
                                 method.mean)$mean
    }
  }
  
  
  ##
  ##
  ## (8) Calculate standard deviation from other information
  ##
  ##
  
  if (!avail.median) {
    median.sd <- mean
    avail.median <- TRUE
    export.median <- FALSE
  }
  else {
    median.sd <- median
    median.sd[is.na(median.sd)] <- mean[is.na(median.sd)]
    export.median <- TRUE
  }
  ##
  if (!avail.approx.sd) {
    ##
    ## (a) Use IQR and range
    ##
    sel.NA <- is.na(sd)
    if (any(sel.NA) &
        avail.q1 & avail.q3 &
        avail.min & avail.max) {
      j <- sel.NA & !is.na(median.sd) & !is.na(q1) & !is.na(q3) &
        !is.na(min) & !is.na(max)
      approx.sd[j] <- "iqr.range"
      ##
      sd[j] <- mean_sd_iqr_range(n[j], median.sd[j], q1[j], q3[j],
                                 min[j], max[j],
                                 method.sd = method.sd)$sd
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA <- is.na(sd)
    if (any(sel.NA) & avail.q1 & avail.q3) {
      j <- sel.NA & !is.na(median.sd) & !is.na(q1) & !is.na(q3)
      approx.sd[j] <- "iqr"
      sd[j] <- mean_sd_iqr(n[j], median.sd[j], q1[j], q3[j])$sd
    }
    ##
    ## (c) Use range
    ##
    sel.NA <- is.na(sd)
    if (any(sel.NA) & avail.min & avail.max) {
      j <- sel.NA & !is.na(median.sd) & !is.na(min) & !is.na(max)
      approx.sd[j] <- "range"
      sd[j] <- mean_sd_range(n[j], median.sd[j], min[j], max[j])$sd
    }
  }
  else {
    j <- 0
    for (i in approx.sd) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        sd[j] <- mean_sd_iqr_range(n[j], median.sd[j], q1[j], q3[j],
                                   min[j], max[j],
                                   method.sd = method.sd)$sd
      else if (i == "iqr")
        sd[j] <- mean_sd_iqr(n[j], median.sd[j], q1[j], q3[j])$sd
      else if (i == "range")
        sd[j] <- mean_sd_range(n[j], median.sd[j], min[j], max[j])$sd
    }
  }
  ##
  if (keepdata) {
    if (!isCol(data, ".subset")) {
      data$.sd <- sd
      data$.mean <- mean
      if (avail.approx.mean)
        data$.approx.mean <- approx.mean
      if (avail.approx.sd)
        data$.approx.sd <- approx.sd
    }
    else {
      data$.sd[data$.subset] <- sd
      data$.mean[data$.subset] <- mean
      if (avail.approx.mean)
        data$.approx.mean[data$.subset] <- approx.mean
      if (avail.approx.sd)
        data$.approx.sd[data$.subset] <- approx.sd
    }
  }
  
  
  ##
  ##
  ## (9) Calculate results for individual studies
  ##
  ##
  
  npn.n <- npn(n)
  ##
  if (any(npn.n) & warn)
    warning("Studies with non-positive sample size get no weight in ",
            "meta-analysis.")
  ##
  if (sm == "MRAW") {
    TE <- ifelse(npn.n, NA, mean)
    ##
    seTE <- ifelse(npn.n, NA, sqrt(sd^2 / n))
    ##
    seTE[is.na(TE)] <- NA
    ##
    if (method.ci == "t")
      ci.study <-
        ci(TE, seTE, level = level, df = n - 1, null.effect = null.effect)
    ##
    transf.null.effect <- null.effect
  }
  ##
  else if (sm == "MLN") {
    npn.mean <- npn(mean)
    ##
    if (any(npn.mean) & warn)
      warning("Studies with negative or zero mean get no weight ",
              "in meta-analysis.")
    ##
    TE <- ifelse(npn.n | npn.mean, NA, log(mean))
    ##
    seTE <- ifelse(npn.n | npn.mean, NA, sqrt(sd^2 / (n * mean^2)))
    ##
    seTE[is.na(TE)] <- NA
    ##
    transf.null.effect <- log(null.effect)
  }
  ##
  ## Studies with non-positive variance get zero weight in meta-analysis
  ##
  sel <- sd <= 0
  ##
  if (any(sel, na.rm = TRUE) & warn)
    warning("Studies with non-positive standard deviation get ",
            "no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  
  
  ##
  ##
  ## (10) Additional checks for three-level model
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
  ## (11) Do meta-analysis
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
               warn = warn,
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
  ## (12) Generate R object
  ##
  ##
  
  res <- list(n = n, mean = mean, sd = sd,
              method.ci = method.ci,
              method.mean = method.mean, method.sd = method.sd)
  ##
  if (export.median)
    res$median <- median
  if (avail.q1)
    res$q1 <- q1
  if (avail.q3)
    res$q3 <- q3
  if (avail.min)
    res$min <- min
  if (avail.max)
    res$max <- max
  ##
  res$approx.sd <- approx.sd
  res$approx.mean <- approx.mean
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
  if (method.ci == "t") {
    res$lower <- ci.study$lower
    res$upper <- ci.study$upper
    res$statistic <- ci.study$statistic
    res$pval <- ci.study$p
    res$df <- ci.study$df
  }
  else if (!is.null(res$df) && all(is.na(res$df)))
    res$df <- NULL
  ##
  if (all(res$approx.mean == "")) {
    res$approx.mean <- NULL
    res$data$.approx.mean <- NULL
  }
  if (all(res$approx.sd == "")) {
    res$approx.sd <- NULL
    res$data$.approx.sd <- NULL
  }
  ##
  ## Remove test statistic and p-value if null effect is missing
  ##
  if (is.na(null.effect)) {
    res$statistic <- NA
    res$pval <- NA
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
  #
  if (is.null(res$approx.mean))
    res$method.mean <- ""
  #
  if (is.null(res$approx.sd))
    res$method.sd <- ""
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")
  
  res
}
