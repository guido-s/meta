#' Meta-analysis of continuous outcome data
#' 
#' @description
#' Calculation of common and random effects estimates for meta-analyses
#' with continuous outcome data; inverse variance weighting is used
#' for pooling.
#' 
#' @param n.e Number of observations in experimental group.
#' @param mean.e Estimated mean in experimental group.
#' @param sd.e Standard deviation in experimental group.
#' @param n.c Number of observations in control group.
#' @param mean.c Estimated mean in control group.
#' @param sd.c Standard deviation in control group.
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
#' @param median.e Median in experimental group (used to estimate the
#'   mean and standard deviation).
#' @param q1.e First quartile in experimental group (used to estimate
#'   the mean and standard deviation).
#' @param q3.e Third quartile in experimental group (used to estimate
#'   the mean and standard deviation).
#' @param min.e Minimum in experimental group (used to estimate the
#'   mean and standard deviation).
#' @param max.e Maximum in experimental group (used to estimate the
#'   mean and standard deviation).
#' @param median.c Median in control group (used to estimate the mean
#'   and standard deviation).
#' @param q1.c First quartile in control group (used to estimate the
#'   mean and standard deviation).
#' @param q3.c Third quartile in control group (used to estimate the
#'   mean and standard deviation).
#' @param min.c Minimum in control group (used to estimate the mean
#'   and standard deviation).
#' @param max.c Maximum in control group (used to estimate the mean
#'   and standard deviation).
#' @param method.mean A character string indicating which method to
#'   use to approximate the mean from the median and other statistics
#'   (see Details).
#' @param method.sd A character string indicating which method to use
#'   to approximate the standard deviation from sample size, median,
#'   interquartile range and range (see Details).
#' @param approx.mean.e Approximation method to estimate means in
#'   experimental group (see Details).
#' @param approx.mean.c Approximation method to estimate means in
#'   control group (see Details).
#' @param approx.sd.e Approximation method to estimate standard
#'   deviations in experimental group (see Details).
#' @param approx.sd.c Approximation method to estimate standard
#'   deviations in control group (see Details).
#' @param method.ci A character string indicating which method is used
#'   to calculate confidence intervals for individual studies (see
#'   Details).
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
#' @param method.bias A character string indicating which test is to
#'   be used. Either \code{"Begg"}, \code{"Egger"}, \code{"Thompson"},
#'   or \code{"Pustejovsky"} (see \code{\link{metabias}}), can be
#'   abbreviated.
#' @param backtransf A logical indicating whether results for ratio of
#'   means (\code{sm="ROM"}) should be back transformed in printouts
#'   and plots. If TRUE (default), results will be presented as ratio
#'   of means; otherwise log ratio of means will be shown.
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
#' @param sm A character string indicating which summary measure
#'   (\code{"MD"}, \code{"SMD"}, or \code{"ROM"}) is to be used for
#'   pooling of studies.
#' @param pooledvar A logical indicating if a pooled variance should
#'   be used for the mean difference (\code{sm="MD"}).
#' @param method.smd A character string indicating which method is
#'   used to estimate the standardised mean difference
#'   (\code{sm="SMD"}). Either \code{"Hedges"} for Hedges' g
#'   (default), \code{"Cohen"} for Cohen's d, or \code{"Glass"} for
#'   Glass' delta, can be abbreviated.
#' @param sd.glass A character string indicating which standard
#'   deviation is used in the denominator for Glass' method to
#'   estimate the standardised mean difference. Either
#'   \code{"control"} using the standard deviation in the control
#'   group (\code{sd.c}) or \code{"experimental"} using the standard
#'   deviation in the experimental group (\code{sd.e}), can be
#'   abbreviated.
#' @param exact.smd A logical indicating whether exact formulae should
#'   be used in estimation of the standardised mean difference and its
#'   standard error (see Details).
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
#' Calculation of common and random effects estimates for meta-analyses
#' with continuous outcome data; inverse variance weighting is used
#' for pooling.
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
#' Three different types of summary measures are available for continuous
#' outcomes:
#' \itemize{
#' \item mean difference (argument \code{sm = "MD"})
#' \item standardised mean difference (\code{sm = "SMD"})
#' \item ratio of means (\code{sm = "ROM"})
#' }
#' 
#' \subsection{Standardised mean difference}{
#' 
#' For the standardised mean difference three methods are implemented:
#' \itemize{
#' \item Hedges' g (default, \code{method.smd = "Hedges"}) - see
#'   Hedges (1981)
#' \item Cohen's d (\code{method.smd = "Cohen"}) - see Cohen (1988)
#' \item Glass' delta (\code{method.smd = "Glass"}) - see Glass (1976)
#' }
#'
#' Hedges (1981) calculated the exact bias in Cohen's d which is a
#' ratio of gamma distributions with the degrees of freedom,
#' i.e. total sample size minus two, as argument. By default (argument
#' \code{exact.smd = FALSE}), an accurate approximation of this bias
#' provided in Hedges (1981) is utilised for Hedges' g as well as its
#' standard error; these approximations are also used in RevMan
#' 5. Following Borenstein et al. (2009) these approximations are not
#' used in the estimation of Cohen's d. White and Thomas (2005) argued
#' that approximations are unnecessary with modern software and
#' accordingly promote to use the exact formulae; this is possible
#' using argument \code{exact.smd = TRUE}. For Hedges' g the exact
#' formulae are used to calculate the standardised mean difference as
#' well as the standard error; for Cohen's d the exact formula is only
#' used to calculate the standard error. In typical applications (with
#' sample sizes above 10), the differences between using the exact
#' formulae and the approximation will be minimal.
#' 
#' For Glass' delta, by default (argument \code{sd.glass =
#' "control"}), the standard deviation in the control group
#' (\code{sd.c}) is used in the denominator of the standard mean
#' difference. The standard deviation in the experimental group
#' (\code{sd.e}) can be used by specifying \code{sd.glass =
#' "experimental"}.
#' }
#' 
#' \subsection{Ratio of means}{
#' 
#' Meta-analysis of ratio of means -- also called response ratios --
#' is described in Hedges et al. (1999) and Friedrich et al. (2008).
#' Calculations are conducted on the log scale and list elements
#' \code{TE}, \code{TE.common}, and \code{TE.random} contain the
#' logarithm of the ratio of means. In printouts and plots these
#' values are back transformed if argument \code{backtransf = TRUE}.
#' }
#' 
#' \subsection{Approximate means from sample sizes, medians and other statistics}{
#' 
#' Missing means in the experimental group (analogously for the
#' control group) can be derived from
#' \enumerate{
#' \item sample size, median, interquartile range and range (arguments
#'   \code{n.e}, \code{median.e}, \code{q1.e}, \code{q3.e},
#'   \code{min.e}, and \code{max.e}),
#' \item sample size, median and interquartile range (arguments
#'   \code{n.e}, \code{median.e}, \code{q1.e}, and \code{q3.e}), or
#' \item sample size, median and range (arguments \code{n.e},
#'   \code{median.e}, \code{min.e}, and \code{max.e}).
#' }
#' 
#' By default, methods described in Luo et al. (2018) are utilized
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
#' ranges (if available) and finally ranges. Arguments
#' \code{approx.mean.e} and \code{approx.mean.c} can be used to
#' overwrite this behaviour for each individual study and treatment
#' arm:
#' \itemize{
#' \item use means directly (entry \code{""} in argument
#'   \code{approx.mean.e} or \code{approx.mean.c});
#' \item median, interquartile range and range (\code{"iqr.range"});
#' \item median and interquartile range (\code{"iqr"});
#' \item median and range (\code{"range"}).
#' }
#' }
#'
#' \subsection{Approximate standard deviations from sample sizes, medians and other statistics}{
#' 
#' Missing standard deviations in the experimental group (analogously
#' for the control group) can be derived from
#' \enumerate{
#' \item sample size, median, interquartile range and range (arguments
#'   \code{n.e}, \code{median.e}, \code{q1.e}, \code{q3.e},
#'   \code{min.e}, and \code{max.e}),
#' \item sample size, median and interquartile range (arguments
#'   \code{n.e}, \code{median.e}, \code{q1.e} and \code{q3.e}), or
#' \item sample size, median and range (arguments \code{n.e},
#'   \code{median.e}, \code{min.e} and \code{max.e}).
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
#' before interquartile ranges before ranges. Arguments
#' \code{approx.sd.e} and \code{approx.sd.c} can be used to overwrite
#' this default for each individual study and treatment arms:
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
#' For the mean difference (argument \code{sm = "MD"}), the confidence
#' interval for individual studies can be based on the
#' \itemize{
#' \item standard normal distribution (\code{method.ci = "z"}, default), or
#' \item t-distribution (\code{method.ci = "t"}).
#' }
#' 
#' Note, this choice does not affect the results of the common effect
#' and random effects meta-analysis.
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
#' calculated regardless of values choosen for arguments \code{common}
#' and \code{random}. Accordingly, the estimate for the random effects
#' model can be extracted from component \code{TE.random} of an object
#' of class \code{"meta"} even if argument \code{random =
#' FALSE}. However, all functions in R package \bold{meta} will
#' adequately consider the values for \code{common} and
#' \code{random}. E.g. function \code{\link{print.meta}} will not
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
#' An object of class \code{c("metacont", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{metabin}}, \code{\link{metagen}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009):
#' \emph{Introduction to Meta-Analysis}.
#' Chichester: Wiley
#'
#' Cai S, Zhou J, Pan J (2021):
#' Estimating the sample mean and standard deviation from order
#' statistics and sample size in meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{30}, 2701--2719
#' 
#' Cohen J (1988):
#' \emph{Statistical Power Analysis for the Behavioral Sciences
#'   (second ed.)}.
#' Lawrence Erlbaum Associates
#' 
#' Friedrich JO, Adhikari NK, Beyene J (2008):
#' The ratio of means method as an alternative to mean differences for
#' analyzing continuous outcome variables in meta-analysis: A
#' simulation study.
#' \emph{BMC Medical Research Methodology},
#' \bold{8}, 32
#' 
#' Glass G (1976):
#' Primary, secondary, and meta-analysis of research.
#' \emph{Educational Researcher},
#' \bold{5}, 3--8
#' 
#' Hedges LV (1981):
#' Distribution theory for Glass's estimator of effect size and
#' related estimators.
#' \emph{Journal of Educational and Behavioral Statistics},
#' \bold{6}, 107--28
#' 
#' Hedges LV, Gurevitch J, Curtis PS (1999):
#' The meta-analysis of response ratios in experimental ecology.
#' \emph{Ecology},
#' \bold{80}, 1150--6
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
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.4.
#' The Cochrane Collaboration, 2020
#' 
#' Shi J, Luo D, Weng H, Zeng XT, Lin L, Chu H, Tong T (2020):
#' Optimally estimating the sample standard deviation from the
#' five-number summary.
#' \emph{Research Synthesis Methods},
#' \bold{11}, 641--54
#'
#' Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca J (2013):
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
#' White IR, Thomas J (2005):
#' Standardized mean differences in individually-randomized and
#' cluster-randomized trials, with applications to meta-analysis.
#' \emph{Clinical Trials},
#' \bold{2}, 141--51
#' 
#' @examples
#' data(Fleiss1993cont)
#'
#' # Meta-analysis with Hedges' g as effect measure
#' #
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
#' m1
#' forest(m1)
#' 
#' # Use Cohen's d instead of Hedges' g as effect measure
#' #
#' update(m1, method.smd = "Cohen")
#' 
#' # Use Glass' delta instead of Hedges' g as effect measure
#' #
#' update(m1, method.smd = "Glass")
#' 
#' # Use Glass' delta based on the standard deviation in the experimental group
#' #
#' update(m1, method.smd = "Glass", sd.glass = "experimental")
#' 
#' # Calculate Hedges' g based on exact formulae
#' #
#' update(m1, exact.smd = TRUE)
#' 
#' data(amlodipine)
#' m2 <- metacont(n.amlo, mean.amlo, sqrt(var.amlo), n.plac, mean.plac, sqrt(var.plac),
#'   data = amlodipine, studlab = study)
#' m2
#' 
#' # Use pooled variance
#' #
#' update(m2, pooledvar = TRUE)
#' 
#' # Meta-analysis of response ratios (Hedges et al., 1999)
#' #
#' data(woodyplants)
#' m3 <- metacont(n.elev, mean.elev, sd.elev, n.amb, mean.amb, sd.amb,
#'   data = woodyplants, sm = "ROM")
#' m3
#' print(m3, backtransf = FALSE)
#' 
#' @export metacont


metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     cluster = NULL,
                     ##
                     median.e, q1.e, q3.e, min.e, max.e,
                     median.c, q1.c, q3.c, min.c, max.c,
                     method.mean = "Luo", method.sd = "Shi",
                     approx.mean.e, approx.mean.c = approx.mean.e,
                     approx.sd.e, approx.sd.c = approx.sd.e,
                     ##
                     sm = gs("smcont"),
                     method.ci = gs("method.ci.cont"),
                     level = gs("level"),
                     ##
                     pooledvar = gs("pooledvar"),
                     method.smd = gs("method.smd"),
                     sd.glass = gs("sd.glass"),
                     exact.smd = gs("exact.smd"),
                     ##
                     common = gs("common"),
                     random = gs("random") | !is.null(tau.preset),
                     overall = common | random,
                     overall.hetstat = common | random,
                     prediction = gs("prediction") | !missing(method.predict),
                     ##
                     method.tau = gs("method.tau"),
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
                     label.e = gs("label.e"), label.c = gs("label.c"),
                     label.left = gs("label.left"),
                     label.right = gs("label.right"),
                     ##
                     subgroup, subgroup.name = NULL,
                     print.subgroup.name = gs("print.subgroup.name"),
                     sep.subgroup = gs("sep.subgroup"),
                     test.subgroup = gs("test.subgroup"),
                     prediction.subgroup = gs("prediction.subgroup"),
                     ##
                     byvar, id, adhoc.hakn,
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
  chknull(sm)
  sm <- setchar(sm, gs("sm4cont"))
  chklevel(level)
  ##
  missing.method.tau <- missing(method.tau)
  method.tau <- setchar(method.tau, gs("meth4tau"))
  ##
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
  adhoc.hakn.pi <- setchar(adhoc.hakn.pi, gs("adhoc4hakn.pi"))
  ## Classic tests + Pustejovsky
  method.bias <-
    setmethodbias(method.bias, c(1:3, if (sm == "SMD") 8))
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
  ## Additional arguments / checks for metacont objects
  ##
  fun <- "metacont"
  if (sm != "MD")
    method.ci <- "z"
  method.ci <- setchar(method.ci, gs("ci4cont"))
  ##
  method.mean <-
    setchar(method.mean, c("Luo", "Wan", "Cai", "QE-McGrath", "BC-McGrath"))
  method.sd <-
    setchar(method.sd, c("Shi", "Wan", "Cai", "QE-McGrath", "BC-McGrath"))
  ##
  if (method.mean %in% c("Cai", "QE-McGrath", "BC-McGrath"))
    is.installed.package("estmeansd", argument = "method.mean",
                         value = method.mean)
  if (method.sd %in% c("Cai", "QE-McGrath", "BC-McGrath"))
    is.installed.package("estmeansd", argument = "method.sd",
                         value = method.sd)
  ##
  chklogical(pooledvar)
  method.smd <- setchar(method.smd, c("Hedges", "Cohen", "Glass"))
  sd.glass <- setchar(sd.glass, c("control", "experimental"))
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
  ## Catch 'n.e', 'mean.e', 'sd.e', 'n.c', 'mean.c', 'sd.c' from data:
  ##
  missing.mean.e <- missing(mean.e)
  missing.sd.e <- missing(sd.e)
  missing.mean.c <- missing(mean.c)
  missing.sd.c <- missing(sd.c)
  ##
  missing.median.e <- missing(median.e)
  missing.q1.e <- missing(q1.e)
  missing.q3.e <- missing(q3.e)
  missing.min.e <- missing(min.e)
  missing.max.e <- missing(max.e)
  ##
  missing.median.c <- missing(median.c)
  missing.q1.c <- missing(q1.c)
  missing.q3.c <- missing(q3.c)
  missing.min.c <- missing(min.c)
  missing.max.c <- missing(max.c)
  ##
  if (missing.mean.e & missing.median.e)
    stop("Provide either argument 'mean.e' or 'median.e'.",
         call. = FALSE)
  ##
  if (missing.mean.c & missing.median.c)
    stop("Provide either argument 'mean.c' or 'median.c'.",
         call. = FALSE)
  ##
  if (missing.sd.e &
      !((!missing.q1.e & !missing.q3.e) |
        (!missing.min.e & !missing.max.e)))
    stop("Provide either argument 'sd.e' and ",
         "arguments 'q1.e' & 'q3.e' or 'min.e & 'max.e'.",
         call. = FALSE)
  ##
  if (missing.sd.c &
      !((!missing.q1.c & !missing.q3.c) |
        (!missing.min.c & !missing.max.c)))
    stop("Provide either argument 'sd.c' and ",
         "arguments 'q1.c' & 'q3.c' or 'min.c & 'max.c'.",
         call. = FALSE)
  ##
  n.e <- catch("n.e", mc, data, sfsp)
  chknull(n.e)
  k.All <- length(n.e)
  ##
  mean.e <- catch("mean.e", mc, data, sfsp)
  if (!missing.mean.e)
    chknull(mean.e)
  else
    mean.e <- rep(NA, k.All)
  ##
  sd.e <- catch("sd.e", mc, data, sfsp)
  if (!missing.sd.e)
    chknull(sd.e)
  else
    sd.e <- rep(NA, k.All)
  ##
  n.c <- catch("n.c", mc, data, sfsp)
  chknull(n.c)
  ##
  mean.c <- catch("mean.c", mc, data, sfsp)
  if (!missing.mean.c)
    chknull(mean.c)
  else
    mean.c <- rep(NA, k.All)
  ##
  sd.c <- catch("sd.c", mc, data, sfsp)
  if (!missing.sd.c)
    chknull(sd.c)
  else
    sd.c <- rep(NA, k.All)
  ##
  ## Catch 'studlab', 'subgroup', 'subset', 'exclude', and
  ## 'cluster' from data:
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
  missing.cluster <- missing(cluster)
  cluster <- catch("cluster", mc, data, sfsp)
  missing.id <- missing(id)
  id <- catch("id", mc, data, sfsp)
  ##
  cluster <- deprecated2(cluster, missing.cluster, id, missing.id,
                         warn.deprecated)
  with.cluster <- !is.null(cluster)
  ##
  if (is.null(method.tau.ci))
    if (with.cluster)
      method.tau.ci <- "PL"
    else if (method.tau == "DL")
      method.tau.ci <- "J"
     else
      method.tau.ci <- "QP"
  method.tau.ci <- setchar(method.tau.ci, gs("meth4tau.ci"))
  ##
  ## Catch 'median.e', 'q1.e', 'q3.e', 'min.e', 'max.e', 'median.c',
  ## 'q1.c', 'q3.c', 'min.c', 'max.c', 'approx.mean.e', 'approx.sd.e',
  ## 'approx.mean.c', and 'approx.sd.c', from data:
  ##
  median.e <- catch("median.e", mc, data, sfsp)
  ##
  q1.e <- catch("q1.e", mc, data, sfsp)
  ##
  q3.e <- catch("q3.e", mc, data, sfsp)
  ##
  min.e <- catch("min.e", mc, data, sfsp)
  ##
  max.e <- catch("max.e", mc, data, sfsp)
  ##
  median.c <- catch("median.c", mc, data, sfsp)
  ##
  q1.c <- catch("q1.c", mc, data, sfsp)
  ##
  q3.c <- catch("q3.c", mc, data, sfsp)
  ##
  min.c <- catch("min.c", mc, data, sfsp)
  ##
  max.c <- catch("max.c", mc, data, sfsp)
  ##
  missing.approx.mean.e <- missing(approx.mean.e)
  approx.mean.e <- catch("approx.mean.e", mc, data, sfsp)
  if (!missing.approx.mean.e)
    missing.approx.mean.c <- FALSE
  else
    missing.approx.mean.c <- missing(approx.mean.c)
  approx.mean.c <- catch("approx.mean.c", mc, data, sfsp)
  ##
  missing.approx.sd.e <- missing(approx.sd.e)
  approx.sd.e <- catch("approx.sd.e", mc, data, sfsp)
  ##
  if (!missing.approx.sd.e)
    missing.approx.sd.c <- FALSE
  else
    missing.approx.sd.c <- missing(approx.sd.c)
  approx.sd.c <- catch("approx.sd.c", mc, data, sfsp)
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
  ## (3) Check length of essential variables
  ##
  ##
  arg <- "n.e"
  chklength(mean.e, k.All, arg)
  chklength(sd.e, k.All, arg)
  chklength(n.c, k.All, arg)
  chklength(mean.c, k.All, arg)
  chklength(sd.c, k.All, arg)
  chklength(studlab, k.All, arg)
  if (with.cluster)
    chklength(cluster, k.All, arg)
  ##
  if (!missing.median.e)
    chklength(median.e, k.All, arg)
  if (!missing.q1.e)
    chklength(q1.e, k.All, arg)
  if (!missing.q3.e)
    chklength(q3.e, k.All, arg)
  if (!missing.min.e)
    chklength(min.e, k.All, arg)
  if (!missing.max.e)
    chklength(max.e, k.All, arg)
  if (!missing.median.c)
    chklength(median.c, k.All, arg)
  if (!missing.q1.c)
    chklength(q1.c, k.All, arg)
  if (!missing.q3.c)
    chklength(q3.c, k.All, arg)
  if (!missing.min.c)
    chklength(min.c, k.All, arg)
  if (!missing.max.c)
    chklength(max.c, k.All, arg)
  ##
  if (!missing.approx.mean.e) {
    if (length(approx.mean.e) == 1)
      rep_len(approx.mean.e, k.All)
    else
      chklength(approx.mean.e, k.All, arg)
    ##
    approx.mean.e <- setchar(approx.mean.e, c("", "iqr.range", "iqr", "range"))
  }
  if (!missing.approx.mean.c) {
    if (length(approx.mean.c) == 1)
      rep_len(approx.mean.c, k.All)
    else
      chklength(approx.mean.c, k.All, arg)
    ##
    approx.mean.c <- setchar(approx.mean.c, c("", "iqr.range", "iqr", "range"))
  }
  ##
  if (!missing.approx.sd.e) {
    if (length(approx.sd.e) == 1)
      rep_len(approx.sd.e, k.All)
    else
      chklength(approx.sd.e, k.All, arg)
    ##
    approx.sd.e <- setchar(approx.sd.e, c("", "iqr.range", "iqr", "range"))
  }
  if (!missing.approx.sd.c) {
    if (length(approx.sd.c) == 1)
      rep_len(approx.sd.c, k.All)
    else
      chklength(approx.sd.c, k.All, arg)
    ##
    approx.sd.c <- setchar(approx.sd.c, c("", "iqr.range", "iqr", "range"))
  }
  ##
  if (by) {
    chklength(subgroup, k.All, arg)
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
      data <- data.frame(.n.e = n.e)
    else
      data$.n.e <- n.e
    ##
    data$.mean.e <- mean.e
    data$.sd.e <- sd.e
    data$.n.c <- n.c
    data$.mean.c <- mean.c
    data$.sd.c <- sd.c
    data$.studlab <- studlab
    ##
    if (!missing.median.e)
      data$.median.e <- median.e
    if (!missing.q1.e)
      data$.q1.e <- q1.e
    if (!missing.q3.e)
      data$.q3.e <- q3.e
    if (!missing.min.e)
      data$.min.e <- min.e
    if (!missing.max.e)
      data$.max.e <- max.e
    if (!missing.median.c)
      data$.median.c <- median.c
    if (!missing.q1.c)
      data$.q1.c <- q1.c
    if (!missing.q3.c)
      data$.q3.c <- q3.c
    if (!missing.min.c)
      data$.min.c <- min.c
    if (!missing.max.c)
      data$.max.c <- max.c
    if (!missing.approx.mean.e)
      data$.approx.mean.e <- approx.mean.e
    if (!missing.approx.mean.c)
      data$.approx.mean.c <- approx.mean.c
    if (!missing.approx.sd.e)
      data$.approx.sd.e <- approx.sd.e
    if (!missing.approx.sd.c)
      data$.approx.sd.c <- approx.sd.c
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
    n.e <- n.e[subset]
    mean.e <- mean.e[subset]
    sd.e <- sd.e[subset]
    n.c <- n.c[subset]
    mean.c <- mean.c[subset]
    sd.c <- sd.c[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (!missing.median.e)
      median.e <- median.e[subset]
    if (!missing.q1.e)
      q1.e <- q1.e[subset]
    if (!missing.q3.e)
      q3.e <- q3.e[subset]
    if (!missing.min.e)
      min.e <- min.e[subset]
    if (!missing.max.e)
      max.e <- max.e[subset]
    if (!missing.median.c)
      median.c <- median.c[subset]
    if (!missing.q1.c)
      q1.c <- q1.c[subset]
    if (!missing.q3.c)
      q3.c <- q3.c[subset]
    if (!missing.min.c)
      min.c <- min.c[subset]
    if (!missing.max.c)
      max.c <- max.c[subset]
    if (!missing.approx.mean.e)
      approx.mean.e <- approx.mean.e[subset]
    if (!missing.approx.mean.c)
      approx.mean.c <- approx.mean.c[subset]
    if (!missing.approx.sd.e)
      approx.sd.e <- approx.sd.e[subset]
    if (!missing.approx.sd.c)
      approx.sd.c <- approx.sd.c[subset]
    ##
    if (by)
      subgroup <- subgroup[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(n.e)
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
  chknumeric(n.e)
  chknumeric(mean.e)
  chknumeric(sd.e)
  chknumeric(n.c)
  chknumeric(mean.c)
  chknumeric(sd.c)
  ##
  if (!missing.median.e)
    chknumeric(median.e)
  if (!missing.q1.e)
    chknumeric(q1.e)
  if (!missing.q3.e)
    chknumeric(q3.e)
  if (!missing.min.e)
    chknumeric(min.e)
  if (!missing.max.e)
    chknumeric(max.e)
  if (!missing.median.c)
    chknumeric(median.c)
  if (!missing.q1.c)
    chknumeric(q1.c)
  if (!missing.q3.c)
    chknumeric(q3.c)
  if (!missing.min.c)
    chknumeric(min.c)
  if (!missing.max.c)
    chknumeric(max.c)
  ##
  ## Recode integer as numeric:
  ##
  n.e    <- int2num(n.e)
  mean.e <- int2num(mean.e)
  sd.e   <- int2num(sd.e)
  n.c    <- int2num(n.c)
  mean.c <- int2num(mean.c)
  sd.c   <- int2num(sd.c)
  ##
  if (!missing.median.e)
    median.e <- int2num(median.e)
  if (!missing.q1.e)
    q1.e <- int2num(q1.e)
  if (!missing.q3.e)
    q3.e <- int2num(q3.e)
  if (!missing.min.e)
    min.e <- int2num(min.e)
  if (!missing.max.e)
    max.e <- int2num(max.e)
  if (!missing.median.c)
    median.c <- int2num(median.c)
  if (!missing.q1.c)
    q1.c <- int2num(q1.c)
  if (!missing.q3.c)
    q3.c <- int2num(q3.c)
  if (!missing.min.c)
    min.c <- int2num(min.c)
  if (!missing.max.c)
    max.c <- int2num(max.c)
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
  if (missing.approx.mean.e) {
    approx.mean.e <- rep_len("", length(n.e))
    ##
    ## (a) Use IQR and range
    ##
    sel.NA.e <- is.na(mean.e)
    if (any(sel.NA.e) & !missing.median.e &
        !missing.q1.e & !missing.q3.e &
        !missing.min.e & !missing.max.e) {
      j <- sel.NA.e & !is.na(median.e) & !is.na(q1.e) & !is.na(q3.e) &
        !is.na(min.e) & !is.na(max.e)
      approx.mean.e[j] <- "iqr.range"
      ##
      mean.e[j] <- mean.sd.iqr.range(n.e[j], median.e[j], q1.e[j], q3.e[j],
                                     min.e[j], max.e[j], method.mean)$mean
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA.e <- is.na(mean.e)
    if (any(sel.NA.e) & !missing.median.e & !missing.q1.e & !missing.q3.e) {
      j <- sel.NA.e & !is.na(median.e) & !is.na(q1.e) & !is.na(q3.e)
      approx.mean.e[j] <- "iqr"
      mean.e[j] <- mean.sd.iqr(n.e[j], median.e[j], q1.e[j], q3.e[j],
                               method.mean)$mean
    }
    ##
    ## (c) Use range
    ##
    sel.NA.e <- is.na(mean.e)
    if (any(sel.NA.e) & !missing.median.e & !missing.min.e & !missing.max.e) {
      j <- sel.NA.e & !is.na(median.e) & !is.na(min.e) & !is.na(max.e)
      approx.mean.e[j] <- "range"
      mean.e[j] <- mean.sd.range(n.e[j], median.e[j], min.e[j], max.e[j],
                                 method.mean)$mean
    }
  }
  else {
    j <- 0
    for (i in approx.mean.e) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        mean.e[j] <- mean.sd.iqr.range(n.e[j], median.e[j], q1.e[j], q3.e[j],
                                     min.e[j], max.e[j], method.mean)$mean
      else if (i == "iqr")
        mean.e[j] <- mean.sd.iqr(n.e[j], median.e[j], q1.e[j], q3.e[j],
                                 method.mean)$mean
      else if (i == "range")
        mean.e[j] <- mean.sd.range(n.e[j], median.e[j], min.e[j], max.e[j],
                                   method.mean)$mean
    }
  }
  ##
  if (missing.approx.mean.c) {
    approx.mean.c <- rep_len("", length(n.c))
    ##
    ## (a) Use IQR and range
    ##
    sel.NA.c <- is.na(mean.c)
    if (any(sel.NA.c) & !missing.median.c &
        !missing.q1.c & !missing.q3.c &
        !missing.min.c & !missing.max.c) {
      j <- sel.NA.c & !is.na(median.c) & !is.na(q1.c) & !is.na(q3.c) &
        !is.na(min.c) & !is.na(max.c)
      approx.mean.c[j] <- "iqr.range"
      ##
      mean.c[j] <- mean.sd.iqr.range(n.c[j], median.c[j], q1.c[j], q3.c[j],
                                     min.c[j], max.c[j], method.mean)$mean
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA.c <- is.na(mean.c)
    if (any(sel.NA.c) & !missing.median.c & !missing.q1.c & !missing.q3.c) {
      j <- sel.NA.c & !is.na(median.c) & !is.na(q1.c) & !is.na(q3.c)
      approx.mean.c[j] <- "iqr"
      mean.c[j] <- mean.sd.iqr(n.c[j], median.c[j], q1.c[j], q3.c[j],
                               method.mean)$mean
    }
    ##
    ## (c) Use range
    ##
    sel.NA.c <- is.na(mean.c)
    if (any(sel.NA.c) & !missing.median.c & !missing.min.c & !missing.max.c) {
      j <- sel.NA.c & !is.na(median.c) & !is.na(min.c) & !is.na(max.c)
      approx.mean.c[j] <- "range"
      mean.c[j] <- mean.sd.range(n.c[j], median.c[j], min.c[j], max.c[j],
                                 method.mean)$mean
    }
  }
  else {
    j <- 0
    for (i in approx.mean.c) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        mean.c[j] <- mean.sd.iqr.range(n.c[j], median.c[j], q1.c[j], q3.c[j],
                                       min.c[j], max.c[j], method.mean)$mean
      else if (i == "iqr")
        mean.c[j] <- mean.sd.iqr(n.c[j], median.c[j], q1.c[j], q3.c[j],
                                 method.mean)$mean
      else if (i == "range")
        mean.c[j] <- mean.sd.range(n.c[j], median.c[j], min.c[j], max.c[j],
                                   method.mean)$mean
    }
  }
  
  
  ##
  ##
  ## (8) Calculate standard deviation from other information
  ##
  ##
  if (missing.median.e) {
    median.e.sd <- mean.e
    missing.median.e <- FALSE
    export.median.e <- FALSE
  }
  else {
    median.e.sd <- median.e
    median.e.sd[is.na(median.e.sd)] <- mean.e[is.na(median.e.sd)]
    export.median.e <- TRUE
  }
  ##
  if (missing.approx.sd.e) {
    approx.sd.e <- rep_len("", length(n.e))
    ##
    ## (a) Use IQR and range
    ##
    sel.NA.e <- is.na(sd.e)
    if (any(sel.NA.e) & !missing.median.e &
        !missing.q1.e & !missing.q3.e &
        !missing.min.e & !missing.max.e) {
      j <- sel.NA.e & !is.na(median.e.sd) & !is.na(q1.e) & !is.na(q3.e) &
        !is.na(min.e) & !is.na(max.e)
      approx.sd.e[j] <- "iqr.range"
      ##
      sd.e[j] <- mean.sd.iqr.range(n.e[j], median.e.sd[j], q1.e[j], q3.e[j],
                                   min.e[j], max.e[j],
                                   method.sd = method.sd)$sd
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA.e <- is.na(sd.e)
    if (any(sel.NA.e) & !missing.median.e & !missing.q1.e & !missing.q3.e) {
      j <- sel.NA.e & !is.na(median.e.sd) & !is.na(q1.e) & !is.na(q3.e)
      approx.sd.e[j] <- "iqr"
      sd.e[j] <- mean.sd.iqr(n.e[j], median.e.sd[j], q1.e[j], q3.e[j])$sd
    }
    ##
    ## (c) Use range
    ##
    sel.NA.e <- is.na(sd.e)
    if (any(sel.NA.e) & !missing.median.e & !missing.min.e & !missing.max.e) {
      j <- sel.NA.e & !is.na(median.e.sd) & !is.na(min.e) & !is.na(max.e)
      approx.sd.e[j] <- "range"
      sd.e[j] <- mean.sd.range(n.e[j], median.e.sd[j], min.e[j], max.e[j])$sd
    }
  }
  else {
    j <- 0
    for (i in approx.sd.e) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        sd.e[j] <- mean.sd.iqr.range(n.e[j], median.e.sd[j], q1.e[j], q3.e[j],
                                     min.e[j], max.e[j],
                                     method.sd = method.sd)$sd
      else if (i == "iqr")
        sd.e[j] <- mean.sd.iqr(n.e[j], median.e.sd[j], q1.e[j], q3.e[j])$sd
      else if (i == "range")
        sd.e[j] <- mean.sd.range(n.e[j], median.e.sd[j], min.e[j], max.e[j])$sd
    }
  }
  ##
  if (missing.median.c) {
    median.c.sd <- mean.c
    missing.median.c <- FALSE
    export.median.c <- FALSE
  }
  else {
    median.c.sd <- median.c
    median.c.sd[is.na(median.c.sd)] <- mean.c[is.na(median.c.sd)]
    export.median.c <- TRUE
  }
  ##
  if (missing.approx.sd.c) {
    approx.sd.c <- rep_len("", length(n.c))
    ##
    ## (a) Use IQR and range
    ##
    sel.NA.c <- is.na(sd.c)
    if (any(sel.NA.c) & !missing.median.c &
        !missing.q1.c & !missing.q3.c &
        !missing.min.c & !missing.max.c) {
      j <- sel.NA.c & !is.na(median.c.sd) & !is.na(q1.c) & !is.na(q3.c) &
        !is.na(min.c) & !is.na(max.c)
      approx.sd.c[j] <- "iqr.range"
      ##
      sd.c[j] <- mean.sd.iqr.range(n.c[j], median.c.sd[j], q1.c[j], q3.c[j],
                                   min.c[j], max.c[j],
                                   method.sd = method.sd)$sd
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA.c <- is.na(sd.c)
    if (any(sel.NA.c) & !missing.median.c & !missing.q1.c & !missing.q3.c) {
      j <- sel.NA.c & !is.na(median.c.sd) & !is.na(q1.c) & !is.na(q3.c)
      approx.sd.c[j] <- "iqr"
      sd.c[j] <- mean.sd.iqr(n.c[j], median.c.sd[j], q1.c[j], q3.c[j])$sd
    }
    ##
    ## (c) Use range
    ##
    sel.NA.c <- is.na(sd.c)
    if (any(sel.NA.c) & !missing.median.c & !missing.min.c & !missing.max.c) {
      j <- sel.NA.c & !is.na(median.c.sd) & !is.na(min.c) & !is.na(max.c)
      approx.sd.c[j] <- "range"
      sd.c[j] <- mean.sd.range(n.c[j], median.c.sd[j], min.c[j], max.c[j])$sd
    }
  }
  else {
    j <- 0
    for (i in approx.sd.c) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        sd.c[j] <- mean.sd.iqr.range(n.c[j], median.c[j], q1.c[j], q3.c[j],
                                     min.c[j], max.c[j],
                                     method.sd = method.sd)$sd
      else if (i == "iqr")
        sd.c[j] <- mean.sd.iqr(n.c[j], median.c.sd[j], q1.c[j], q3.c[j])$sd
      else if (i == "range")
        sd.c[j] <- mean.sd.range(n.c[j], median.c.sd[j], min.c[j], max.c[j])$sd
    }
  }
  ##
  if (keepdata) {
    if (!isCol(data, ".subset")) {
      data$.sd.e <- sd.e
      data$.mean.e <- mean.e
      data$.sd.c <- sd.c
      data$.mean.c <- mean.c
      if (!missing.approx.sd.e)
        data$.approx.sd.e <- approx.sd.e
      if (!missing.approx.sd.c)
        data$.approx.sd.c <- approx.sd.c
      if (!missing.approx.mean.e)
        data$.approx.mean.e <- approx.mean.e
      if (!missing.approx.mean.c)
        data$.approx.mean.c <- approx.mean.c
    }
    else {
      data$.sd.e[data$.subset] <- sd.e
      data$.mean.e[data$.subset] <- mean.e
      data$.sd.c[data$.subset] <- sd.c
      data$.mean.c[data$.subset] <- mean.c
      if (!missing.approx.sd.e)
        data$.approx.sd.e[data$.subset] <- approx.sd.e
      if (!missing.approx.sd.c)
        data$.approx.sd.c[data$.subset] <- approx.sd.c
      if (!missing.approx.mean.e)
        data$.approx.mean.e[data$.subset] <- approx.mean.e
      if (!missing.approx.mean.c)
        data$.approx.mean.c[data$.subset] <- approx.mean.c
    }
  }
  
  
  ##
  ##
  ## (9) Calculate results for individual studies
  ##
  ##
  npn.n <- npn(n.e) | npn(n.c)
  ##
  N <- n.e + n.c
  if (sm == "MD" | sm == "ROM")
    var.pooled <- ((n.e - 1) * sd.e^2 + (n.c - 1) * sd.c^2) / (N - 2)
  ##
  if (any(npn.n) & warn)
    warning("Note, studies with non-positive values for n.e and / or n.c get no weight in meta-analysis.")
  ##
  if (sm == "MD") {
    TE <- ifelse(npn.n, NA, mean.e - mean.c)
    ##
    if (pooledvar)
      seTE <-
        ifelse(npn.n, NA, sqrt(var.pooled * (1 / n.e + 1 / n.c)))
    else
      seTE <-
        ifelse(npn.n, NA, sqrt(sd.e^2 / n.e + sd.c^2 / n.c))
    ##
    seTE[is.na(TE)] <- NA
    ##
    if (method.ci == "t")
      ci.study <- ci(TE, seTE, level = level, df = n.e + n.c - 2)
  }
  else if (sm == "SMD") {
    J <- function(x)
      exp(lgamma(x / 2) - log(sqrt(x / 2)) - lgamma((x - 1) / 2))
    ##
    K <- function(x) 1 - (x - 2) / (x * J(x)^2)
    ##
    if (method.smd %in% c("Hedges", "Cohen"))
      S.within <- sqrt(((n.e - 1) * sd.e^2 + (n.c - 1) * sd.c^2) / (N - 2))
    else
      S.within <- if (sd.glass == "control") sd.c else sd.e
    ##
    smd <- ifelse(npn.n, NA, (mean.e - mean.c) / S.within)
    ##
    if (method.smd == "Cohen") {
      ##
      ## Borenstein et al. (2009), p. 26-27;
      ## White and Thomas (2005), p. 143
      ##
      TE <- smd
      if (exact.smd) {
        J <- function(x)
          exp(lgamma(x / 2) - log(sqrt(x / 2)) - lgamma((x - 1) / 2))
        ##
        K <- function(x) 1 - (x - 2) / (x * J(x)^2)
        ##
        seTE <-
          ifelse(npn.n, NA,
                 sqrt(1 / n.e + 1 / n.c + (J(N - 2) * smd)^2 * K(N - 2)))
      }
      else
        seTE <-
          ifelse(npn.n, NA, sqrt(1 / n.e + 1 / n.c + TE^2 / (2 * N)))
    }
    else if (method.smd == "Hedges") {
      ##
      ## Hedges and Olkin (1985); White and Thomas (2005), p. 143;
      ## formulae used in RevMan 5 (exact.smd = FALSE)
      ##
      if (exact.smd) {
        J <- function(x)
          exp(lgamma(x / 2) - log(sqrt(x / 2)) - lgamma((x - 1) / 2))
        ##
        K <- function(x) 1 - (x - 2) / (x * J(x)^2)
      }
      else {
        J <- function(x) 1 - 3 / (4 * x - 1)
        K <- function(x) 1 / (2 * (x - 1.94))
      }
      ##
      TE   <- J(N - 2) * smd
      seTE <-
        ifelse(npn.n, NA, sqrt(1 / n.e + 1 / n.c + TE^2 * K(N - 2)))
    }
    else if (method.smd == "Glass") {
      ##
      ## see Cooper & Hedges (1994), p. 238
      ##
      n.g  <- if (sd.glass == "control") n.c else n.e
      ##
      TE <- smd
      seTE <-
        ifelse(npn.n, NA, sqrt(1 / n.e + 1 / n.c + TE^2 / (2 * (n.g - 1))))
    }
    ##
    seTE[is.na(TE)] <- NA
  }
  ##
  else if (sm == "ROM") {
    npn.mean <- npn(mean.e) | npn(mean.c)
    ##
    if (any(npn.mean) & warn)
      warning("Note, studies with negative or zero means get no weight in meta-analysis.")
    
    TE <- ifelse(npn.n | npn.mean, NA, log(mean.e / mean.c))
    ##
    if (pooledvar)
      seTE <- ifelse(npn.n, NA,
                     sqrt(var.pooled * (1 / (n.e * mean.e^2) + 1 / (n.c * mean.c^2))))
    else
      seTE <- ifelse(npn.n | npn.mean, NA,
                     sqrt(sd.e^2 / (n.e * mean.e^2) + sd.c^2 / (n.c * mean.c^2)))
    ##
    seTE[is.na(TE)] <- NA
  }
  ##
  ## Studies with non-positive variance get zero weight in meta-analysis
  ##
  sel <- sd.e <= 0 | sd.c <= 0
  ##
  if (any(sel, na.rm = TRUE) & warn)
    warning("Note, studies with non-positive values for sd.e or sd.c get no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  ##
  if (sm == "SMD")
    TE[sel] <- NA
  
  
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
    common <- FALSE
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
  ## (11) Do meta-analysis
  ##
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
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               ##
               text.common = text.common, text.random = text.random,
               text.predict = text.predict,
               text.w.common = text.w.common, text.w.random = text.w.random,
               ##
               title = title, complab = complab, outclab = outclab,
               label.e = label.e, label.c = label.c,
               label.left = label.left, label.right = label.right,
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (by & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   TE.tau, level.ma, subgroup, control)
  }
  
  
  ##
  ##
  ## (11) Generate R object
  ##
  ##
  res <- list(n.e = n.e, mean.e = mean.e, sd.e = sd.e,
              n.c = n.c, mean.c = mean.c, sd.c = sd.c,
              pooledvar = pooledvar,
              method.smd = method.smd, sd.glass = sd.glass,
              exact.smd = exact.smd,
              method.ci = method.ci,
              method.mean = method.mean, method.sd = method.sd,
              n.e.pooled = sum(n.e, na.rm = TRUE),
              n.c.pooled = sum(n.c, na.rm = TRUE))
  ##
  if (export.median.e)
    res$median.e <- median.e
  if (!missing.q1.e)
    res$q1.e <- q1.e
  if (!missing.q3.e)
    res$q3.e <- q3.e
  if (!missing.min.e)
    res$min.e <- min.e
  if (!missing.max.e)
    res$max.e <- max.e
  ##
  if (export.median.c)
    res$median.c <- median.c
  if (!missing.q1.c)
    res$q1.c <- q1.c
  if (!missing.q3.c)
    res$q3.c <- q3.c
  if (!missing.min.c)
    res$min.c <- min.c
  if (!missing.max.c)
    res$max.c <- max.c
  ##
  res$approx.sd.e <- approx.sd.e
  res$approx.sd.c <- approx.sd.c
  res$approx.mean.e <- approx.mean.e
  res$approx.mean.c <- approx.mean.c
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
  res <- c(res, m)
  ##
  ## Add data
  ##
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
  if (all(res$approx.mean.e == "")) {
    res$approx.mean.e <- NULL
    res$data$.approx.mean.e <- NULL
  }
  if (all(res$approx.sd.e == "")) {
    res$approx.sd.e <- NULL
    res$data$.approx.sd.e <- NULL
  }
  if (all(res$approx.mean.c == "")) {
    res$approx.mean.c <- NULL
    res$data$.approx.mean.c <- NULL
  }
  if (all(res$approx.sd.c == "")) {
    res$approx.sd.c <- NULL
    res$data$.approx.sd.c <- NULL
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
      if (res$three.level)
        res <- c(res, subgroup(res, NULL,
                               factor(res$subgroup, bylevs(res$subgroup))))
      else
        res <- c(res, subgroup(res, hcc$tau.resid))
    }
    ##
    if (tau.common && is.null(tau.preset))
      res <- addHet(res, hcc)
    ##
    res$n.w <- NULL
    res$event.w <- NULL
    ##
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
