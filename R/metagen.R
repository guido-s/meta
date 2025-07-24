#' Generic inverse variance meta-analysis
#' 
#' @description
#' Common effect and random effects meta-analysis based on estimates
#' (e.g. log hazard ratios) and their standard errors. The inverse
#' variance method is used for pooling.
#'
#' Three-level random effects meta-analysis (Van den Noortgate et al., 2013) is
#' available by internally calling \code{\link[metafor]{rma.mv}} function from
#' R package \bold{metafor} (Viechtbauer, 2010).
#' 
#' @param TE Estimate of treatment effect, e.g., log hazard ratio or
#'   risk difference or an R object created with \code{\link{pairwise}}.
#' @param seTE Standard error of treatment estimate or standard deviation of
#'   n-of-1 trials.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used (see Details).
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots (see Details).
#' @param cluster An optional vector specifying which estimates come
#'   from the same cluster resulting in the use of a three-level
#'   meta-analysis model.
#' @param rho Assumed correlation of estimates within a cluster.
#' @param cycles A numeric vector with the number of cycles per patient / study
#'   in n-of-1 trials.
#' @param weights A single numeric or vector with user-specified weights.
#' @param weights.common User-specified weights (common effect model).
#' @param weights.random User-specified weights (random effects model).
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param method.ci A character string indicating which method is used
#'   to calculate confidence intervals for individual studies, see
#'   Details.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
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
#' @param n.e Number of observations in experimental group (or total
#'   sample size in study).
#' @param n.c Number of observations in control group.
#' @param pval P-value (used to estimate the standard error).
#' @param df Degrees of freedom (used in test or to construct
#'   confidence intervals).
#' @param lower Lower limit of confidence interval (used to estimate
#'   the standard error).
#' @param upper Upper limit of confidence interval (used to estimate
#'   the standard error).
#' @param level.ci Level of confidence interval.
#' @param median Median (used to estimate the treatment effect and
#'   standard error).
#' @param q1 First quartile (used to estimate the treatment effect and
#'   standard error).
#' @param q3 Third quartile (used to estimate the treatment effect and
#'   standard error).
#' @param min Minimum (used to estimate the treatment effect and
#'   standard error).
#' @param max Maximum (used to estimate the treatment effect and
#'   standard error).
#' @param method.mean A character string indicating which method to
#'   use to approximate the mean from the median and other statistics
#'   (see Details).
#' @param method.sd A character string indicating which method to use
#'   to approximate the standard deviation from sample size, median,
#'   interquartile range and range (see Details).
#' @param approx.TE Approximation method to estimate treatment
#'   estimate (see Details).
#' @param approx.seTE Approximation method to estimate standard error
#'   (see Details).
#' @param transf A logical indicating whether inputs for arguments
#'   \code{TE}, \code{lower} and \code{upper} are already
#'   appropriately transformed to conduct the meta-analysis or on the
#'   original scale. If \code{transf = TRUE} (default), inputs are
#'   expected to be log odds ratios instead of odds ratios for
#'   \code{sm = "OR"} and Fisher's z transformed correlations instead
#'   of correlations for \code{sm = "ZCOR"}, for example.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE} (default), results for \code{sm = "OR"} are printed as odds
#'   ratios rather than log odds ratios and results for \code{sm =
#'   "ZCOR"} are printed as correlations rather than Fisher's z
#'   transformed correlations, for example.
#' @param func.transf A function used to transform inputs for
#'   arguments \code{TE}, \code{lower} and \code{upper}.
#' @param func.backtransf A function used to back-transform results.
#' @param args.transf An optional list to provide additional arguments
#'   to \code{func.transf}.
#' @param args.backtransf An optional list to provide additional
#'   arguments to \code{func.backtransf}.
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
#' @param adhoc.hakn Deprecated argument (replaced by
#'   'adhoc.hakn.ci').
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param keeprma A logical indicating whether \code{\link[metafor]{rma.mv}}
#'   object from three-level meta-analysis should be stored.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.mv}}.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' This function provides the \emph{generic inverse variance method}
#' for meta-analysis which requires treatment estimates and their
#' standard errors (Borenstein et al., 2010). The method is useful,
#' e.g., for pooling of survival data (using log hazard ratio and
#' standard errors as input). Arguments \code{TE} and \code{seTE} can
#' be used to provide treatment estimates and standard errors
#' directly. However, it is possible to derive these quantities from
#' other information.
#' 
#' Argument \code{cycles} can be used to conduct a meta-analysis of n-of-1
#' trials according to Senn (2024). In this case, argument \code{seTE} does
#' not contain the standard error but standard deviation for individual
#' trials / patients. Trial-specific standard errors are calculated from an
#' average standard deviation multiplied by the number of cycles minus 1, i.e.,
#' the degrees of freedom. Details of the meta-analysis method are provided in
#' Senn (2024). Note, arguments used in the approximation of means or
#' standard errors, like \code{lower} and \code{upper}, or \code{df}, are
#' ignored for the meta-analysis of n-of-1 trials.
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
#' \subsection{Approximate treatment estimates}{
#' 
#' Missing treatment estimates can be derived from
#' \enumerate{
#' \item confidence limits provided by arguments \code{lower} and
#'   \code{upper};
#' \item median, interquartile range and range (arguments
#'   \code{median}, \code{q1}, \code{q3}, \code{min}, and \code{max});
#' \item median and interquartile range (arguments \code{median},
#'   \code{q1} and \code{q3});
#' \item median and range (arguments \code{median}, \code{min} and
#'   \code{max}).
#' }
#' For confidence limits, the treatment estimate is defined as the
#' center of the confidence interval (on the log scale for relative
#' effect measures like the odds ratio or hazard ratio).
#'
#' If the treatment effect is a mean it can be approximated from
#' sample size, median, interquartile range and range.
#'
#' By default, methods described in Luo et al. (2018) are utilised
#' (argument \code{method.mean = "Luo"}):
#' \itemize{
#' \item equation (7) if sample size, median and range are available,
#' \item equation (11) if sample size, median and interquartile range
#'   are available,
#' \item equation (15) if sample size, median, range and interquartile
#'   range are available.
#' }
#' 
#' Instead the methods described in Wan et al. (2014) are used if
#' argument \code{method.mean = "Wan"}:
#' \itemize{
#' \item equation (2) if sample size, median and range are available,
#' \item equation (14) if sample size, median and interquartile range
#'   are available,
#' \item equation (10) if sample size, median, range and interquartile
#'   range are available.
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
#' By default, missing treatment estimates are replaced successively
#' using these method, i.e., confidence limits are utilised before
#' interquartile ranges. Argument \code{approx.TE} can be used to
#' overwrite this default for each individual study:
#' \itemize{
#' \item Use treatment estimate directly (entry \code{""} in argument
#'   \code{approx.TE});
#' \item confidence limits (\code{"ci"} in argument \code{approx.TE});
#' \item median, interquartile range and range (\code{"iqr.range"});
#' \item median and interquartile range (\code{"iqr"});
#' \item median and range (\code{"range"}).
#' }
#' }
#'
#' \subsection{Approximate standard errors}{
#' 
#' Missing standard errors can be derived from
#' \enumerate{
#' \item p-value provided by arguments \code{pval} and (optional)
#'   \code{df};
#' \item confidence limits (arguments \code{lower}, \code{upper}, and
#'   (optional) \code{df});
#' \item sample size, median, interquartile range and range (arguments
#'   \code{n.e} and / or \code{n.c}, \code{median}, \code{q1},
#'   \code{q3}, \code{min}, and \code{max});
#' \item sample size, median and interquartile range (arguments
#'   \code{n.e} and / or \code{n.c}, \code{median}, \code{q1} and
#'   \code{q3});
#' \item sample size, median and range (arguments \code{n.e} and / or
#'   \code{n.c}, \code{median}, \code{min} and \code{max}).
#' }
#' For p-values and confidence limits, calculations are either based
#' on the standard normal or \emph{t}-distribution if argument
#' \code{df} is provided. Furthermore, argument \code{level.ci} can be
#' used to provide the level of the confidence interval.
#'
#' Wan et al. (2014) describe methods to estimate the standard
#' deviation (and thus the standard error by deviding the standard
#' deviation with the square root of the sample size) from the sample
#' size, median and additional statistics. Shi et al. (2020) provide
#' an improved estimate of the standard deviation if the interquartile
#' range and range are available in addition to the sample size and
#' median. Accordingly, equation (11) in Shi et al. (2020) is the
#' default (argument \code{method.sd = "Shi"}), if the median,
#' interquartile range and range are provided (arguments
#' \code{median}, \code{q1}, \code{q3}, \code{min} and
#' \code{max}). The method by Wan et al. (2014) is used if argument
#' \code{method.sd = "Wan"} and, depending on the sample size, either
#' equation (12) or (13) is used. If only the interquartile range or
#' range is available, equations (15) / (16) and (7) / (9) in Wan et
#' al. (2014) are used, respectively. The sample size of individual
#' studies must be provided with arguments \code{n.e} and / or
#' \code{n.c}. The total sample size is calculated as \code{n.e} +
#' \code{n.c} if both arguments are provided.
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
#' By default, missing standard errors are replaced successively using
#' these method, e.g., p-value before confidence limits before
#' interquartile range and range. Argument \code{approx.seTE} can be
#' used to overwrite this default for each individual study:
#' 
#' \itemize{
#' \item Use standard error directly (entry \code{""} in argument
#'   \code{approx.seTE});
#' \item p-value (\code{"pval"} in argument \code{approx.seTE});
#' \item confidence limits (\code{"ci"});
#' \item median, interquartile range and range (\code{"iqr.range"});
#' \item median and interquartile range (\code{"iqr"});
#' \item median and range (\code{"range"}).
#' }
#' }
#'
#' \subsection{Confidence intervals for individual studies}{
#' 
#' For the mean difference (argument \code{sm = "MD"}), the confidence
#' interval for individual studies can be based on the
#' \itemize{
#' \item standard normal distribution (\code{method.ci = "z"}), or
#' \item \emph{t}-distribution (\code{method.ci = "t"}).
#' }
#'
#' By default, the first method is used if argument \code{df} is
#' missing and the second method otherwise.
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
#' \subsection{Specify the null hypothesis of test for an overall effect}{
#'
#' Argument \code{null.effect} can be used to specify the (treatment)
#' effect under the null hypothesis in a test for an overall
#' effect.
#'
#' By default (\code{null.effect = 0}), the null hypothesis
#' corresponds to "no difference" (which is obvious for absolute
#' effect measures like the mean difference (\code{sm = "MD"}) or
#' standardised mean difference (\code{sm = "SMD"})). For relative
#' effect measures, e.g., risk ratio (\code{sm = "RR"}) or odds ratio
#' (\code{sm = "OR"}), the null effect is defined on the log scale,
#' i.e., \emph{log}(RR) = 0 or \emph{log}(OR) = 0 which is equivalent
#' to testing RR = 1 or OR = 1.
#'
#' Use of argument \code{null.effect} is especially useful for summary
#' measures without a "natural" null effect, i.e., in situations
#' without a second (treatment) group. For example, an overall
#' proportion of 50\% could be tested in the meta-analysis of single
#' proportions with argument \code{null.effect = 0.5}.
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
#' \code{exclude} (see Examples).
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
#' \code{random}. For example, functions \code{\link{print.meta}} and
#' \code{\link{forest.meta}} will not show results for the random
#' effects model if \code{random = FALSE}.
#'
#' A prediction interval will only be shown if \code{prediction =
#' TRUE}.
#' 
#' Argument \code{pscale} can be used to rescale single proportions or
#' risk differences, e.g. \code{pscale = 1000} means that proportions
#' are expressed as events per 1000 observations. This is useful in
#' situations with (very) low event probabilities.
#' 
#' Argument \code{irscale} can be used to rescale single rates or rate
#' differences, e.g. \code{irscale = 1000} means that rates are
#' expressed as events per 1000 time units, e.g. person-years. This is
#' useful in situations with (very) low rates. Argument \code{irunit}
#' can be used to specify the time unit used in individual studies
#' (default: "person-years"). This information is printed in summaries
#' and forest plots if argument \code{irscale} is not equal to 1.
#'
#' Default settings for \code{common}, \code{random},
#' \code{pscale}, \code{irscale}, \code{irunit} and several other
#' arguments can be set for the whole R session using
#' \code{\link{settings.meta}}.
#' }
#'
#' @note
#' R function \code{\link[metafor]{rma.uni}} from R package
#' \pkg{metafor} (Viechtbauer 2010) is called internally to estimate
#' the between-study variance \eqn{\tau^2}.
#' 
#' @return
#' An object of class \code{c("metagen", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{pairwise}},
#'   \code{\link{print.meta}}, \code{\link{settings.meta}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
#' A basic introduction to fixed-effect and random-effects models for
#' meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{1}, 97--111
#'
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
#' Senn S (2024):
#' The analysis of continuous data from n-of-1 trials using paired cycles:
#' a simple tutorial.
#' \emph{Trials},
#' \bold{25}.
#' 
#' Shi J, Luo D, Weng H, Zeng X-T, Lin L, Chu H, et al. (2020):
#' Optimally estimating the sample standard deviation from the
#' five-number summary.
#' \emph{Research Synthesis Methods}.
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
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
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, study,
#'   data = Fleiss1993bin, sm = "RR", method = "I")
#' m1
#' 
#' # Identical results using the generic inverse variance method with
#' # log risk ratio and its standard error:
#' # Note, argument 'n.e' in metagen() is used to provide the total
#' # sample size which is calculated from the group sample sizes n.e
#' # and n.c in meta-analysis m1.
#' m1.gen <- metagen(TE, seTE, studlab, n.e = n.e + n.c, data = m1, sm = "RR")
#' m1.gen
#' forest(m1.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))
#' 
#' # Meta-analysis with prespecified between-study variance
#' #
#' metagen(m1$TE, m1$seTE, sm = "RR", tau.preset = sqrt(0.1))
#' 
#' # Meta-analysis of survival data:
#' #
#' logHR <- log(c(0.95, 1.5))
#' selogHR <- c(0.25, 0.35)
#' metagen(logHR, selogHR, sm = "HR")
#' 
#' # Paule-Mandel method to estimate between-study variance for data
#' # from Paule & Mandel (1982)
#' #
#' average <- c(27.044, 26.022, 26.340, 26.787, 26.796)
#' variance <- c(0.003, 0.076, 0.464, 0.003, 0.014)
#' #
#' metagen(average, sqrt(variance), sm = "MD", method.tau = "PM")
#' 
#' # Conduct meta-analysis using hazard ratios and 95% confidence intervals
#' #
#' # Data from Steurer et al. (2006), Analysis 1.1 Overall survival
#' # https://doi.org/10.1002/14651858.CD004270.pub2
#' #
#' study <- c("FCG on CLL 1996", "Leporrier 2001", "Rai 2000", "Robak 2000")
#' HR <- c(0.55, 0.92, 0.79, 1.18)
#' lower.HR <- c(0.28, 0.79, 0.59, 0.64)
#' upper.HR <- c(1.09, 1.08, 1.05, 2.17)
#' #
#' # Hazard ratios and confidence intervals as input
#' #
#' summary(metagen(HR, lower = lower.HR, upper = upper.HR,
#'   studlab = study, sm = "HR", transf = FALSE))
#' #
#' # Same result with log hazard ratios as input
#' #
#' summary(metagen(log(HR), lower = log(lower.HR), upper = log(upper.HR),
#'   studlab = study, sm = "HR"))
#' #
#' # Again, same result using an unknown summary measure and
#' # arguments 'func.transf' and 'func.backtransf'
#' #
#' summary(metagen(HR, lower = lower.HR, upper = upper.HR,
#'   studlab = study, sm = "Hazard ratio",
#'   func.transf = log, func.backtransf = exp))
#' #
#' # Finally, same result only providing argument 'func.transf' as the
#' # back-transformation for the logarithm is known
#' #
#' summary(metagen(HR, lower = lower.HR, upper = upper.HR,
#'   studlab = study, sm = "Hazard ratio",
#'   func.transf = log))
#'
#' # Exclude MRC-1 and MRC-2 studies from meta-analysis, however,
#' # show them in printouts and forest plots
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'   data = Fleiss1993bin, sm = "RR", method = "I",
#'   exclude = study %in% c("MRC-1", "MRC-2"))
#' #
#' # Exclude MRC-1 and MRC-2 studies completely from meta-analysis
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'   data = Fleiss1993bin, sm = "RR", method = "I",
#'   subset = !(study %in% c("MRC-1", "MRC-2")))
#'
#' # Exclude studies with total sample size above 1500
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'   data = Fleiss1993bin, sm = "RR", method = "I",
#'   exclude = (n.asp + n.plac) > 1500)
#'
#' # Exclude studies containing "MRC" in study name
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'   data = Fleiss1993bin, sm = "RR", method = "I",
#'   exclude = grep("MRC", study))
#'
#' # Use both arguments 'subset' and 'exclude'
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'   data = Fleiss1993bin, sm = "RR", method = "I",
#'   subset = (n.asp + n.plac) > 1500,
#'   exclude = grep("MRC", study))
#'
#' \dontrun{
#' # Three-level model: effects of modified school calendars on
#' # student achievement
#' data(dat.konstantopoulos2011, package = "metadat")
#' metagen(yi, sqrt(vi), studlab = study, data = dat.konstantopoulos2011,
#'   sm = "SMD",
#'   cluster = district, detail.tau = c("district", "district/school"))
#' }
#' 
#' @export metagen


metagen <- function(TE, seTE, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    cluster = NULL, rho = 0,
                    #
                    cycles = NULL,
                    #
                    weights = NULL,
                    weights.common = weights, weights.random = weights,
                    #
                    sm = "",
                    ##
                    method.ci = if (missing(df)) "z" else "t",
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
                    n.e = NULL, n.c = NULL,
                    ##
                    pval, df, lower, upper, level.ci = 0.95,
                    median, q1, q3, min, max,
                    method.mean = "Luo",
                    method.sd = "Shi",
                    ##
                    approx.TE, approx.seTE,
                    ##
                    transf = gs("transf") & missing(func.transf),
                    backtransf = gs("backtransf") | !missing(func.backtransf),
                    func.transf,
                    func.backtransf,
                    args.transf,
                    args.backtransf,
                    pscale = 1,
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
                    byvar, id, adhoc.hakn,
                    ##
                    keepdata = gs("keepdata"),
                    keeprma = gs("keeprma"),
                    #
                    warn = gs("warn"),
                    warn.deprecated = gs("warn.deprecated"),
                    ##
                    control = NULL,
                    ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  missing.studlab <- missing(studlab)
  missing.sm <- missing(sm)
  missing.subgroup <- missing(subgroup)
  missing.overall <- missing(overall)
  missing.overall.hetstat <- missing(overall.hetstat)
  missing.test.subgroup <- missing(test.subgroup)
  missing.label.e <- missing(label.e)
  missing.label.c <- missing(label.c)
  missing.complab <- missing(complab)
  #
  missing.adhoc.hakn.pi <- missing(adhoc.hakn.pi)
  #
  missing.method.tau <- missing(method.tau)
  missing.tau.common <- missing(tau.common)
  avail.detail.tau <- !missing(detail.tau) & !is.null(detail.tau)
  missing.method.predict <- missing(method.predict)
  #
  missing.func.transf <- missing(func.transf)
  missing.args.transf <- missing(args.transf)
  missing.func.backtransf <- missing(func.backtransf)
  missing.args.backtransf <- missing(args.backtransf)
  #
  sm <- replaceNULL(sm, "")
  sm <- setchar(sm,
                unique(c(gs("sm4bin"), gs("sm4cont"), gs("sm4cor"),
                         gs("sm4inc"), gs("sm4mean"),
                         gs("sm4prop"), gs("sm4rate"), "")),
                stop.at.error = FALSE, return.NULL = FALSE,
                nchar.equal = TRUE)
  ##
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
  chklevel(level)
  #
  if (length(method.common.ci) != 1)
    stop("Argument 'method.common.ci' must be of length 1.",
         call. = FALSE)
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  method.tau <- setchar(method.tau, gs("meth4tau"))
  #
  tau.common <- replaceNULL(tau.common, FALSE)
  chklogical(tau.common)
  #
  method.I2 <- setchar(method.I2, gs("meth4i2"))
  #
  chklogical(prediction)
  chklevel(level.predict)
  ##
  chknumeric(rho, min = -1, max = 1)
  #
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
  if (!is.null(seed.predict))
    chknumeric(seed.predict, length = 1)
  ##
  chknumeric(null.effect, length = 1)
  #
  chklogical(transf)
  chklogical(backtransf)
  ##
  avail.func.transf <- !missing.func.transf && !is.null(func.transf)
  avail.args.transf <- !missing.args.transf && !is.null(args.transf)
  avail.func.backtransf <-
    !missing.func.backtransf && !is.null(func.backtransf)
  avail.args.backtransf <-
    !missing.args.backtransf && !is.null(args.backtransf)
  ##
  if (avail.func.transf) {
    chkfunc(func.transf)
    if (is.function(func.transf))
      func.transf <- deparse(substitute(func.transf))
  }
  else
    func.transf <- NULL
  ##
  if (avail.args.transf)
    chklist(args.transf)
  else
    args.transf <- NULL
  ##
  if (avail.func.backtransf) {
    chkfunc(func.backtransf)    
    if (is.function(func.backtransf))
      func.backtransf <- deparse(substitute(func.backtransf))
  }
  else
    func.backtransf <- NULL
  ##
  if (avail.args.backtransf)
    chklist(args.backtransf)    
  else
    args.backtransf <- NULL
  ##
  if (is.null(func.transf) & !is.null(args.transf)) {
    warning("Argument 'args.transf' ignored as argument ",
            "'func.transf' is ",
            if (!avail.func.transf) "missing." else "NULL.",
            call. = FALSE)
    args.transf <- NULL
  }
  ##
  if (is.null(func.backtransf) & !is.null(args.backtransf)) {
    warning("Argument 'args.backtransf' ignored as argument ",
            "'func.backtransf' is ",
            if (!avail.func.backtransf) "missing." else "NULL.",
            call. = FALSE)
    args.backtransf <- NULL
  }
  ##
  if (!is.null(func.transf) & is.null(func.backtransf)) {
    if (func.transf == "log" & is.null(args.transf))
      func.backtransf <- "exp"
    else if (func.transf == "cor2z" & is.null(args.transf))
      func.backtransf <- "z2cor"
    else if (func.transf == "p2logit" & is.null(args.transf))
      func.backtransf <- "logit2p"
    else if (func.transf == "p2asin" & is.null(args.transf))
      func.backtransf <- "asin2p"
    else if (func.transf == "VE2logVR" & is.null(args.transf))
      func.backtransf <- "logVR2VE"
    else
      stop("Argument 'func.backtransf' must be specified.",
           call. = FALSE)
  }
  ##
  missing.pscale <- missing(pscale)
  chknumeric(pscale, length = 1)
  if (!backtransf & pscale != 1) {
    if (!missing.pscale)
      warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.",
              call. = FALSE)
    pscale <- 1
  }
  missing.irscale <- missing(irscale)
  chknumeric(irscale, length = 1)
  if (!backtransf & irscale != 1) {
    if (!missing.irscale)
      warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.",
              call. = FALSE)
    irscale <- 1
  }
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
  #
  chkchar(label.e, length = 1)
  chkchar(label.c, length = 1)
  #
  chklogical(keepdata)
  chklogical(keeprma)
  ##
  ## Additional arguments / checks
  ##
  fun <- "metagen"
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
  missing.adhoc.hakn.ci <- missing(adhoc.hakn.ci)
  adhoc.hakn.ci <-
    deprecated2(adhoc.hakn.ci, missing.adhoc.hakn.ci,
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
  missing.byvar <- missing(byvar)
  byvar <- catch("byvar", mc, data, sfsp)
  ##
  ## Catch 'TE', 'seTE', 'median', 'lower', 'upper', 'n.e', 'n.c', and
  ## 'cluster' from data:
  ##
  missing.TE <- missing(TE)
  missing.seTE <- missing(seTE)
  missing.median <- missing(median)
  missing.lower <- missing(lower)
  missing.upper <- missing(upper)
  ##
  if (missing.TE & missing.median & (missing.lower | missing.upper))
    stop("Treatment estimates missing. ",
         "Provide either argument 'TE' or 'median', ",
         "or arguments 'lower' and 'upper'.",
         call. = FALSE)
  #
  TE <- catch("TE", mc, data, sfsp)
  avail.TE <- !(missing.TE || is.null(TE))
  #
  if (inherits(TE, "pairwise")) {
    is.pairwise <- TRUE
    #
    txt.ignore <- "as first argument is a pairwise object"
    #
    warn_ignore_input(sm, !missing.sm, txt.ignore)
    warn_ignore_input(seTE, !missing.seTE, txt.ignore)
    warn_ignore_input(studlab, !missing.studlab, txt.ignore)
    #
    sm <- attr(TE, "sm")
    reference.group <- attr(TE, "reference.group")
    #
    if (is.null(attr(TE, "varnames")))
      seTE <- TE$seTE
    else
      seTE <- TE[[attr(TE, "varnames")[2]]]
    #
    missing.sm <- FALSE
    missing.seTE <- FALSE
    missing.studlab <- FALSE
    #
    studlab <- TE$studlab
    #
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    #
    if (!is.null(TE$n1))
      n.e <- TE$n1
    else
      n.e <- NULL
    #
    if (!is.null(TE$n2))
      n.c <- TE$n2
    else
      n.c <- NULL
    #
    pairdata <- TE
    data <- TE
    nulldata <- FALSE
    #
    if (is.null(attr(TE, "varnames")))
      TE <- TE$TE
    else
      TE <- TE[[attr(TE, "varnames")[1]]]
    #
    avail.TE <- !is.null(TE)
    #
    wo <- treat1 == reference.group
    #
    if (any(wo)) {
      TE[wo] <- -TE[wo]
      #
      ttreat1 <- treat1
      treat1[wo] <- treat2[wo]
      treat2[wo] <- ttreat1[wo]
      #
      if (!is.null(n.e) & !is.null(n.c)) {
        tn.e <- n.e
        n.e[wo] <- n.c[wo]
        n.c[wo] <- tn.e[wo]
      }
    }
    #
    if (missing.subgroup) {
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
    else
      subgroup <- catch("subgroup", mc, data, sfsp)
  }
  else {
    is.pairwise <- FALSE
    #
    if (missing.sm)
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
    else
      sm <- ""
    ##
    seTE <- catch("seTE", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    subgroup <- catch("subgroup", mc, data, sfsp)
    subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                            warn.deprecated)
    #
    if (!missing(n.e))
      n.e <- catch("n.e", mc, data, sfsp)
    if (!missing(n.c))
      n.c <- catch("n.c", mc, data, sfsp)
  }
  #
  method.bias <- setmethodbias(method.bias)
  #
  by <- !is.null(subgroup)
  #
  median <- catch("median", mc, data, sfsp)
  lower <- catch("lower", mc, data, sfsp)
  upper <- catch("upper", mc, data, sfsp)
  ##
  avail.median <- !(missing.median || is.null(median))
  avail.lower <- !(missing.lower || is.null(lower))
  avail.upper <- !(missing.upper || is.null(upper))
  ##
  if (!avail.TE & !avail.median & (!avail.lower | !avail.upper))
    stop("Treatment estimates missing. ",
         "Provide either argument 'TE' or 'median', ",
         "or arguments 'lower' and 'upper'.",
         call. = FALSE)
  ##
  TE.orig <- NULL
  lower.orig <- NULL
  upper.orig <- NULL
  ##
  if (!transf) {
    if (avail.TE) {
      TE.orig <- TE
      TE <- transf(TE, sm, func.transf, args.transf)
    }
    if (avail.lower) {
      lower.orig <- lower
      lower <- transf(lower, sm, func.transf, args.transf)
    }
    if (avail.upper) {
      upper.orig <- upper
      upper <- transf(upper, sm, func.transf, args.transf)
    }
    if (sm == "VE" && avail.lower & avail.upper) {
      tmp.l <- lower
      lower <- upper
      upper <- tmp.l
      ##
      tmp.l <- lower.orig
      lower.orig <- upper.orig
      upper.orig <- tmp.l
    }   
  }
  ##
  missing.cluster <- missing(cluster)
  cluster <- catch("cluster", mc, data, sfsp)
  missing.id <- missing(id)
  id <- catch("id", mc, data, sfsp)
  ##
  cluster <- deprecated2(cluster, missing.cluster, id, missing.id,
                         warn.deprecated)
  with.cluster <- !is.null(cluster)
  #
  if (with.cluster)
    idx <- seq_along(cluster)
  #
  missing.cycles <- missing(cycles)
  cycles <- catch("cycles", mc, data, sfsp)
  with.cycles <- !is.null(cycles)
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
  if (with.cycles) {
    chknumeric(cycles, min = 1)
    #
    if (method.ci != "z")
      method.ci <- "z"
  }
  #
  if (with.cluster & with.cycles)
    stop("Arguments 'cluster' (multi-level model) and 'cycles' (n-of-1 trials)",
         " cannot be used together.",
         call. = FALSE)
  #
  missing.method.tau.ci <- missing(method.tau.ci)
  ##
  k.All <- if (avail.TE)
             length(TE)
           else if (avail.median)
             length(median)
           else if (avail.lower)
             length(lower)
           else if (avail.upper)
             length(upper)
           else
             NA
  ##
  if (!avail.TE)
    TE <- rep_len(NA, k.All)
  ##
  if (missing.seTE)
    seTE <- rep_len(NA, k.All)
  ##
  ## Catch 'studlab', 'subset', and 'exclude' from data:
  ##
  studlab <- catch("studlab", mc, data, sfsp)
  studlab <- setstudlab(studlab, k.All)
  ##
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  ##
  exclude <- catch("exclude", mc, data, sfsp)
  missing.exclude <- is.null(exclude)
  ##
  ## Catch 'pval', 'df', 'level.ci', 'q1', 'q3', 'min', 'max',
  ## 'approx.TE' and 'approx.seTE', from data:
  ##
  missing.pval <- missing(pval)
  pval <- catch("pval", mc, data, sfsp)
  avail.pval <- !(missing.pval || is.null(pval))
  ##
  missing.df <- missing(df)
  df <- catch("df", mc, data, sfsp)
  avail.df <- !(missing.df || is.null(df))
  ##
  if (!missing(level.ci))
    level.ci <- catch("level.ci", mc, data, sfsp)
  ##
  missing.q1 <- missing(q1)
  q1 <- catch("q1", mc, data, sfsp)
  avail.q1 <- !(missing.q1 || is.null(q1))
  ##
  missing.q3 <- missing(q3)
  q3 <- catch("q3", mc, data, sfsp)
  avail.q3 <- !(missing.q3 || is.null(q3))
  ##
  missing.min <- missing(min)
  min <- catch("min", mc, data, sfsp)
  avail.min <- !(missing.min || is.null(min))
  ##
  missing.max <- missing(max)
  max <- catch("max", mc, data, sfsp)
  avail.max <- !(missing.max || is.null(max))
  ##
  missing.approx.TE <- missing(approx.TE)
  approx.TE <- catch("approx.TE", mc, data, sfsp)
  avail.approx.TE <- !(missing.approx.TE || is.null(approx.TE))
  ##
  missing.approx.seTE <- missing(approx.seTE)
  approx.seTE <- catch("approx.seTE", mc, data, sfsp)
  avail.approx.seTE <- !(missing.approx.seTE || is.null(approx.seTE))
  #
  missing.text.random <- missing(text.random)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  
  arg <- if (avail.TE) "TE" else "median"
  chklength(seTE, k.All, arg)
  chklength(studlab, k.All, arg)
  ##
  if (by) {
    chklength(subgroup, k.All, arg)
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
  #
  if (usw.common) {
    if (length(weights.common) == 1)
      weights.common <- rep(weights.common, k.All)
    else
      chklength(weights.common, k.All, arg)
  }
  #
  if (usw.random) {
    if (length(weights.random) == 1)
      weights.random <- rep(weights.random, k.All)
    else
      chklength(weights.random, k.All, arg)
  }
  #
  if (!is.null(n.e))
    chklength(n.e, k.All, arg)
  if (!is.null(n.c))
    chklength(n.c, k.All, arg)
  if (with.cluster)
    chklength(cluster, k.All, arg)
  ##
  if (avail.approx.TE) {
    if (length(approx.TE) == 1)
      rep_len(approx.TE, k.All)
    else
      chklength(approx.TE, k.All, arg)
    ##
    approx.TE <- setchar(approx.TE, c("", "ci", "iqr.range", "iqr", "range"))
  }
  ##
  if (avail.approx.seTE) {
    if (length(approx.seTE) == 1)
      rep_len(approx.seTE, k.All)
    else
      chklength(approx.seTE, k.All, arg)
    ##
    approx.seTE <- setchar(approx.seTE,
                           c("", "pval", "ci", "iqr.range", "iqr", "range"))
  }
  ##
  if (avail.pval)
    chklength(pval, k.All, arg)
  if (avail.df)
    chklength(df, k.All, arg)
  if (avail.lower)
    chklength(lower, k.All, arg)
  if (avail.upper)
    chklength(upper, k.All, arg)
  if (length(level.ci) == 1)
    level.ci <- rep_len(level.ci, k.All)
  else if (!is.null(level.ci))
    chklength(level.ci, k.All, arg)
  if (avail.median)
    chklength(median, k.All, arg)
  if (avail.q1)
    chklength(q1, k.All, arg)
  if (avail.q3)
    chklength(q3, k.All, arg)
  if (avail.min)
    chklength(min, k.All, arg)
  if (avail.max)
    chklength(max, k.All, arg)
  if (with.cycles)
    chklength(cycles, k.All, arg)
  
  
  ##
  ##
  ## (4) Subset, exclude studies, and subgroups
  ##
  ##
  
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of argument 'subset' is larger than number of studies.")
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
    if (inherits(data, "meta")) {
      data <- data$data
      if (isCol(data, ".subset"))
        data <- data[data$.subset, ]
    }
    else if (nulldata & !is.pairwise)
      data <- data.frame(.studlab = studlab)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    #
    data$.TE <- TE
    data$.seTE <- seTE
    #
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
    if (with.cluster) {
      data$.id <- data$.cluster <- cluster
      data$.idx <- idx
    }
    #
    if (with.cycles)
      data$.cycles <- cycles
    #
    if (usw.common)
      data$.weights.common <- weights.common
    #
    if (usw.random)
      data$.weights.random <- weights.random
    #
    if (avail.pval)
      data$.pval <- pval
    if (avail.df)
      data$.df <- df
    ##
    if (avail.lower)
      data$.lower <- lower
    if (avail.upper)
      data$.upper <- upper
    if (avail.lower | avail.upper)
      data$.level.ci <- level.ci
    ##
    if (avail.median)
      data$.median <- median
    if (avail.q1)
      data$.q1 <- q1
    if (avail.q3)
      data$.q3 <- q3
    ##
    if (avail.min)
      data$.min <- min
    if (avail.max)
      data$.max <- max
    ##
    if (avail.approx.TE)
      data$.approx.TE <- approx.TE
    if (avail.approx.seTE)
      data$.approx.seTE <- approx.seTE
    ##
    if (!is.null(n.e))
      data$.n.e <- n.e
    if (!is.null(n.c))
      data$.n.c <- n.c
    ##
    if (!is.null(TE.orig))
      data$.TE.orig <- TE.orig
    ##
    if (!is.null(lower.orig))
      data$.lower.orig <- lower.orig
    ##
    if (!is.null(upper.orig))
      data$.upper.orig <- upper.orig
  }
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  
  if (!missing.subset) {
    TE <- TE[subset]
    seTE <- seTE[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (with.cluster) {
      cluster <- cluster[subset]
      idx <- idx[subset]
    }
    #
    if (with.cycles)
      cycles <- cycles[subset]
    #
    if (by)
      subgroup <- subgroup[subset]
    ##
    if (!is.null(n.e))
      n.e <- n.e[subset]
    if (!is.null(n.c))
      n.c <- n.c[subset]
    ##
    if (avail.pval)
      pval <- pval[subset]
    if (avail.df)
      df <- df[subset]
    if (avail.lower)
      lower <- lower[subset]
    if (avail.upper)
      upper <- upper[subset]
    level.ci <- level.ci[subset]
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
    if (avail.approx.TE)
      approx.TE <- approx.TE[subset]
    if (avail.approx.seTE)
      approx.seTE <- approx.seTE[subset]
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
  k.all <- length(TE)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1) {
    common <- FALSE
    ##
    random <- FALSE
    method.random.ci <- "classic"
    adhoc.hakn.ci <- ""
    ##
    prediction <- FALSE
    method.predict <- "V"
    adhoc.hakn.pi <- ""
    ##
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(TE)
  chknumeric(seTE, 0)
  
  
  ##
  ##
  ## (7) Calculate standard error from other information
  ##
  ##
  
  if (!with.cycles) {
    if (!avail.approx.seTE) {
      approx.seTE <- rep_len("", length(TE))
      ##
      ## Use confidence limits
      ##
      sel.NA <- is.na(seTE)
      if (any(sel.NA) & avail.lower & avail.upper) {
        j <- sel.NA & !is.na(lower) & !is.na(upper)
        approx.seTE[j] <- "ci"
        if (!avail.df)
          seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$seTE
        else
          seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j], df[j])$seTE
      }
      ##
      ## Use p-values
      ##
      sel.NA <- is.na(seTE)
      if (any(sel.NA) & avail.pval) {
        j <- sel.NA & !is.na(TE) & !is.na(pval)
        approx.seTE[j] <- "pval"
        if (!avail.df)
          seTE[j] <- seTE.pval(TE[j], pval[j])$seTE
        else
          seTE[j] <- seTE.pval(TE[j], pval[j], df[j])$seTE
      }
      ##
      ## Use IQR and range
      ##
      sel.NA <- is.na(seTE)
      if (any(sel.NA) &
          avail.median & avail.q1 & avail.q3 & avail.min & avail.max &
          !(is.null(n.e) & is.null(n.c))) {
        j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3) &
          !is.na(min) & !is.na(max)
        approx.seTE[j] <- "iqr.range"
        if (is.null(n.c))
          seTE[j] <- mean_sd_iqr_range(n.e[j], median[j], q1[j], q3[j],
                                       min[j], max[j],
                                       method.sd = method.sd)$se
        else if (is.null(n.e))
          seTE[j] <- mean_sd_iqr_range(n.c[j], median[j], q1[j], q3[j],
                                       min[j], max[j],
                                       method.sd = method.sd)$se
        else
          seTE[j] <- mean_sd_iqr_range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                       min[j], max[j],
                                       method.sd = method.sd)$se
      }
      ##
      ## Use IQR
      ##
      sel.NA <- is.na(seTE)
      if (any(sel.NA) &
          avail.median & avail.q1 & avail.q3 &
          !(is.null(n.e) & is.null(n.c))) {
        j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3)
        approx.seTE[j] <- "iqr"
        if (is.null(n.c))
          seTE[j] <- mean_sd_iqr(n.e[j], median[j], q1[j], q3[j])$se
        else if (is.null(n.e))
          seTE[j] <- mean_sd_iqr(n.c[j], median[j], q1[j], q3[j])$se
        else
          seTE[j] <- mean_sd_iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$se
      }
      ##
      ## Use range
      ##
      sel.NA <- is.na(seTE)
      if (any(sel.NA) &
          avail.median & avail.min & avail.max &
          !(is.null(n.e) & is.null(n.c))) {
        j <- sel.NA & !is.na(median) & !is.na(min) & !is.na(max)
        approx.seTE[j] <- "range"
        if (is.null(n.c))
          seTE[j] <- mean_sd_range(n.e[j], median[j], min[j], max[j])$se
        else if (is.null(n.e))
          seTE[j] <- mean_sd_range(n.c[j], median[j], min[j], max[j])$se
        else
          seTE[j] <- mean_sd_range(n.e[j] + n.c[j], median[j],
                                   min[j], max[j])$se
      }
    }
    else {
      j <- 0
      for (i in approx.seTE) {
        j <- j + 1
        ##
        if (i == "ci") {
          if (!avail.df)
            seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$seTE
          else
            seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j], df[j])$seTE
        }
        else if (i == "pval") {
          if (!avail.df)
            seTE[j] <- seTE.pval(TE[j], pval[j])$seTE
          else
            seTE[j] <- seTE.pval(TE[j], pval[j], df[j])$seTE
        }
        else if (i == "iqr.range") {
          if (is.null(n.e) & is.null(n.c))
            stop("Sample size needed if argument 'approx.seTE' = \"iqr\".",
                 call. = FALSE)
          else if (is.null(n.c))
            seTE[j] <- mean_sd_iqr_range(n.e[j], median[j], q1[j], q3[j],
                                         min[j], max[j],
                                         method.sd = method.sd)$se
          else if (is.null(n.e))
            seTE[j] <- mean_sd_iqr_range(n.c[j], median[j], q1[j], q3[j],
                                         min[j], max[j],
                                         method.sd = method.sd)$se
          else
            seTE[j] <- mean_sd_iqr_range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                         min[j], max[j],
                                         method.sd = method.sd)$se
        }
        else if (i == "iqr") {
          if (is.null(n.e) & is.null(n.c))
            stop("Sample size needed if argument 'approx.seTE' = \"iqr\".",
                 call. = FALSE)
          else if (is.null(n.c))
            seTE[j] <- mean_sd_iqr(n.e[j], median[j], q1[j], q3[j])$se
          else if (is.null(n.e))
            seTE[j] <- mean_sd_iqr(n.c[j], median[j], q1[j], q3[j])$se
          else
            seTE[j] <- mean_sd_iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$se
        }
        else if (i == "range") {
          if (is.null(n.e) & is.null(n.c))
            stop("Sample size needed if argument 'approx.seTE' = \"range\".",
                 call. = FALSE)
          else if (is.null(n.c))
            seTE[j] <- mean_sd_range(n.e[j], median[j], min[j], max[j])$se
          else if (is.null(n.e))
            seTE[j] <- mean_sd_range(n.c[j], median[j], min[j], max[j])$se
          else
            seTE[j] <- mean_sd_range(n.e[j] + n.c[j], median[j],
                                     min[j], max[j])$se
        }
      }
    }
  }
  
  
  ##
  ##
  ## (8) Calculate treatment estimate from other information
  ##
  ##
  
  if (!with.cycles) {
    if (!avail.approx.TE) {
      approx.TE <- rep_len("", length(TE))
      ##
      ## Use confidence limits
      ##
      sel.NA <- is.na(TE)
      if (any(sel.NA) & avail.lower & avail.upper) {
        j <- sel.NA & !is.na(lower) & !is.na(upper)
        approx.TE[j] <- "ci"
        TE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$TE
      }
      ##
      ## Use IQR and range
      ##
      sel.NA <- is.na(TE)
      if (any(sel.NA) &
          avail.median & avail.q1 & avail.q3 & avail.min & avail.max &
          !(is.null(n.e) & is.null(n.c))) {
        j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3) &
          !is.na(min) & !is.na(max)
        approx.TE[j] <- "iqr.range"
        if (is.null(n.c))
          TE[j] <- mean_sd_iqr_range(n.e[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
        else if (is.null(n.e))
          TE[j] <- mean_sd_iqr_range(n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
        else
          TE[j] <- mean_sd_iqr_range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
      }
      ##
      ## Use IQR
      ##
      sel.NA <- is.na(TE)
      if (any(sel.NA) &
          avail.median & avail.q1 & avail.q3 &
          !(is.null(n.e) & is.null(n.c))) {
        j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3)
        approx.TE[j] <- "iqr"
        if (is.null(n.c))
          TE[j] <- mean_sd_iqr(n.e[j], median[j], q1[j], q3[j], method.mean)$mean
        else if (is.null(n.e))
          TE[j] <- mean_sd_iqr(n.c[j], median[j], q1[j], q3[j], method.mean)$mean
        else
          TE[j] <- mean_sd_iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                               method.mean)$mean
      }
      ##
      ## Use range
      ##
      sel.NA <- is.na(TE)
      if (any(sel.NA) &
          avail.median & avail.min & avail.max &
          !(is.null(n.e) & is.null(n.c))) {
        j <- sel.NA & !is.na(median) & !is.na(min) & !is.na(max)
        approx.TE[j] <- "range"
        if (is.null(n.c))
          TE[j] <- mean_sd_range(n.e[j], median[j], min[j], max[j],
                                 method.mean)$mean
        else if (is.null(n.e))
          TE[j] <- mean_sd_range(n.c[j], median[j], min[j], max[j],
                                 method.mean)$mean
        else
          TE[j] <- mean_sd_range(n.e[j] + n.c[j], median[j], min[j], max[j],
                                 method.mean)$mean
      }
    }
    else {
      j <- 0
      for (i in approx.TE) {
        j <- j + 1
        ##
        if (i == "ci")
          TE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$TE
        else if (i == "iqr.range")
          if (is.null(n.e) & is.null(n.c))
            stop("Sample size needed if argument 'approx.TE' = \"iqr.range\".",
                 call. = FALSE)
        else if (is.null(n.c))
          TE[j] <- mean_sd_iqr_range(n.e[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
        else if (is.null(n.e))
          TE[j] <- mean_sd_iqr_range(n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
        else
          TE[j] <- mean_sd_iqr_range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
        else if (i == "iqr") {
          if (is.null(n.e) & is.null(n.c))
            stop("Sample size needed if argument 'approx.TE' = \"iqr\".",
                 call. = FALSE)
          else if (is.null(n.c))
            TE[j] <-
              mean_sd_iqr(n.e[j], median[j], q1[j], q3[j], method.mean)$mean
          else if (is.null(n.e))
            TE[j] <-
              mean_sd_iqr(n.c[j], median[j], q1[j], q3[j], method.mean)$mean
          else
            TE[j] <-
              mean_sd_iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                          method.mean)$mean
        }
        else if (i == "range") {
          cat(paste0("Use 'range' for study", j, "\n"))
          if (is.null(n.e) & is.null(n.c))
            stop("Sample size needed if argument 'approx.TE' = \"range\".",
                 call. = FALSE)
          else if (is.null(n.c))
            TE[j] <- mean_sd_range(n.e[j], median[j], min[j], max[j],
                                   method.mean)$mean
          else if (is.null(n.e))
            TE[j] <- mean_sd_range(n.c[j], median[j], min[j], max[j],
                                   method.mean)$mean
          else
            TE[j] <- mean_sd_range(n.e[j] + n.c[j], median[j], min[j], max[j],
                                   method.mean)$mean
        }
      }
    }
  }
  ##
  if (keepdata) {
    if (!isCol(data, ".subset")) {
      data$.TE <- TE
      data$.seTE <- seTE
      #
      if (!avail.approx.TE & any(approx.TE != ""))
        data$.approx.TE <- approx.TE
      if (!avail.approx.seTE & any(approx.seTE != ""))
        data$.approx.seTE <- approx.seTE
    }
    else {
      data$.TE[data$.subset] <- TE
      data$.seTE[data$.subset] <- seTE
      #
      if (!avail.approx.TE & any(approx.TE != "")) {
        data$.approx.TE <- ""
        data$.approx.TE[data$.subset] <- approx.TE
      }
      if (!avail.approx.seTE & any(approx.seTE != "")) {
        data$.approx.seTE <- ""
        data$.approx.seTE[data$.subset] <- approx.seTE
      }
    }
  }
  
  
  ##
  ##
  ## (9) Check standard errors
  ##
  ##
  
  TE   <- int2num(TE)
  seTE <- int2num(seTE)
  ##
  if (any(seTE[!is.na(seTE)] <= 0)) {
    if (warn & !with.cycles)
      warning("Zero values in seTE replaced by NAs.", call. = FALSE)
    seTE[!is.na(seTE) & seTE == 0] <- NA
  }
  ##
  tau2.calc <- NA
  
  
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
  if (is.null(method.tau.ci))
    if (three.level)
      method.tau.ci <- "PL"
    else if (method.tau == "DL")
      method.tau.ci <- "J"
    else
      method.tau.ci <- "QP"
  #
  method.tau.ci <- setchar(method.tau.ci, gs("meth4tau.ci"))
  ##
  if (!three.level & method.tau.ci == "PL") {
    if (method.tau == "DL")
      method.tau.ci <- "J"
    else
      method.tau.ci <- "QP"
  }
  #
  chklevel(level.hetstat)
  #
  if (three.level) {
    chkmlm(method.tau, missing.method.tau, method.predict)
    ##
    common <- FALSE
    ##
    if (!(method.tau %in% c("REML", "ML")))
      method.tau <- "REML"
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
  ## (11) Additional checks and calculations for n-of-1 trials
  ##
  ##
  
  if (with.cycles) {
    df.n_of_1 <- cycles - 1
    sd.n_of_1 <-
      sqrt(sum(df.n_of_1 * seTE^2, na.rm = TRUE) / sum(df.n_of_1, na.rm = TRUE))
    #
    seTE <- sd.n_of_1 / sqrt(cycles)
    #
    if (keepdata) {
      data$.seTE.orig <- data$.seTE
      #
      if (!isCol(data, ".subset"))
        data$.seTE <- seTE
      else
        data$.seTE[data$.subset] <- seTE
    }
  }
  
  
  ##
  ##
  ## (12) Do meta-analysis
  ##
  ##
  
  k <- sum(!is.na(TE[!exclude]) & !is.na(seTE[!exclude]))
  ##
  if (three.level) {
    cluster.incl <- cluster[!exclude]
    k.study <-
      length(unique(
        cluster.incl[!is.na(TE[!exclude]) & !is.na(seTE[!exclude])]))
  }
  else
    k.study <- k
  ##
  seTE.hakn.ci <- seTE.hakn.adhoc.ci <-
    seTE.hakn.pi <- seTE.hakn.adhoc.pi <-
      seTE.kero <- NA
  ##
  pi <- list(seTE = NA, lower = NA, upper = NA, df = NA)
  ##
  df.random <- df.predict <-
    df.hakn <- df.hakn.ci <- df.hakn.pi <- df.kero <- NA
  ##
  if (k == 0) {
    TE.common <- seTE.common <-
      lower.common <- upper.common <-
        statistic.common <- pval.common <- NA
    w.common <- rep(0, k.all)
    ##
    TE.random <- seTE.random <-
      lower.random <- upper.random <-
        statistic.random <- pval.random <- NA
    seTE.classic <- seTE.random
    w.random <- rep(0, k.all)
    ##
    seTE.predict <- df.predict <-
      lower.predict <- upper.predict <-
        seTE.hakn.pi <- seTE.hakn.adhoc.pi <- NA
    ##
    hc <- list(tau2 = NA, se.tau2 = NA, lower.tau2 = NA, upper.tau2 = NA,
               tau = NA, lower.tau = NA, upper.tau = NA,
               method.tau.ci = "", level = NA,
               sign.lower.tau = "", sign.upper.tau = "",
               ##
               Q = NA, df.Q = NA, pval.Q = NA,
               H = NA, lower.H = NA, upper.H = NA,
               I2 = NA, lower.I2 = NA, upper.I2 = NA,
               ##
               Q.resid = NA, df.Q.resid = NA, pval.Q.resid = NA,
               H.resid = NA, lower.H.resid = NA, upper.H.resid = NA,
               I2.resid = NA, lower.I2.resid = NA, upper.I2.resid = NA)
  }
  else {
    ## At least two studies to perform Hartung-Knapp or Kenward-Roger
    ## method
    if (k == 1) {
      method.random.ci <-
        ifelse(method.random.ci %in% c("HK", "KR"), "classic",
               method.random.ci)
      ##
      method.predict <- 
        ifelse(method.predict %in% c("HK", "KR"), "V", method.predict)
      ##
      adhoc.hakn.ci <- ""
      adhoc.hakn.pi <- ""
    }
    ##
    ## Estimate tau-squared
    ##
    hc <- hetcalc(TE[!exclude], seTE[!exclude],
                  method.tau, method.tau.ci, TE.tau,
                  method.I2, level.hetstat, control = control,
                  cluster = cluster[!exclude], rho = rho)
    #
    # Estimate common tau-squared across subgroups
    #
    if (by & tau.common)
      hcc <- hetcalc(TE[!exclude], seTE[!exclude],
                     method.tau, method.tau.ci, TE.tau,
                     method.I2, level.hetstat,
                     subgroup = subgroup, control = control,
                     cluster = cluster[!exclude], rho = rho)
    ##
    ## Different calculations for three-level models
    ##
    if (!three.level) {
      #
      # Between-study heterogeneity
      #
      if (is.null(tau.preset))
        tau2.calc <- if (is.na(sum(hc$tau2))) 0 else sum(hc$tau2)
      else {
        tau2.calc <- tau.preset^2
        ##
        hc$tau2 <- tau.preset^2
        hc$se.tau2 <- hc$lower.tau2 <- hc$upper.tau2 <- NA
        hc$tau <- tau.preset
        hc$lower.tau <- hc$upper.tau <- NA
        hc$method.tau.ci <- ""
        hc$sign.lower.tau <- ""
        hc$sign.upper.tau <- ""
      }
      #
      # Common effect model
      #
      if (!usw.common) {
        #
        # Classic meta-analysis (Cooper & Hedges, 1994, p. 265-6)
        #
        w.common <- 1 / seTE^2
        w.common[is.na(w.common) | is.na(TE) | exclude] <- 0
        ##
        TE.common <- weighted.mean(TE, w.common, na.rm = TRUE)
        #
        if (method.common.ci == "classic") {
          seTE.common <- sqrt(1 / sum(w.common, na.rm = TRUE))
        }
        else if (method.common.ci == "IVhet") {
          seTE.common <-
            sqrt(sum((w.common / (sum(w.common)))^2 * (seTE^2 + tau2.calc)))
        }
        else
          seTE.common <- NA
        #
        ci.c <- ci(TE.common, seTE.common, level = level.ma,
                   null.effect = null.effect)
        statistic.common <- ci.c$statistic
        pval.common <- ci.c$p
        lower.common <- ci.c$lower
        upper.common <- ci.c$upper
      }
      else {
        #
        # Conduct common effect meta-analysis with user-specified weights
        #
        sel.4 <- !is.na(TE) & !is.na(seTE) & !exclude
        #
        m4.usw.c <-
          runUSW(list(yi = TE[sel.4], sei = seTE[sel.4],
                      weights = weights.common[sel.4]),
                 method.tau = "FE",
                 method.random.ci = NULL,
                 level = level.ma,
                 control = control)
        #
        res.usw.c <- extrUSW(m4.usw.c, "FE", null.effect)
        #
        w.common <- res.usw.c$w.common
        TE.common <- res.usw.c$TE.common
        seTE.common <- res.usw.c$seTE.common
        lower.common <- res.usw.c$lower.common
        upper.common <- res.usw.c$upper.common
        statistic.common <- res.usw.c$statistic.common
        pval.common <- res.usw.c$pval.common
        #
        # Drop list for single random effects meta-analysis
        #
        if (length(m4.usw.c) == 1)
          m4.usw.c <- m4.usw.c[[1]]
      }
      #
      # Random effects model
      #
      if (!usw.random) {
        #
        # Classic meta-analysis (Cooper & Hedges, 1994, p. 265, 274-5)
        #
        w.random <- 1 / (seTE^2 + tau2.calc)
        w.random[is.na(w.random) | is.na(TE) | exclude] <- 0
        ##
        TE.random   <- weighted.mean(TE, w.random, na.rm = TRUE)
        seTE.random <- sqrt(1 / sum(w.random, na.rm = TRUE))
        ##
        seTE.classic <- seTE.random
        ##
        ## Kenward-Roger method for confidence or prediction interval
        ##
        kr <- kenwardroger(w.random)
        seTE.kero <- kr$se
        df.kero <- kr$df
        #
        if (any(method.random.ci == "KR") | any(method.predict == "KR")) {
          # Fallback: classic random effects meta-analysis
          if (is.nan(seTE.kero) | is.nan(df.kero)) {
            method.random.ci[method.random.ci == "KR"] <- "classic-KR"
            method.predict[method.predict == "KR"] <- "V-KR"
          }
        }
        ##
        ## Hartung-Knapp method for confidence or prediction interval
        ##
        df.hakn <- k - 1
        q <- 1 / (k - 1) * sum(w.random * (TE - TE.random)^2, na.rm = TRUE)
        ##
        seTE.hakn.ci <- seTE.hakn.adhoc.ci <-
          seTE.hakn.pi <- seTE.hakn.adhoc.pi <-
          sqrt(q / sum(w.random))
        ##
        ## Confidence interval for random effects model
        ##
        ci.r <- as.data.frame(ci(1, NA, level = 0.99999))
        ##
        if (length(adhoc.hakn.ci) == 1) {
          adhoc.hakn.ci <- ifelse(method.random.ci == "HK", adhoc.hakn.ci, "")
        }
        else if (length(adhoc.hakn.ci) == sum(method.random.ci == "HK")) {
          adho <- rep("", length(method.random.ci))
          adho[method.random.ci == "HK"] <- adhoc.hakn.ci
          adhoc.hakn.ci <- adho
        }
        else if (length(method.random.ci) == 1 && method.random.ci == "HK") {
          method.random.ci <- rep("HK", length(adhoc.hakn.ci))
        }
        else if (length(adhoc.hakn.ci) != length(method.random.ci))
          stop("Argument 'adhoc.hakn.ci' must be of same length as ",
               "'method.random.ci' or number of meta-analyses with ",
               "Hartung-Knapp method",
               call. = FALSE)
        ##
        seTE.hakn.adhoc.ci <- rep(seTE.hakn.adhoc.ci, length(method.random.ci))
        df.hakn.ci <- rep(df.hakn, length(method.random.ci))
        ##
        for (i in seq_along(method.random.ci)) {
          if (method.random.ci[i] %in% c("classic", "classic-KR")) {
            ci.r.i <- ci(TE.random, seTE.classic, level = level.ma,
                         null.effect = null.effect)
            ##
            seTE.hakn.adhoc.ci[i] <- NA
          }
          else if (method.random.ci[i] == "HK") {
            if (adhoc.hakn.ci[i] == "se") {
              ##
              ## Variance correction if SE_HK < SE_notHK
              ## (Knapp and Hartung, 2003), i.e., if q < 1
              ##
              if (q < 1)
                seTE.hakn.adhoc.ci[i] <- seTE.classic
            }
            else if (adhoc.hakn.ci[i] == "ci") {
              ##
              ## Use wider confidence interval, i.e., confidence interval
              ## from classic random effects meta-analysis if CI_HK is
              ## smaller
              ## (Wiksten et al., 2016; Jackson et al., 2017, hybrid 2)
              ##
              ci.hk <-
                ci(TE.random, seTE.hakn.ci, level = level.ma, df = df.hakn.ci[i])
              ci.re <-
                ci(TE.random, seTE.classic, level = level.ma)
              ##
              width.hk <- ci.hk$upper - ci.hk$lower
              width.re <- ci.re$upper - ci.re$lower
              ##
              if (width.hk < width.re) {
                seTE.hakn.adhoc.ci[i] <- seTE.classic
                df.hakn.ci[i] <- NA
              }
            }
            else if (adhoc.hakn.ci[i] == "IQWiG6") {
              ##
              ## Variance correction if CI_HK < CI_DL (IQWiG, 2020)
              ##
              ci.hk <-
                ci(TE.random, seTE.hakn.ci, level = level.ma, df = df.hakn.ci[i])
              ##
              m.dl <- metagen(TE, seTE, method.tau = "DL", method.tau.ci = "",
                              method.random.ci = "classic")
              ci.dl <- ci(m.dl$TE.random, m.dl$seTE.classic, level = level.ma)
              ##
              width.hk <- ci.hk$upper - ci.hk$lower
              width.dl <- ci.dl$upper - ci.dl$lower
              ##
              if (width.hk < width.dl)
                seTE.hakn.adhoc.ci[i] <- seTE.classic
            }
            ##
            ci.r.i <- ci(TE.random, seTE.hakn.adhoc.ci[i],
                         level = level.ma, df = df.hakn.ci[i],
                         null.effect = null.effect)
          }
          else if (method.random.ci[i] == "KR") {
            ci.r.i <- ci(TE.random, seTE.kero, level = level.ma, df = df.kero,
                         null.effect = null.effect)
          }
          ##
          ci.r <- rbind(ci.r, as.data.frame(ci.r.i))
        }
        ##
        ci.r <- ci.r[-1, ]
        ##
        seTE.random <- ci.r$seTE
        lower.random <- ci.r$lower
        upper.random <- ci.r$upper
        df.random <- ci.r$df
        statistic.random <- ci.r$statistic
        pval.random <- ci.r$p
        ##
        if (missing.text.random ||
            (length(text.random) == 1 & length(method.random.ci) > 1)) {
          text.random <-
            ifelse(method.random.ci == "classic",
                   text.random,
                   ifelse(method.random.ci %in% c("KR", "classic-KR"),
                          paste0(text.random, " (", method.random.ci, ")"),
                          paste0(text.random, " (HK")))
          text.random <-
            paste0(text.random,
                   ifelse(method.random.ci != "HK",
                          "",
                          ifelse(adhoc.hakn.ci == "",
                                 ")",
                                 paste0("-", toupper(substring(adhoc.hakn.ci, 1, 2)),
                                        ")"))))
        }
        ##
        ## Prediction interval
        ##
        pi <- data.frame()
        #
        if (length(adhoc.hakn.pi) == 1) {
          adhoc.hakn.pi <- ifelse(method.predict == "HK", adhoc.hakn.pi, "")
        }
        else if (length(adhoc.hakn.pi) == sum(method.predict == "HK")) {
          adho <- rep("", length(method.predict))
          adho[method.predict == "HK"] <- adhoc.hakn.pi
          adhoc.hakn.pi <- adho
        }
        else if (length(method.predict) == 1 && method.predict == "HK") {
          method.predict <- rep("HK", length(adhoc.hakn.pi))
        }
        else if (length(adhoc.hakn.pi) != length(method.predict)) {
          stop("Argument 'adhoc.hakn.pi' must be of same length as ",
               "'method.predict' or number of prediction intervals using ",
               "Hartung-Knapp method",
               call. = FALSE)
        }
        #
        seTE.hakn.adhoc.pi <- rep(seTE.hakn.adhoc.pi, length(method.predict))
        df.hakn.pi <- rep(df.hakn, length(method.predict))
        #
        for (i in seq_along(method.predict)) {
          if (method.predict[i] == "HK" && df.hakn.pi[i] > 1) {
            if (adhoc.hakn.pi[i] == "se") {
              #
              # Variance correction if SE_HK < SE_notHK (Knapp and
              # Hartung, 2003), i.e., if q < 1
              #
              if (q < 1)
                seTE.hakn.adhoc.pi[i] <- seTE.classic
            }
            #
            pi.i <- ci(TE.random, sqrt(seTE.hakn.adhoc.pi[i]^2 + tau2.calc),
                       level = level.predict, df = df.hakn.pi[i])
          }
          else if (method.predict[i] == "HK-PR" && df.hakn.pi[i] > 2) {
            if (adhoc.hakn.pi[i] == "se") {
              #
              # Variance correction if SE_HK < SE_notHK (Knapp and
              # Hartung, 2003), i.e., if q < 1
              #
              if (q < 1)
                seTE.hakn.adhoc.pi[i] <- seTE.classic
            }
            #
            pi.i <- ci(TE.random, sqrt(seTE.hakn.adhoc.pi[i]^2 + tau2.calc),
                       level = level.predict, df = df.hakn.pi[i] - 1)
          }
          else if (method.predict[i] %in% c("V", "V-KR") & k > 1) {
            pi.i <- ci(TE.random, sqrt(seTE.classic^2 + tau2.calc),
                       level.predict, k - 1)
          }
          else if (method.predict[i] %in% c("HTS", "HTS-KR") & k > 2) {
            pi.i <- ci(TE.random, sqrt(seTE.classic^2 + tau2.calc),
                       level.predict, k - 2)
          }
          else if (method.predict[i] == "KR" & df.kero > 0) {
            pi.i <- ci(TE.random, sqrt(seTE.kero^2 + tau2.calc),
                       level.predict, df.kero)
          }
          else if (method.predict[i] == "KR-PR" & df.kero > 1) {
            pi.i <- ci(TE.random, sqrt(seTE.kero^2 + tau2.calc),
                       level.predict, df.kero - 1)
          }
          else if (method.predict[i] == "PR" & df.kero > 0) {
            pi.i <- ci(TE.random, sqrt(seTE.kero^2 + tau2.calc),
                       level.predict, df.kero - 1)
          }
          else if (method.predict[i] == "NNF") {
            res.pima <- pimeta::pima(TE[!exclude], seTE[!exclude],
                                     method = "boot",
                                     alpha = 1 - level.predict,
                                     seed = seed.predict)
            #
            pi.i <- as.data.frame(ci(1, NA, level = level.predict))
            pi.i$seTE <- NA
            pi.i$lower <- res.pima$lpi
            pi.i$upper <- res.pima$upi
            pi.i$df <- res.pima$nup
          }
          else if (method.predict[i] == "S")
            pi.i <- ci(TE.random, sqrt(seTE.classic^2 + tau2.calc), level.predict)
          else if (method.predict[i] == "KR")
            pi.i <- ci(TE.random, NA, level.predict, df.kero - 1)
          else if (method.predict[i] == "V")
            pi.i <- ci(TE.random, NA, level.predict, k - 1)
          else if (method.predict[i] == "HTS")
            pi.i <- ci(TE.random, NA, level.predict, k - 2)
          else
            pi.i <- ci(TE.random, NA, level.predict)
          ##
          pi <- rbind(pi, as.data.frame(pi.i))
        }
        #
        seTE.predict <- pi$seTE
        lower.predict <- pi$lower
        upper.predict <- pi$upper
        df.predict <- pi$df
      }
      else {
        #
        # No adhoc method for meta-analysis with user-specified weights
        #
        if ((!missing.adhoc.hakn.ci && any(adhoc.hakn.ci != "")) |
            (!missing.adhoc.hakn.pi && any(adhoc.hakn.pi != ""))) {
          warning("Ad hoc variance correction not implemented ",
                  "for user-specified weights.",
                  call. = FALSE)
          adhoc.hakn.ci[adhoc.hakn.ci != ""] <- ""
          adhoc.hakn.pi[adhoc.hakn.pi != ""] <- ""
        }
        #
        # Conduct random effect meta-analysis with user-specified weights
        #
        sel.4 <- !is.na(TE) & !is.na(seTE) & !exclude
        #
        m4.usw.r <-
          runUSW(list(yi = TE[sel.4], sei = seTE[sel.4],
                      weights = weights.random[sel.4]),
                 method.tau = method.tau,
                 method.random.ci = method.random.ci,
                 level = level.ma,
                 control = control)
        #
        res.usw.r <-
          extrUSW(m4.usw.r, method.tau, null.effect,
                  k, length(TE), sel.4,
                  method.random.ci, method.predict,
                  level.ma, level.predict)
        #
        w.random <- res.usw.r$w.random
        tau2.calc <- sum(res.usw.r$tau2)
        if (is.na(tau2.calc))
          tau2.calc <- 0
        #
        TE.random <- res.usw.r$TE.random
        seTE.random <- res.usw.r$seTE.random
        lower.random <- res.usw.r$lower.random
        upper.random <- res.usw.r$upper.random
        statistic.random <- res.usw.r$statistic.random
        pval.random <- res.usw.r$pval.random
        #
        seTE.classic <- m4.usw.r[[1]]$se
        #
        df.random <- df.hakn <- ifelse(method.random.ci == "HK", k - 1, NA)
        #
        if (missing(text.random) ||
            (length(text.random) == 1 & length(method.random.ci) > 1))
          text.random <-
          ifelse(method.random.ci == "classic",
                 text.random,
                 ifelse(method.random.ci == "HK",
                        paste0(text.random, " (T)"),
                        ""))
        #
        # Prediction interval
        #
        seTE.predict <- res.usw.r$seTE.predict
        lower.predict <- res.usw.r$lower.predict
        upper.predict <- res.usw.r$upper.predict
        df.predict <- res.usw.r$df.predict
        #
        # Drop list for single random effects meta-analysis
        #
        if (length(m4.usw.r) == 1)
          m4.usw.r <- m4.usw.r[[1]]
      }
    }
    else {
      ##
      ## Conduct three-level meta-analysis
      ##
      if (common) {
        ##
        ## No common effect method for three-level model
        ##
        if (!missing.common & common)
          warning(gs("text.common"), " not calculated for three-level model.",
                  call. = FALSE)
        common <- FALSE
      }
      ##
      w.common <- rep_len(NA, length(seTE))
      ##
      TE.common <- seTE.common <- lower.common <- upper.common <-
        statistic.common <- pval.common <- NA
      ##
      ## No adhoc method for three-level models
      ##
      if ((!missing.adhoc.hakn.ci && any(adhoc.hakn.ci != "")) |
          (!missing.adhoc.hakn.pi && any(adhoc.hakn.pi != ""))) {
        warning("Ad hoc variance correction not implemented ",
                "for three-level model.",
                call. = FALSE)
        adhoc.hakn.ci[adhoc.hakn.ci != ""] <- ""
        adhoc.hakn.pi[adhoc.hakn.pi != ""] <- ""
      }
      ##
      ## Conduct three-level meta-analysis
      ##
      sel.4 <- !is.na(TE) & !is.na(seTE) & !exclude
      ##
      list.mlm <- list(yi = TE[sel.4],
                       V = vcalc(vi = seTE[sel.4]^2,
                                 cluster = cluster[sel.4],
                                 obs = idx[sel.4], rho = rho))
             
      ##
      m4 <- runMLM(c(list.mlm,
                     list(data = data.frame(cluster = cluster[sel.4],
                                            idx = idx[sel.4]))),
                   method.tau = method.tau,
                   method.random.ci = method.random.ci,
                   level = level.ma,
                   control = control)
      ##
      res.mlm <-
        extrMLM(m4, k, length(TE), sel.4,
                method.random.ci, method.predict,
                level.ma, level.predict, null.effect)
      ##
      w.random <- res.mlm$w.random
      tau2.calc <- sum(res.mlm$tau2)
      if (is.na(tau2.calc))
        tau2.calc <- 0
      ##
      TE.random <- res.mlm$TE.random
      seTE.random <- res.mlm$seTE.random
      lower.random <- res.mlm$lower.random
      upper.random <- res.mlm$upper.random
      statistic.random <- res.mlm$statistic.random
      pval.random <- res.mlm$pval.random
      ##
      seTE.classic <- m4[[1]]$se
      ##
      df.random <- df.hakn <- ifelse(method.random.ci == "HK", k - 1, NA)
      ##
      if (missing(text.random) ||
          (length(text.random) == 1 & length(method.random.ci) > 1))
        text.random <-
          ifelse(method.random.ci == "classic",
                 text.random,
          ifelse(method.random.ci == "HK",
                 paste0(text.random, " (T)"),
                 ""))
      ##
      ## Prediction interval
      ##
      seTE.predict <- res.mlm$seTE.predict
      lower.predict <- res.mlm$lower.predict
      upper.predict <- res.mlm$upper.predict
      df.predict <- res.mlm$df.predict
      ##
      ## Drop list for single random effects meta-analysis
      ##
      if (length(m4) == 1)
        m4 <- m4[[1]]
    }
  }
  ##
  if (missing(text.predict) ||
      (length(text.predict) == 1 & length(method.predict) > 1)) {
    if (length(method.predict) > 1) {
      text.predict <- paste0(text.predict, " (", method.predict)
      text.predict <-
        paste0(text.predict,
               ifelse(method.predict != "HK",
                      ")",
               ifelse(adhoc.hakn.pi == "",
                      ")",
                      paste0("-", toupper(adhoc.hakn.pi), ")"))))
    }
  }
  
  
  ##
  ##
  ## (13) Heterogeneity measures
  ##
  ##
  
  ##
  ## Calculate Rb (but not for three-level model)
  ##
  if (length(tau2.calc) == 1)
    Rbres <- Rb(seTE[!is.na(seTE)], seTE.classic,
                tau2.calc, hc$Q, hc$df.Q, level.ma)
  else
    Rbres <- list(TE = NA, lower = NA, upper = NA)
  
  
  ##
  ##
  ## (14) Generate R object
  ##
  ##
  if (!avail.detail.tau) {
    if (k != k.study)
      detail.tau <- c("between cluster", "within cluster")
    else
      detail.tau <- ""
  }
  ##
  ci.study <- ci(TE, seTE, level = level,
                 df =
                   if (!is.null(method.ci) && method.ci == "t")
                     df
                   else
                     NULL,
                 null.effect = null.effect)
  ##
  ## Keep original confidence limits
  ##
  if (avail.lower)
    ci.study$lower[!is.na(lower)] <- lower[!is.na(lower)]
  if (avail.upper)
    ci.study$upper[!is.na(upper)] <- upper[!is.na(upper)]
  ##
  if (length(lower.random) > 1) {
    methci <- paste(method.random.ci,
                    toupper(substring(adhoc.hakn.ci, 1, 2)),
                    sep = "-")
    methci <- gsub("-$", "", methci)
    ##
    if (length(TE.random) == 1)
      TE.random <- rep(TE.random, length(lower.random))
    ##
    if (length(seTE.random) == 1)
      seTE.random <- rep(seTE.random, length(lower.random))
    ##
    names(TE.random) <- names(seTE.random) <-
      names(statistic.random) <- names(pval.random) <-
      names(df.random) <- names(lower.random) <- names(upper.random) <-
      methci
    ##
    if (!three.level & !usw.random)
      names(adhoc.hakn.ci) <- names(df.hakn.ci) <-
        names(seTE.hakn.adhoc.ci) <-
        methci
  }
  ##
  if (length(lower.predict) > 1) {
    methpi <- paste(method.predict,
                    toupper(substring(adhoc.hakn.pi, 1, 2)),
                    sep = "-")
    methpi <- gsub("-$", "", methpi)
    ##
    names(seTE.predict) <- names(df.predict) <-
      names(lower.predict) <- names(upper.predict) <-
      methpi
    ##
    if (!three.level & !usw.random)
      names(adhoc.hakn.pi) <- names(df.hakn.pi) <-
        names(seTE.hakn.adhoc.pi) <-
        methpi
  }
  ##
  res <- list(studlab = studlab,
              ##
              sm = sm,
              null.effect = null.effect,
              ##
              TE = TE, seTE = seTE,
              statistic = ci.study$statistic,
              pval = ci.study$p,
              df =
                if (!is.null(method.ci) && method.ci == "t")
                  df
                else
                  rep_len(NA, length(TE)),
              level = level,
              lower = ci.study$lower, upper = ci.study$upper,
              ##
              three.level = three.level,
              cluster = cluster, rho = rho,
              ##
              k = k, k.study = k.study, k.all = k.all, k.TE = sum(!is.na(TE)),
              #
              cycles = cycles,
              sd.n_of_1 = if (with.cycles) sd.n_of_1 else NULL,
              #
              overall = overall,
              overall.hetstat = overall.hetstat,
              common = common,
              random = random,
              prediction = prediction,
              transf = transf,
              backtransf = backtransf,
              func.transf = func.transf,
              func.backtransf = func.backtransf,
              args.transf = args.transf,
              args.backtransf = args.backtransf,
              ##
              method = "Inverse",
              method.random = "Inverse",
              ##
              w.common = w.common,
              TE.common = TE.common,
              seTE.common = seTE.common,
              statistic.common = statistic.common,
              pval.common = pval.common,
              level.ma = level.ma,
              method.common.ci = method.common.ci,
              lower.common = lower.common,
              upper.common = upper.common,
              ##
              w.random = w.random,
              TE.random = TE.random,
              seTE.random = seTE.random,
              statistic.random = statistic.random,
              pval.random = pval.random,
              method.random.ci = method.random.ci,
              df.random = df.random,
              lower.random = lower.random,
              upper.random = upper.random,
              ##
              seTE.classic = seTE.classic,
              #
              weights.common = weights.common,
              weights.random = weights.random,
              #
              adhoc.hakn.ci = adhoc.hakn.ci,
              df.hakn.ci =
                if (any(method.random.ci == "HK")) df.hakn.ci else NA,
              seTE.hakn.ci = seTE.hakn.ci,
              seTE.hakn.adhoc.ci = seTE.hakn.adhoc.ci,
              ##
              df.kero = if (any(method.random.ci == "KR") |
                            any(method.predict == "KR")) df.kero else NA,
              seTE.kero = seTE.kero,
              ##
              method.predict = method.predict,
              adhoc.hakn.pi = adhoc.hakn.pi,
              df.hakn.pi = if (any(method.predict == "HK")) df.hakn.pi else NA,
              ##
              seTE.predict = seTE.predict,
              df.predict = df.predict,
              level.predict = level.predict,
              lower.predict = lower.predict,
              upper.predict = upper.predict,
              seTE.hakn.pi = seTE.hakn.pi,
              seTE.hakn.adhoc.pi = seTE.hakn.adhoc.pi,
              ##
              Q = hc$Q, df.Q = hc$df.Q, pval.Q = hc$pval.Q,
              ##
              method.tau = method.tau,
              control = control,
              method.tau.ci = hc$method.tau.ci,
              level.hetstat = level.hetstat,
              tau2 = hc$tau2,
              se.tau2 = hc$se.tau2,
              lower.tau2 = hc$lower.tau2, upper.tau2 = hc$upper.tau2,
              tau = hc$tau,
              lower.tau = hc$lower.tau, upper.tau = hc$upper.tau,
              tau.preset = tau.preset,
              TE.tau =
                if (!missing(TE.tau) & method.tau == "DL") TE.tau else NULL,
              detail.tau = detail.tau,
              sign.lower.tau = hc$sign.lower.tau,
              sign.upper.tau = hc$sign.upper.tau,
              #
              method.I2 = method.I2,
              #
              H = hc$H, lower.H = hc$lower.H, upper.H = hc$upper.H,
              ##
              I2 = hc$I2, lower.I2 = hc$lower.I2, upper.I2 = hc$upper.I2,
              ##
              Rb = Rbres$TE, lower.Rb = Rbres$lower, upper.Rb = Rbres$upper,
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
              label.left = label.left, label.right = label.right,
              col.label.left = replaceNULL(col.label.left, "black"),
              col.label.right = replaceNULL(col.label.right, "black"),
              #
              keepdata = keepdata,
              data = if (keepdata) data else NULL,
              subset = if (keepdata) subset else NULL,
              exclude = if (!missing.exclude) exclude else NULL,
              ## No general list elements
              n.e = n.e,
              n.c = n.c,
              pscale = pscale,
              irscale = irscale, irunit = irunit,
              method.ci = method.ci,              
              method.mean = method.mean,
              approx.TE = if (all(approx.TE == "")) NULL else approx.TE,
              approx.seTE = if (all(approx.seTE == "")) NULL else approx.seTE,
              ##
              seed.predict = seed.predict,
              #
              pairwise = is.pairwise,
              #
              warn = warn,
              call = match.call(),
              version = packageDescription("meta")$Version,
              # Keep debug information
              tau2.calc = tau2.calc,
              rma.usw.c = if (usw.common & keeprma) m4.usw.c else NULL,
              rma.usw.r = if (usw.random & keeprma) m4.usw.r else NULL,
              rma.three.level = if (three.level & keeprma) m4 else NULL,
              # Deprecated list elements
              zval = ci.study$statistic,
              hakn = any(method.random.ci == "HK"),
              zval.common = statistic.common,
              zval.random = statistic.random
              )
  #
  if (is.pairwise | data.pairwise) {
    res$pairwise <- TRUE
    res$k.study <- length(unique(res$studlab[!is.na(res$TE)]))
  }
  #
  class(res) <- c(fun, "meta")
  ##
  if (avail.lower | avail.upper)
    res$level.ci <- level.ci
  ##
  ## Add results from subgroup analysis
  ##
  if (by) {
    res$subgroup <- subgroup
    res$subgroup.name <- subgroup.name
    ##
    res$tau.common <- tau.common
    res$print.subgroup.name <- print.subgroup.name
    res$sep.subgroup <- sep.subgroup
    res$test.subgroup <- test.subgroup
    res$prediction.subgroup <- prediction.subgroup
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
      if (three.level)
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
    else {
      res$tau2.resid <- res$lower.tau2.resid <- res$upper.tau2.resid <- NA
      res$tau.resid <- res$lower.tau.resid <- res$upper.tau.resid <- NA
      ##
      res$Q.resid <- res$df.Q.resid <- res$pval.Q.resid <- NA
      res$H.resid <- res$lower.H.resid <- res$upper.H.resid <- NA
      res$I2.resid <- res$lower.I2.resid <- res$upper.I2.resid <- NA
    }
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$event.w   <- NULL
    res$n.w       <- NULL
    res$time.e.w  <- NULL
    res$time.c.w  <- NULL
    ##
    res <- setNAwithin(res, res$three.level)
  }
  ##
  ## Add names to tau2 & rest (if necessary)
  ##
  if (length(res$tau2) > 1)
    names(res$tau2) <- res$detail.tau
  ##
  if (length(res$tau) > 1)
    names(res$tau) <- res$detail.tau
  ##
  if (length(res$tau2.resid) > 1)
    names(res$tau2.resid) <- res$detail.tau
  ##
  if (length(res$tau.resid) > 1)
    names(res$tau.resid) <- res$detail.tau
  ##
  ## Unset variables for prediction intervals
  ##
  res$method.predict <-
    ifelse(is.na(res$lower.predict) & is.na(res$upper.predict),
           "", res$method.predict)
  res$df.predict <-
    ifelse(is.na(res$lower.predict) & is.na(res$upper.predict),
           NA, res$df.predict)
  res$adhoc.hakn.pi <-
    ifelse(is.na(res$lower.predict) & is.na(res$upper.predict),
           "", res$adhoc.hakn.pi)
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
