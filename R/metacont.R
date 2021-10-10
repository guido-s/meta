#' Meta-analysis of continuous outcome data
#' 
#' @description
#' Calculation of fixed and random effects estimates for meta-analyses
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
#' @param id An optional vector specifying which estimates come from
#'   the same study resulting in the use of a three-level
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
#'   to calculate confidence intervals for individual studies, see
#'   Details.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
#' @param fixed A logical indicating whether a fixed effect / common
#'   effect meta-analysis should be conducted.
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
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param adhoc.hakn A character string indicating whether an \emph{ad
#'   hoc} variance correction should be applied in the case of an
#'   arbitrarily small Hartung-Knapp variance estimate, see Details.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"PM"},
#'   \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
#'   \code{"HE"}, or \code{"EB"}, can be abbreviated.
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau}. Either \code{"QP"}, \code{"BJ"}, or \code{"J"}, or
#'   \code{""}, can be abbreviated.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param method.bias A character string indicating which test is to
#'   be used. Either \code{"Begg"}, \code{"Egger"}, \code{"Thompson"},
#'   or \code{"Pustejovsky"}, can be abbreviated. See function
#'   \code{\link{metabias}}.
#' @param backtransf A logical indicating whether results for ratio of
#'   means (\code{sm="ROM"}) should be back transformed in printouts
#'   and plots. If TRUE (default), results will be presented as ratio
#'   of means; otherwise log ratio of means will be shown.
#' @param text.fixed A character string used in printouts and forest
#'   plot to label the pooled fixed effect estimate.
#' @param text.random A character string used in printouts and forest
#'   plot to label the pooled random effects estimate.
#' @param text.predict A character string used in printouts and forest
#'   plot to label the prediction interval.
#' @param text.w.fixed A character string used to label weights of
#'   fixed effect model.
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
#' @param byvar Deprecated argument (replaced by 'subgroup').
##' @param keepdata A logical indicating whether original data (set)
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
#' Calculation of fixed and random effects estimates for meta-analyses
#' with continuous outcome data; inverse variance weighting is used
#' for pooling.
#' 
#' Three different types of summary measures are available for continuous
#' outcomes:
#' \itemize{
#' \item mean difference (argument \code{sm = "MD"})
#' \item standardised mean difference (\code{sm = "SMD"})
#' \item ratio of means (\code{sm = "ROM"})
#' }
#' 
#' Default settings are utilised for several arguments (assignments
#' using \code{\link{gs}} function). These defaults can be changed for
#' the current R session using the \code{\link{settings.meta}}
#' function.
#' 
#' Furthermore, R function \code{\link{update.meta}} can be used to
#' rerun a meta-analysis with different settings.
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
#' \code{TE}, \code{TE.fixed}, and \code{TE.random} contain the
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
#' argument \code{method.mean = "Wan"}):
#' \itemize{
#' \item equation (10) if sample size, median, interquartile range and 
#'   range are available,
#' \item equation (14) if sample size, median and interquartile range
#'   are available,
#' \item equation (2) if sample size, median and range are available.
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
#' Note, this choice does not affect the results of the fixed effect
#' and random effects meta-analysis.
#' }
#' 
#' \subsection{Estimation of between-study variance}{
#' 
#' The following methods to estimate the between-study variance
#' \eqn{\tau^2} are available:
#' \itemize{
#' \item DerSimonian-Laird estimator (\code{method.tau = "DL"})
#' \item Paule-Mandel estimator (\code{method.tau = "PM"})
#' \item Restricted maximum-likelihood estimator (\code{method.tau =
#'   "REML"})
#' \item Maximum-likelihood estimator (\code{method.tau = "ML"})
#' \item Hunter-Schmidt estimator (\code{method.tau = "HS"})
#' \item Sidik-Jonkman estimator (\code{method.tau = "SJ"})
#' \item Hedges estimator (\code{method.tau = "HE"})
#' \item Empirical Bayes estimator (\code{method.tau = "EB"})
#' }
#' See \code{\link{metagen}} for more information on these
#' estimators.
#' }
#' 
#' \subsection{Confidence interval for the between-study variance}{
#'
#' The following methods to calculate a confidence interval for
#' \eqn{\tau^2} and \eqn{\tau} are available.
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{Method} \cr 
#' \code{method.tau.ci = "J"}\tab Method by Jackson (2013) \cr
#' \code{method.tau.ci = "BJ"}\tab Method by Biggerstaff and Jackson (2008) \cr
#' \code{method.tau.ci = "QP"}\tab Q-Profile method (Viechtbauer, 2007) \cr
#' \code{method.tau.ci = "PL"}\tab Profile-Likelihood method for
#'  three-level meta-analysis model \cr
#' \tab (Van den Noortgate et al., 2013)
#' }
#' See \code{\link{metagen}} for more information on these methods. No
#' confidence intervals for \eqn{\tau^2} and \eqn{\tau} are calculated
#' if \code{method.tau.ci = ""}.
#' }
#' 
#' \subsection{Hartung-Knapp method}{
#' 
#' Hartung and Knapp (2001) proposed an alternative method for random
#' effects meta-analysis based on a refined variance estimator for the
#' treatment estimate. Simulation studies (Hartung and Knapp, 2001;
#' IntHout et al., 2014; Langan et al., 2019) show improved coverage
#' probabilities compared to the classic random effects
#' method.
#'
#' In rare settings with very homogeneous treatment estimates, the
#' Hartung-Knapp variance estimate can be arbitrarily small resulting
#' in a very narrow confidence interval (Knapp and Hartung, 2003;
#' Wiksten et al., 2016). In such cases, an
#' \emph{ad hoc} variance correction has been proposed by utilising
#' the variance estimate from the classic random effects model with
#' the HK method (Knapp and Hartung, 2003; IQWiQ, 2020). An
#' alternative approach is to use the wider confidence interval of
#' classic fixed or random effects meta-analysis and the HK method
#' (Wiksten et al., 2016; Jackson et al., 2017).
#'
#' Argument \code{adhoc.hakn} can be used to choose the \emph{ad hoc}
#' method:
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{\emph{Ad hoc} method} \cr
#' \code{adhoc.hakn = ""}\tab not used \cr
#' \code{adhoc.hakn = "se"}\tab use variance correction if HK standard
#'  error is smaller \cr
#'  \tab than standard error from classic random effects
#'  \cr
#'  \tab meta-analysis (Knapp and Hartung, 2003) \cr
#' \code{adhoc.hakn = "iqwig6"}\tab use variance correction if HK
#'  confidence interval \cr
#'  \tab is narrower than CI from classic random effects model \cr
#'  \tab with DerSimonian-Laird estimator (IQWiG, 2020) \cr
#' \code{adhoc.hakn = "ci"}\tab use wider confidence interval of
#'  classic random effects \cr
#'  \tab and HK meta-analysis \cr
#'  \tab (Hybrid method 2 in Jackson et al., 2017)
#' }
#' }
#' 
#' \subsection{Prediction interval}{
#' 
#' A prediction interval for the proportion in a new study (Higgins et
#' al., 2009) is calculated if arguments \code{prediction} and
#' \code{random} are \code{TRUE}. Note, the definition of
#' prediction intervals varies in the literature. This function
#' implements equation (12) of Higgins et al., (2009) which proposed a
#' \emph{t} distribution with \emph{K-2} degrees of freedom where
#' \emph{K} corresponds to the number of studies in the meta-analysis.
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
#' Internally, both fixed effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{fixed} and \code{random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{fixed} and \code{random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{random = FALSE}.
#' }
#' 
#' @note
#' The function \code{\link{metagen}} is called internally to
#' calculate individual and overall treatment estimates and standard
#' errors.
#' 
#' @return
#' An object of class \code{c("metacont", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#'
#' \item{n.e, mean.e, sd.e,}{As defined above.}
#' \item{n.c, mean.c, sd.c,}{As defined above.}
#' \item{studlab, exclude, sm, method.ci,}{As defined above.}
#' \item{median.e, q1.e, q3.e, min.e, max.e,}{As defined above.}
#' \item{median.c, q1.c, q3.c, min.c, max.c,}{As defined above.}
#' \item{method.mean, method.sd,}{As defined above.}
#' \item{approx.mean.e, approx.sd.e, approx.mean.c, approx.sd.c,}{As defined above.}
#' \item{level, level.ma,}{As defined above.}
#' \item{fixed, random,}{As defined above.}
#' \item{overall, overall.hetstat,}{As defined above.}
#' \item{pooledvar, method.smd, sd.glass,}{As defined above.}
#' \item{hakn, adhoc.hakn, method.tau, method.tau.ci,}{As defined above.}
#' \item{tau.preset, TE.tau, method.bias,}{As defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{label.e, label.c, label.left, label.right,}{As defined
#'   above.}
#' \item{subgroup, subgroup.name,}{As defined above.}
#' \item{print.subgroup.name, sep.subgroup, warn,}{As defined above.}
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{statistic, pval}{Statistic and p-value for test of treatment
#'   effect for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
#'   standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{statistic.fixed, pval.fixed}{Statistic and p-value for test of
#'   overall treatment effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall treatment effect
#'   and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{statistic.random, pval.random}{Statistic and p-value for test
#'   of overall treatment effect (random effects model).}
#' \item{prediction, level.predict}{As defined above.}
#' \item{seTE.predict}{Standard error utilised for prediction
#'   interval.}
#' \item{lower.predict, upper.predict}{Lower and upper limits of
#'   prediction interval.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic Q.}
#' \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
#' \item{pval.Q}{P-value of heterogeneity test.}
#' \item{tau2}{Between-study variance \eqn{\tau^2}.}
#' \item{se.tau2}{Standard error of \eqn{\tau^2}.}
#' \item{lower.tau2, upper.tau2}{Lower and upper limit of confidence
#'   interval for \eqn{\tau^2}.}
#' \item{tau}{Square-root of between-study variance \eqn{\tau}.}
#' \item{lower.tau, upper.tau}{Lower and upper limit of confidence
#'   interval for \eqn{\tau}.}
#' \item{H}{Heterogeneity statistic H.}
#' \item{lower.H, upper.H}{Lower and upper confidence limit for
#'  heterogeneity statistic H.}
#' \item{I2}{Heterogeneity statistic I\eqn{^2}.}
#' \item{lower.I2, upper.I2}{Lower and upper confidence limit for
#'   heterogeneity statistic I\eqn{^2}.}
#' \item{Rb}{Heterogeneity statistic R\eqn{_b}.}
#' \item{lower.Rb, upper.Rb}{Lower and upper confidence limit for
#'   heterogeneity statistic R\eqn{_b}.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{method}{Pooling method: \code{"Inverse"}.}
#' \item{bylevs}{Levels of grouping variable - if \code{subgroup} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
#'   standard error in subgroups (fixed effect model) - if
#'   \code{subgroup} is not missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{subgroup} is not missing.}
#' \item{statistic.fixed.w, pval.fixed.w}{Statistics and p-values for
#'   test of treatment effect in subgroups (fixed effect model) - if
#'   \code{subgroup} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
#'   standard error in subgroups (random effects model) - if
#'   \code{subgroup} is not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{subgroup} is not missing.}
#' \item{statistic.random.w, pval.random.w}{Statistics and p-values
#'   for test of treatment effect in subgroups (random effects model)
#'   - if \code{subgroup} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{subgroup} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of treatment effect
#'   for Hartung-Knapp method in subgroups - if \code{subgroup} is not
#'   missing and \code{hakn = TRUE}.}
#' \item{n.e.w}{Number of observations in experimental group in
#'   subgroups - if \code{subgroup} is not missing.}
#' \item{n.c.w}{Number of observations in control group in subgroups -
#'   if \code{subgroup} is not missing.}
#' \item{k.w}{Number of studies combined within subgroups - if
#'   \code{subgroup} is not missing.}
#' \item{k.all.w}{Number of all studies in subgroups - if \code{subgroup}
#'   is not missing.}
#' \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
#'   (based on fixed effect model) - if \code{subgroup} is not missing.}
#' \item{Q.w.random}{Overall within subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{subgroup} is not
#'   missing (only calculated if argument \code{tau.common} is TRUE).}
#' \item{df.Q.w}{Degrees of freedom for test of overall within
#'   subgroups heterogeneity - if \code{subgroup} is not missing.}
#' \item{pval.Q.w.fixed}{P-value of within subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{subgroup} is
#'   not missing.}
#' \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{subgroup} is
#'   not missing.}
#' \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic
#'   Q (based on fixed effect model) - if \code{subgroup} is not
#'   missing.}
#' \item{Q.b.random}{Overall between subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{subgroup} is not
#'   missing.}
#' \item{df.Q.b}{Degrees of freedom for test of overall between
#'   subgroups heterogeneity - if \code{subgroup} is not missing.}
#' \item{pval.Q.b.fixed}{P-value of between subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{subgroup} is
#'   not missing.}
#' \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{subgroup} is
#'   not missing.}
#' \item{tau.w}{Square-root of between-study variance within subgroups
#'   - if \code{subgroup} is not missing.}
#' \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{subgroup} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limit for
#'   heterogeneity statistic H within subgroups - if \code{subgroup} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I\eqn{^2} within subgroups - if
#'   \code{subgroup} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limit for
#'   heterogeneity statistic I\eqn{^2} within subgroups - if \code{subgroup} is
#'   not missing.}
#' \item{keepdata}{As defined above.}
#' \item{data}{Original data (set) used in function call (if
#'   \code{keepdata = TRUE}).}
#' \item{subset}{Information on subset of original data used in
#'   meta-analysis (if \code{keepdata = TRUE}).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{update.meta}}, \code{\link{metabin}},
#'   \code{\link{metagen}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009):
#' \emph{Introduction to Meta-Analysis}.
#' Chichester: Wiley
#' 
#' Cohen J (1988):
#' \emph{Statistical Power Analysis for the Behavioral Sciences
#'   (second ed.)}.
#' Lawrence Erlbaum Associates
#' 
#' Cooper H & Hedges LV (1994):
#' \emph{The Handbook of Research Synthesis}.
#' Newbury Park, CA: Russell Sage Foundation
#' 
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
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
#' Hartung J & Knapp G (2001):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
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
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#'
#' IntHout J, Ioannidis JPA, Borm GF (2014):
#' The Hartung-Knapp-Sidik-Jonkman method for random effects
#' meta-analysis is straightforward and considerably outperforms the
#' standard DerSimonian-Laird method.
#' \emph{BMC Medical Research Methodology},
#' \bold{14}, 25
#' 
#' IQWiG (2020):
#' General Methods: Version 6.0.
#' \url{https://www.iqwig.de/en/about-us/methods/methods-paper/}
#' 
#' Jackson D, Law M, Rücker G, Schwarzer G (2017): 
#' The Hartung-Knapp modification for random-effects meta-analysis: A
#' useful refinement but are there any residual concerns?
#' \emph{Statistics in Medicine},
#' \bold{36}, 3923--34
#' 
#' Knapp G & Hartung J (2003):
#' Improved tests for a random effects meta-regression with a single
#' covariate.
#' \emph{Statistics in Medicine},
#' \bold{22}, 2693--710
#'
#' Langan D, Higgins JPT, Jackson D, Bowden J, Veroniki AA,
#' Kontopantelis E, et al. (2019):
#' A comparison of heterogeneity variance estimators in simulated
#' random-effects meta-analyses.
#' \emph{Research Synthesis Methods},
#' \bold{10}, 83--98
#' 
#' Luo D, Wan X, Liu J, Tong T (2018):
#' Optimally estimating the sample mean from the sample size, median,
#' mid-range, and/or mid-quartile range.
#' \emph{Statistical Methods in Medical Research},
#' \bold{27}, 1785--805
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
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
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
#' Wiksten A, Rücker G, Schwarzer G (2016):
#' Hartung-Knapp method is not always conservative compared with
#' fixed-effect meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 2503--15
#' 
#' @examples
#' data(Fleiss1993cont)
#'
#' # Meta-analysis with Hedges' g as effect measure
#' #
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                data = Fleiss1993cont, sm = "SMD")
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
#' m2 <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
#'                n.plac, mean.plac, sqrt(var.plac),
#'                data = amlodipine, studlab = study)
#' m2
#' 
#' # Use pooled variance
#' #
#' update(m2, pooledvar = TRUE)
#' 
#' # Meta-analysis of response ratios (Hedges et al., 1999)
#' #
#' data(woodyplants)
#' m3 <- metacont(n.elev, mean.elev, sd.elev,
#' 		  n.amb, mean.amb, sd.amb,
#'                data = woodyplants, sm = "ROM")
#' m3
#' print(m3, backtransf = FALSE)
#' 
#' @export metacont


metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL, id = NULL,
                     ##
                     median.e, q1.e, q3.e, min.e, max.e,
                     median.c, q1.c, q3.c, min.c, max.c,
                     method.mean = "Luo", method.sd = "Shi",
                     approx.mean.e, approx.mean.c = approx.mean.e,
                     approx.sd.e, approx.sd.c = approx.sd.e,
                     ##
                     sm = gs("smcont"),
                     ##
                     pooledvar = gs("pooledvar"),
                     method.smd = gs("method.smd"),
                     sd.glass = gs("sd.glass"),
                     exact.smd = gs("exact.smd"),
                     ##
                     method.ci = gs("method.ci.cont"),
                     level = gs("level"), level.ma = gs("level.ma"),
                     fixed = gs("fixed"),
                     random = gs("random") | !is.null(tau.preset),
                     overall = fixed | random,
                     overall.hetstat = fixed | random,
                     ##
                     hakn = gs("hakn"), adhoc.hakn = gs("adhoc.hakn"),
                     method.tau = gs("method.tau"),
                     method.tau.ci = gs("method.tau.ci"),
                     tau.preset = NULL, TE.tau = NULL,
                     tau.common = gs("tau.common"),
                     ##
                     prediction = gs("prediction"),
                     level.predict = gs("level.predict"),
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     ##
                     text.fixed = gs("text.fixed"),
                     text.random = gs("text.random"),
                     text.predict = gs("text.predict"),
                     text.w.fixed = gs("text.w.fixed"),
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
                     ##
                     byvar,
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
  chklevel(level)
  ##
  chklogical(hakn)
  adhoc.hakn <- setchar(adhoc.hakn, .settings$adhoc4hakn)
  method.tau <- setchar(method.tau, .settings$meth4tau)
  ##
  missing.id <- missing(id)
  if (is.null(method.tau.ci))
    if (method.tau == "DL")
      method.tau.ci <- "J"
    else if (!missing.id)
      method.tau.ci <- "PL"
    else
      method.tau.ci <- "QP"
  method.tau.ci <- setchar(method.tau.ci, .settings$meth4tau.ci)
  ##
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  method.bias <- setmethodbias(method.bias)
  ##
  if (!is.null(text.fixed))
    chkchar(text.fixed, length = 1)
  if (!is.null(text.random))
    chkchar(text.random, length = 1)
  if (!is.null(text.predict))
    chkchar(text.predict, length = 1)
  if (!is.null(text.w.fixed))
    chkchar(text.w.fixed, length = 1)
  if (!is.null(text.w.random))
    chkchar(text.w.random, length = 1)
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  fun <- "metacont"
  sm <- setchar(sm, .settings$sm4cont)
  if (sm != "MD")
    method.ci <- "z"
  method.ci <- setchar(method.ci, .settings$ci4cont)
  ##
  method.mean <- setchar(method.mean, c("Luo", "Wan"))
  method.sd <- setchar(method.sd, c("Shi", "Wan"))
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
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
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
  print.subgroup.name <- replaceNULL(print.subgroup.name, FALSE)
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
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'n.e', 'mean.e', 'sd.e', 'n.c', 'mean.c', 'sd.c', and 'id'
  ## from data:
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
  if (!missing.id & is.null(id))
    missing.id <- TRUE
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
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  k.All <- length(n.e)
  ##
  mean.e <- eval(mf[[match("mean.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  if (!missing.mean.e)
    chknull(mean.e)
  else
    mean.e <- rep(NA, k.All)
  ##
  sd.e <- eval(mf[[match("sd.e", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  if (!missing.sd.e)
    chknull(sd.e)
  else
    sd.e <- rep(NA, k.All)
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
  ##
  mean.c <- eval(mf[[match("mean.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  if (!missing.mean.c)
    chknull(mean.c)
  else
    mean.c <- rep(NA, k.All)
  ##
  sd.c <- eval(mf[[match("sd.c", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  if (!missing.sd.c)
    chknull(sd.c)
  else
    sd.c <- rep(NA, k.All)
  ##
  id <- eval(mf[[match("id", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  ## Catch 'studlab', 'subgroup', 'subset' and 'exclude' from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  missing.subgroup <- missing(subgroup)
  subgroup <- eval(mf[[match("subgroup", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  missing.byvar <- missing(byvar)
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                          warn.deprecated)
  by <- !is.null(subgroup)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  missing.exclude <- is.null(exclude)
  ##
  ## Catch 'median.e', 'q1.e', 'q3.e', 'min.e', 'max.e', 'median.c',
  ## 'q1.c', 'q3.c', 'min.c', 'max.c', 'approx.mean.e', 'approx.sd.e',
  ## 'approx.mean.c', and 'approx.sd.c', from data:
  ##
  median.e <- eval(mf[[match("median.e", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  ##
  q1.e <- eval(mf[[match("q1.e", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  q3.e <- eval(mf[[match("q3.e", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  min.e <- eval(mf[[match("min.e", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  max.e <- eval(mf[[match("max.e", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  median.c <- eval(mf[[match("median.c", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  ##
  q1.c <- eval(mf[[match("q1.c", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  q3.c <- eval(mf[[match("q3.c", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  min.c <- eval(mf[[match("min.c", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  max.c <- eval(mf[[match("max.c", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  missing.approx.mean.e <- missing(approx.mean.e)
  approx.mean.e <- eval(mf[[match("approx.mean.e", names(mf))]],
                        data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing.approx.mean.e)
    missing.approx.mean.c <- FALSE
  else
    missing.approx.mean.c <- missing(approx.mean.c)
  approx.mean.c <- eval(mf[[match("approx.mean.c", names(mf))]],
                        data, enclos = sys.frame(sys.parent()))
  ##
  missing.approx.sd.e <- missing(approx.sd.e)
  approx.sd.e <- eval(mf[[match("approx.sd.e", names(mf))]],
                      data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing.approx.sd.e)
    missing.approx.sd.c <- FALSE
  else
    missing.approx.sd.c <- missing(approx.sd.c)
  approx.sd.c <- eval(mf[[match("approx.sd.c", names(mf))]],
                      data, enclos = sys.frame(sys.parent()))
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
  if (!missing.id)
    chklength(id, k.All, arg)
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
    if (!missing.id)
      data$.id <- id
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
    if (!missing.id)
      id <- id[subset]
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
    fixed  <- FALSE
    random <- FALSE
    prediction  <- FALSE
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
        subgroup.name <- byvarname(mf[[match("subgroup", names(mf))]])
      else if (!missing.byvar)
        subgroup.name <- byvarname(mf[[match("byvar", names(mf))]])
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
      seTE <- ifelse(npn.n, NA,
                     sqrt(var.pooled * (1 / n.e + 1 / n.c)))
    else
      seTE <- ifelse(npn.n, NA,
                     sqrt(sd.e^2 / n.e + sd.c^2 / n.c))
    ##
    seTE[is.na(TE)] <- NA
    ##
    if (method.ci == "t")
      ci.study <- ci(TE, seTE, df = n.e + n.c - 2)
  }
  else if (sm == "SMD") {
    J <- function(x) gamma(x / 2) / (sqrt(x / 2) * gamma((x - 1) / 2))
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
        J <- function(x) gamma(x / 2) / (sqrt(x / 2) * gamma((x - 1) / 2))
        K <- function(x) 1 - (x - 2) / (x * J(x)^2)
        seTE <- ifelse(npn.n, NA,
                       sqrt(N / (n.e * n.c) + (J(N - 2) * smd)^2 * K(N - 2)))
      }
      else
        seTE <- ifelse(npn.n, NA,
                       sqrt(N / (n.e * n.c) + TE^2 / (2 * N)))
    }
    else if (method.smd == "Hedges") {
      ##
      ## Hedges and Olkin (1985); White and Thomas (2005), p. 143;
      ## formulae used in RevMan 5 (exact.smd = FALSE)
      ##
      if (exact.smd) {
        J <- function(x) gamma(x / 2) / (sqrt(x / 2) * gamma((x - 1) / 2))
        K <- function(x) 1 - (x - 2) / (x * J(x)^2)
      }
      else {
        J <- function(x) 1 - 3 / (4 * x - 1)
        K <- function(x) 1 / (2 * (x - 1.94))
      }
      ##
      TE   <- J(N - 2) * smd
      seTE <- ifelse(npn.n, NA,
                     sqrt(N / (n.e * n.c) + TE^2 * K(N - 2)))
    }
    else if (method.smd == "Glass") {
      ##
      ## see Cooper & Hedges (1994), p. 238
      ##
      n.g  <- if (sd.glass == "control") n.c else n.e
      ##
      TE <- smd
      seTE <- ifelse(npn.n, NA,
                     sqrt(N / (n.e * n.c) + TE^2 / (2 * (n.g - 1))))
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
  ## No three-level meta-analysis conducted if variable 'id' contains
  ## different values for each estimate
  ##
  multi.level <- FALSE
  ##
  sel.ni <- !is.infinite(TE) & !is.infinite(seTE)
  if (!missing.id && length(unique(id[sel.ni])) != length(id[sel.ni]))
    multi.level <- TRUE
  ##
  if (multi.level) {
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
  ## (10) Do meta-analysis
  ##
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               id = id,
               ##
               sm = sm,
               level = level,
               level.ma = level.ma,
               fixed = fixed,
               random = random,
               overall = overall,
               overall.hetstat = overall.hetstat,
               ##
               hakn = hakn, adhoc.hakn = adhoc.hakn,
               method.tau = method.tau, method.tau.ci = method.tau.ci,
               tau.preset = tau.preset,
               TE.tau = TE.tau,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               ##
               text.fixed = text.fixed, text.random = text.random,
               text.predict = text.predict,
               text.w.fixed = text.w.fixed, text.w.random = text.w.random,
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
              exact.smd = exact.smd, method.ci = method.ci)
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
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$n.e.pooled <- sum(res$n.e, na.rm = TRUE)
  res$n.c.pooled <- sum(res$n.c, na.rm = TRUE)
  ##
  res$method.mean <- method.mean
  res$method.sd <- method.sd
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
    res$tau.common <- tau.common
    ##
    if (!tau.common)
      res <- c(res, subgroup(res))
    else if (!is.null(tau.preset))
      res <- c(res, subgroup(res, tau.preset))
    else
      res <- c(res, subgroup(res, hcc$tau.resid))
    ##
    if (!tau.common || !is.null(tau.preset)) {
      res$tau2.resid <- res$lower.tau2.resid <- res$upper.tau2.resid <- NA
      res$tau.resid <- res$lower.tau.resid <- res$upper.tau.resid <- NA
      ##
      res$Q.resid <- res$df.Q.resid <- res$pval.Q.resid <- NA
      res$H.resid <- res$lower.H.resid <- res$upper.H.resid <- NA
      res$I2.resid <- res$lower.I2.resid <- res$upper.I2.resid <- NA
    }
    else {
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
    res$event.w <- NULL
    res$n.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  ## Backward compatibility
  ##
  res$comb.fixed <- fixed
  res$comb.random <- random
  res$level.comb <- level.ma
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
