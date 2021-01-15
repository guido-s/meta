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
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
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
#'   be used.  Either \code{"rank"}, \code{"linreg"}, or \code{"mm"},
#'   can be abbreviated.  See function \code{\link{metabias}}
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots for \code{sm = "MLN"}. If
#'   TRUE (default), results will be presented as means; otherwise
#'   logarithm of means will be shown.
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
#' @param sm A character string indicating which summary measure
#'   (\code{"MRAW"} or \code{"MLN"}) is to be used for pooling of
#'   studies.
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{n}).
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard deviations).
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}}.
#' 
#' @details
#' Fixed effect and random effects meta-analysis of single means to
#' calculate an overall mean; inverse variance weighting is used for
#' pooling. The following transformations of means are implemented to
#' calculate an overall mean:
#' \itemize{
#' \item Raw, i.e. untransformed, means (\code{sm = "MRAW"}, default)
#' \item Log transformed means (\code{sm = "MLN"})
#' }
#' 
#' Note, you should use R function \code{\link{metacont}} to compare
#' means of pairwise comparisons instead of using \code{metamean} for
#' each treatment arm separately which will break randomisation in
#' randomised controlled trials.
#' 
#' Calculations are conducted on the log scale if \code{sm =
#' "ROM"}. Accordingly, list elements \code{TE}, \code{TE.fixed}, and
#' \code{TE.random} contain the logarithm of means. In printouts and
#' plots these values are back transformed if argument
#' \code{backtransf = TRUE}.
#' 
#' Default settings are utilised for several arguments (assignments
#' using \code{\link{gs}} function). These defaults can be changed for
#' the current R session using the \code{\link{settings.meta}}
#' function.
#' 
#' Furthermore, R function \code{\link{update.meta}} can be used to
#' rerun a meta-analysis with different settings.
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
#' \code{method.tau.ci = "J"}\tab Method by Jackson \cr
#' \code{method.tau.ci = "BJ"}\tab Method by Biggerstaff and Jackson \cr
#' \code{method.tau.ci = "QP"}\tab Q-Profile method
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
#' \code{comb.random} are \code{TRUE}. Note, the definition of
#' prediction intervals varies in the literature. This function
#' implements equation (12) of Higgins et al., (2009) which proposed a
#' \emph{t} distribution with \emph{K-2} degrees of freedom where
#' \emph{K} corresponds to the number of studies in the meta-analysis.
#' }
#'
#' \subsection{Subgroup analysis}{
#' 
#' Argument \code{byvar} can be used to conduct subgroup analysis for
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
#' \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{comb.random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{comb.fixed} and \code{comb.random}. E.g. functions
#' \code{\link{print.meta}} and \code{\link{forest.meta}} will not
#' print results for the random effects model if \code{comb.random =
#' FALSE}.
#' }
#' 
#' @note
#' The function \code{\link{metagen}} is called internally to
#' calculate individual and overall treatment estimates and standard
#' errors.
#' 
#' @return
#' An object of class \code{c("metamean", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{n, mean, sd,}{As defined above.}
#' \item{studlab, exclude, sm, method.ci,}{As defined above.}
#' \item{median, q1, q3, min, max,}{As defined above.}
#' \item{method.mean, method.sd,}{As defined above.}
#' \item{approx.mean, approx.sd,}{As defined above.}
#' \item{level, level.comb,}{As defined above.}
#' \item{comb.fixed, comb.random,}{As defined above.}
#' \item{overall, overall.hetstat,}{As defined above.}
#' \item{hakn, adhoc.hakn, method.tau, method.tau.ci,}{As defined above.}
#' \item{tau.preset, TE.tau, method.bias,}{As defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{byvar, bylab, print.byvar, byseparator, warn}{As defined
#'   above.}
#' \item{TE, seTE}{Estimated effect (mean or log mean) and standard
#'   error of individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{statistic, pval}{Statistic and p-value for test of treatment
#'   effect for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall effect (mean or log
#'   mean) and standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{statistic.fixed, pval.fixed}{Statistic and p-value for test of
#'   overall treatment effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall effect (mean or log
#'   mean) and standard error (random effects model).}
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
#' \item{Q}{Heterogeneity statistic.}
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
#' \item{method}{Pooling method: \code{"Inverse"}.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated effect and standard error
#'   in subgroups (fixed effect model) - if \code{byvar} is not
#'   missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{statistic.fixed.w, pval.fixed.w}{Statistics and p-values for
#'   test of treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated effect and standard
#'   error in subgroups (random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{statistic.random.w, pval.random.w}{Statistics and p-values
#'   for test of treatment effect in subgroups (random effects model)
#'   - if \code{byvar} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{byvar} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of effect for
#'   Hartung-Knapp method in subgroups - if \code{byvar} is not
#'   missing and \code{hakn = TRUE}.}
#' \item{n.e.w}{Number of observations in experimental group in
#'   subgroups - if \code{byvar} is not missing.}
#' \item{n.c.w}{Number of observations in control group in subgroups -
#'   if \code{byvar} is not missing.}
#' \item{k.w}{Number of studies combined within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{k.all.w}{Number of all studies in subgroups - if \code{byvar}
#'   is not missing.}
#' \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
#'   (based on fixed effect model) - if \code{byvar} is not missing.}
#' \item{Q.w.random}{Overall within subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{byvar} is not
#'   missing (only calculated if argument \code{tau.common} is TRUE).}
#' \item{df.Q.w}{Degrees of freedom for test of overall within
#'   subgroups heterogeneity - if \code{byvar} is not missing.}
#' \item{pval.Q.w.fixed}{P-value of within subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{byvar} is
#'   not missing.}
#' \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic
#'   Q (based on fixed effect model) - if \code{byvar} is not
#'   missing.}
#' \item{Q.b.random}{Overall between subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{byvar} is not
#'   missing.}
#' \item{df.Q.b}{Degrees of freedom for test of overall between
#'   subgroups heterogeneity - if \code{byvar} is not missing.}
#' \item{pval.Q.b.fixed}{P-value of between subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{byvar} is
#'   not missing.}
#' \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{tau.w}{Square-root of between-study variance within subgroups
#'   - if \code{byvar} is not missing.}
#' \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limit for
#'   heterogeneity statistic H within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I\eqn{^2} within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limit for
#'   heterogeneity statistic I\eqn{^2} within subgroups - if \code{byvar} is
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
#' @seealso \code{\link{update.meta}}, \code{\link{metamean}},
#'   \code{\link{metagen}}
#' 
#' @references
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Hartung J & Knapp G (2001):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
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
#' Langan D, Higgins JPT, Jackson D, Bowden J, Veroniki AA,
#' Kontopantelis E, et al. (2019):
#' A comparison of heterogeneity variance estimators in simulated
#' random-effects meta-analyses.
#' \emph{Research Synthesis Methods},
#' \bold{10}, 83--98
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' Wiksten A, Rücker G, Schwarzer G (2016):
#' Hartung-Knapp method is not always conservative compared with
#' fixed-effect meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 2503--15
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
                     ##
                     median, q1, q3, min, max,
                     method.mean = "Luo", method.sd = "Shi",
                     approx.mean, approx.sd,
                     ##
                     sm = gs("smmean"),
                     ##
                     method.ci = gs("method.ci.cont"),
                     level = gs("level"), level.comb = gs("level.comb"),
                     comb.fixed = gs("comb.fixed"),
                     comb.random = gs("comb.random"),
                     overall = comb.fixed | comb.random,
                     overall.hetstat = comb.fixed | comb.random,
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
                     null.effect = NA,
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
                     ##
                     byvar, bylab, print.byvar = gs("print.byvar"),
                     byseparator = gs("byseparator"),
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"),
                     ##
                     control = NULL
                     ) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chknull(sm)
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  chklogical(overall.hetstat)
  ##
  chklogical(hakn)
  adhoc.hakn <- setchar(adhoc.hakn, .settings$adhoc4hakn)
  method.tau <- setchar(method.tau, .settings$meth4tau)
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
  method.tau.ci <- setchar(method.tau.ci, .settings$meth4tau.ci)
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  chknumeric(null.effect, length = 1)
  ##
  method.bias <- setchar(method.bias, .settings$meth4bias)
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
  ## Additional arguments / checks for metamean objects
  ##
  fun <- "metamean"
  sm <- setchar(sm, .settings$sm4mean)
  if (sm != "MRAW")
    method.ci <- "z"
  method.ci <- setchar(method.ci, .settings$ci4cont)
  ##
  method.mean <- setchar(method.mean, c("Luo", "Wan"))
  method.sd <- setchar(method.sd, c("Shi", "Wan"))
  ##
  chklogical(warn)
  
  
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
  ## Catch 'n', 'mean', and 'sd' from data:
  ##
  missing.mean <- missing(mean)
  missing.sd <- missing(sd)
  ##
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
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(n)
  k.All <- length(n)
  ##
  mean <- eval(mf[[match("mean", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  if (!missing.mean)
    chknull(mean)
  else
    mean <- rep(NA, k.All)
  ##
  sd <- eval(mf[[match("sd", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  if (!missing.sd)
    chknull(sd)
  else
    sd <- rep(NA, k.All)
  ##
  ## Catch 'studlab', 'byvar', 'subset' and 'exclude' from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  by <- !is.null(byvar)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  missing.exclude <- is.null(exclude)
  ##
  ## Catch 'median', 'q1', 'q3', 'min', 'max', 'approx.mean', and
  ## 'approx.sd', from data:
  ##
  median <- eval(mf[[match("median", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  q1 <- eval(mf[[match("q1", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  q3 <- eval(mf[[match("q3", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  min <- eval(mf[[match("min", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  max <- eval(mf[[match("max", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  missing.approx.mean <- missing(approx.mean)
  approx.mean <- eval(mf[[match("approx.mean", names(mf))]],
                      data, enclos = sys.frame(sys.parent()))
  ##
  missing.approx.sd <- missing(approx.sd)
  approx.sd <- eval(mf[[match("approx.sd", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(mean, k.All, fun)
  chklength(sd, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (!missing.median)
    chklength(median, k.All, fun)
  if (!missing.q1)
    chklength(q1, k.All, fun)
  if (!missing.q3)
    chklength(q3, k.All, fun)
  if (!missing.min)
    chklength(min, k.All, fun)
  if (!missing.max)
    chklength(max, k.All, fun)
  ##
  if (!missing.approx.mean) {
    if (length(approx.mean) == 1)
      rep_len(approx.mean, k.All)
    else
      chklength(approx.mean, k.All, fun)
    ##
    approx.mean <- setchar(approx.mean, c("", "iqr.range", "iqr", "range"))
  }
  ##
  if (!missing.approx.sd) {
    if (length(approx.sd) == 1)
      rep_len(approx.sd, k.All)
    else
      chklength(approx.sd, k.All, fun)
    ##
    approx.sd <- setchar(approx.sd, c("", "iqr.range", "iqr", "range"))
  }
  ##
  if (by)
    chklength(byvar, k.All, fun)
  ##
  ## Additional checks
  ##
  if (!by & tau.common) {
    warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
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
  if (by) {
    chkmiss(byvar)
    byvar.name <- byvarname(mf[[match("byvar", names(mf))]])
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  
  
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
    if (!missing.median)
      data$.median <- median
    if (!missing.q1)
      data$.q1 <- q1
    if (!missing.q3)
      data$.q3 <- q3
    if (!missing.min)
      data$.min <- min
    if (!missing.max)
      data$.max <- max
    if (!missing.approx.mean)
      data$.approx.mean <- approx.mean
    if (!missing.approx.sd)
      data$.approx.sd <- approx.sd
    ##
    if (by)
      data$.byvar <- byvar
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
    exclude <- exclude[subset]
    ##
    if (!missing.median)
      median <- median[subset]
    if (!missing.q1)
      q1 <- q1[subset]
    if (!missing.q3)
      q3 <- q3[subset]
    if (!missing.min)
      min <- min[subset]
    if (!missing.max)
      max <- max[subset]
    if (!missing.approx.mean)
      approx.mean <- approx.mean[subset]
    if (!missing.approx.sd)
      approx.sd <- approx.sd[subset]
    ##
    if (by)
      byvar <- byvar[subset]
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
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
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
  if (!missing.median)
    chknumeric(median)
  if (!missing.q1)
    chknumeric(q1)
  if (!missing.q3)
    chknumeric(q3)
  if (!missing.min)
    chknumeric(min)
  if (!missing.max)
    chknumeric(max)
  ##
  ## Recode integer as numeric:
  ##
  n    <- int2num(n)
  mean <- int2num(mean)
  sd   <- int2num(sd)
  ##
  if (!missing.median)
    median <- int2num(median)
  if (!missing.q1)
    q1 <- int2num(q1)
  if (!missing.q3)
    q3 <- int2num(q3)
  if (!missing.min)
    min <- int2num(min)
  if (!missing.max)
    max <- int2num(max)
  
  
  ##
  ##
  ## (7) Calculate means from other information
  ##
  ##
  if (missing.approx.mean) {
    approx.mean <- rep_len("", length(n))
    ##
    ## (a) Use IQR and range
    ##
    sel.NA <- is.na(mean)
    if (any(sel.NA) & !missing.median &
        !missing.q1 & !missing.q3 &
        !missing.min & !missing.max) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3) &
        !is.na(min) & !is.na(max)
      approx.mean[j] <- "iqr.range"
      ##
      mean[j] <- mean.sd.iqr.range(n[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA <- is.na(mean)
    if (any(sel.NA) & !missing.median & !missing.q1 & !missing.q3) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3)
      approx.mean[j] <- "iqr"
      mean[j] <- mean.sd.iqr(n[j], median[j], q1[j], q3[j],
                               method.mean)$mean
    }
    ##
    ## (c) Use range
    ##
    sel.NA <- is.na(mean)
    if (any(sel.NA) & !missing.median & !missing.min & !missing.max) {
      j <- sel.NA & !is.na(median) & !is.na(min) & !is.na(max)
      approx.mean[j] <- "range"
      mean[j] <- mean.sd.range(n[j], median[j], min[j], max[j],
                                 method.mean)$mean
    }
  }
  else {
    j <- 0
    for (i in approx.mean) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        mean[j] <- mean.sd.iqr.range(n[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
      else if (i == "iqr")
        mean[j] <- mean.sd.iqr(n[j], median[j], q1[j], q3[j],
                                 method.mean)$mean
      else if (i == "range")
        mean[j] <- mean.sd.range(n[j], median[j], min[j], max[j],
                                   method.mean)$mean
    }
  }
  
  
  ##
  ##
  ## (8) Calculate standard deviation from other information
  ##
  ##
  if (missing.median) {
    median.sd <- mean
    missing.median <- FALSE
    export.median <- FALSE
  }
  else {
    median.sd <- median
    median.sd[is.na(median.sd)] <- mean[is.na(median.sd)]
    export.median <- TRUE
  }
  ##
  if (missing.approx.sd) {
    approx.sd <- rep_len("", length(n))
    ##
    ## (a) Use IQR and range
    ##
    sel.NA <- is.na(sd)
    if (any(sel.NA) & !missing.median &
        !missing.q1 & !missing.q3 &
        !missing.min & !missing.max) {
      j <- sel.NA & !is.na(median.sd) & !is.na(q1) & !is.na(q3) &
        !is.na(min) & !is.na(max)
      approx.sd[j] <- "iqr.range"
      ##
      sd[j] <- mean.sd.iqr.range(n[j], median.sd[j], q1[j], q3[j],
                                   min[j], max[j],
                                   method.sd = method.sd)$sd
    }
    ##
    ## (b) Use IQR
    ##
    sel.NA <- is.na(sd)
    if (any(sel.NA) & !missing.median & !missing.q1 & !missing.q3) {
      j <- sel.NA & !is.na(median.sd) & !is.na(q1) & !is.na(q3)
      approx.sd[j] <- "iqr"
      sd[j] <- mean.sd.iqr(n[j], median.sd[j], q1[j], q3[j])$sd
    }
    ##
    ## (c) Use range
    ##
    sel.NA <- is.na(sd)
    if (any(sel.NA) & !missing.median & !missing.min & !missing.max) {
      j <- sel.NA & !is.na(median.sd) & !is.na(min) & !is.na(max)
      approx.sd[j] <- "range"
      sd[j] <- mean.sd.range(n[j], median.sd[j], min[j], max[j])$sd
    }
  }
  else {
    j <- 0
    for (i in approx.sd) {
      j <- j + 1
      ##
      if (i == "iqr.range")
        sd[j] <- mean.sd.iqr.range(n[j], median.sd[j], q1[j], q3[j],
                                     min[j], max[j],
                                     method.sd = method.sd)$sd
      else if (i == "iqr")
        sd[j] <- mean.sd.iqr(n[j], median.sd[j], q1[j], q3[j])$sd
      else if (i == "range")
        sd[j] <- mean.sd.range(n[j], median.sd[j], min[j], max[j])$sd
    }
  }
  ##
  if (keepdata) {
    if (!isCol(data, ".subset")) {
      data$.sd <- sd
      data$.mean <- mean
      data$.approx.sd <- approx.sd
      data$.approx.mean <- approx.mean
    }
    else {
      data$.sd[data$.subset] <- sd
      data$.mean[data$.subset] <- mean
      data$.approx.sd[data$.subset] <- approx.sd
      data$.approx.mean[data$.subset] <- approx.mean
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
    warning("Studies with non-positive sample size get no weight in meta-analysis.")
  ##
  if (sm == "MRAW") {
    TE <- ifelse(npn.n, NA, mean)
    ##
    seTE <- ifelse(npn.n, NA, sqrt(sd^2 / n))
    ##
    seTE[is.na(TE)] <- NA
    ##
    if (method.ci == "t")
      ci.study <- ci(TE, seTE, df = n - 1)
    ##
    transf.null.effect <- null.effect
  }
  ##
  else if (sm == "MLN") {
    npn.mean <- npn(mean)
    ##
    if (any(npn.mean) & warn)
      warning("Studies with negative or zero mean get no weight in meta-analysis.")

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
    warning("Studies with non-positive standard deviation get no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  
  
  ##
  ##
  ## (10) Do meta-analysis
  ##
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               ##
               sm = sm,
               level = level,
               level.comb = level.comb,
               comb.fixed = comb.fixed,
               comb.random = comb.random,
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
               null.effect = transf.null.effect,
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
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (by & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   TE.tau, level.comb, byvar, control)
  }
  
  
  ##
  ##
  ## (11) Generate R object
  ##
  ##
  res <- list(n = n, mean = mean, sd = sd, method.ci = method.ci)
  ##
  if (export.median)
    res$median <- median
  if (!missing.q1)
    res$q1 <- q1
  if (!missing.q3)
    res$q3 <- q3
  if (!missing.min)
    res$min <- min
  if (!missing.max)
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
  m$label.e <- ""
  m$label.c <- ""
  m$label.left <- ""
  m$label.right <- ""
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Add data
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
  if (all(res$approx.mean == "")) {
    res$approx.mean <- NULL
    res$data$.approx.mean <- NULL
  }
  if (all(res$approx.sd == "")) {
    res$approx.sd <- NULL
    res$data$.approx.sd <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (by) {
    res$byvar <- byvar
    res$bylab <- bylab
    res$print.byvar <- print.byvar
    res$byseparator <- byseparator
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
      res$Q.w.random <- hcc$Q.resid
      res$df.Q.w.random <- hcc$df.Q.resid
      res$pval.Q.w.random <- hcc$pval.Q.resid
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
  class(res) <- c(fun, "meta")
  
  
  res
}
