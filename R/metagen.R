#' Generic inverse variance meta-analysis
#' 
#' @description
#' 
#' Fixed effect and random effects meta-analysis based on estimates
#' (e.g. log hazard ratios) and their standard errors. The inverse
#' variance method is used for pooling.
#' 
#' @param TE Estimate of treatment effect, e.g., log hazard ratio or
#'   risk difference.
#' @param seTE Standard error of treatment estimate.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used (see Details).
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots (see Details).
#' @param id An optional vector specifying which estimates come from
#'   the same study resulting in the use of a three-level
#'   meta-analysis model.
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
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
#' @param n.e Number of observations in experimental group (or total
#'   sample size in study).
#' @param n.c Number of observations in control group.
#' @param pval P-value (used to estimate the standard error).
#' @param df Degrees of freedom (used in test or to construct
#'   confidence interval).
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
#' @param hakn A logical indicating whether method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param adhoc.hakn A character string indicating whether an \emph{ad
#'   hoc} variance correction should be applied in the case of an
#'   arbitrarily small Hartung-Knapp variance estimate. Either
#'   \code{""}, \code{"se"}, \code{"ci"}, or \code{"iqwig6"} (see
#'   Details), can be abbreviated.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"PM"},
#'   \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
#'   \code{"HE"}, or \code{"EB"}, can be abbreviated.
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau}. Either \code{"QP"}, \code{"BJ"}, \code{"J"},
#'   \code{"PL"}, or \code{""}, can be abbreviated.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param detail.tau Detail on between-study variance estimate.
#' @param method.bias A character string indicating which test is to
#'   be used.  Either \code{"rank"}, \code{"linreg"}, or \code{"mm"},
#'   can be abbreviated.  See function \code{\link{metabias}}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE} (default), results for \code{sm = "OR"} are printed as odds
#'   ratios rather than log odds ratios and results for \code{sm =
#'   "ZCOR"} are printed as correlations rather than Fisher's z
#'   transformed correlations, for example.
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
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{TE}).
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
#'   standard errors).
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.mv}}.
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
#' sample size, median, interquartile range and range. By default,
#' methods described in Luo et al. (2018) are utilized (argument
#' \code{method.mean = "Luo"}):
#' \itemize{
#' \item equation (7) if sample size, median and range are available,
#' \item equation (11) if sample size, median and interquartile range
#'   are available,
#' \item equation (15) if sample size, median, range and interquartile
#'   range are available.
#' }
#' 
#' Instead the methods described in Wan et al. (2014) are used if
#' argument \code{method.mean = "Wan"}):
#' \itemize{
#' \item equation (2) if sample size, median and range are available,
#' \item equation (14) if sample size, median and interquartile range
#'   are available,
#' \item equation (10) if sample size, median, range and interquartile
#'   range are available.
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
#' on the standard normal or \emph{t} distribution if argument
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
#' \item t-distribution (\code{method.ci = "t"}).
#' }
#'
#' By default, the first method is used if argument \code{df} is
#' missing and the second method otherwise.
#' 
#' Note, this choice does not affect the results of the fixed effect
#' and random effects meta-analysis.
#' }
#' 
#' \subsection{Estimation of between-study variance}{
#' 
#' The following methods are available to estimate the between-study
#' variance \eqn{\tau^2}.
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{Method} \cr 
#' \code{method.tau = "DL"}
#'  \tab DerSimonian-Laird estimator (DerSimonian and Laird, 1986) \cr
#' \code{method.tau = "PM"}
#'  \tab Paule-Mandel estimator (Paule and Mandel, 1982) \cr
#' \code{method.tau = "REML"}
#'  \tab Restricted maximum-likelihood estimator (Viechtbauer, 2005) \cr
#' \code{method.tau = "ML"}
#'  \tab Maximum-likelihood estimator (Viechtbauer, 2005) \cr
#' \code{method.tau = "HS"}
#'  \tab Hunter-Schmidt estimator (Hunter and Schmidt, 2015) \cr
#' \code{method.tau = "SJ"}
#'  \tab Sidik-Jonkman estimator (Sidik and Jonkman, 2005) \cr
#' \code{method.tau = "HE"}
#'  \tab Hedges estimator (Hedges and Olkin, 1985) \cr
#' \code{method.tau = "EB"}
#'  \tab Empirical Bayes estimator (Morris, 1983)
#' }
#'
#' Historically, the DerSimonian-Laird method was the de facto
#' standard to estimate the between-study variance \eqn{\tau^2} and is
#' still the default in many software packages including Review
#' Manager 5 (RevMan 5) and R package \pkg{meta}. However, its role
#' has been challenged and especially the Paule-Mandel and REML
#' estimators have been recommended (Veroniki et al.,
#' 2016). Accordingly, the following R command can be used to use the
#' Paule-Mandel estimator in all meta-analyses of the R session:
#' \code{settings.meta(method.tau = "PM")}
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
#' The first three methods have been recommended by Veroniki et
#' al. (2016). By default, the Jackson method is used for the
#' DerSimonian-Laird estimator of \eqn{\tau^2} and the Q-profile
#' method for all other estimators of \eqn{\tau^2}. The
#' Profile-Likelihood method is the only method available for the
#' three-level meta-analysis model. No confidence intervals for
#' \eqn{\tau^2} and \eqn{\tau} are calculated if \code{method.tau.ci =
#' ""}.
#' }
#' 
#' \subsection{Hartung-Knapp method}{
#' 
#' Hartung and Knapp (2001a,b) proposed an alternative method for
#' random effects meta-analysis based on a refined variance estimator
#' for the treatment estimate. Simulation studies (Hartung and Knapp,
#' 2001a,b; IntHout et al., 2014; Langan et al., 2019) show improved
#' coverage probabilities compared to the classic random effects
#' method. However, in rare settings with very homogeneous treatment
#' estimates, the Hartung-Knapp (HK) variance estimate can be
#' arbitrarily small resulting in a very narrow confidence interval
#' (Knapp and Hartung, 2003; Wiksten et al., 2016). In such cases, an
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
#' A prediction interval for the treatment effect of a new study
#' (Higgins et al., 2009) is calculated if arguments \code{prediction}
#' and \code{comb.random} are \code{TRUE}. Note, the definition of
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
#' i.e., \emph{ln}(RR) = 0 or \emph{ln}(OR) = 0 which is equivalent to
#' testing RR = 1 or OR = 1.
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
#' Internally, both fixed effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{comb.random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{comb.fixed} and \code{comb.random}. For example, functions
#' \code{\link{print.meta}} and \code{\link{forest.meta}} will not
#' show results for the random effects model if \code{comb.random =
#' FALSE}.
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
#' Default settings for \code{comb.fixed}, \code{comb.random},
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
#' 
#' @return
#' An object of class \code{c("metagen", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{TE, seTE, studlab, exclude, n.e, n.c}{As defined above.}
#' \item{id, sm, method.ci, level, level.comb,}{As defined above.}
#' \item{comb.fixed, comb.random,}{As defined above.}
#' \item{overall, overall.hetstat,}{As defined above.}
#' \item{hakn, adhoc.hakn, method.tau, method.tau.ci,}{As defined above.}
#' \item{tau.preset, TE.tau, method.bias,}{As defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{label.e, label.c, label.left, label.right,}{As defined
#'   above.}
#' \item{byvar, bylab, print.byvar, byseparator, warn}{As defined
#'   above.}
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
#' \item{null.effect}{As defined above.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic.}
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
#' \item{approx.TE, approx.seTE}{As defined above.}
#' \item{method}{Pooling method: \code{"Inverse"}.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
#'   standard error in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{statistic.fixed.w, pval.fixed.w}{Statistics and p-values for
#'   test of treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
#'   standard error in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{statistic.random.w, pval.random.w}{Statistics and p-values
#'   for test of treatment effect in subgroups (random effects model)
#'   - if \code{byvar} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{byvar} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of treatment effect
#'   for Hartung-Knapp method in subgroups - if \code{byvar} is not
#'   missing and \code{hakn = TRUE}.}
#' \item{n.harmonic.mean.w}{Harmonic mean of number of observations in
#'   subgroups (for back transformation of Freeman-Tukey Double
#'   arcsine transformation) - if \code{byvar} is not missing.}
#' \item{n.e.w}{Number of observations in experimental group in
#'   subgroups - if \code{byvar} is not missing.}
#' \item{n.c.w}{Number of observations in control group in subgroups -
#'   if \code{byvar} is not missing.}
#' \item{k.w}{Number of studies combined within
#'   subgroups - if \code{byvar} is not missing.}
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
#'   not missing.}  \item{pval.Q.b.random}{P-value of between
#'   subgroups heterogeneity statistic Q (based on random effects
#'   model) - if \code{byvar} is not missing.}
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
#' @seealso \code{\link{update.meta}}, \code{\link{metabin}},
#'   \code{\link{metacont}}, \code{\link{print.meta}},
#'   \code{\link{settings.meta}}
#' 
#' @references
#' Biggerstaff BJ, Jackson D (2008):
#' The exact distribution of Cochran’s heterogeneity statistic in
#' one-way random effects meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{27}, 6093--110
#'
#' Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
#' A basic introduction to fixed-effect and random-effects models for
#' meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{1}, 97--111
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
#' Hedges LV & Olkin I (1985):
#' \emph{Statistical methods for meta-analysis}.
#' San Diego, CA: Academic Press
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' Hunter JE & Schmidt FL (2015):
#' \emph{Methods of Meta-Analysis: Correcting Error and Bias in
#' Research Findings} (Third edition).
#' Thousand Oaks, CA: Sage
#' 
#' Hartung J, Knapp G (2001a):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
#' 
#' Hartung J, Knapp G (2001b):
#' A refined method for the meta-analysis of controlled clinical
#' trials with binary outcome.
#' \emph{Statistics in Medicine},
#' \bold{20}, 3875--89
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
#' Jackson D (2013):
#' Confidence intervals for the between-study variance in random
#' effects meta-analysis using generalised Cochran heterogeneity
#' statistics.
#' \emph{Research Synthesis Methods},
#' \bold{4}, 220--229
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
#' Morris CN (1983):
#' Parametric empirical Bayes inference: Theory and applications (with
#' discussion).
#' \emph{Journal of the American Statistical Association}
#' \bold{78}, 47--65
#' 
#' Paule RC & Mandel J (1982):
#' Consensus values and weighting factors.
#' \emph{Journal of Research of the National Bureau of Standards},
#' \bold{87}, 377--85
#' 
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.4.
#' The Cochrane Collaboration, 2020
#'
#' Shi J, Luo D, Weng H, Zeng X-T, Lin L, Chu H, et al. (2020):
#' Optimally estimating the sample standard deviation from the
#' five-number summary.
#' \emph{Research Synthesis Methods}.
#'
#' Sidik K & Jonkman JN (2005):
#' Simple heterogeneity variance estimation for meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)},
#' \bold{54}, 367--84
#'
#' Veroniki AA, Jackson D, Viechtbauer W, Bender R, Bowden J, Knapp G,
#' et al. (2016):
#' Methods to estimate the between-study variance and its uncertainty
#' in meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{7}, 55--79 
#' 
#' Viechtbauer W (2005):
#' Bias and efficiency of meta-analytic variance estimators in the
#' random-effects model.
#' \emph{Journal of Educational and Behavioral Statistics},
#' \bold{30}, 261--93
#' 
#' Viechtbauer W (2007):
#' Confidence intervals for the amount of heterogeneity in
#' meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{26}, 37--52
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
#' Wiksten A, Rücker G, Schwarzer G (2016):
#' Hartung-Knapp method is not always conservative compared with
#' fixed-effect meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 2503--15
#' 
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, study,
#'               data = Fleiss1993bin, sm = "RR", method = "I")
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
#' 
#' # Meta-analysis with prespecified between-study variance
#' #
#' summary(metagen(m1$TE, m1$seTE, sm = "RR", tau.preset = sqrt(0.1)))
#'
#' 
#' # Meta-analysis of survival data:
#' #
#' logHR <- log(c(0.95, 1.5))
#' selogHR <- c(0.25, 0.35)
#' metagen(logHR, selogHR, sm = "HR")
#'
#' 
#' # Paule-Mandel method to estimate between-study variance for data
#' # from Paule & Mandel (1982)
#' #
#' average <- c(27.044, 26.022, 26.340, 26.787, 26.796)
#' variance <- c(0.003, 0.076, 0.464, 0.003, 0.014)
#' #
#' summary(metagen(average, sqrt(variance), sm = "MD", method.tau = "PM"))
#'
#' 
#' # Conduct meta-analysis using hazard ratios and 95% confidence intervals
#' #
#' # Data from Steurer et al. (2006), Analysis 1.1 Overall survival
#' # https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD004270.pub2/abstract
#' #
#' study <- c("FCG on CLL 1996", "Leporrier 2001", "Rai 2000", "Robak 2000")
#' HR <- c(0.55, 0.92, 0.79, 1.18)
#' lower.HR <- c(0.28, 0.79, 0.59, 0.64)
#' upper.HR <- c(1.09, 1.08, 1.05, 2.17)
#' #
#' # Input must be log hazard ratios, not hazard ratios
#' #
#' metagen(log(HR), lower = log(lower.HR), upper = log(upper.HR),
#'         studlab = study, sm = "HR")
#'
#' 
#' # Exclude MRC-1 and MRC-2 studies from meta-analysis, however,
#' # show them in printouts and forest plots
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'         data = Fleiss1993bin, sm = "RR", method = "I",
#'         exclude = study %in% c("MRC-1", "MRC-2"))
#' #
#' # Exclude MRC-1 and MRC-2 studies completely from meta-analysis
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'         data = Fleiss1993bin, sm = "RR", method = "I",
#'         subset = !(study %in% c("MRC-1", "MRC-2")))
#'
#' 
#' # Exclude studies with total sample size above 1500
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'         data = Fleiss1993bin, sm = "RR", method = "I",
#'         exclude = (n.asp + n.plac) > 1500)
#'
#' # Exclude studies containing "MRC" in study name
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'         data = Fleiss1993bin, sm = "RR", method = "I",
#'         exclude = grep("MRC", study))
#'
#' # Use both arguments 'subset' and 'exclude'
#' #
#' metabin(d.asp, n.asp, d.plac, n.plac, study,
#'         data = Fleiss1993bin, sm = "RR", method = "I",
#'         subset = (n.asp + n.plac) > 1500,
#'         exclude = grep("MRC", study))
#' 
#' @export metagen


metagen <- function(TE, seTE, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL, id = NULL,
                    ##
                    sm = "",
                    ##
                    method.ci = if (missing(df)) "z" else "t",
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
                    detail.tau = "",
                    ##
                    prediction = gs("prediction"),
                    level.predict = gs("level.predict"),
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
                    backtransf = gs("backtransf"),
                    pscale = 1,
                    irscale = 1, irunit = "person-years",
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
  ##
  method.ci <- setchar(method.ci, .settings$ci4cont)
  ##
  method.mean <- setchar(method.mean, c("Luo", "Wan"))
  method.sd <- setchar(method.sd, c("Shi", "Wan"))
  ##
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  chklogical(overall.hetstat)
  ##
  chklogical(hakn)
  adhoc.hakn <- setchar(adhoc.hakn, .settings$adhoc4hakn)
  missing.method.tau <- missing(method.tau)
  method.tau <- setchar(method.tau, .settings$meth4tau)
  ##
  missing.id <- missing(id)
  ##
  missing.method.tau.ci <- missing(method.tau.ci)
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
  chknumeric(null.effect, length = 1)
  ##
  method.bias <- setchar(method.bias, .settings$meth4bias)
  ##
  chklogical(backtransf)
  if (!is.prop(sm))
    pscale <- 1
  chknumeric(pscale, length = 1)
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate(sm))
    irscale <- 1
  chknumeric(irscale, length = 1)
  if (!backtransf & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
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
  ## Additional arguments / checks
  ##
  fun <- "metagen"
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
  ## Catch 'TE', 'seTE', 'median', 'lower', 'upper', 'n.e', 'n.c', and
  ## 'id' from data:
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
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  seTE <- eval(mf[[match("seTE", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  median <- eval(mf[[match("median", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  lower <- eval(mf[[match("lower", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  upper <- eval(mf[[match("upper", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  id <- eval(mf[[match("id", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing.id & is.null(id))
    missing.id <- TRUE
  ##
  k.All <- if (!missing.TE)
             length(TE)
           else if (!missing.median)
             length(median)
           else if (!missing.lower)
             length(lower)
           else
             length(upper)
  ##
  if (!missing.TE)
    chknull(TE)
  else
    TE <- rep_len(NA, k.All)
  ##
  if (!missing.seTE)
    chknull(seTE)
  else
    seTE <- rep_len(NA, k.All)
  ##
  if (!missing.median)
    chknull(median)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  ## Catch 'studlab', 'byvar', 'subset', and 'exclude' from data:
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
  ## Catch 'pval', 'df', 'level.ci', 'q1', 'q3', 'min', 'max',
  ## 'approx.TE' and 'approx.seTE', from data:
  ##
  missing.pval <- missing(pval)
  pval <- eval(mf[[match("pval", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  missing.df <- missing(df)
  df <- eval(mf[[match("df", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing(level.ci))
    level.ci <- eval(mf[[match("level.ci", names(mf))]],
                     data, enclos = sys.frame(sys.parent()))
  ##
  missing.q1 <- missing(q1)
  q1 <- eval(mf[[match("q1", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  missing.q3 <- missing(q3)
  q3 <- eval(mf[[match("q3", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  missing.min <- missing(min)
  min <- eval(mf[[match("min", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  missing.max <- missing(max)
  max <- eval(mf[[match("max", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  missing.approx.TE <- missing(approx.TE)
  approx.TE <- eval(mf[[match("approx.TE", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
  ##
  missing.approx.seTE <- missing(approx.seTE)
  approx.seTE <- eval(mf[[match("approx.seTE", names(mf))]],
                      data, enclos = sys.frame(sys.parent()))
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  arg <- if (!missing.TE) "TE" else "median"
  chklength(seTE, k.All, arg)
  chklength(studlab, k.All, arg)
  ##
  if (by)
    chklength(byvar, k.All, arg)
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
  if (!is.null(n.e))
    chklength(n.e, k.All, arg)
  if (!is.null(n.c))
    chklength(n.c, k.All, arg)
  if (!missing.id)
    chklength(id, k.All, arg)
  if (by) {
    chklogical(print.byvar)
    chkchar(byseparator)
  }
  ##
  if (!missing.approx.TE) {
    if (length(approx.TE) == 1)
      rep_len(approx.TE, k.All)
    else
      chklength(approx.TE, k.All, arg)
    ##
    approx.TE <- setchar(approx.TE, c("", "ci", "iqr.range", "iqr", "range"))
  }
  ##
  if (!missing.approx.seTE) {
    if (length(approx.seTE) == 1)
      rep_len(approx.seTE, k.All)
    else
      chklength(approx.seTE, k.All, arg)
    ##
    approx.seTE <- setchar(approx.seTE,
                           c("", "pval", "ci", "iqr.range", "iqr", "range"))
  }
  ##
  if (!missing.pval)
    chklength(pval, k.All, arg)
  if (!missing.df)
    chklength(df, k.All, arg)
  if (!missing.lower)
    chklength(lower, k.All, arg)
  if (!missing.upper)
    chklength(upper, k.All, arg)
  if (length(level.ci) == 1)
    level.ci <- rep_len(level.ci, k.All)
  else
    chklength(level.ci, k.All, arg)
  if (!missing.median)
    chklength(median, k.All, arg)
  if (!missing.q1)
    chklength(q1, k.All, arg)
  if (!missing.q3)
    chklength(q3, k.All, arg)
  if (!missing.min)
    chklength(min, k.All, arg)
  if (!missing.max)
    chklength(max, k.All, arg)
  
  
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
    if (inherits(data, "meta")) {
      data <- data$data
      if (isCol(data, ".subset"))
        data <- data[data$.subset, ]
    }
    ##
    if (nulldata)
      data <- data.frame(.TE = TE)
    else
      data$.TE <- TE
    ##
    data$.seTE <- seTE
    data$.studlab <- studlab
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
    ##
    if (!missing.id)
      data$.id <- id
    ##
    if (!missing.pval)
      data$.pval <- pval
    if (!missing.df)
      data$.df <- df
    if (!missing.lower)
      data$.lower <- lower
    if (!missing.upper)
      data$.upper <- upper
    if (!missing.lower | !missing.upper)
      data$.level.ci <- level.ci
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
    if (!missing.approx.TE)
      data$.approx.TE <- approx.TE
    if (!missing.approx.seTE)
      data$.approx.seTE <- approx.seTE
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
    if (by)
      byvar <- byvar[subset]
    ##
    if (!is.null(n.e))
      n.e <- n.e[subset]
    if (!is.null(n.c))
      n.c <- n.c[subset]
    ##
    if (!missing.id)
      id <- id[subset]
    ##
    if (!missing.pval)
      pval <- pval[subset]
    if (!missing.df)
      df <- df[subset]
    if (!missing.lower)
      lower <- lower[subset]
    if (!missing.upper)
      upper <- upper[subset]
    level.ci <- level.ci[subset]
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
    if (!missing.approx.TE)
      approx.TE <- approx.TE[subset]
    if (!missing.approx.seTE)
      approx.seTE <- approx.seTE[subset]
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
    comb.fixed <- FALSE
    comb.random <- FALSE
    prediction <- FALSE
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(TE)
  chknumeric(seTE, 0)
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
  if (!multi.level & method.tau.ci == "PL") {
    if (method.tau == "DL")
      method.tau.ci <- "J"
    else
      method.tau.ci <- "QP"
  }
  ##
  if (multi.level) {
    if (!(method.tau %in% c("REML", "ML"))) {
      if (!missing.method.tau)
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
  ## (7) Calculate standard error from other information
  ##
  ##
  if (missing.approx.seTE) {
    approx.seTE <- rep_len("", length(TE))
    ##
    ## Use confidence limits
    ##
    sel.NA <- is.na(seTE)
    if (any(sel.NA) & !missing.lower & !missing.upper) {
      j <- sel.NA & !is.na(lower) & !is.na(upper)
      approx.seTE[j] <- "ci"
      if (missing.df)
        seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$seTE
      else
        seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j], df[j])$seTE
    }
    ##
    ## Use p-values
    ##
    sel.NA <- is.na(seTE)
    if (any(sel.NA) & !missing.pval) {
      j <- sel.NA & !is.na(TE) & !is.na(pval)
      approx.seTE[j] <- "pval"
      if (missing.df)
        seTE[j] <- seTE.pval(TE[j], pval[j])$seTE
      else
        seTE[j] <- seTE.pval(TE[j], pval[j], df[j])$seTE
    }
    ##
    ## Use IQR and range
    ##
    sel.NA <- is.na(seTE)
    if (any(sel.NA) & !missing.median &
        !missing.q1 & !missing.q3 &
        !missing.min & !missing.max &
        !(is.null(n.e) & is.null(n.c))) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3) &
        !is.na(min) & !is.na(max)
      approx.seTE[j] <- "iqr.range"
      if (is.null(n.c))
        seTE[j] <- mean.sd.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                     min[j], max[j],
                                     method.sd = method.sd)$se
      else if (is.null(n.e))
        seTE[j] <- mean.sd.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j],
                                     method.sd = method.sd)$se
      else
        seTE[j] <- mean.sd.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j],
                                     method.sd = method.sd)$se
    }
    ##
    ## Use IQR
    ##
    sel.NA <- is.na(seTE)
    if (any(sel.NA) & !missing.median &
        !missing.q1 & !missing.q3 &
        !(is.null(n.e) & is.null(n.c))) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3)
      approx.seTE[j] <- "iqr"
      if (is.null(n.c))
        seTE[j] <- mean.sd.iqr(n.e[j], median[j], q1[j], q3[j])$se
      else if (is.null(n.e))
        seTE[j] <- mean.sd.iqr(n.c[j], median[j], q1[j], q3[j])$se
      else
        seTE[j] <- mean.sd.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$se
    }
    ##
    ## Use range
    ##
    sel.NA <- is.na(seTE)
    if (any(sel.NA) & !missing.median &
        !missing.min & !missing.max &
        !(is.null(n.e) & is.null(n.c))) {
      j <- sel.NA & !is.na(median) & !is.na(min) & !is.na(max)
      approx.seTE[j] <- "range"
      if (is.null(n.c))
        seTE[j] <- mean.sd.range(n.e[j], median[j], min[j], max[j])$se
      else if (is.null(n.e))
        seTE[j] <- mean.sd.range(n.c[j], median[j], min[j], max[j])$se
      else
        seTE[j] <- mean.sd.range(n.e[j] + n.c[j], median[j],
                                 min[j], max[j])$se
    }
  }
  else {
    j <- 0
    for (i in approx.seTE) {
      j <- j + 1
      ##
      if (i == "ci") {
        if (missing.df)
          seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$seTE
        else
          seTE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j], df[j])$seTE
      }
      else if (i == "pval") {
        if (missing.df)
          seTE[j] <- seTE.pval(TE[j], pval[j])$seTE
        else
          seTE[j] <- seTE.pval(TE[j], pval[j], df[j])$seTE
      }
      else if (i == "iqr.range") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.seTE' = \"iqr\".",
               call. = FALSE)
        else if (is.null(n.c))
          seTE[j] <- mean.sd.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                       min[j], max[j],
                                       method.sd = method.sd)$se
        else if (is.null(n.e))
          seTE[j] <- mean.sd.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                       min[j], max[j],
                                       method.sd = method.sd)$se
        else
          seTE[j] <- mean.sd.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                       min[j], max[j],
                                       method.sd = method.sd)$se
      }
      else if (i == "iqr") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.seTE' = \"iqr\".",
               call. = FALSE)
        else if (is.null(n.c))
          seTE[j] <- mean.sd.iqr(n.e[j], median[j], q1[j], q3[j])$se
        else if (is.null(n.e))
          seTE[j] <- mean.sd.iqr(n.c[j], median[j], q1[j], q3[j])$se
        else
          seTE[j] <- mean.sd.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$se
      }
      else if (i == "range") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.seTE' = \"range\".",
               call. = FALSE)
        else if (is.null(n.c))
          seTE[j] <- mean.sd.range(n.e[j], median[j], min[j], max[j])$se
        else if (is.null(n.e))
          seTE[j] <- mean.sd.range(n.c[j], median[j], min[j], max[j])$se
        else
          seTE[j] <- mean.sd.range(n.e[j] + n.c[j], median[j],
                                   min[j], max[j])$se
      }
    }
  }
  
  
  ##
  ##
  ## (8) Calculate treatment estimate from other information
  ##
  ##
  if (missing.approx.TE) {
    approx.TE <- rep_len("", length(TE))
    ##
    ## Use confidence limits
    ##
    sel.NA <- is.na(TE)
    if (any(sel.NA) & !missing.lower & !missing.upper) {
      j <- sel.NA & !is.na(lower) & !is.na(upper)
      approx.TE[j] <- "ci"
      TE[j] <- TE.seTE.ci(lower[j], upper[j], level.ci[j])$TE
    }
    ##
    ## Use IQR and range
    ##
    sel.NA <- is.na(TE)
    if (any(sel.NA) & !missing.median &
        !missing.q1 & !missing.q3 &
        !missing.min & !missing.max &
        !(is.null(n.e) & is.null(n.c))) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3) &
        !is.na(min) & !is.na(max)
      approx.TE[j] <- "iqr.range"
      if (is.null(n.c))
        TE[j] <- mean.sd.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                   min[j], max[j], method.mean)$mean
      else if (is.null(n.e))
        TE[j] <- mean.sd.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j], method.mean)$mean
      else
        TE[j] <- mean.sd.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j], method.mean)$mean
    }
    ##
    ## Use IQR
    ##
    sel.NA <- is.na(TE)
    if (any(sel.NA) & !missing.median &
        !missing.q1 & !missing.q3 &
        !(is.null(n.e) & is.null(n.c))) {
      j <- sel.NA & !is.na(median) & !is.na(q1) & !is.na(q3)
      approx.TE[j] <- "iqr"
      if (is.null(n.c))
        TE[j] <- mean.sd.iqr(n.e[j], median[j], q1[j], q3[j], method.mean)$mean
      else if (is.null(n.e))
        TE[j] <- mean.sd.iqr(n.c[j], median[j], q1[j], q3[j], method.mean)$mean
      else
        TE[j] <- mean.sd.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                             method.mean)$mean
    }
    ##
    ## Use range
    ##
    sel.NA <- is.na(TE)
    if (any(sel.NA) & !missing.median &
        !missing.min & !missing.max &
        !(is.null(n.e) & is.null(n.c))) {
      j <- sel.NA & !is.na(median) & !is.na(min) & !is.na(max)
      approx.TE[j] <- "range"
      if (is.null(n.c))
        TE[j] <- mean.sd.range(n.e[j], median[j], min[j], max[j],
                               method.mean)$mean
      else if (is.null(n.e))
        TE[j] <- mean.sd.range(n.c[j], median[j], min[j], max[j],
                               method.mean)$mean
      else
        TE[j] <- mean.sd.range(n.e[j] + n.c[j], median[j], min[j], max[j],
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
          TE[j] <- mean.sd.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                     min[j], max[j], method.mean)$mean
      else if (is.null(n.e))
        TE[j] <- mean.sd.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j], method.mean)$mean
      else
        TE[j] <- mean.sd.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j], method.mean)$mean
      else if (i == "iqr") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.TE' = \"iqr\".",
               call. = FALSE)
        else if (is.null(n.c))
          TE[j] <- mean.sd.iqr(n.e[j], median[j], q1[j], q3[j], method.mean)$mean
        else if (is.null(n.e))
          TE[j] <- mean.sd.iqr(n.c[j], median[j], q1[j], q3[j], method.mean)$mean
        else
          TE[j] <- mean.sd.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                               method.mean)$mean
      }
      else if (i == "range") {
        cat(paste0("Use 'range' for study", j, "\n"))
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.TE' = \"range\".",
               call. = FALSE)
        else if (is.null(n.c))
          TE[j] <- mean.sd.range(n.e[j], median[j], min[j], max[j],
                                 method.mean)$mean
        else if (is.null(n.e))
          TE[j] <- mean.sd.range(n.c[j], median[j], min[j], max[j],
                                 method.mean)$mean
        else
          TE[j] <- mean.sd.range(n.e[j] + n.c[j], median[j], min[j], max[j],
                                 method.mean)$mean
      }
    }
  }
  ##
  if (keepdata) {
    if (!isCol(data, ".subset")) {
      data$.TE <- TE
      data$.seTE <- seTE
    }
    else {
      data$.TE[data$.subset] <- TE
      data$.seTE[data$.subset] <- seTE
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
    if (warn)
      warning("Zero values in seTE replaced by NAs.")
    seTE[!is.na(seTE) & seTE == 0] <- NA
  }
  ##
  tau2.calc <- NA
  
  
  ##
  ##
  ## (10) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(seTE[!exclude]))
  ##
  if (!missing.id) {
    id.incl <- id[!exclude]
    k.study <- length(unique(id.incl[!is.na(seTE[!exclude])]))
  }
  else
    k.study <- k
  ##
  seTE.random.hakn.orig <- NULL
  ##
  if (k == 0) {
    TE.fixed <- NA
    seTE.fixed <- NA
    statistic.fixed <- NA
    pval.fixed <- NA
    lower.fixed <- NA
    upper.fixed <- NA
    w.fixed <- rep(0, k.all)
    ##
    TE.random <- NA
    seTE.random <- NA
    statistic.random <- NA
    pval.random <- NA
    lower.random <- NA
    upper.random <- NA
    w.random <- rep(0, k.all)
    if (hakn)
      df.hakn <- NA
    ##
    hc <- list(tau2 = NA, se.tau2 = NA, lower.tau2 = NA, upper.tau2 = NA,
               tau = NA, lower.tau = NA, upper.tau = NA,
               method.tau.ci = "", sign.lower.tau = "", sign.upper.tau = "",
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
    ## At least two studies to perform Hartung-Knapp method
    if (k == 1 & hakn)
      hakn <- FALSE
    ##
    ## Estimate tau-squared
    ##
    hc <- hetcalc(TE[!exclude], seTE[!exclude],
                  method.tau, method.tau.ci,
                  TE.tau, level.comb,
                  control = control, id = id)
    ##
    if (by & tau.common) {
      ## Estimate common tau-squared across subgroups
      hcc <- hetcalc(TE[!exclude], seTE[!exclude],
                     method.tau, method.tau.ci,
                     TE.tau, level.comb,
                     byvar = byvar,
                     control = control, id = id)
    }
    ##
    ## Different calculations for three-level models
    ##
    if (!multi.level) {
      ##
      ## Classic meta-analysis
      ##
      if (is.null(tau.preset))
        tau2.calc <- if (is.na(hc$tau2)) 0 else hc$tau2
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
      ##
      ## Fixed effect estimate (Cooper & Hedges, 1994, p. 265-6)
      ##
      w.fixed <- 1 / seTE^2
      w.fixed[is.na(w.fixed) | is.na(TE) | exclude] <- 0
      ##
      TE.fixed   <- weighted.mean(TE, w.fixed, na.rm = TRUE)
      seTE.fixed <- sqrt(1 / sum(w.fixed, na.rm = TRUE))
      ##
      ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb,
                 null.effect = null.effect)
      statistic.fixed <- ci.f$statistic
      pval.fixed <- ci.f$p
      lower.fixed <- ci.f$lower
      upper.fixed <- ci.f$upper
      ##
      ## Random effects estimate (Cooper & Hedges, 1994, p. 265, 274-5
      ##
      w.random <- 1 / (seTE^2 + sum(tau2.calc, na.rm = TRUE))
      w.random[is.na(w.random) | is.na(TE) | exclude] <- 0
      ##
      TE.random   <- weighted.mean(TE, w.random, na.rm = TRUE)
      seTE.random <- sqrt(1 / sum(w.random, na.rm = TRUE))
    }
    else {
      ##
      ## Three-level meta-analysis
      ##
      id.TE <- seq_along(TE)
      ##
      ## No fixed effect method for three-level model
      ##
      if (comb.fixed) {
        if (!missing(comb.fixed))
          warning("Fixed effect model not calculated for three-level model.",
                  call. = FALSE)
        comb.fixed <- FALSE
      }
      ##
      w.fixed <- rep_len(NA, length(seTE))
      ##
      TE.fixed <- seTE.fixed <- lower.fixed <- upper.fixed <-
        statistic.fixed <- pval.fixed <- NA
      ##
      ## Conduct three-level meta-analysis
      ##
      if (by & !tau.common) {
        if (!missing(tau.common))
          warning("Argument 'tau.common' set to TRUE for three-level model.",
                  call. = FALSE)
        tau.common <- TRUE
      }
      ##
      m4 <- rma.mv(TE, seTE^2,
                   method = method.tau,
                   random = ~ 1 | id / id.TE,
                   control = control)
      ##
      w.random <- rep_len(NA, length(w.fixed))
      ##
      TE.random <- m4$b
      seTE.random <- m4$se
      ##
      tau2.calc <- sum(m4$tau2)
    }
    ##
    ## Hartung-Knapp adjustment
    ##
    if (hakn) {
      alpha <- (1 - level.comb) / 2
      df.hakn <- k - 1
      q <- 1 / (k - 1) * sum(w.random * (TE - TE.random)^2, na.rm = TRUE)
      ##
      seTE.random <- sqrt(q / sum(w.random))
      ##
      if (adhoc.hakn == "se") {
        ##
        ## Variance correction if SE_HK < SE_notHK (Knapp and Hartung, 2003)
        ##
        if (q < 1) {
          seTE.random.hakn.orig <- seTE.random
          seTE.random <- sqrt(1 /  sum(w.random))
        }
      }
      else if (adhoc.hakn == "ci") {
        ##
        ## Use wider confidence interval, i.e., confidence interval
        ## from classic random effects meta-analysis if HK CI is
        ## smaller
        ## (Wiksten et al., 2016; Jackson et al., 2017, hybrid 2)
        ##
        if (q > qnorm(alpha) / qt(alpha, df = df.hakn)) {
          seTE.random.hakn.orig <- seTE.random
          seTE.random <- sqrt(1 /  sum(w.random))
          df.hakn <- NULL
        }
      }
      else if (adhoc.hakn == "iqwig6") {
        ##
        ## Variance correction if CI_HK < CI_DL (IQWiG, 2020)
        ##
        ci.hk <- ci(TE.random, seTE.random, level = level.comb, df = df.hakn)
        ##
        m.dl <- metagen(TE, seTE, method.tau = "DL", method.tau.ci = "",
                        hakn = FALSE)
        ci.dl <- ci(m.dl$TE.random, m.dl$seTE.random, level = level.comb)
        ##
        width.hk <- ci.hk$upper - ci.hk$lower
        width.dl <- ci.dl$upper - ci.dl$lower
        ##
        if (width.hk < width.dl) {
          seTE.random.hakn.orig <- seTE.random
          seTE.random <- sqrt(max(q, 1) /  sum(w.random))
        }
      }
      ##
      ci.r <- ci(TE.random, seTE.random, level = level.comb, df = df.hakn,
                 null.effect = null.effect)
    }
    else
      ci.r <- ci(TE.random, seTE.random, level = level.comb,
                 null.effect = null.effect)
    ##
    statistic.random <- ci.r$statistic
    pval.random <- ci.r$p
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }
  ##
  ## Individual study results
  ##
  ci.study <- ci(TE, seTE, level = level,
                 df = if (method.ci == "t") df else NULL,
                 null.effect = null.effect)
  ##
  ## Prediction interval
  ##
  if (k >= 3) {
    seTE.predict <- sqrt(seTE.random^2 + tau2.calc)
    ci.p <- ci(TE.random, seTE.predict, level.predict, k - 2)
    p.lower <- ci.p$lower
    p.upper <- ci.p$upper
  }
  else {
    seTE.predict <- NA
    p.lower <- NA
    p.upper <- NA
  }
  ##
  ## Calculate Rb (but not for three-level model)
  ##
  if (length(tau2.calc) == 1)
    Rbres <- Rb(seTE[!is.na(seTE)], seTE.random,
                tau2.calc, hc$Q, hc$df.Q, level.comb)
  else
    Rbres <- list(TE = NA, lower = NA, upper = NA)
  
  
  ##
  ##
  ## (11) Generate R object
  ##
  ##
  if (missing(detail.tau) && k != k.study)
    detail.tau <- c("level 2", "level 1")
  ##
  res <- list(studlab = studlab,
              ##
              TE = TE, seTE = seTE,
              lower = ci.study$lower, upper = ci.study$upper,
              zval = ci.study$statistic,
              statistic = ci.study$statistic,
              pval = ci.study$p,
              df = if (method.ci == "t") df else rep_len(NA, length(TE)),
              w.fixed = w.fixed, w.random = w.random,
              id = id,
              ##
              TE.fixed = TE.fixed, seTE.fixed = seTE.fixed,
              lower.fixed = lower.fixed, upper.fixed = upper.fixed,
              zval.fixed = statistic.fixed,
              statistic.fixed = statistic.fixed, pval.fixed = pval.fixed,
              ##
              TE.random = TE.random, seTE.random = seTE.random,
              lower.random = lower.random, upper.random = upper.random,
              zval.random = statistic.random,
              statistic.random = statistic.random, pval.random = pval.random,
              ##
              null.effect = null.effect,
              ##
              seTE.predict = seTE.predict,
              lower.predict = p.lower, upper.predict = p.upper,
              level.predict = level.predict,
              ##
              k = k, k.study = k.study,
              Q = hc$Q, df.Q = hc$df.Q, pval.Q = hc$pval.Q,
              tau2 = hc$tau2, se.tau2 = hc$se.tau2,
              lower.tau2 = hc$lower.tau2, upper.tau2 = hc$upper.tau2,
              tau = hc$tau, lower.tau = hc$lower.tau, upper.tau = hc$upper.tau,
              method.tau.ci = hc$method.tau.ci,
              sign.lower.tau = hc$sign.lower.tau,
              sign.upper.tau = hc$sign.upper.tau,
              ##
              H = hc$H, lower.H = hc$lower.H, upper.H = hc$upper.H,
              ##
              I2 = hc$I2, lower.I2 = hc$lower.I2, upper.I2 = hc$upper.I2,
              ##
              Rb = Rbres$TE, lower.Rb = Rbres$lower, upper.Rb = Rbres$upper,
              ##
              method.mean = method.mean,
              approx.TE = approx.TE,
              approx.seTE = approx.seTE,
              ##
              sm = sm, method = "Inverse",
              level = level,
              level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              overall = overall,
              overall.hetstat = overall.hetstat,
              hakn = hakn, adhoc.hakn = adhoc.hakn,
              df.hakn = if (hakn) df.hakn else NA,
              seTE.random.hakn.orig = seTE.random.hakn.orig,
              method.tau = method.tau, method.tau.ci = hc$method.tau.ci,
              tau.preset = tau.preset,
              TE.tau =
                if (!missing(TE.tau) & method.tau == "DL") TE.tau else NULL,
              tau.common = tau.common,
              detail.tau = detail.tau,
              prediction = prediction,
              method.bias = method.bias,
              n.e = n.e,
              n.c = n.c,
              ##
              text.fixed = text.fixed, text.random = text.random,
              text.predict = text.predict,
              text.w.fixed = text.w.fixed, text.w.random = text.w.random,
              ##
              title = title, complab = complab, outclab = outclab,
              label.e = label.e,
              label.c = label.c,
              label.left = label.left,
              label.right = label.right,
              data = if (keepdata) data else NULL,
              subset = if (keepdata) subset else NULL,
              exclude = if (!missing.exclude) exclude else NULL,
              print.byvar = print.byvar,
              byseparator = byseparator,
              warn = warn,
              call = match.call(),
              backtransf = backtransf,
              pscale = pscale,
              irscale = irscale, irunit = irunit,
              control = control,
              version = packageDescription("meta")$Version)
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (by) {
    res$byvar <- byvar
    res$bylab <- bylab
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
      ##
      res$tau2.resid <- hcc$tau2.resid
      res$lower.tau2.resid <- hcc$lower.tau2.resid
      res$upper.tau2.resid <- hcc$upper.tau2.resid
      ##
      res$tau.resid <- hcc$tau.resid
      res$lower.tau.resid <- hcc$lower.tau.resid
      res$upper.tau.resid <- hcc$upper.tau.resid
      res$sign.lower.tau.resid <- hcc$sign.lower.tau
      res$sign.upper.tau.resid <- hcc$sign.upper.tau
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
    res$event.w   <- NULL
    res$n.w       <- NULL
    res$time.e.w  <- NULL
    res$time.c.w  <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
