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
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param n.e Number of observations in experimental group.
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
#' @param approx.TE Approximation method to estimate treatment
#'   estimate (see Details).
#' @param approx.seTE Approximation method to estimate standard error
#'   (see Details).
#' @param hakn A logical indicating whether method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2}. Either
#'   \code{"DL"}, \code{"PM"}, \code{"REML"}, \code{"ML"},
#'   \code{"HS"}, \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be
#'   abbreviated.
#' @param tau.preset Prespecified value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param method.bias A character string indicating which test is to
#'   be used.  Either \code{"rank"}, \code{"linreg"}, or \code{"mm"},
#'   can be abbreviated.  See function \code{\link{metabias}}
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
#'   estimate the between-study variance tau^2. This argument is
#'   passed on to \code{\link[metafor]{rma.uni}}.
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
#' For several arguments defaults settings are utilised (see
#' assignments with \code{\link{gs}} under \bold{Usage}). These
#' defaults can be changed using \code{\link{settings.meta}}.
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
#' effect measures like the odds ratio or hazard ratio). For median,
#' interquartile range and range, equation (10) in Wan et al. (2014)
#' is used to approximate the treatment effect (i.e.,
#' mean). Similarly, equations (14) and (2) in Wan et al. (2014) are
#' used if median and interquartile range or range, respectively, are
#' provided.
#'
#' By default, missing treatment estimates are replaced successively
#' using these method, e.g., confidence limits are utilised before
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
#' used to provide the level of the confidence interval. For median,
#' interquartile range and range, depending on the sample size,
#' equation (12) or (13) in Wan et al. (2014) is used to approximate
#' the standard error. Similarly, equations (15) / (16) and (7) / (9)
#' in Wan et al. (2014) are used if median and interquartile range or
#' range, respectively, are provided. The sample size of individual
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
#' \subsection{Estimation of between-study variance}{
#' 
#' The following methods are available to estimate the between-study
#' variance \eqn{\tau^2}.
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{Method} \cr 
#' \code{method.tau = "DL"}\tab DerSimonian-Laird estimator (DerSimonian and Laird, 1986) \cr
#' \code{method.tau = "PM"}\tab Paule-Mandel estimator (Paule and Mandel, 1982) \cr
#' \code{method.tau = "REML"}\tab Restricted maximum-likelihood estimator (Viechtbauer, 2005) \cr
#' \code{method.tau = "ML"}\tab Maximum-likelihood estimator (Viechtbauer, 2005) \cr
#' \code{method.tau = "HS"}\tab Hunter-Schmidt estimator (Hunter and Schmidt, 2015) \cr
#' \code{method.tau = "SJ"}\tab Sidik-Jonkman estimator (Sidik and Jonkman, 2005) \cr
#' \code{method.tau = "HE"}\tab Hedges estimator (Hedges and Olkin, 1985) \cr
#' \code{method.tau = "EB"}\tab Empirical Bayes estimator (Morris, 1983)
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
#' 
#' The DerSimonian-Laird and Paule-Mandel estimators are implemented
#' in R package \pkg{meta}. The other estimators are available if R
#' package \pkg{metafor} (Viechtbauer 2010) is installed by internally
#' calling R function \code{\link[metafor]{rma.uni}}.
#' }
#' 
#' \subsection{Hartung-Knapp method}{
#' 
#' Hartung and Knapp (2001a,b) proposed an alternative method for
#' random effects meta-analysis based on a refined variance estimator
#' for the treatment estimate. Simulation studies (Hartung and Knapp,
#' 2001a,b; IntHout et al., 2014; Langen et al., 2018) show improved
#' coverage probabilities compared to the classic random effects
#' method. However, in rare settings with very homogeneous treatment
#' estimates, the Hartung-Knapp method can be anti-conservative
#' (Wiksten et al., 2016). The Hartung-Knapp method is used if
#' argument \code{hakn = TRUE}.
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
#' @return
#' An object of class \code{c("metagen", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{TE, seTE, studlab, exclude, n.e, n.c}{As defined above.}
#' \item{sm, level, level.comb,}{As defined above.}
#' \item{comb.fixed, comb.random,}{As defined above.}
#' \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{As
#'   defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{label.e, label.c, label.left, label.right,}{As defined
#'   above.}
#' \item{byvar, bylab, print.byvar, byseparator, warn}{As defined
#'   above.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{zval, pval}{z-value and p-value for test of treatment effect
#'   for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
#'   standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{zval.fixed, pval.fixed}{z-value and p-value for test of
#'   overall treatment effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall treatment effect
#'   and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{zval.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall treatment effect
#'   (random effects model).}
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
#' \item{tau}{Square-root of between-study variance.}
#' \item{se.tau}{Standard error of square-root of between-study
#'   variance.}
#' \item{C}{Scaling factor utilised internally to
#'   calculate common tau-squared across subgroups.}
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
#' \item{zval.fixed.w, pval.fixed.w}{z-value and p-value for test of
#'   treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
#'   standard error in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of treatment effect in subgroups
#'   (random effects model) - if \code{byvar} is not missing.}
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
#'   - if \code{byvar} is not missing.}  \item{C.w}{Scaling factor
#'   utilised internally to calculate common tau-squared across
#'   subgroups - if \code{byvar} is not missing.}
#'   \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limti for
#'   heterogeneity statistic H within subgroups - if \code{byvar} is
#'   not missing.}  \item{I2.w}{Heterogeneity statistic I2 within
#'   subgroups - if \code{byvar} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limti for
#'   heterogeneity statistic I2 within subgroups - if \code{byvar} is
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
#' Langan D, Higgins JPT, Jackson D, Bowden J, Veroniki AA,
#' Kontopantelis E, et al. (2018):
#' A comparison of heterogeneity variance estimators in simulated
#' random-effects meta-analyses.
#' \emph{Research Synthesis Methods}
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
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.3.
#' Copenhagen: The Nordic Cochrane Centre, The Cochrane Collaboration, 2014
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
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
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
#' data(Fleiss93)
#' m1 <- metabin(event.e, n.e, event.c, n.c,
#'               data = Fleiss93, sm = "RR", method = "I")
#' m1
#' # Identical results by using the generic inverse variance method
#' metagen(m1$TE, m1$seTE, sm = "RR")
#' #
#' forest(metagen(m1$TE, m1$seTE, sm = "RR"))
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
#' @export metagen


metagen <- function(TE, seTE, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    ##
                    sm = "",
                    ##
                    level = gs("level"), level.comb = gs("level.comb"),
                    comb.fixed = gs("comb.fixed"),
                    comb.random = gs("comb.random"),
                    ##
                    hakn = gs("hakn"),
                    method.tau = gs("method.tau"),
                    tau.preset = NULL, TE.tau = NULL,
                    tau.common = gs("tau.common"),
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
                    ##
                    approx.TE, approx.seTE,
                    ##
                    backtransf = gs("backtransf"),
                    pscale = 1,
                    irscale = 1, irunit = "person-years",
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
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(hakn)
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  chknumeric(null.effect, single = TRUE)
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  if (!is.prop(sm))
    pscale <- 1
  chknumeric(pscale, single = TRUE)
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate(sm))
    irscale <- 1
  chknumeric(irscale, single = TRUE)
  if (!backtransf & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks
  ##
  fun <- "metagen"
  chklogical(warn)
  ##
  if (tau.common & method.tau == "PM") {
    warning("Argument 'method.tau' set to \"DL\" as argument tau.common = TRUE.")
    method.tau <- "DL"
  }
  ##
  if (!is.null(tau.preset) & method.tau == "PM") {
    warning("Argument 'tau.preset' not considered as",
            "argument method.tau = \"PM\".")
    tau.preset <- NULL
  }
  
  
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
  ## Catch 'TE', 'seTE', 'median', 'lower', 'upper', 'n.e', and 'n.c'
  ## from data:
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
      if (!is.null(data$.subset))
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
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
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
        seTE[j] <- TE.seTE.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                     min[j], max[j])$seTE
      else if (is.null(n.e))
        seTE[j] <- TE.seTE.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j])$seTE
      else
        seTE[j] <- TE.seTE.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                     min[j], max[j])$seTE
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
        seTE[j] <- TE.seTE.iqr(n.e[j], median[j], q1[j], q3[j])$seTE
      else if (is.null(n.e))
        seTE[j] <- TE.seTE.iqr(n.c[j], median[j], q1[j], q3[j])$seTE
      else
        seTE[j] <- TE.seTE.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$seTE
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
        seTE[j] <- TE.seTE.range(n.e[j], median[j], min[j], max[j])$seTE
      else if (is.null(n.e))
        seTE[j] <- TE.seTE.range(n.c[j], median[j], min[j], max[j])$seTE
      else
        seTE[j] <- TE.seTE.range(n.e[j] + n.c[j], median[j],
                                 min[j], max[j])$seTE
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
          seTE[j] <- TE.seTE.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                       min[j], max[j])$seTE
        else if (is.null(n.e))
          seTE[j] <- TE.seTE.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                       min[j], max[j])$seTE
        else
          seTE[j] <- TE.seTE.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                       min[j], max[j])$seTE
      }
      else if (i == "iqr") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.seTE' = \"iqr\".",
               call. = FALSE)
        else if (is.null(n.c))
          seTE[j] <- TE.seTE.iqr(n.e[j], median[j], q1[j], q3[j])$seTE
        else if (is.null(n.e))
          seTE[j] <- TE.seTE.iqr(n.c[j], median[j], q1[j], q3[j])$seTE
        else
          seTE[j] <- TE.seTE.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$seTE
      }
      else if (i == "range") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.seTE' = \"range\".",
               call. = FALSE)
        else if (is.null(n.c))
          seTE[j] <- TE.seTE.range(n.e[j], median[j], min[j], max[j])$seTE
        else if (is.null(n.e))
          seTE[j] <- TE.seTE.range(n.c[j], median[j], min[j], max[j])$seTE
        else
          seTE[j] <- TE.seTE.range(n.e[j] + n.c[j], median[j],
                                   min[j], max[j])$seTE
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
        TE[j] <- TE.seTE.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                   min[j], max[j])$TE
      else if (is.null(n.e))
        TE[j] <- TE.seTE.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j])$TE
      else
        TE[j] <- TE.seTE.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j])$TE
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
        TE[j] <- TE.seTE.iqr(n.e[j], median[j], q1[j], q3[j])$TE
      else if (is.null(n.e))
        TE[j] <- TE.seTE.iqr(n.c[j], median[j], q1[j], q3[j])$TE
      else
        TE[j] <- TE.seTE.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$TE
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
        TE[j] <- TE.seTE.range(n.e[j], median[j], min[j], max[j])$TE
      else if (is.null(n.e))
        TE[j] <- TE.seTE.range(n.c[j], median[j], min[j], max[j])$TE
      else
        TE[j] <- TE.seTE.range(n.e[j] + n.c[j], median[j], min[j], max[j])$TE
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
          TE[j] <- TE.seTE.iqr.range(n.e[j], median[j], q1[j], q3[j],
                                     min[j], max[j])$TE
      else if (is.null(n.e))
        TE[j] <- TE.seTE.iqr.range(n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j])$TE
      else
        TE[j] <- TE.seTE.iqr.range(n.e[j] + n.c[j], median[j], q1[j], q3[j],
                                   min[j], max[j])$TE
      else if (i == "iqr") {
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.TE' = \"iqr\".",
               call. = FALSE)
        else if (is.null(n.c))
          TE[j] <- TE.seTE.iqr(n.e[j], median[j], q1[j], q3[j])$TE
        else if (is.null(n.e))
          TE[j] <- TE.seTE.iqr(n.c[j], median[j], q1[j], q3[j])$TE
        else
          TE[j] <- TE.seTE.iqr(n.e[j] + n.c[j], median[j], q1[j], q3[j])$TE
      }
      else if (i == "range") {
        cat(paste0("Use 'range' for study", j, "\n"))
        if (is.null(n.e) & is.null(n.c))
          stop("Sample size needed if argument 'approx.TE' = \"range\".",
               call. = FALSE)
        else if (is.null(n.c))
          TE[j] <- TE.seTE.range(n.e[j], median[j], min[j], max[j])$TE
        else if (is.null(n.e))
          TE[j] <- TE.seTE.range(n.c[j], median[j], min[j], max[j])$TE
        else
          TE[j] <- TE.seTE.range(n.e[j] + n.c[j], median[j], min[j], max[j])$TE
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
    if (warn)
      warning("Zero values in seTE replaced by NAs.")
    seTE[!is.na(seTE) & seTE == 0] <- NA
  }
  
  
  ##
  ##
  ## (10) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(seTE[!exclude]))
  ##
  tau2 <- NA
  ##
  if (k == 0) {
    TE.fixed <- NA
    seTE.fixed <- NA
    zval.fixed <- NA
    pval.fixed <- NA
    lower.fixed <- NA
    upper.fixed <- NA
    w.fixed <- rep(0, k.all)
    ##
    TE.random <- NA
    seTE.random <- NA
    zval.random <- NA
    pval.random <- NA
    lower.random <- NA
    upper.random <- NA
    w.random <- rep(0, k.all)
    ##
    Q <- NA
    df.Q <- NA
    se.tau2 <- NA
    ##
    if (hakn)
      df.hakn <- NA
    ##
    Cval <- NA
    ##
    hc <- list(H = NA, lower.H = NA, upper.H = NA,
               I2 = NA, lower.I2 = NA, upper.I2 = NA)
  }
  else {
    ## At least two studies to perform Hartung-Knapp method
    if (k == 1 & hakn)
      hakn <- FALSE
    ## Estimate tau-squared
    hc <- hetcalc(TE[!exclude], seTE[!exclude], method.tau, TE.tau,
                  level.comb, control = control)
    ##
    if (by & tau.common) {
      ## Estimate common tau-squared across subgroups
      hcc <- hetcalc(TE[!exclude], seTE[!exclude], method.tau, TE.tau,
                     level.comb, byvar, control)
    }
    ##
    if (is.null(tau.preset)) {
      if (k > 1)
        tau2 <- hc$tau^2
      tau2.calc <- if (is.na(tau2)) 0 else tau2
      se.tau2 <- hc$se.tau2
    }
    else {
      tau2 <- tau2.calc <- tau.preset^2
      se.tau2 <- NULL
    }
    ##
    Q    <- hc$Q
    df.Q <- hc$df.Q
    Cval <- hc$Cval
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
    zval.fixed <- ci.f$z
    pval.fixed <- ci.f$p
    lower.fixed <- ci.f$lower
    upper.fixed <- ci.f$upper
    ##
    ## Random effects estimate
    ##
    if (method.tau == "PM") {
      if (Q < k - 1) {
        TE.random <- TE.fixed
        seTE.random <- seTE.fixed
        w.random <- w.fixed
        ##
        if (k > 1)
          tau2 <- 0
        tau2.calc <- 0
      }
      else {
        pm <- paulemandel(TE[!exclude], seTE[!exclude])
        TE.random <- pm$TE.random
        seTE.random <- pm$seTE.random
        w.random <- pm$w.random
        ##
        tau2 <- tau2.calc <- pm$tau^2
      }
    }
    else {
      ##
      ## Cooper & Hedges (1994), p. 265, 274-5
      ##
      w.random <- 1 / (seTE^2 + tau2.calc)
      w.random[is.na(w.random) | is.na(TE) | exclude] <- 0
      ##
      TE.random   <- weighted.mean(TE, w.random, na.rm = TRUE)
      seTE.random <- sqrt(1 / sum(w.random, na.rm = TRUE))
    }
    ##
    ## Hartung-Knapp adjustment
    ##
    if (hakn) {
      seTE.random <- sqrt(1 / (k - 1) * sum(w.random * (TE - TE.random)^2 /
                                              sum(w.random), na.rm = TRUE))
      df.hakn <- k - 1
      ci.r <- ci(TE.random, seTE.random, level = level.comb, df = df.hakn,
                 null.effect = null.effect)
    }
    else
      ci.r <- ci(TE.random, seTE.random, level = level.comb,
                 null.effect = null.effect)
    ##
    zval.random <- ci.r$z
    pval.random <- ci.r$p
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }
  ##
  ## Individual study results
  ##
  ci.study <- ci(TE, seTE, level = level, null.effect = null.effect)
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
  ## Calculate Rb
  ##
  Rbres <- Rb(seTE[!is.na(seTE)], seTE.random, tau2.calc, Q, df.Q, level.comb)
  
  
  ##
  ##
  ## (11) Generate R object
  ##
  ##
  res <- list(studlab = studlab,
              ##
              TE = TE, seTE = seTE,
              lower = ci.study$lower, upper = ci.study$upper,
              zval = ci.study$z, pval = ci.study$p,
              w.fixed = w.fixed, w.random = w.random,
              ##
              TE.fixed = TE.fixed, seTE.fixed = seTE.fixed,
              lower.fixed = lower.fixed, upper.fixed = upper.fixed,
              zval.fixed = zval.fixed, pval.fixed = pval.fixed,
              ##
              TE.random = TE.random, seTE.random = seTE.random,
              lower.random = lower.random, upper.random = upper.random,
              zval.random = zval.random, pval.random = pval.random,
              ##
              null.effect = null.effect,
              ##
              seTE.predict = seTE.predict,
              lower.predict = p.lower, upper.predict = p.upper,
              level.predict = level.predict,
              ##
              k = k, Q = Q, df.Q = df.Q, pval.Q = pvalQ(Q, df.Q),
              tau = sqrt(tau2), se.tau2 = se.tau2,
              C = Cval,
              ##
              H = hc$H,
              lower.H = hc$lower.H,
              upper.H = hc$upper.H,
              ##
              I2 = hc$I2,
              lower.I2 = hc$lower.I2,
              upper.I2 = hc$upper.I2,
              ##
              Rb = Rbres$TE,
              lower.Rb = Rbres$lower,
              upper.Rb = Rbres$upper,
              ##
              approx.TE = approx.TE,
              approx.seTE = approx.seTE,
              ##
              sm = sm, method = "Inverse",
              level = level,
              level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              hakn = hakn,
              df.hakn = if (hakn) df.hakn else NULL,
              method.tau = method.tau,
              tau.preset = tau.preset,
              TE.tau = if (!missing(TE.tau) & method.tau == "DL") TE.tau else NULL,
              tau.common = tau.common,
              prediction = prediction,
              method.bias = method.bias,
              n.e = n.e,
              n.c = n.c,
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
    if (!tau.common) {
      res <- c(res, subgroup(res))
      res$tau.resid <- NA
    }
    else if (!is.null(tau.preset)) {
      res <- c(res, subgroup(res, tau.preset))
      res$tau.resid <- NA
    }
    else {
      res <- c(res, subgroup(res, hcc$tau))
      res$Q.w.random <- hcc$Q
      res$df.Q.w.random <- hcc$df.Q
      res$tau.resid <- hcc$tau
    }
    ##
    if (!tau.common || method.tau == "DL") {
      ci.H.resid <- calcH(res$Q.w.fixed, res$df.Q.w, level.comb)
      ##
      res$H.resid <- ci.H.resid$TE
      res$lower.H.resid <- ci.H.resid$lower
      res$upper.H.resid <- ci.H.resid$upper
    }
    else {
      res$H.resid <- hcc$H.resid
      res$lower.H.resid <- hcc$lower.H.resid
      res$upper.H.resid <- hcc$upper.H.resid
    }
    ##
    if (!tau.common || method.tau == "DL") {
      ci.I2.resid <- isquared(res$Q.w.fixed, res$df.Q.w, level.comb)
      ##
      res$I2.resid <- ci.I2.resid$TE
      res$lower.I2.resid <- ci.I2.resid$lower
      res$upper.I2.resid <- ci.I2.resid$upper
    }
    else {
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
