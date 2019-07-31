#' Meta-analysis of single proportions
#' 
#' @description
#' Calculation of an overall proportion from studies reporting a
#' single proportion. Inverse variance method and generalised linear
#' mixed model (GLMM) are available for pooling. For GLMMs, the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} (Viechtbauer 2010) is called internally.
#' 
#' @param event Number of events.
#' @param n Number of observations.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., event and n.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"} and
#'   \code{"GLMM"}, can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"PFT"}, \code{"PAS"}, \code{"PRAW"}, \code{"PLN"}, or
#'   \code{"PLOGIT"}) is to be used for pooling of studies, see
#'   Details.
#' @param incr A numeric which is added to event number and sample
#'   size of studies with zero or all events, i.e., studies with an
#'   event probability of either 0 or 1.
#' @param allincr A logical indicating if \code{incr} is considered
#'   for all studies if at least one study has either zero or all
#'   events. If FALSE (default), \code{incr} is considered only in
#'   studies with zero or all events.
#' @param addincr A logical indicating if \code{incr} is used for all
#'   studies irrespective of number of events.
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
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2}, see
#'   Details.
#' @param tau.preset Prespecified value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param method.bias A character string indicating which test is to
#'   be used.  Either \code{"rank"}, \code{"linreg"}, or \code{"mm"},
#'   can be abbreviated.  See function \code{\link{metabias}}.
#' @param backtransf A logical indicating whether results for
#'   transformed proportions (argument \code{sm != "PRAW"}) should be
#'   back transformed in printouts and plots. If TRUE (default),
#'   results will be presented as proportions; otherwise transformed
#'   proportions will be shown. See Details for presentation of
#'   confidence intervals.
#' @param pscale A numeric defining a scaling factor for printing of
#'   single event probabilities.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{event}).
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether the addition of
#'   \code{incr} to studies with zero or all events should result in a
#'   warning.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance tau^2. This argument is
#'   passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}, respectively.
#' @param \dots Additional arguments passed on to
#'   \code{\link[metafor]{rma.glmm}} function.
#' 
#' @details
#' This function provides methods for fixed effect and random effects
#' meta-analysis of single proportions to calculate an overall
#' proportion. Note, you should use R function \code{\link{metabin}}
#' to compare proportions of pairwise comparisons instead of using
#' \code{metaprop} for each treatment arm separately which will break
#' randomisation in randomised controlled trials.
#' 
#' The following transformations of proportions are
#' implemented to calculate an overall proportion:
#' 
#' \itemize{
#' \item Logit transformation (\code{sm = "PLOGIT"}, default)
#' \item Log transformation (\code{sm = "PLN"})
#' \item Freeman-Tukey Double arcsine transformation (\code{sm = "PFT"})
#' \item Arcsine transformation (\code{sm = "PAS"})
#' \item Raw, i.e. untransformed, proportions (\code{sm = "PRAW"})
#' }
#'
#' Classic meta-analysis (Borenstein et al., 2010) utilises the
#' transformed proportions and corresponding standard errors in the
#' generic inverse variance method. A distinctive and frequently
#' overlooked advantage of binary data is that individual patient data
#' can be extracted. Accordingly, a generalised linear mixed model
#' (GLMM) - more specific, a random intercept logistic regression
#' model - can be utilised for the meta-analysis of proportions
#' (Stijnen et al., 2010). This method - implicitly using the logit
#' transformation - is available (argument \code{method = "GLMM"}) by
#' calling the \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} internally.
#'
#' For the logit transformation, a random intercept logistic
#' regression model is used by default, i.e., argument \code{method =
#' "GLMM"}. The classic meta-analysis model based on the inverse
#' variance method can be used instead by setting argument
#' \code{method} equal to \code{"Inverse"}.
#'
#' Contradictory recommendations on the use of transformations of
#' proportions have been published in the literature. For example,
#' Barendregt et al. (2013) recommend the use of the Freeman-Tukey
#' double arcsine transformation instead of the logit transformation
#' whereas Warton & Hui (2011) strongly advise to use generalised
#' linear mixed models with the logit transformation instead of the
#' arcsine transformation. Schwarzer et al. (2019) describe seriously
#' misleading results in a meta-analysis with very different sample
#' sizes due to problems with the back-transformation of the
#' Freeman-Tukey transformation which requires a single sample
#' size. Accordingly, Schwarzer et al. (2019) also recommend to use
#' GLMMs for the meta-analysis of single proportions, however, admit
#' that individual study weights are not available with this
#' method. Meta-analysts which require individual study weights should
#' consider the arcsine or logit transformation.
#'
#' In order to prevent misleading conclusions for the Freeman-Tukey
#' double arcsine transformation, sensitivity analyses using other
#' transformations or using a range of sample sizes should be
#' conducted (Schwarzer et al., 2019).
#' 
#' Various methods are available to calculate confidence intervals for
#' individual study results (see Agresti & Coull 1998 and Newcombe
#' 1988):
#' \itemize{
#' \item Clopper-Pearson interval also called 'exact' binomial
#'   interval (\code{method.ci = "CP"}, default)
#' \item Wilson Score interval (\code{method.ci = "WS"})
#' \item Wilson Score interval with continuity correction
#'   (\code{method.ci = "WSCC"})
#' \item Agresti-Coull interval (\code{method.ci = "AC"})
#' \item Simple approximation interval (\code{method.ci = "SA"})
#' \item Simple approximation interval with continuity correction
#'   (\code{method.ci = "SACC"})
#' \item Normal approximation interval based on summary measure,
#'   i.e. defined by argument \code{sm} (\code{method.ci = "NAsm"})
#' }
#' 
#' Note, with exception of the normal approximation based on the
#' summary measure, i.e. \code{method.ci = "NAsm"}, the same
#' confidence interval is calculated for individual studies for any
#' summary measure (argument \code{sm}) as only number of events and
#' observations are used in the calculation disregarding the chosen
#' summary measure. Results will be presented for transformed
#' proportions if argument \code{backtransf = FALSE} in the
#' \code{\link{print.meta}}, \code{\link{print.summary.meta}}, or
#' \code{\link{forest.meta}} function. In this case, argument
#' \code{method.ci = "NAsm"} is used, i.e. confidence intervals based
#' on the normal approximation based on the summary measure.
#' 
#' Argument \code{pscale} can be used to rescale proportions, e.g.
#' \code{pscale = 1000} means that proportions are expressed as events
#' per 1000 observations. This is useful in situations with (very) low
#' event probabilities.
#' 
#' For several arguments defaults settings are utilised (assignments
#' using \code{\link{gs}} function). These defaults can be changed
#' using the \code{\link{settings.meta}} function.
#' 
#' Internally, both fixed effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{comb.random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{comb.fixed} and \code{comb.random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{comb.random = FALSE}.
#' 
#' If the summary measure is equal to "PRAW", "PLN", or "PLOGIT", a
#' continuity correction is applied if any study has either zero or
#' all events, i.e., an event probability of either 0 or 1. By
#' default, 0.5 is used as continuity correction (argument
#' \code{incr}). This continuity correction is used both to calculate
#' individual study results with confidence limits and to conduct
#' meta-analysis based on the inverse variance method. For GLMMs no
#' continuity correction is used.
#' 
#' Argument \code{byvar} can be used to conduct subgroup analysis for
#' all methods but GLMMs. Instead use the \code{\link{metareg}}
#' function for GLMMs which can also be used for continuous
#' covariates.
#' 
#' A prediction interval for the treatment effect of a new study is
#' calculated (Higgins et al., 2009) if arguments \code{prediction}
#' and \code{comb.random} are \code{TRUE}.
#' 
#' R function \code{\link{update.meta}} can be used to redo the
#' meta-analysis of an existing metaprop object by only specifying
#' arguments which should be changed.
#' 
#' For the random effects, the method by Hartung and Knapp (2003) is
#' used to adjust test statistics and confidence intervals if argument
#' \code{hakn = TRUE}.
#' 
#' The DerSimonian-Laird estimate (1986) is used in the random effects
#' model if \code{method.tau = "DL"}. The iterative Paule-Mandel
#' method (1982) to estimate the between-study variance is used if
#' argument \code{method.tau = "PM"}.  Internally, R function
#' \code{paulemandel} is called which is based on R function
#' mpaule.default from R package \bold{metRology} from S.L.R. Ellison
#' <s.ellison at lgc.co.uk>.
#' 
#' If R package \bold{metafor} (Viechtbauer 2010) is installed, the
#' following methods to estimate the between-study variance
#' \eqn{\tau^2} (argument \code{method.tau}) are also available:
#' \itemize{
#' \item Restricted maximum-likelihood estimator (\code{method.tau =
#'   "REML"})
#' \item Maximum-likelihood estimator (\code{method.tau = "ML"})
#' \item Hunter-Schmidt estimator (\code{method.tau = "HS"})
#' \item Sidik-Jonkman estimator (\code{method.tau = "SJ"})
#' \item Hedges estimator (\code{method.tau = "HE"})
#' \item Empirical Bayes estimator (\code{method.tau = "EB"})
#' }
#' For these methods the R function \code{rma.uni} of R package
#' \bold{metafor} is called internally. See help page of R function
#' \code{rma.uni} for more details on these methods to estimate
#' between-study variance.
#' 
#' @return
#' An object of class \code{c("metaprop", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{event, n, studlab, exclude,}{As defined above.}
#' \item{sm, incr, allincr, addincr, method.ci,}{As defined above.}
#' \item{level, level.comb,}{As defined above.}
#' \item{comb.fixed, comb.random,}{As defined above.}
#' \item{hakn, method.tau, tau.preset, TE.tau, null.hypothesis,}{As
#'   defined above.}
#' \item{method.bias, tau.common, title, complab, outclab,}{As defined
#'   above.}
#' \item{byvar, bylab, print.byvar, byseparator, warn}{As defined
#'   above.}
#' \item{TE, seTE}{Estimated (un)transformed proportion and its
#'   standard error for individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{zval, pval}{z-value and p-value for test of treatment effect
#'   for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall (un)transformed
#'   proportion and standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{zval.fixed, pval.fixed}{z-value and p-value for test of
#'   overall effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall (un)transformed
#'   proportion and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{zval.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall effect (random effects
#'   model).}
#' \item{prediction, level.predict}{As defined above.}
#' \item{seTE.predict}{Standard error utilised for prediction
#'   interval.}
#' \item{lower.predict, upper.predict}{Lower and upper limits of
#'   prediction interval.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic Q.}
#' \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
#' \item{pval.Q}{P-value of heterogeneity test.}
#' \item{Q.LRT}{Heterogeneity statistic for likelihood-ratio test
#'   (only if \code{method = "GLMM"}).}
#' \item{df.Q.LRT}{Degrees of freedom for likelihood-ratio test}
#' \item{pval.Q.LRT}{P-value of likelihood-ratio test.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{se.tau}{Standard error of square-root of between-study
#'   variance.}
#' \item{C}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups.}
#' \item{method}{A character string indicating method used for
#'   pooling: \code{"Inverse"}}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn=TRUE}).}
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
#'   missing and \code{hakn=TRUE}.}
#' \item{n.harmonic.mean.w}{Harmonic mean of number of observations in
#'   subgroups (for back transformation of Freeman-Tukey Double
#'   arcsine transformation) - if \code{byvar} is not missing.}
#' \item{event.w}{Number of events in subgroups - if \code{byvar} is
#'   not missing.}
#' \item{n.w}{Number of observations in subgroups - if \code{byvar} is
#'   not missing.}
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
#' \item{C.w}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups - if \code{byvar} is not missing.}
#' \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limti for
#'   heterogeneity statistic H within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I2 within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limti for
#'   heterogeneity statistic I2 within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{incr.event}{Increment added to number of events.}
#' \item{keepdata}{As defined above.}
#' \item{data}{Original data (set) used in function call (if
#'   \code{keepdata=TRUE}).}
#' \item{subset}{Information on subset of original data used in
#'   meta-analysis (if \code{keepdata=TRUE}).}
#' \item{.glmm.fixed}{GLMM object generated by call of
#'   \code{\link[metafor]{rma.glmm}} function (fixed effect model).}
#' \item{.glmm.random}{GLMM object generated by call of
#'   \code{\link[metafor]{rma.glmm}} function (random effects model).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' \item{version.metafor}{Version of R package \bold{metafor} used for
#'   GLMMs.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{update.meta}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{print.meta}}
#' 
#' @references
#' Agresti A & Coull BA (1998):
#' Approximate is better than "exact" for interval estimation of
#' binomial proportions.
#' \emph{The American Statistician},
#' \bold{52}, 119--26
#' 
#' Barendregt JJ, Doi SA, Lee YY, Norman RE, Vos T (2013):
#' Meta-analysis of prevalence.
#' \emph{Journal of Epidemiology and Community Health},
#' \bold{67}, 974--8
#' 
#' Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
#' A basic introduction to fixed-effect and random-effects models for
#' meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{1}, 97--111
#' 
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Edward JM et al. (2006):
#' Adherence to antiretroviral therapy in sub-saharan
#' Africa and North America - a meta-analysis.
#' \emph{Journal of the American Medical Association},
#' \bold{296}, 679--90
#' 
#' Freeman MF & Tukey JW (1950):
#' Transformations related to the angular and the square root.
#' \emph{Annals of Mathematical Statistics},
#' \bold{21}, 607--11
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' Knapp G & Hartung J (2003):
#' Improved tests for a random effects meta-regression with a single
#' covariate.
#' \emph{Statistics in Medicine},
#' \bold{22}, 2693--710
#' 
#' Miller JJ (1978):
#' The inverse of the Freeman-Tukey double arcsine transformation.
#' \emph{The American Statistician},
#' \bold{32}, 138
#' 
#' Newcombe RG (1998):
#' Two-sided confidence intervals for the single proportion:
#' comparison of seven methods.
#' \emph{Statistics in Medicine},
#' \bold{17}, 857--72
#' 
#' Paule RC & Mandel J (1982):
#' Consensus values and weighting factors.
#' \emph{Journal of Research of the National Bureau of Standards},
#' \bold{87}, 377--85
#' 
#' Pettigrew HM, Gart JJ, Thomas DG (1986):
#' The bias and higher cumulants of the logarithm of a binomial
#' variate.
#' \emph{Biometrika},
#' \bold{73}, 425--35
#' 
#' Schwarzer G, Chemaitelly H, Abu-Raddad LJ, Rücker G (2019):
#' Seriously misleading results using inverse of Freeman-Tukey double
#' arcsine transformation in meta-analysis of single proportions.
#' \emph{Research Synthesis Methods}, 1--8.
#' https://doi.org/10.1002/jrsm.1348
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67
#' 
#' Viechtbauer W (2010):
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#'
#' Warton DI, Hui FKC (2011):
#' The arcsine is asinine: the analysis of proportions in ecology.
#' \emph{Ecology},
#' \bold{92}, 3--10
#' 
#' @examples
#' # Meta-analysis using generalised linear mixed model
#' #
#' metaprop(4:1, 10 * 1:4)
#' 
#' # Apply various classic meta-analysis methods to estimate
#' # proportions
#' #
#' m1 <- metaprop(4:1, 10 * 1:4, method = "Inverse")
#' m2 <- update(m1, sm = "PAS")
#' m3 <- update(m1, sm = "PRAW")
#' m4 <- update(m1, sm = "PLN")
#' m5 <- update(m1, sm = "PFT")
#' #
#' m1
#' m2
#' m3
#' m4
#' m5
#' #
#' forest(m1)
#' \dontrun{
#' forest(m2)
#' forest(m3)
#' forest(m3, pscale = 100)
#' forest(m4)
#' forest(m5)
#' }
#' 
#' # Do not back transform results, e.g. print logit transformed
#' # proportions if sm = "PLOGIT" and store old settings
#' #
#' oldset <- settings.meta(backtransf = FALSE)
#' #
#' m6  <- metaprop(4:1, c(10, 20, 30, 40), method = "Inverse")
#' m7  <- update(m6, sm = "PAS")
#' m8  <- update(m6, sm = "PRAW")
#' m9  <- update(m6, sm = "PLN")
#' m10 <- update(m6, sm = "PFT")
#' #
#' forest(m6)
#' \dontrun{
#' forest(m7)
#' forest(m8)
#' forest(m8, pscale = 100)
#' forest(m9)
#' forest(m10)
#' }
#' 
#' # Use old settings
#' #
#' settings.meta(oldset)
#' 
#' # Examples with zero events
#' #
#' m1 <- metaprop(c(0, 0, 10, 10), rep(100, 4), method = "Inverse")
#' m2 <- metaprop(c(0, 0, 10, 10), rep(100, 4), incr = 0.1, method = "Inverse")
#' #
#' summary(m1)
#' summary(m2)
#' #
#' \dontrun{
#' forest(m1)
#' forest(m2)
#' }
#' 
#' # Example from Miller (1978):
#' #
#' death <- c(3, 6, 10, 1)
#' animals <- c(11, 17, 21, 6)
#' #
#' m3 <- metaprop(death, animals, sm = "PFT")
#' forest(m3)
#' 
#' # Data examples from Newcombe (1998)
#' # - apply various methods to estimate confidence intervals for
#' #   individual studies
#' #
#' event <- c(81, 15, 0, 1)
#' n <- c(263, 148, 20, 29)
#' #
#' m1 <- metaprop(event, n, method.ci = "SA", method = "Inverse")
#' m2 <- update(m1, method.ci = "SACC")
#' m3 <- update(m1, method.ci = "WS")
#' m4 <- update(m1, method.ci = "WSCC")
#' m5 <- update(m1, method.ci = "CP")
#' #
#' lower <- round(rbind(NA, m1$lower, m2$lower, NA, m3$lower,
#'                      m4$lower, NA, m5$lower), 4)
#' upper <- round(rbind(NA, m1$upper, m2$upper, NA, m3$upper,
#'                      m4$upper, NA, m5$upper), 4)
#' #
#' tab1 <- data.frame(
#'   scen1 = meta:::formatCI(lower[, 1], upper[, 1]),
#'   scen2 = meta:::formatCI(lower[, 2], upper[, 2]),
#'   scen3 = meta:::formatCI(lower[, 3], upper[, 3]),
#'   scen4 = meta:::formatCI(lower[, 4], upper[, 4]),
#'   stringsAsFactors = FALSE
#'   )
#' names(tab1) <- c("r=81, n=263", "r=15, n=148",
#'                  "r=0, n=20", "r=1, n=29")
#' row.names(tab1) <- c("Simple", "- SA", "- SACC",
#'                      "Score", "- WS", "- WSCC",
#'                      "Binomial", "- CP")
#' tab1[is.na(tab1)] <- ""
#' # Newcombe (1998), Table I, methods 1-5:
#' tab1
#' 
#' # Same confidence interval, i.e. unaffected by choice of summary
#' # measure
#' #
#' print(metaprop(event, n, method.ci = "WS", method = "Inverse"), ma = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "WS"), ma = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "WS"), ma = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "WS"), ma = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "WS"), ma = FALSE)
#' 
#' # Different confidence intervals as argument sm = "NAsm"
#' #
#' print(metaprop(event, n, method.ci = "NAsm", method = "Inverse"), ma = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "NAsm"), ma = FALSE)
#' 
#' # Different confidence intervals as argument backtransf = FALSE.
#' # Accordingly, method.ci = "NAsm" used internally.
#' #
#' print(metaprop(event, n, method.ci = "WS", method = "Inverse"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' 
#' # Same results (printed on original and log scale, respectively)
#' #
#' print(metaprop(event, n, sm = "PLN", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PLN"), ma = FALSE, backtransf = FALSE)
#' # Results for first study (on log scale)
#' round(log(c(0.3079848, 0.2569522, 0.3691529)), 4)
#' 
#' # Print results as events per 1000 observations
#' #
#' print(metaprop(6:8, c(100, 1200, 1000), method = "Inverse"),
#'       pscale = 1000, digits = 1)
#' 
#' @export metaprop


metaprop <- function(event, n, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     method,
                     ##
                     sm = gs("smprop"),
                     ##
                     incr = gs("incr"), allincr = gs("allincr"),
                     addincr = gs("addincr"),
                     method.ci = gs("method.ci"),
                     ##
                     level = gs("level"), level.comb = gs("level.comb"),
                     comb.fixed = gs("comb.fixed"),
                     comb.random = gs("comb.random"),
                     ##
                     hakn = gs("hakn"),
                     method.tau,
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
                     pscale = 1,
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
                     ##
                     byvar, bylab, print.byvar = gs("print.byvar"),
                     byseparator = gs("byseparator"),
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"),
                     ##
                     control = NULL,
                     ...
                     ) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  if (missing(method))
    method <- if (sm == "PLOGIT") "GLMM" else "Inverse"
  else
    method <- setchar(method, c("Inverse", "GLMM"))
  is.glmm <- method == "GLMM"
  ##
  chknull(sm)
  sm <- setchar(sm, .settings$sm4prop)
  ##
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(hakn)
  ##
  if (missing(method.tau))
    method.tau <- if (method == "GLMM") "ML" else gs("method.tau")
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  if (!anyNA(null.effect) | length(null.effect) != 1)
    chknumeric(null.effect, min = 0, max = 1, single = TRUE)
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  ##
  chknumeric(pscale, single = TRUE)
  if (!backtransf & pscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metaprop"
  ##
  chklogical(allincr)
  chklogical(addincr)
  method.ci <- setchar(method.ci, .settings$ci4prop)
  chklogical(warn)
  ##
  if (is.glmm & sm != "PLOGIT")
    stop("Generalised linear mixed models only possible with argument 'sm = \"PLOGIT\"'.")
  ##
  if (is.glmm & method.tau != "ML")
    stop("Generalised linear mixed models only possible with argument 'method.tau = \"ML\"'.")
  
  
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
  ## Catch 'event' and 'n' from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  chknull(event)
  k.All <- length(event)
  ##
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(n)
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
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
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(n, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (by)
    chklength(byvar, k.All, fun)
  ##
  ## Additional checks
  ##
  if (is.glmm) {
    if (!is.null(TE.tau)) {
      if (warn)
        warning("Argument 'TE.tau' not considered for GLMM.")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)) {
      if (warn)
        warning("Argument 'tau.preset' not considered for GLMM.")
      tau.preset <- NULL
    }
  }
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
      data <- data.frame(.event = event)
    else
      data$.event <- event
    ##
    data$.n <- n
    data$.studlab <- studlab
    ##
    data$.incr <- incr
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
    event <- event[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (by)
      byvar <- byvar[subset]
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
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(event, 0)
  chknumeric(n, 0, zero = TRUE)
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must not be larger than number of observations")
  ##
  ## Recode integer as numeric:
  ##
  event <- int2num(event)
  n     <- int2num(n)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                PFT = rep(FALSE, length(event)),
                PAS = rep(FALSE, length(event)),
                PRAW =   event == 0 | (n - event) == 0,
                PLN =    event == 0 | (n - event) == 0,
                PLOGIT = event == 0 | (n - event) == 0)
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
  if (sm == "PFT") {
    TE <- 0.5 * (asin(sqrt(event / (n + 1))) + asin(sqrt((event + 1) / (n + 1))))
    seTE <- sqrt(1 / (4 * n + 2))
    transf.null.effect <- asin(sqrt(null.effect))
  }
  else if (sm == "PAS") {
    TE <- asin(sqrt(event / n))
    seTE <- sqrt(0.25 * (1 / n))
    transf.null.effect <- asin(sqrt(null.effect))
  }
  else if (sm == "PRAW") {
    TE <- event / n
    seTE <- sqrt((event + incr.event) * (n - event + incr.event) /
                 (n + 2 * incr.event)^3)
    transf.null.effect <- null.effect
  }
  else if (sm == "PLN") {
    TE <- log((event + incr.event) / (n + incr.event))
    ## Hartung, Knapp (2001), p. 3880, formula (18):
    seTE <- ifelse(event == n,
                   sqrt(1 / event                - 1 / (n + incr.event)),
                   sqrt(1 / (event + incr.event) - 1 / (n + incr.event))
                   )
    transf.null.effect <- log(null.effect)
  }
  else if (sm == "PLOGIT") {
    TE <- log((event + incr.event) / (n - event + incr.event))
    seTE <- sqrt(1 / (event + incr.event) +
                 1 / ((n - event + incr.event)))
    transf.null.effect <- log(null.effect / (1 - null.effect))
  }
  ##
  ## Calculate confidence intervals
  ##
  NAs <- rep(NA, k.all)
  ##
  if (method.ci == "CP") {
    lower.study <- upper.study <- NAs
    for (i in 1:k.all) {
      if (!is.na(event[i] & !is.na(n[i]))) {
        cint <- binom.test(event[i], n[i], conf.level = level)
        ##
        lower.study[i] <- cint$conf.int[[1]]
        upper.study[i] <- cint$conf.int[[2]]
      }
      else {
        lower.study[i] <- NA
        upper.study[i] <- NA
      }
    }
  }
  ##
  else {
    if (method.ci == "WS")
      ci.study <- ciWilsonScore(event, n, level = level)
    ##
    else if (method.ci == "WSCC")
      ci.study <- ciWilsonScore(event, n, level = level, correct = TRUE)
    ##
    else if (method.ci == "AC")
      ci.study <- ciAgrestiCoull(event, n, level = level)
    ##
    else if (method.ci == "SA")
      ci.study <- ciSimpleAsymptotic(event, n, level = level)
    ##
    else if (method.ci == "SACC")
      ci.study <- ciSimpleAsymptotic(event, n, level = level, correct = TRUE)
    ##
    else if (method.ci == "NAsm")
      ci.study <- ci(TE, seTE, level = level)
    ##
    lower.study <- ci.study$lower
    upper.study <- ci.study$upper
  }
  ##  
  if (method.ci == "NAsm") {
    if (sm == "PLN") {
      lower.study <- exp(lower.study)
      upper.study <- exp(upper.study)
    }
    ##
    else if (sm == "PLOGIT") {
      lower.study <- logit2p(lower.study)
      upper.study <- logit2p(upper.study)
    }
    ##
    else if (sm == "PAS") {
      lower.study <- asin2p(lower.study, value = "lower", warn = FALSE)
      upper.study <- asin2p(upper.study, value = "upper", warn = FALSE)
    }
    ##
    else if (sm == "PFT") {
      lower.study <- asin2p(lower.study, n, value = "lower", warn = FALSE)
      upper.study <- asin2p(upper.study, n, value = "upper", warn = FALSE)
    }
    ##
    lower.study[lower.study < 0] <- 0
    upper.study[upper.study > 1] <- 1
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event[!exclude]) & !is.na(n[!exclude]))
  ##
  if (is.glmm & k > 0) {
    glmm.fixed <- rma.glmm(xi = event[!exclude], ni = n[!exclude],
                           method = "FE", test = ifelse(hakn, "t", "z"),
                           level = 100 * level.comb,
                           measure = "PLO", control = control,
                           ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(event))
  }
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               ##
               sm = sm,
               level = level,
               level.comb = level.comb,
               comb.fixed = comb.fixed,
               comb.random = comb.random,
               ##
               hakn = hakn,
               method.tau = method.tau,
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
               title = title, complab = complab, outclab = outclab,
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (method != "GLMM" & by & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, TE.tau,
                   level.comb, byvar, control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event = event, n = n,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              method.ci = method.ci,
              incr.event = incr.event)
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
  if (is.glmm & k > 0) {
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb,
               null.effect = transf.null.effect)
    ##
    res$method <- "GLMM"
    ##
    res$TE.fixed <- TE.fixed
    res$seTE.fixed <- seTE.fixed
    res$w.fixed <- w.fixed
    res$lower.fixed <- ci.f$lower
    res$upper.fixed <- ci.f$upper
    res$zval.fixed <- ci.f$z
    res$pval.fixed <- ci.f$p
    ##
    if (sum(!exclude) > 1)
      glmm.random <- rma.glmm(xi = event[!exclude], ni = n[!exclude],
                              method = method.tau,
                              test = ifelse(hakn, "t", "z"),
                              level = 100 * level.comb,
                              measure = "PLO", control = control,
                              ...)
    else {
      ##
      ## Fallback to fixed effect model due to small number of studies
      ##
      glmm.random <- glmm.fixed
    }
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.comb,
               null.effect = transf.null.effect)
    ##
    res$w.random <- rep(NA, length(event))
    ##
    res$TE.random <- TE.random
    res$seTE.random <- seTE.random
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$zval.random <- ci.r$z
    res$pval.random <- ci.r$p
    ##
    res$se.tau2 <- NA
    ci.p <- predict.rma(glmm.random, level = 100 * level.predict)
    res$seTE.predict <- NA
    res$lower.predict <- ci.p$cr.lb
    res$upper.predict <- ci.p$cr.ub
    if (is.null(res$lower.predict))
      res$lower.predict <- NA
    if (is.null(res$upper.predict))
      res$upper.predict <- NA
    ##
    res$Q <- glmm.random$QE.Wld
    res$df.Q <- glmm.random$QE.df
    res$pval.Q <- pvalQ(res$Q, res$df.Q)
    ##
    res$Q.LRT      <- glmm.random$QE.LRT
    res$df.Q.LRT   <- res$df.Q
    res$pval.Q.LRT <- pvalQ(res$Q.LRT, res$df.Q.LRT)
    ##
    if (k > 1)
      res$tau <- sqrt(glmm.random$tau2)
    ##
    res$H <- sqrt(glmm.random$H2)
    res$lower.H <- NA
    res$upper.H <- NA
    ##
    res$I2 <- glmm.random$I2 / 100
    res$lower.I2 <- NA
    res$upper.I2 <- NA
    ##
    res$.glmm.fixed  <- glmm.fixed
    res$.glmm.random <- glmm.random
    res$version.metafor <- packageDescription("metafor")$Version
  }
  ##
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$pscale <- pscale
  ##
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
    res$byvar <- byvar
    res$bylab <- bylab
    res$print.byvar <- print.byvar
    res$byseparator <- byseparator
    res$tau.common <- tau.common
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
      if (is.glmm) {
        res <- c(res, subgroup(res, NULL,
                               factor(res$byvar, bylevs(res$byvar)), ...))
        res$tau.resid <- NA
      }
      else {
        res <- c(res, subgroup(res, hcc$tau))
        res$Q.w.random <- hcc$Q
        res$df.Q.w.random <- hcc$df.Q
        res$tau.resid <- hcc$tau
      }
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
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")

  
  res
}
