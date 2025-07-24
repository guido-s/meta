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
#' @param cluster An optional vector specifying which estimates come
#'   from the same cluster resulting in the use of a three-level
#'   meta-analysis model.
#' @param rho Assumed correlation of estimates within a cluster.
#' @param weights A single numeric or vector with user-specified weights.
#' @param weights.common User-specified weights (common effect model).
#' @param weights.random User-specified weights (random effects model).
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"} and
#'   \code{"GLMM"}, can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"PLOGIT"}, \code{"PAS"}, \code{"PFT"}, \code{"PLN"}, or
#'   \code{"PRAW"}) is to be used for pooling of studies, see Details.
#' @param incr A numeric which is added to event number and sample
#'   size of studies with zero or all events, i.e., studies with an
#'   event probability of either 0 or 1. Or a numeric vector with the
#'   continuity correction for each study.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, or \code{"all"}), see Details.
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
#' @param backtransf A logical indicating whether results for
#'   transformed proportions (argument \code{sm != "PRAW"}) should be
#'   back transformed in printouts and plots. If TRUE (default),
#'   results will be presented as proportions; otherwise transformed
#'   proportions will be shown. See Details for presentation of
#'   confidence intervals.
#' @param pscale A numeric defining a scaling factor for printing of
#'   single event probabilities.
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
#' @param hakn Deprecated argument (replaced by 'method.random.ci').
#' @param adhoc.hakn Deprecated argument (replaced by 'adhoc.hakn.ci').
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if estimation problems exist in fitting a GLMM).
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
#' This function provides methods for common effect and random effects
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
#' \item Arcsine transformation (\code{sm = "PAS"})
#' \item Freeman-Tukey Double arcsine transformation (\code{sm = "PFT"})
#' \item Log transformation (\code{sm = "PLN"})
#' \item No transformation (\code{sm = "PRAW"})
#' }
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain the transformed proportions. In printouts and plots these
#' values are back transformed if argument \code{backtransf = TRUE}
#' (default).
#'
#' A generalised linear mixed model (GLMM) - more specific, a random
#' intercept logistic regression model - can be utilised for the
#' meta-analysis of proportions (Stijnen et al., 2010). This is the
#' default method for the logit transformation (argument \code{sm =
#' "PLOGIT"}). Internally, the \code{\link[metafor]{rma.glmm}}
#' function from R package \bold{metafor} is called to fit a GLMM.
#'
#' Classic meta-analysis (Borenstein et al., 2010) utilising the
#' (un)transformed proportions and corresponding standard errors in
#' the inverse variance method is conducted by calling the
#' \code{\link{metagen}} function internally. This is the only
#' available method for all transformations but the logit
#' transformation. The classic meta-analysis model with logit
#' transformed proportions is used by setting argument \code{method =
#' "Inverse"}.
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
#' \subsection{Choice of transformation / meta-analysis method}{
#' 
#' Contradictory recommendations on the use of transformations of
#' proportions have been published in the literature. For example,
#' Barendregt et al. (2013) recommend the use of the Freeman-Tukey
#' double arcsine transformation instead of the logit transformation
#' whereas Warton & Hui (2011) strongly advise to use generalised
#' linear mixed models with the logit transformation instead of the
#' arcsine transformation.
#'
#' Schwarzer et al. (2019) describe seriously misleading results in a
#' meta-analysis with very different sample sizes due to problems with
#' the back-transformation of the Freeman-Tukey transformation which
#' requires a single sample size (Miller, 1978). Accordingly,
#' Schwarzer et al. (2019) also recommend to use GLMMs for the
#' meta-analysis of single proportions, however, admit that individual
#' study weights are not available with this method. Meta-analysts
#' which require individual study weights should consider the inverse
#' variance method with the arcsine or logit transformation.
#'
#' In order to prevent misleading conclusions for the Freeman-Tukey
#' double arcsine transformation, sensitivity analyses using other
#' transformations or using a range of sample sizes should be
#' conducted (Schwarzer et al., 2019).
#' }
#' 
#' \subsection{Continuity correction}{
#'
#' Three approaches are available to apply a continuity correction:
#' \itemize{
#' \item Only studies with a zero cell count (\code{method.incr =
#'   "only0"})
#' \item All studies if at least one study has a zero cell count
#'   (\code{method.incr = "if0all"})
#' \item All studies irrespective of zero cell counts
#'   (\code{method.incr = "all"})
#' }
#' 
#' If the summary measure is equal to "PLOGIT", "PLN", or "PRAW", the
#' continuity correction is applied if a study has either zero or all
#' events, i.e., an event probability of either 0 or 1.
#'
#' By default, 0.5 is used as continuity correction (argument
#' \code{incr}). This continuity correction is used both to calculate
#' individual study results with confidence limits and to conduct
#' meta-analysis based on the inverse variance method. For GLMMs no
#' continuity correction is used. Furthermore, the value of \code{incr} is
#' only considered in the calculation of confidence intervals for individual
#' studies if \code{method.ci = "NAsm"} (see next subsection).
#' }
#' 
#' \subsection{Confidence intervals for individual studies}{
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
#' transformation. Furthermore, the continuity correction 
#'
#' Results will be presented for transformed proportions if argument
#' \code{backtransf = FALSE}. In this case, argument \code{method.ci =
#' "NAsm"} is used, i.e. confidence intervals based on the normal
#' approximation based on the summary measure.
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
#' \subsection{Specify the null hypothesis of test for an overall proportion}{
#'
#' Argument \code{null.effect} can be used to specify the proportion
#' used under the null hypothesis in a test for an overall effect.
#'
#' By default (\code{null.effect = NA}), no hypothesis test is
#' conducted as it is unclear which value is a sensible choice for the
#' data at hand.  An overall proportion of 50\%, for example, could be
#' tested by setting argument \code{null.effect = 0.5}.
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
#' Argument \code{pscale} can be used to rescale proportions, e.g.
#' \code{pscale = 1000} means that proportions are expressed as events
#' per 1000 observations. This is useful in situations with (very) low
#' event probabilities.
#'
#' A prediction interval will only be shown if \code{prediction =
#' TRUE}.
#' }
#' 
#' @return
#' An object of class \code{c("metaprop", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{metacont}}, \code{\link{metagen}},
#'   \code{\link{print.meta}}, \code{\link{forest.meta}}
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
#' Freeman MF & Tukey JW (1950):
#' Transformations related to the angular and the square root.
#' \emph{Annals of Mathematical Statistics},
#' \bold{21}, 607--11
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
#' Pettigrew HM, Gart JJ, Thomas DG (1986):
#' The bias and higher cumulants of the logarithm of a binomial
#' variate.
#' \emph{Biometrika},
#' \bold{73}, 425--35
#' 
#' Schwarzer G, Chemaitelly H, Abu-Raddad LJ, Rücker G (2019):
#' Seriously misleading results using inverse of Freeman-Tukey double
#' arcsine transformation in meta-analysis of single proportions.
#' \emph{Research Synthesis Methods},
#' \bold{10}, 476--83
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
#' m1
#' m2
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
#' lower <- round(logit2p(rbind(NA, m1$lower, m2$lower, NA, m3$lower,
#'   m4$lower, NA, m5$lower)), 4)
#' upper <- round(logit2p(rbind(NA, m1$upper, m2$upper, NA, m3$upper,
#'   m4$upper, NA, m5$upper)), 4)
#' #
#' tab1 <- data.frame(
#'   scen1 = meta:::formatCI(lower[, 1], upper[, 1]),
#'   scen2 = meta:::formatCI(lower[, 2], upper[, 2]),
#'   scen3 = meta:::formatCI(lower[, 3], upper[, 3]),
#'   scen4 = meta:::formatCI(lower[, 4], upper[, 4])
#'   )
#' names(tab1) <- c("r=81, n=263", "r=15, n=148",
#'   "r=0, n=20", "r=1, n=29")
#' row.names(tab1) <- c("Simple", "- SA", "- SACC",
#'   "Score", "- WS", "- WSCC", "Binomial", "- CP")
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
#'   ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "WS"),
#'   ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "WS"),
#'   ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "WS"),
#'   ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "WS"),
#'   ma = FALSE, backtransf = FALSE)
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
#'   pscale = 1000, digits = 1)
#' 
#' @export metaprop


metaprop <- function(event, n, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     cluster = NULL, rho = 0,
                     #
                     weights = NULL,
                     weights.common = weights, weights.random = weights,
                     #
                     method,
                     ##
                     sm = gs("smprop"),
                     ##
                     incr = gs("incr"), method.incr = gs("method.incr"),
                     ##
                     method.ci = gs("method.ci.prop"), level = gs("level"),
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
                     null.effect = NA,
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     pscale = 1,
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
  ## (1) Check and set arguments
  ##
  ##
  
  chknumeric(rho, min = -1, max = 1)
  #
  chknull(sm)
  sm <- setchar(sm, gs("sm4prop"))
  #
  missing.method <- missing(method)
  if (missing.method)
    method <- if (sm == "PLOGIT") "GLMM" else "Inverse"
  else
    method <- setchar(method, gs("meth4prop"))
  is.glmm <- method == "GLMM"
  #
  missing.method.common.ci <- missing(method.common.ci)
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  if (method != "Inverse" & method.common.ci == "IVhet") {
    if (!missing.method.common.ci)
      warning("Argument 'method.common.ci = \"IVhet\"' only available ",
              "if 'method = \"Inverse\".",
              call. = FALSE)
    method.common.ci <- "classic"
  }
  ##
  missing.method.incr <- missing(method.incr)
  method.incr <- setchar(method.incr, gs("meth4incr"))
  ##
  chklevel(level)
  ##
  missing.method.tau <- missing(method.tau)
  if (missing.method.tau)
    method.tau <- if (method == "GLMM") "ML" else gs("method.tau")
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
  if (!anyNA(null.effect) | length(null.effect) != 1)
    chknumeric(null.effect, min = 0, max = 1, length = 1)
  ##
  method.bias <- setmethodbias(method.bias)
  ##
  chklogical(backtransf)
  ##
  chknumeric(pscale, length = 1)
  if (!backtransf & pscale != 1 & !is_untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.",
            call. = FALSE)
    pscale <- 1
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
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metaprop"
  ##
  method.ci <- setchar(method.ci, gs("ci4prop"))
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
    deprecated2(method.random.ci, missing(method.random.ci),
                hakn, missing(hakn),
                warn.deprecated)
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
  addincr <-
    deprecated(method.incr, missing.method.incr, args, "addincr",
               warn.deprecated)
  allincr <-
    deprecated(method.incr, missing.method.incr, args, "allincr",
               warn.deprecated)
  if (missing.method.incr) {
    method.incr <- gs("method.incr")
    ##
    if (is.logical(addincr) && addincr)
      method.incr <- "all"
    else if (is.logical(allincr) && allincr)
      method.incr <- "if0all"
  }
  ##
  addincr <- allincr <- FALSE
  if (method.incr == "all")
    addincr <- TRUE
  else if (method.incr == "if0all")
    allincr <- TRUE
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
  ## Catch 'event' and 'n' from data:
  ##
  event <- catch("event", mc, data, sfsp)
  chknull(event)
  k.All <- length(event)
  ##
  n <- catch("n", mc, data, sfsp)
  chknull(n)
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- catch("incr", mc, data, sfsp)
  chknumeric(incr, min = 0)
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
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  
  chklength(n, k.All, fun)
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
  ##
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
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as ",
            "argument tau.preset is not NULL.",
            call. = FALSE)
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
      data <- data.frame(.event = event)
    else
      data$.event <- event
    ##
    data$.n <- n
    data$.studlab <- studlab
    ##
    data$.incr <- NA
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
    event <- event[subset]
    n <- n[subset]
    studlab <- studlab[subset]
    ##
    cluster <- cluster[subset]
    exclude <- exclude[subset]
    #
    weights.common <- weights.common[subset]
    weights.random <- weights.random[subset]
    #
    if (length(incr) > 1)
      incr <- incr[subset]
    #
    if (by)
      subgroup <- subgroup[subset]
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
    common <- FALSE
    random <- FALSE
    prediction <- FALSE
    overall <- FALSE
    overall.hetstat <- FALSE
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
  ## Check for whole numbers
  ##
  if (method.ci != "NAsm") {
    if (any(!is_wholenumber(event), na.rm = TRUE)) {
      warning("Normal approximation confidence interval ",
              "(argument method.ci = \"NAsm\") used as\n",
              "at least one number of events contains a non-integer value.",
              call. = FALSE)
      method.ci <- "NAsm"
    }
    else if (any(!is_wholenumber(n), na.rm = TRUE)) {
      warning("Normal approximation confidence interval ",
              "(argument method.ci = \"NAsm\") used as\n",
              "at least one sample size contains a non-integer value.",
              call. = FALSE)
      method.ci <- "NAsm"
    }
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
                PLOGIT = event == 0 | (n - event) == 0,
                PAS = rep(FALSE, length(event)),
                PFT = rep(FALSE, length(event)),
                PLN =    event == 0 | (n - event) == 0,
                PRAW =   event == 0 | (n - event) == 0)
  #
  sparse <- any(sel, na.rm = TRUE)
  #
  # No need to add anything to cell counts for arcsine transformation
  #
  if (addincr | method.incr == "user")
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
  if (keepdata) {
    if (missing.subset) {
      data$.incr <- incr.event
    }
    else {
      data$.incr <- NA
      #
      data$.incr[subset] <- incr.event
    }
  }
  
  
  ##
  ##
  ## (8) Calculate results for individual studies
  ##
  ##
  
  if (sm == "PLOGIT") {
    TE <- log((event + incr.event) / (n - event + incr.event))
    seTE <- sqrt(1 / (event + incr.event) +
                 1 / ((n - event + incr.event)))
    transf.null.effect <- p2logit(null.effect)
  }
  else if (sm == "PAS") {
    TE <- asin(sqrt(event / n))
    seTE <- sqrt(1 / (4 * n))
    transf.null.effect <- p2asin(null.effect)
  }
  else if (sm == "PFT") {
    TE <-
      0.5 * (asin(sqrt(event / (n + 1))) + asin(sqrt((event + 1) / (n + 1))))
    seTE <- sqrt(1 / (4 * n + 2))
    transf.null.effect <- p2asin(null.effect)
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
  else if (sm == "PRAW") {
    TE <- event / n
    seTE <- sqrt((event + incr.event) * (n - event + incr.event) /
                 (n + 2 * incr.event)^3)
    transf.null.effect <- null.effect
  }
  ##
  ## Calculate confidence intervals
  ##
  if (method.ci == "CP")
    ci.study <- ciClopperPearson(event, n, level, null.effect)
  else if (method.ci == "WS")
    ci.study <- ciWilsonScore(event, n, level)
  else if (method.ci == "WSCC")
    ci.study <- ciWilsonScore(event, n, level, correct = TRUE)
  else if (method.ci == "AC")
    ci.study <- ciAgrestiCoull(event, n, level)
  else if (method.ci == "SA")
    ci.study <- ciSimpleAsymptotic(event, n, level)
  else if (method.ci == "SACC")
    ci.study <- ciSimpleAsymptotic(event, n, level, correct = TRUE)
  else if (method.ci == "NAsm")
    ci.study <- ci(TE, seTE, level, null.effect = null.effect)
  ##
  lower.study <- ci.study$lower
  upper.study <- ci.study$upper
  ##
  if (method.ci != "NAsm") {
    if (sm == "PLOGIT") {
      lower.study <- p2logit(lower.study)
      upper.study <- p2logit(upper.study)
    }
    ##
    else if (sm == "PAS") {
      lower.study <- p2asin(lower.study)
      upper.study <- p2asin(upper.study)
    }
    ##
    else if (sm == "PFT") {
      lower.ev <- n * lower.study 
      upper.ev <- n * upper.study 
      ##
      lower.study <-
        0.5 * (asin(sqrt(lower.ev / (n + 1))) +
               asin(sqrt((lower.ev + 1) / (n + 1))))
      upper.study <-
        0.5 * (asin(sqrt(upper.ev / (n + 1))) +
               asin(sqrt((upper.ev + 1) / (n + 1))))
    }
    ##
    else if (sm == "PLN") {
      lower.study <- log(lower.study)
      upper.study <- log(upper.study)
    }
  }
  
  
  ##
  ##
  ## (9) Additional checks for three-level model
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
  ## (10) Additional checks for GLMM
  ##
  ##
  
  if (is.glmm) {
    chkglmm(sm, method.tau, method.random.ci, method.predict,
            adhoc.hakn.ci, adhoc.hakn.pi,
            "PLOGIT")
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
  ## (11) Do meta-analysis
  ##
  ##
  
  k <- sum(!is.na(event[!exclude]) & !is.na(n[!exclude]))
  ##
  for (i in seq_along(method.random.ci))
    if (k == 1 & method.random.ci[i] == "HK")
      method.random.ci[i] <- "classic"
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
  if (by & tau.common & !is.glmm)
    hcc <- hetcalc(TE, seTE, method.tau, "", TE.tau,
                   method.I2, level.hetstat, subgroup, control)
  
  
  ##
  ##
  ## (12) Generate R object
  ##
  ##
  
  res <- list(event = event, n = n,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              method.incr = method.incr,
              sparse = sparse,
              method.ci = method.ci,
              incr.event = incr.event)
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
  if (method.ci == "CP") {
    m$statistic <- rep(NA, length(m$statistic))
    m$pval <- ci.study$p
  }    
  else if (method.ci != "NAsm") {
    m$statistic <- rep(NA, length(m$statistic))
    m$pval <- rep(NA, length(m$pval))
  }
  ##
  if (is.glmm) 
    m$method.tau <- method.tau
  #
  if (is.glmm | three.level) {
    m$seTE.hakn.ci <- m$seTE.hakn.adhoc.ci <-
      m$seTE.hakn.pi <- m$seTE.hakn.adhoc.pi <-
        m$seTE.kero <- NA
    ##
    m$text.random <- gsub("(HK)", "(T)", m$text.random, fixed = TRUE)
  }
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Run GLMM and add data
  ##
  if (is.glmm & k > 0) {
    res$method <- "GLMM"
    res$method.random <- "GLMM"
    ##
    list.prop <- list(xi = event[!exclude], ni = n[!exclude], measure = "PLO")
    ##
    use.random <-
      sum(!exclude) > 1 &
      sum(event[!exclude], na.rm = TRUE) > 0 &
      any(event[!exclude] != n[!exclude])
    ##
    res.glmm <-
      runGLMM(list.prop,
              method.tau = method.tau,
              method.random.ci = method.random.ci,
              level = level.ma,
              control = control, use.random = use.random,
              warn = warn)
    ##
    res <- addGLMM(res, res.glmm, method.I2, transf.null.effect)
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
          runGLMM(list.prop,
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
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$pscale <- pscale
  #
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
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    ##
    res$time.e.w <- NULL
    res$time.c.w <- NULL
    res$t.harmonic.mean.w <- NULL
    ##
    res <- setNAwithin(res, res$three.level | is.glmm)
  }
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
